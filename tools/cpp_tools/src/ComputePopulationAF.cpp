/**
 * Variant statistics pipeline for VCF/BCF files.
 * 
 * Computes per-variant statistics including NS, AN, AC, AF, population-specific
 * allele frequencies, heterozygous/homozygous counts, with optional INFO field
 * carryover, region subsetting, and multi-threaded processing support.
 */

#include <cstdint>
#include <cstdlib>

#include <algorithm>
#include <iostream>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

#include <boost/program_options.hpp>
#include "hts.h"
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "bgzf.h"
#include "kstring.h"

#include "aux.h"
#include "TypeSwitcher.h"

namespace po = boost::program_options;
using namespace std;
using namespace aux;
 
/** Command-line argument settings. */
struct Arguments {
    string input_file = "";
    string output_file = "";
    string samples_file = "";
    string region = "";
    vector<string> info_fields;
    int threads = 1; // Number of threads for multi-threaded BGZF/HTSlib
    
    bool is_valid() const {
        return !input_file.empty() && !output_file.empty();
    }
};
 
/** Population lookup and allele counts. */
struct PopulationData {
    string samples;
    vector<int> sample_to_pop_id;  // sample index -> population ID (-1 if invalid)
    vector<string> pop_id_to_name; // population ID -> population name
    vector<vector<unsigned int>> pop_ac_by_id;  // population ID -> allele counts
    vector<unsigned int> pop_an_by_id;          // population ID -> allele number
};
 
/** Per-variant allele statistics. */
struct VariantStats {
    unsigned int ns = 0;
    unsigned int an = 0;
    vector<unsigned int> ac;
    vector<unsigned int> het;
    vector<unsigned int> hom;
    
    explicit VariantStats(int n_alleles) 
        : ac(n_alleles, 0), het(n_alleles, 0), hom(n_alleles, 0) {}
};
 
/** Parse and validate CLI arguments. */
Arguments parse_command_line_arguments(int argc, char* argv[]) {
    Arguments arguments;
    
    po::options_description desc("Options");
    desc.add_options()
        ("help,h", "Produce help message")
        ("in,i", po::value<string>(&arguments.input_file)->required(), 
         "Input VCF/BCF file. Must be indexed using tabix.")
        ("samples,s", po::value<string>(&arguments.samples_file), 
         "Input file with samples. One sample per line.")
        ("region,r", po::value<string>(&arguments.region), 
         "Region to be processed. Must follow <CHR>:<START_BP>-<END_BP> format.")
        ("fields,f", po::value<vector<string>>(&arguments.info_fields)->multitoken(),
         "Whitespace delimited list of INFO column fields to carry over to the output file.")
        ("threads,t", po::value<int>(&arguments.threads)->default_value(1),
         "Number of threads for BGZF (de)compression. Use 1 to disable multithreading.")
        ("out,o", po::value<string>(&arguments.output_file)->required(), 
         "Output file. Compressed using gzip.");

    po::variables_map vm;
    
    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            cout << "For each variant this program computes: NS, AN, AC, AF, population-specific AF, Het, Hom." << endl << endl;
            cout << desc << endl;
            exit(0);
        }
        po::notify(vm);
    } catch (po::error &e) {
        cerr << "Error in command line:" << endl;
        cerr << e.what() << endl;
        exit(1);
    }
    
    return arguments;
}

/** Validate population label. */
bool is_valid_population_label(string pop_label) {
    // Change to lowercase
    transform(pop_label.begin(), pop_label.end(), pop_label.begin(), [](unsigned char c){ return tolower(c); });

    // Check for na, NA, nA, Na, empty and 0
    if (pop_label.empty() || pop_label == "na" || pop_label == "0")
        return false;
    return true;
}

/** Build population lookup tables. */
void initialize_population_data(const string& samples_file, bcf_hdr_t* header, PopulationData& pop_data) {
    if (samples_file.empty()) {
        return;
    }

    cout << "Reading samples file and building lookup... " << flush;

    unordered_map<string, string> sample_to_pop_map;
    pop_data.samples = read_samples(samples_file.c_str(), sample_to_pop_map);

    // Discover unique populations and assign IDs
    set<string> unique_pops; 
    for (const auto& pair : sample_to_pop_map) {
        if (is_valid_population_label(pair.second)) {
            unique_pops.insert(pair.second);
        }
    }

    unordered_map<string, int> pop_name_to_id;
    int pop_id_counter = 0;
    for (const auto& pop_name : unique_pops) {
        pop_data.pop_id_to_name.push_back(pop_name);
        pop_name_to_id[pop_name] = pop_id_counter++;
    }

    // Resize allele count/number vectors based on the number of unique populations
    pop_data.pop_ac_by_id.resize(pop_data.pop_id_to_name.size());
    pop_data.pop_an_by_id.resize(pop_data.pop_id_to_name.size(), 0u);

    // Build sample -> population ID mapping
    pop_data.sample_to_pop_id.resize(header->n[BCF_DT_SAMPLE], -1);
    for (int i = 0; i < header->n[BCF_DT_SAMPLE]; ++i) {
        string sample_name = header->samples[i];
        auto it = sample_to_pop_map.find(sample_name);
        if (it != sample_to_pop_map.end()) {
            auto pop_it = pop_name_to_id.find(it->second);
            if (pop_it != pop_name_to_id.end()) {
                pop_data.sample_to_pop_id[i] = pop_it->second;
            }
        }
    }

    cout << "Done." << endl;
    cout << "Found " << count(pop_data.samples.begin(), pop_data.samples.end(), ',') + 1 << " sample(s) across " << pop_data.pop_id_to_name.size() << " populations." << endl;
}
 
/** Validate requested INFO fields against header. */
vector<string> validate_info_fields(bcf_hdr_t* header, const vector<string>& requested_fields) {
    vector<string> valid_fields;
    
    for (const auto& field : requested_fields) {
        bool found = false;
        for (int i = 0; i < header->nhrec; ++i) {
            if (strcmp(header->hrec[i]->key, "INFO") == 0) {
                if (header->hrec[i]->nkeys > 0 && header->hrec[i]->vals && 
                    strcmp(header->hrec[i]->vals[0], field.c_str()) == 0) {
                    // Skip our computed fields to avoid conflicts
                    if ((strcmp(field.c_str(), "NS") == 0) || (strcmp(field.c_str(), "AN") == 0) ||
                        (strcmp(field.c_str(), "AC") == 0) || (strcmp(field.c_str(), "AF") == 0) ||
                        (strcmp(field.c_str(), "Het") == 0) || (strcmp(field.c_str(), "Hom") == 0)) {
                        continue;
                    }
                    valid_fields.push_back(field);
                    found = true;
                    break;
                }
            }
        }
        if (!found) {
            throw runtime_error("Field '" + field + "' was not found in INFO meta-information lines!");
        }
    }
    
    return valid_fields;
}
 
/** Set up BCF/VCF reader with filters. */
pair<bcf_srs_t*, bcf_hdr_t*> initialize_bcf_reader(const Arguments& arguments, const PopulationData& pop_data) {
    bcf_srs_t *sr = bcf_sr_init();
    if (!sr) {
        throw runtime_error("Error while initializing BCF reader!");
    }
    
    // Region handling
    if ((arguments.input_file.compare("-") != 0) && (!arguments.region.empty())) {
        bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
    }
    if (!arguments.region.empty()) {
        cout << "Setting region... " << flush;
        if (bcf_sr_set_regions(sr, arguments.region.c_str(), 0) < 0) {
            throw runtime_error("Error while subsetting region!");
        }
        cout << "Done." << endl;
        cout << "Set " << arguments.region.c_str() << " region." << endl;
    }
    
    // Add reader
    if (bcf_sr_add_reader(sr, arguments.input_file.c_str()) <= 0) {
        throw runtime_error("Error while initializing VCF/BCF reader!");
    }
    
    bcf_hdr_t* header = bcf_sr_get_header(sr, 0);
    
    // Sample subsetting
    if (!pop_data.samples.empty()) {
        cout << "Setting sample(s)... " << flush;
        if (bcf_hdr_set_samples(header, pop_data.samples.c_str(), 0) != 0) {
            throw runtime_error("Error while subsetting samples!");
        }
        cout << "Done." << endl;
    }
    return make_pair(sr, header);
}
 
/** Write trimmed VCF header with computed INFO tags. */
void write_vcf_header(BGZF* ofp, bcf_hdr_t* header, const PopulationData& pop_data, const vector<string>& info_fields) {
    // Set of INFO IDs computed by this program to skip existing definitions
    unordered_set<string> computed_info_tags = {"NS", "AN", "AC", "AF", "Het", "Hom"};
    for (const auto& pop_name : pop_data.pop_id_to_name) {
        computed_info_tags.insert(pop_name + "_AF");
    }

    // Helper function to write a header record
    auto write_hrec = [&](bcf_hrec_t *hrec) {
        if (!hrec || !hrec->key) {
            cerr << "[warning] Null or invalid header record encountered, skipping." << endl;
            return;
        }

        kstring_t str = {0, 0, 0};
        if (bcf_hrec_format(hrec, &str) != 0) {
            cerr << "[warning] Failed to format header record: ##" << hrec->key << endl;
            if (str.s) free(str.s);
            return;
        }

        string formatted(str.s, str.l);
        if (str.s) free(str.s);

        // Remove internal ,IDX=... attribute if present
        size_t idx_pos = formatted.find(",IDX=");
        if (idx_pos != string::npos) {
            size_t end_pos = formatted.find_first_of(",>", idx_pos + 1);
            if (end_pos != string::npos) {
                formatted.erase(idx_pos, end_pos - idx_pos);
            }
        }

        write(ofp, "%s", formatted.c_str());
    };

    // Emit ##fileformat=VCFv4.2
    write(ofp, "##fileformat=VCFv4.2\n");
    
    // Loop through all hrec in the original header
    for (int i = 0; i < header->nhrec; ++i) {
        bcf_hrec_t *hrec = header->hrec[i];

        // Skip fileformat - we already wrote our own
        if (strcmp(hrec->key, "fileformat") == 0) {
            continue;
        }
        
        // Keep essential VCF metadata lines (contig, FILTER, source, reference, ALT)
        if (strcmp(hrec->key, "contig") == 0 ||
            strcmp(hrec->key, "FILTER") == 0 ||
            strcmp(hrec->key, "source") == 0 ||
            strcmp(hrec->key, "reference") == 0 ||
            strcmp(hrec->key, "ALT") == 0) {
            write_hrec(hrec);
            continue;
        }

        // Handle INFO records: keep only user-requested fields
        if (strcmp(hrec->key, "INFO") == 0) {
            // Extract INFO record ID from first key-value pair
            string id;
            if (hrec->nkeys > 0 && hrec->keys[0] && hrec->vals[0] && 
                strcmp(hrec->keys[0], "ID") == 0) {
                id = hrec->vals[0];
            } else {
                cerr << "[warning] INFO record missing ID field" << endl;
                continue;
            }

            // Skip INFO records that will be recomputed
            if (computed_info_tags.count(id)) {
                continue;
            }

            // Only keep INFO tags explicitly requested by user (default: none)
            if (info_fields.empty() ||
                find(info_fields.begin(), info_fields.end(), id) == info_fields.end()) {
                continue;
            }
            
            cout << "Carrying over INFO field: " << id << endl;
            write_hrec(hrec);
            continue;
        }

        // Skip all other header records (only keep those explicitly allowed above)
    }
    
    // Add INFO definitions for computed statistics
    write(ofp, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage\">\n");
    write(ofp, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage\">\n");
    write(ofp, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
    write(ofp, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
    
    // Add population-specific allele frequency INFO headers 
    for (const auto& pop_name : pop_data.pop_id_to_name) {
        string info_str = "##INFO=<ID=" + pop_name + "_AF,Number=A,Type=Float,Description=\"Alternate Allele Frequencies in " + pop_name + "\">\n";
        write(ofp, info_str.c_str());
    }
    
    write(ofp, "##INFO=<ID=Het,Number=A,Type=Integer,Description=\"Heterozygous Counts\">\n");
    write(ofp, "##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Homozygous Alternate Counts\">\n");
    
    // Output the columns line
    write(ofp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");
}
 
/** Compute per-variant statistics from genotypes. */
VariantStats calculate_variant_stats(bcf1_t* record, bcf_hdr_t* header, PopulationData& pop_data) {
    VariantStats stats(record->n_allele);
    
    // Resize and reset fast population structures
    for (auto& pop_ac : pop_data.pop_ac_by_id) {
        if (pop_ac.size() < (size_t)record->n_allele) {
            pop_ac.resize(record->n_allele, 0u);
        }
        fill(pop_ac.begin(), pop_ac.end(), 0u);
    }
    fill(pop_data.pop_an_by_id.begin(), pop_data.pop_an_by_id.end(), 0u);

    // Find GT format field
    int gt_id = bcf_hdr_id2int(header, BCF_DT_ID, "GT");
    int gt_index = -1;
    
    for (int i = 0; i < record->n_fmt; ++i) {
        if (record->d.fmt[i].id == gt_id) {
            gt_index = i;
            break;
        }
    }
    
    if (gt_index == -1) {
        throw runtime_error("GT field was not found in FORMAT!");
    }
    
    // Initialize TypeSwitcher for GT field
    TypeSwitcher gt_switcher;
    gt_switcher.init(&record->d.fmt[gt_index]);
    
    // Pre-allocate gt_values vector to avoid repeated allocations
    vector<int32_t> gt_values;
    gt_values.reserve(record->n_sample * 2); // Most genotypes are diploid

    // Process each sample using TypeSwitcher
    for (int i = 0; i < record->n_sample; ++i) {
        //  Use a vector to handle any number of alleles
        vector<unsigned int> allele_count(record->n_allele, 0u);
        unsigned int ac_sample = 0;
        
        // Read genotype values for this sample
        (gt_switcher.*(gt_switcher.read))(gt_values);
        
        // Extract alleles for this sample 
        for (auto&& v : gt_values) {
            if (!bcf_gt_is_missing(v)) {
                int allele = bcf_gt_allele(v);
                if (allele >= 0 && allele < record->n_allele) {
                    // Count this allele
                    allele_count[allele]++;
                    stats.ac[allele]++;
                    ac_sample++;
                    
                    // Population-specific counting 
                    int pop_id = pop_data.sample_to_pop_id[i];
                    if (pop_id >= 0 && (size_t)allele < pop_data.pop_ac_by_id[pop_id].size()) {
                        pop_data.pop_ac_by_id[pop_id][allele]++;
                    }
                }
            }
        }
        
        if (ac_sample == 0) {
            gt_values.clear();
            continue;
        }
        
        stats.ns++;
        stats.an += ac_sample;
        
        // Population allele number
        int pop_id = pop_data.sample_to_pop_id[i];
        if (pop_id >= 0) {
            pop_data.pop_an_by_id[pop_id] += ac_sample;
        }
        
        // Het/Hom determination
        int total_distinct_alleles = 0;
        for (size_t k = 0; k < allele_count.size(); ++k) {
            if (allele_count[k] > 0) {
                total_distinct_alleles++;
            }
        }
        
        if (total_distinct_alleles == 1) {
            // Homozygous (ref or alt)
            for (int a = 1; a < record->n_allele; ++a) {
                if (allele_count[a] == ac_sample) { 
                    stats.hom[a]++;
                    break; 
                }
            }
        } else if (total_distinct_alleles > 1) {
            // Heterozygous: increment count for every alternate allele present
            for (int a = 1; a < record->n_allele; ++a) {
                if (allele_count[a] > 0) {
                    stats.het[a]++;
                }
            }
        }
        
        gt_values.clear();
    }
    
    return stats;
}

void write_carryover_info_fields(BGZF* ofp, bcf1_t* record, bcf_hdr_t* header, const vector<string>& info_fields) {
    if (info_fields.empty()) return;
    
    // Unpack INFO if needed
    if ((record->unpacked & BCF_UN_INFO) == 0) {
        bcf_unpack(record, BCF_UN_INFO);
    }
    
    for (const auto& field : info_fields) {
        bcf_info_t *info = bcf_get_info(header, record, field.c_str());
        if (!info) continue;
        
        int header_type = bcf_hdr_id2type(header, BCF_HL_INFO, info->key);
        
        // Special case for flags (no value to extract)
        if (header_type == BCF_HT_FLAG) {
            write(ofp, ";%s", field.c_str());
            continue;
        }
        
        // Handle types with values (int, float, string)
        void *dst = nullptr;
        int ndst = 0;
        int nret = bcf_get_info_values(header, record, field.c_str(), &dst, &ndst, header_type);
        
        if (nret > 0) {
            // Write field name and first value with appropriate format
            write(ofp, ";%s=", field.c_str());
            
            if (header_type == BCF_HT_INT) {
                int32_t *values = (int32_t*)dst;
                write(ofp, "%d", values[0]);
                for (int i = 1; i < nret; ++i) write(ofp, ",%d", values[i]);
            } 
            else if (header_type == BCF_HT_REAL) {
                float *values = (float*)dst;
                write(ofp, "%g", values[0]);
                for (int i = 1; i < nret; ++i) write(ofp, ",%g", values[i]);
            }
            else if (header_type == BCF_HT_STR) {
                write(ofp, "%s", (char*)dst);
            }
            else {
                cerr << "[warning] Unknown INFO field type for: " << field << endl;
            }
        }
        
        if (dst) free(dst);
    }
}

/** Write a complete VCF record with computed stats. */
void write_variant_line(BGZF* ofp, bcf1_t* record, bcf_hdr_t* header, const VariantStats& stats, const PopulationData& pop_data, const vector<string>& info_fields) {
    // Skip monomorphic sites
    unsigned int ac_total = 0;
    for (int i = 1; i < record->n_allele; ++i) {
        ac_total += stats.ac[i];
    }
    if (ac_total == 0) return;
    
    // Basic variant info following existing format
    write(ofp, "%s\t%lu\t%s\t%s", bcf_seqname(header, record), record->pos + 1, record->d.id, record->d.allele[0]);
    
    if (record->n_allele > 1) {
        write(ofp, "\t%s", record->d.allele[1]);
        for (int i = 2; i < record->n_allele; ++i) {
            write(ofp, ",%s", record->d.allele[i]);
        }
    } else {
        write(ofp, "\t.");
    }
    
    // QUAL and FILTER
    if (isnan(record->qual)) {
        write(ofp, "\t.");
    } else {
        write(ofp, "\t%g", record->qual);
    }
    
    if (record->d.n_flt <= 0) {
        write(ofp, "\t.");
    } else {
        write(ofp, "\t%s", bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[0]));
        for (int i = 1; i < record->d.n_flt; ++i) {
            write(ofp, ";%s", bcf_hdr_int2id(header, BCF_DT_ID, record->d.flt[i]));
        }
    }
    
    // INFO fields
    write(ofp, "\tNS=%d;AN=%d", stats.ns, stats.an);
    if (record->n_allele > 1) {
        write(ofp, ";AC=%d", stats.ac[1]);
        for (int i = 2; i < record->n_allele; ++i) {
            write(ofp, ",%d", stats.ac[i]);
        }
        write(ofp, ";AF=%g", stats.an > 0 ? (stats.ac[1] / (double)stats.an) : 0.0);
        for (int i = 2; i < record->n_allele; ++i) {
            write(ofp, ",%g", stats.an > 0 ? (stats.ac[i] / (double)stats.an) : 0.0);
        }
    }
    
    // Population AF
    for (size_t pop_id = 0; pop_id < pop_data.pop_id_to_name.size(); ++pop_id) {
        const string& pop_name = pop_data.pop_id_to_name[pop_id];
        unsigned int pop_an = pop_data.pop_an_by_id[pop_id];
        const auto& pop_ac = pop_data.pop_ac_by_id[pop_id];
        
        if (record->n_allele > 1) {
            string afp_str = ";" + pop_name + "_AF=%g";
            write(ofp, afp_str.c_str(), pop_an > 0 ? (pop_ac[1] / (double)pop_an) : 0.0);
            for (int i = 2; i < record->n_allele; ++i) {
                write(ofp, ",%g", pop_an > 0 ? (pop_ac[i] / (double)pop_an) : 0.0);
            }
        }
    }
    
    // Het/Hom counts
    if (record->n_allele > 1) {
        write(ofp, ";Het=%d", stats.het[1]);
        for (int i = 2; i < record->n_allele; ++i) {
            write(ofp, ",%d", stats.het[i]);
        }
        write(ofp, ";Hom=%d", stats.hom[1]);
        for (int i = 2; i < record->n_allele; ++i) {
            write(ofp, ",%d", stats.hom[i]);
        }
    }
    
    write_carryover_info_fields(ofp, record, header, info_fields);
    
    write(ofp, "\n");
}
 

int main(int argc, char* argv[]) {
    try {
        Arguments arguments = parse_command_line_arguments(argc, argv);
        
        if (!arguments.is_valid()) {
            cerr << "Error: Input and output files required. Use -h for help." << endl;
            return 1;
        }
        
        // Load population data
        PopulationData pop_data;
        if (!arguments.samples_file.empty()) {
            unordered_map<string, string> temp_populations;
            pop_data.samples = read_samples(arguments.samples_file.c_str(), temp_populations);
        }
        
        auto reader_pair = initialize_bcf_reader(arguments, pop_data);
        bcf_srs_t* sr = reader_pair.first;
        bcf_hdr_t* header = reader_pair.second;
        
        initialize_population_data(arguments.samples_file, header, pop_data);

        // Enable multi-threaded decompression if requested
        if (arguments.threads > 1) {
            if (bcf_sr_set_threads(sr, arguments.threads) != 0) {
                cerr << "[warning] Could not enable multi-threaded reading (" << arguments.threads << " threads). Falling back to single-threaded mode." << endl;
            }
        }
        
        // Validate INFO fields to carry over
        vector<string> valid_info_fields = validate_info_fields(header, arguments.info_fields);
        
        // Setup output file
        BGZF *ofp = bgzf_open(arguments.output_file.c_str(), "w");
        if (!ofp) {
            throw runtime_error("Error while opening output file!");
        }
        
        // Enable multi-threaded compression on output if requested
        if (arguments.threads > 1) {
            if (bgzf_mt(ofp, arguments.threads, 256) != 0) { 
                cerr << "[warning] Could not enable multi-threaded compression (" << arguments.threads << " threads). Continuing single-threaded." << endl;
            }
        }
        
        write_vcf_header(ofp, header, pop_data, valid_info_fields);
        
        cout << "Processing... " << flush;
        while (bcf_sr_next_line(sr) > 0) {
            bcf1_t* record = bcf_sr_get_line(sr, 0);
            
            bcf_unpack(record, BCF_UN_ALL);
            VariantStats stats = calculate_variant_stats(record, header, pop_data);
            write_variant_line(ofp, record, header, stats, pop_data, valid_info_fields);
        }
        cout << "Done." << endl;
        
        // Cleanup
        bcf_sr_destroy(sr);
        if (bgzf_close(ofp) != 0) {
            throw runtime_error("Error while closing output file!");
        }
        
    } catch (exception &e) {
        cerr << "Error: " << endl;
        cerr << e.what() << endl;
        return 1;
    }
    
    return 0;
}