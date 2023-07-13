#include <iostream>
#include <fstream>
#include <cstdint>
#include <cstdlib>
#include <functional>
#include <algorithm>
#include <vector>
#include <set>
#include "hts.h"
#include "vcf.h"
#include "synced_bcf_reader.h"
#include "Histogram.h"
#include "TypeSwitcher.h"
#include "aux.h"

#include <boost/program_options.hpp>

#include <numeric> //HX: to use accumulate()

namespace po = boost::program_options;

using namespace std;
using namespace aux;

int main(int argc, char* argv[]) {
    string input_file("");
    string output_file("");
    string samples_file("");
    string samples("");
    string region("");
    vector<string> info_fields;

    vector<double> hist_borders = { 0, 6, 11, 16, 21, 26, 31, 36, 41, 46, 51, 56, 61, 66, 71, 76, 81, 86, 91, 96, numeric_limits<double>::max() };

    po::options_description desc("Options");

    desc.add_options()
            ("help,h", "Produce help message")
            ("in,i", po::value<string>(&input_file)->required(), "Input VCF/BCF file. Must be indexed using tabix.")
            ("samples,s", po::value<string>(&samples_file), "Input file with samples. One sample per line.")
            ("region,r", po::value<string>(&region), "Region to be processed. Must follow <CHR>:<START_BP>-<END_BP> format.")
            ("fields,f", po::value<vector<string>>(&info_fields)->multitoken(), "Whitespace delimited list of INFO column fields to carry over to the output file.")
            ("out,o", po::value<string>(&output_file)->required(), "Output file. Compressed using gzip.")
            ;

    po::variables_map vm;

    try {
        po::store(po::parse_command_line(argc, argv, desc), vm);
        if (vm.count("help")) {
            cout << "For each variant this program computes: NS, AN, AC, AF, Hom, Het, AVGDP, AVGGQ, histograms of depth (DP) and genotype qualities (GQ)." << endl << endl;
            cout << desc << endl;
            return 0;
        }
        po::notify(vm);
    } catch (po::error &e) {
        cout << "Error in command line:" << endl;
        cout << e.what() << endl;
        return 1;
    }

    try {
        unordered_map<string, string> populations; // HX
        unordered_map<string, vector<unsigned int> > populations_ac; // HX: create map that uses populations as keys and stores vectors of counts as values
        unordered_map<string, unsigned int> populations_an; // HX: create map that uses populations as keys and stores allele numbers as values
        if (!samples_file.empty()) {
            cout << "Reading samples file... " << flush;
            //samples = read_samples(samples_file.c_str()); // HX
            samples = read_samples(samples_file.c_str(), populations); // HX
            //cout << "sample: " << samples << endl;  // HX --> samples is verified to be a string of all sample ids
            //cout << "Population unordered map: " << endl;
            for (const auto& pair : populations) { // HX
                //cout << "Sample: " << pair.first << ", Population: " << pair.second << endl; 
                populations_ac[pair.second] = vector<unsigned int>() ;// pair.second is the populations
                populations_an[pair.second] = 0u ;
            }
            //HX
            cout << "Found these population(s): " << endl;
            for (const auto& pair : populations_ac) { // HX
                cout << "Population: " << pair.first; // << ", Count: ";
                for(const auto& num : pair.second) {
                cout << num << ' ';
                }
                cout << '\n';
            }
            cout << "Done." << endl;
            cout << "Found " << count(samples.begin(), samples.end(), ',')  + 1 << " sample(s)." << endl;
        }

        BGZF *ofp = bgzf_open(output_file.c_str(), "w");
        if (!ofp) {
            throw runtime_error("Error while opening output file!");
        }

        bcf_srs_t *sr = bcf_sr_init();

        if ((input_file.compare("-") != 0) && (!region.empty())) {
            bcf_sr_set_opt(sr, BCF_SR_REQUIRE_IDX);
        }

        if (!region.empty()) {
            cout << "Setting region... " << flush;
            if (bcf_sr_set_regions(sr, region.c_str(), 0) < 0) {
                throw runtime_error("Error while subsetting region!");
            }
            cout << "Done." << endl;
            cout << "Set " << region.c_str() << " region." << endl;
        }

        if (bcf_sr_add_reader(sr, input_file.c_str()) <= 0) {
            throw runtime_error("Error while initializing VCF/BCF reader!");
        }

        bcf_hdr_t* header = bcf_sr_get_header(sr, 0);

        if (!samples.empty()) {
            cout << "Setting sample(s)... " << flush;
            //HX
            //cout << "Number of samples before:" << header->n[2] << endl;

            if (bcf_hdr_set_samples(header, samples.c_str(), 0) != 0) {
                throw runtime_error("Error while subsetting samples!");
            }

            //HX
            // cout << "Number of samples after:" << header->n[2] << endl;
            // for (int i = 0; i < header->n[2]; ++i) {
            // const char* sample_name = bcf_hdr_int2id(header, BCF_DT_SAMPLE, i);
            // std::cout << "Sample " << i+1 << ": " << sample_name << std::endl;
            // }
            
            cout << "Done." << endl;
        }

        int gt_id = bcf_hdr_id2int(header, BCF_DT_ID, "GT");
        int dp_id = bcf_hdr_id2int(header, BCF_DT_ID, "DP");
        int gq_id = bcf_hdr_id2int(header, BCF_DT_ID, "GQ");
        int ad_id = bcf_hdr_id2int(header, BCF_DT_ID, "AD"); //HX

        if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, gt_id) == 0) {
            throw runtime_error("GT field was not found!");
        }


        //HX

        /*
        if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, dp_id) == 0) {
            throw runtime_error("DP field was not found!");
        }
        */
        if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, dp_id) == 0) {
            //HX: if cannot find dp, try to use info. of ad
            cerr << "DP field was not found, try to find the AD field" << endl;
            if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, ad_id) == 0) {
                throw runtime_error("Neither DP nor AD field were found!");
            } else {
                cout << "AD field has been found and will be used to calculate DP" << endl;
            }
        }



        if (bcf_hdr_idinfo_exists(header, BCF_HL_FMT, gq_id) == 0) {
            throw runtime_error("GQ field was not found!");
        }

        vector<string> input_info_fields;

        write(ofp, "##fileformat=VCFv4.2\n");
        for (int i = 0; i < header->nhrec; ++i) {
            if (strcmp(header->hrec[i]->key, "FILTER") == 0) {
                write(ofp, "##%s=<%s=%s", header->hrec[i]->key, header->hrec[i]->keys[0], header->hrec[i]->vals[0]);
                for (int j = 1; j < header->hrec[i]->nkeys; ++j) {
                    if (strcmp(header->hrec[i]->keys[j], "IDX") == 0) {
                        continue;
                    }
                    write(ofp, ",%s=%s", header->hrec[i]->keys[j], header->hrec[i]->vals[j]);
                }
                write(ofp, ">\n");
            }
            if (strcmp(header->hrec[i]->key, "INFO") == 0) {
                if ((strcmp(header->hrec[i]->vals[0], "NS") == 0) || (strcmp(header->hrec[i]->vals[0], "AN") == 0) ||
                    (strcmp(header->hrec[i]->vals[0], "AC") == 0) || (strcmp(header->hrec[i]->vals[0], "AF") == 0) ||
                    (strcmp(header->hrec[i]->vals[0], "Het") == 0) || (strcmp(header->hrec[i]->vals[0], "Hom") == 0) ||
                    (strcmp(header->hrec[i]->vals[0], "DP") == 0) || (strcmp(header->hrec[i]->vals[0], "AVGDP") == 0) ||
                    (strcmp(header->hrec[i]->vals[0], "AVGDP_R") == 0) || (strcmp(header->hrec[i]->vals[0], "AVGGQ") == 0) ||
                    (strcmp(header->hrec[i]->vals[0], "AVGGQ_R") == 0) || (strcmp(header->hrec[i]->vals[0], "DP_HIST") == 0) ||
                    (strcmp(header->hrec[i]->vals[0], "DP_HIST_R") == 0) || (strcmp(header->hrec[i]->vals[0], "GQ_HIST") == 0) ||
                    (strcmp(header->hrec[i]->vals[0], "GQ_HIST_R") == 0) || (strcmp(header->hrec[i]->vals[0], "AD") == 0)) { //HX: do i need to add AD here?
                    continue;
                }
                if (find(info_fields.begin(), info_fields.end(), header->hrec[i]->vals[0]) != info_fields.end()) {
                    input_info_fields.emplace_back(header->hrec[i]->vals[0]);
                    cout << "Carrying over '" << header->hrec[i]->vals[0] << "' INFO field." << endl;
                    write(ofp, "##%s=<%s=%s", header->hrec[i]->key, header->hrec[i]->keys[0], header->hrec[i]->vals[0]);
                    for (int j = 1; j < header->hrec[i]->nkeys; ++j) {
                        if (strcmp(header->hrec[i]->keys[j], "IDX") == 0) {
                            continue;
                        }
                        write(ofp, ",%s=%s", header->hrec[i]->keys[j], header->hrec[i]->vals[j]);
                    }
                    write(ofp, ">\n");
                }
            }
        }

        for (auto&& field : info_fields) {
            if (find(input_info_fields.begin(), input_info_fields.end(), field) == input_info_fields.end()) {
                throw runtime_error("Field '" + string(field) + "' was not found in INFO meta-information lines!");
            }
        }

        write(ofp, "##INFO=<ID=NS,Number=1,Type=Integer,Description=\"Number of Samples With Coverage\">\n");
        write(ofp, "##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Number of Alleles in Samples with Coverage\">\n");
        write(ofp, "##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Alternate Allele Counts in Samples with Coverage\">\n");
        write(ofp, "##INFO=<ID=AF,Number=A,Type=Float,Description=\"Alternate Allele Frequencies\">\n");
        
        //HX loop through your population_ac and
        for (const auto& pair : populations_ac) { 
            if (pair.first != "NA") {
                string info_str = "##INFO=<ID=" + pair.first + "_AF,Number=A,Type=Float,Description=\"Alternate Allele Frequencies in " + pair.first + "\">\n";
                write(ofp, info_str.c_str());
            }    
        }

        write(ofp, "##INFO=<ID=Het,Number=A,Type=Integer,Description=\"Heterozygous Counts\">\n");
        write(ofp, "##INFO=<ID=Hom,Number=A,Type=Integer,Description=\"Homozygous Alternate Counts\">\n");
        write(ofp, "##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Total depth at site calculated by summing the Allelic depth\">\n");
        write(ofp, "##INFO=<ID=AVGDP,Number=1,Type=Float,Description=\"Average depth per sample\">\n");
        write(ofp, "##INFO=<ID=AVGDP_R,Number=R,Type=Float,Description=\"Average depth per sample carrying allele\">\n");
        write(ofp, "##INFO=<ID=AVGGQ,Number=1,Type=Float,Description=\"Average genotype quality per sample\">\n");
        write(ofp, "##INFO=<ID=AVGGQ_R,Number=R,Type=Float,Description=\"Average genotype quality per sample carrying allele\">\n");
        write(ofp, "##INFO=<ID=DP_HIST,Number=1,Type=String,Description=\"Histogram of DP across all samples; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
        write(ofp, "##INFO=<ID=DP_HIST_R,Number=R,Type=String,Description=\"Histograms of DP across samples carrying allele; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
        write(ofp, "##INFO=<ID=GQ_HIST,Number=1,Type=String,Description=\"Histogram of GQ across all samples; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
        write(ofp, "##INFO=<ID=GQ_HIST_R,Number=R,Type=String,Description=\"Histograms of GQ across samples carrying allele; Mids: 2.5|7.5|12.5|17.5|22.5|27.5|32.5|37.5|42.5|47.5|52.5|57.5|62.5|67.5|72.5|77.5|82.5|87.5|92.5|97.5\">\n");
        write(ofp, "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n");

        int gt_index, dp_index, gq_index, ad_index; //HX
        TypeSwitcher gt_switcher, dp_switcher, gq_switcher, ad_switcher; //HX
        unsigned int ns, an, ac_sample, ac_total; //HX
        set<unsigned int> hom_sample;
        vector<unsigned int> ac, hom, het; // HX: here is where ac is defined
        int allele = 0;
        vector<int32_t> gt_values, dp_values, gq_values, ad_values; //HX
        vector<Histogram> dp_histograms, gq_histograms;
        set<int> unique_alleles;

        cout << "Processing... " << flush;
        while (bcf_sr_next_line(sr) > 0) {
            bcf1_t* rec = bcf_sr_get_line(sr, 0);

            if ((sr->streaming == 0) && (rec->pos < sr->regions->start)) {
                continue;
            }

            if ((rec->unpacked & BCF_UN_FMT) == 0) {
                bcf_unpack(rec, BCF_UN_FMT);
            }

            if ((rec->unpacked & BCF_UN_FLT) == 0) {
                bcf_unpack(rec, BCF_UN_FLT);
            }

            if (!input_info_fields.empty() && ((rec->unpacked & BCF_UN_INFO) == 0)) {
                bcf_unpack(rec, BCF_UN_INFO);
            }

            gt_index = -1;
            dp_index = -1;
            gq_index = -1;
            ad_index = -1; //HX

            for (int i = 0; i < rec->n_fmt; ++i) {
                if (rec->d.fmt[i].id == gt_id) {
                    gt_index = i;
                } else if (rec->d.fmt[i].id == dp_id) {
                    dp_index = i;
                } else if (rec->d.fmt[i].id == gq_id) {
                    gq_index = i;
                    //cout << "gq_index: " << gq_index << endl;
                } else if (rec->d.fmt[i].id == ad_id) {
                    ad_index = i;
                    //cout << "ad_index: " << ad_index << endl; //HX: ad_idex was found to be 2
                }
            }

            if (gt_index == -1) {
                throw runtime_error("GT field was not found in FORMAT!");
            }
            gt_switcher.init(&rec->d.fmt[gt_index]);


            //HX
            
            // if (dp_index == -1) {
            //     cerr << "[warning] DP field missing for " << bcf_seqname(header, rec) << ":" << rec->pos + 1 << ":" << rec->d.allele[0] << "/" << rec->d.allele[1];
            //     for (int i = 2; i < rec->n_allele; ++i) {
            //         cerr << "/" << rec->d.allele[i];
            //     }
            //     cerr << endl;
            // } else {
            //     dp_switcher.init(&rec->d.fmt[dp_index]);
            // }
            
            if (dp_index == -1)
            {
                cerr << "DP field missing, try to use AD field to calculate DP field..." << endl;
                if (ad_index == -1) {
                    cerr << "[warning] DP and AD field missing for " << bcf_seqname(header, rec) << ":" << rec->pos + 1 << ":" << rec->d.allele[0] << "/" << rec->d.allele[1];
                    for (int i = 2; i < rec->n_allele; ++i) {
                        cerr << "/" << rec->d.allele[i];
                    }
                    cerr << endl;
                } else {
                    ad_switcher.init(&rec->d.fmt[ad_index]);
                }
            } else {
                dp_switcher.init(&rec->d.fmt[dp_index]);
            }
            


        //    if (ad_index == -1) {
        //         cerr << "[warning] AD field missing for " << bcf_seqname(header, rec) << ":" << rec->pos + 1 << ":" << rec->d.allele[0] << "/" << rec->d.allele[1];
        //         for (int i = 2; i < rec->n_allele; ++i) {
        //             cerr << "/" << rec->d.allele[i];
        //         }
        //         cerr << endl;
        //     } else {
        //         ad_switcher.init(&rec->d.fmt[ad_index]);
        //     }


            if (gq_index == -1) {
                cerr << "[warning] GQ field missing for " << bcf_seqname(header, rec) << ":" << rec->pos + 1 << ":" << rec->d.allele[0] << "/" << rec->d.allele[1];
                for (int i = 2; i < rec->n_allele; ++i) {
                    cerr << "/" << rec->d.allele[i];
                }
                cerr << endl;
            } else {
                gq_switcher.init(&rec->d.fmt[gq_index]);
            }

            ns = 0u;
            an = 0u;
            fill(ac.begin(), ac.end(), 0u); // cleanup allele, hom and het counts

            // HX: here you loop over populations_ac
            // for &pop1, &pop_ac in population_ac {
            //     fill(pop_ac.begin(), pop_ac.end(), 0u);
            // }
            for(auto& pair : populations_ac) {
                fill(pair.second.begin(), pair.second.end(), 0u);
                populations_an[pair.first] = 0u ;
            }
            

            

            fill(hom.begin(), hom.end(), 0u);
            fill(het.begin(), het.end(), 0u);
            if (rec->n_allele > ac.size()) { // append additional allele, hom and het counts if needed
                ac.resize(rec->n_allele, 0u);
                // HX: here in similar way resize you pop ac
                for(auto& pair : populations_ac) {
                    pair.second.resize(rec->n_allele, 0u);
                }

                hom.resize(rec->n_allele, 0u);
                het.resize(rec->n_allele, 0u);
            }

            for (auto&& h : dp_histograms) { // cleanup DP histograms and reuse them
                h.clear();
            }
            for (int i = dp_histograms.size(); i < rec->n_allele; ++i) { // append additional DP histograms if needed
                dp_histograms.emplace_back(Histogram(hist_borders));
            }

            for (auto &&h: gq_histograms) { // cleanup DP histohgrams and reuse them
                h.clear();
            }
            for (int i = gq_histograms.size(); i < rec->n_allele; ++i) { // append additional GQ histograms if needed
                gq_histograms.emplace_back(Histogram(hist_borders));
            }

            Histogram dp_histogram(hist_borders), gq_histogram(hist_borders);

            for (int i = 0; i < rec->n_sample; ++i) {
                //HX
                // cout << "Iteration " << i << " --- ";
                // const char* sample_name_test = ;
                // cout << "Sample " << i+1 << ": " << sample_name_test << endl;
                // cout << "Population: " << const char* pop = populations[header->samples[i]] << endl;
                // we can call the population sequencially using the index i:
                const char* pop = populations[header->samples[i]].c_str(); // set pop as the population of current sample, here populations[header->samples[i]] is a std::string
                // cout << populations[header->samples[i]] << endl;



                (gt_switcher.*(gt_switcher.read))(gt_values);
                if (dp_index != -1) (dp_switcher.*(dp_switcher.read))(dp_values); 
                if (gq_index != -1) {
                    (gq_switcher.*(gq_switcher.read))(gq_values);
                    // cout << "gq_values is: ";
                    // for (const auto& gqs : gq_values){
                    //     cout << gqs << " ";
                    // }
                    // cout << endl;
                    
                    }
                //HX: try this to output the sum of ad_values
                if (ad_index != -1){
                    (ad_switcher.*(ad_switcher.read))(ad_values);
                    //HX: try to cout the elements in the ad_values vector if not missing
                    if (ad_values[0] != bcf_int32_missing){
                        // cout << "i is: " << i << " | ";
                        // cout << "ad_values is: ";
                        // for (const auto& ads : ad_values){
                        //     cout << ads << " ";
                        // }
                        // cout << "| sum of ad_values is: " << accumulate(ad_values.begin(), ad_values.end(), 0) << endl;
                    } else {
                        cout << "ad_values are missing" << endl;
                    }
                }

                ac_sample = 0u;
                hom_sample.clear();
                unique_alleles.clear();
                for (auto&& v : gt_values) {
                    if (!bcf_gt_is_missing(v)) {
                        allele = bcf_gt_allele(v);
                        unique_alleles.insert(allele);
                        ac[allele] += 1u;
                        // HX
                        // populations_ac[populations[header->samples[i]]][allele] += 1u;
                        populations_ac[pop][allele] += 1u;
                        // pop_ac[pop][allele] += 1u; // something like that


                        ++ac_sample;
                        hom_sample.insert(allele);
                    }
                }
                // HX:
                

                if (ac_sample == 0u) { // all alleles missing for this sample
                    continue;
                }
 
                ++ns;
                an += ac_sample;
                populations_an[pop] += ac_sample;

                if (hom_sample.size() == 1) { // homozygous for some allele
                    hom[*hom_sample.begin()] += 1;
                } else if (hom_sample.size() > 1) {
                    if (hom_sample.erase(0) != 0) { // heterozygous with 1 REF and 1 ALT
                        het[*hom_sample.begin()] += 1;
                    }
                }

                if (dp_index != -1 && dp_values[0] != bcf_int32_missing) {
                    dp_histogram.add((int)dp_values[0]);
                } else {
                    if (ad_index != -1 && ad_values[0] != bcf_int32_missing) {
                        //cout << "| sum of ad_values is: " << accumulate(ad_values.begin(), ad_values.end(), 0) << endl;
                        dp_histogram.add((int)accumulate(ad_values.begin(), ad_values.end(), 0)); //HX
                    }
                    
                }
                

                if (gq_index != -1 && gq_values[0] != bcf_int32_missing) {
                    gq_histogram.add((int)gq_values[0]);
                }

                for (auto&& allele: unique_alleles) {
                    if (dp_index != -1 && dp_values[0] != bcf_int32_missing) {
                        dp_histograms[allele].add((int)dp_values[0]);
                    } else {
                        if (ad_index != -1 && ad_values[0] != bcf_int32_missing) {
                            dp_histograms[allele].add((int)accumulate(ad_values.begin(), ad_values.end(), 0)); //HX
                        }
                    }


                    if (gq_index != -1 && gq_values[0] != bcf_int32_missing) {
                        gq_histograms[allele].add((int)gq_values[0]);
                    }
                }

                gt_values.clear();
                dp_values.clear();
                gq_values.clear();
                ad_values.clear(); //HX
            }
            //HX
            
            cout << "Summary for each variant:" << endl;
            for (const auto& pair : populations_ac) { // HX
                cout << "Population: " << pair.first << ", Count: ";
                for(const auto& num : pair.second) {
                cout << num << ' ';
                }
                cout << " AN: " << populations_an[pair.first];
                cout << '\n';
            }
            cout << "Allele Count: ";
            for(const auto& num : ac) {
                cout << num << ' ';
                }
                cout << '\n';
            cout << "Allele Number: " << an << endl;





            //HX
            // cout << "dp_histogram is: ";
            // for (const auto& dps : dp_histograms){
            //                 cout << dps. << " ";
            //             }
            // cout << endl;


            ac_total = ac[1];
            for (int i = 2; i < rec->n_allele; ++i) {
                ac_total += ac[i];
            }
            if (ac_total == 0) {
                continue;
            }

            write(ofp, "%s\t%lu\t%s\t%s\t%s", bcf_seqname(header, rec), rec->pos + 1, rec->d.id, rec->d.allele[0], rec->d.allele[1]);
            for (int i = 2; i < rec->n_allele; ++i) {
                write(ofp, ",%s", rec->d.allele[i]);
            }
            if (isnan(rec->qual)) {
                write(ofp, "\t.");
            } else {
                write(ofp, "\t%g", rec->qual);
            }
            if (rec->d.n_flt <= 0) {
                write(ofp, "\t.");
            } else {
                write(ofp, "\t%s", bcf_hdr_int2id(header, BCF_DT_ID, rec->d.flt[0]));
                for (int i = 1; i < rec->d.n_flt; ++i) {
                    write(ofp, ";%s", bcf_hdr_int2id(header, BCF_DT_ID, rec->d.flt[i]));
                }
            }
            write(ofp, "\tNS=%d;AN=%d", ns, an);
            write(ofp, ";AC=%d", ac[1]);
            for (int i = 2; i < rec->n_allele; ++i) {
                write(ofp, ",%d", ac[i]);
            }
            write(ofp, ";AF=%g", an > 0 ? (ac[1] / (double)an) : 0.0);
            for (int i = 2; i < rec->n_allele; ++i) {
                write(ofp, ",%g", an > 0 ? (ac[i] / (double)an) : 0.0);
            }
            // {POP}_AF; but when POP == NA 
            for (const auto& pair : populations_ac) { 
                if (pair.first != "NA") {
                    string afp_str = ";" + pair.first + "_AF=%g";
                    write(ofp, afp_str.c_str(),  an > 0 ? (populations_ac[pair.first][1] / (double)populations_an[pair.first]) : 0.0);
                    for (int i = 2; i < rec->n_allele; ++i) {
                        write(ofp, ",%g", an > 0 ? (populations_ac[pair.first][i] / (double)populations_an[pair.first]) : 0.0);
                    }
                }    
            }



            write(ofp, ";Het=%d", het[1]);
            for (int i = 2; i < rec->n_allele; ++i) {
                write(ofp, ",%d", het[i]);
            }
            write(ofp, ";Hom=%d", hom[1]);
            for (int i = 2; i < rec->n_allele; ++i) {
                write(ofp, ",%d", hom[i]);
            }
            if (dp_index != -1 || ad_index != -1) {
                write(ofp, ";DP=%d", (long long int)dp_histogram.get_total());
                write(ofp, ";AVGDP=%g", dp_histogram.get_average());
                write(ofp, ";AVGDP_R=%g", dp_histograms[0].get_average());
                for (int i = 1; i < rec->n_allele; ++i) {
                    write(ofp, ",%g", dp_histograms[i].get_average());
                }
            }
            if (gq_index != -1) {
                write(ofp, ";AVGGQ=%g", gq_histogram.get_average());
                write(ofp, ";AVGGQ_R=%g", gq_histograms[0].get_average());
                for (int i = 1; i < rec->n_allele; ++i) {
                    write(ofp, ",%g", gq_histograms[i].get_average());
                }
            }
            if (dp_index != -1 || ad_index != -1) {
                write(ofp, ";DP_HIST=%s", dp_histogram.get_text());
                write(ofp, ";DP_HIST_R=%s", dp_histograms[0].get_text());
                for (int i = 1; i < rec->n_allele; ++i) {
                    write(ofp, ",%s", dp_histograms[i].get_text());
                }
            }
            if (gq_index != -1) {
                write(ofp, ";GQ_HIST=%s", gq_histogram.get_text());
                write(ofp, ";GQ_HIST_R=%s", gq_histograms[0].get_text());
                for (int i = 1; i < rec->n_allele; ++i) {
                    write(ofp, ",%s", gq_histograms[i].get_text());
                }
            }
            if (!input_info_fields.empty()) {
                for (auto&& field : input_info_fields) {
                    bcf_info_t *info = bcf_get_info(header, rec, field.c_str());
                    int header_type = bcf_hdr_id2type(header, BCF_HL_INFO, info->key);
                    if ((header_type != BCF_HT_REAL) || (info->type != BCF_BT_FLOAT)) {
                        throw runtime_error("This version can carry over only float INFO fields.");
                    }
                    void *dst = nullptr;
                    int ndst = 0;
                    if (bcf_get_info_values(header, rec, field.c_str(), &dst, &ndst, header_type) <= 0) {
                        throw runtime_error("Error while writing '" + string(field) + "' value!");
                    }
                    write(ofp, ";%s=%g", field.c_str(), ((float*)dst)[0]);
                    for (int i = 1; i < ndst; ++i) {
                        write(ofp, ",%g", ((float*)dst)[i]);
                    }
                    if (dst != nullptr) {
                        free(dst);
                    }
                }
            }
            write(ofp, "\n");
        }
        cout << "Done." << endl;

        bcf_sr_destroy(sr);

        if (bgzf_close(ofp) != 0) {
            throw runtime_error("Error while closing output file!");
        }
    } catch (exception &e) {
        cout << "Error: " << endl;
        cout << e.what() << endl;
        return 1;
    }

    return 0;
}
