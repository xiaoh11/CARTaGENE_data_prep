#include "aux.h"

/* to read samples from list saved in the text file - one sample per line */
string aux::read_samples(const char* samples_file) throw (runtime_error) {
   string line;
   stringstream buffer;
   ifstream if_samples;
   if_samples.exceptions(ifstream::failbit | ifstream::badbit);
   try {
      if_samples.open(samples_file);
      if (getline(if_samples, line)) {
         buffer << line;
      } 
      while (getline(if_samples, line)) {
         buffer << "," << line;
      }
      if_samples.close();
   } catch (exception& e) {
      if (!if_samples.eof()) {
         throw runtime_error("Error while reading samples file! Check if file exists and has read permissions.");
      }
      if_samples.close();
   }
   buffer.flush();
   return buffer.str();
}

/* to read samples from tab-delimited (without header) table saved in the text file */
string aux::read_samples(const char* samples_file, unsigned int col_idx) throw (runtime_error) {
    string line;
    vector<string> tokens;
    auto separator = regex("[ \t]");
    stringstream buffer;
    ifstream if_samples;
    if_samples.exceptions(ifstream::failbit | ifstream::badbit);
    try {
        if_samples.open(samples_file);
        while (getline(if_samples, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
            if (tokens.size() < col_idx) {
                throw runtime_error("Samples file has less columns than required.");
            }
            if (buffer.tellp() > 0) {
                buffer << ",";
            }
            buffer << tokens.at(col_idx - 1);
            tokens.clear();
        }
        if_samples.close();
    } catch (exception& e) {
        if (!if_samples.eof()) {
            throw runtime_error("Error while reading samples file! Check if file exists and has read permissions.");
        }
        if_samples.close();
    }
    buffer.flush();
    return buffer.str();
}

/* HX: to read samples from list saved in the text file - one sample per line. Two space- or tab-delimited columns (no header):  sample name, population label */
string aux::read_samples(const char* samples_file, unordered_map<string, string>& populations) throw (runtime_error) {
    string line;
    vector<string> tokens;
    auto separator = regex("[ \t]");
    stringstream buffer;
    ifstream if_samples;
    if_samples.exceptions(ifstream::failbit | ifstream::badbit);
    try {
        if_samples.open(samples_file);
        while (getline(if_samples, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
            
            // make sure there are exactly two columns per line
            if (tokens.size() != 2) {
                throw runtime_error("Samples file has number of columns other than 2.");
            }

            // check if the key exists before
            if(populations.find(tokens.at(0)) == populations.end()) {
                // if not, create pairing of key and value
                populations[tokens.at(0)] = tokens.at(1);
            } else {
                // if exists, throw an error
                throw runtime_error("Samples file does not satisfied one sample per line.");
            }
            
            if (buffer.tellp() > 0) {
                buffer << ",";
            }
            buffer << tokens.at(0);

            tokens.clear();
        }
        if_samples.close();
    } catch (exception& e) {
        if (!if_samples.eof()) {
            throw runtime_error("Error while reading samples file! Check if file exists and has read permissions.");
        }
        if_samples.close();
    }
    buffer.flush();
    return buffer.str();
}


/* to read sample new and current names from tab-delimited (without header) table saved in the text file */
map<string, string> aux::read_sample_names_map(const char* samples_file, unsigned int name_col_idx, unsigned int new_name_col_idx) throw (runtime_error) {
    string line;
    vector<string> tokens;
    map<string, string> names;
    auto separator = regex("[ \t]");
    ifstream if_samples;
    if_samples.exceptions(ifstream::failbit | ifstream::badbit);
    try {
        if_samples.open(samples_file);
        while (getline(if_samples, line)) {
            copy(sregex_token_iterator(line.begin(), line.end(), separator, -1), sregex_token_iterator(), back_inserter(tokens));
            if ((tokens.size() < name_col_idx) || (tokens.size() < new_name_col_idx)) {
                throw runtime_error("Samples file has less columns than required.");
            }
            names.emplace(make_pair(tokens.at(name_col_idx - 1), tokens.at(new_name_col_idx - 1)));
            tokens.clear();
        }
        if_samples.close();
    } catch (exception& e) {
        if (!if_samples.eof()) {
            throw runtime_error("Error while reading samples file! Check if file exists and has read permissions.");
        }
        if_samples.close();
    }
    return names;
}


void aux::write(BGZF* f, const char* format, ...) throw (runtime_error) {
   va_list arguments;
   long int n = 0;
   const unsigned int BUFFER_SIZE = 32768u;
   unsigned int max_string_length = BUFFER_SIZE - 1u;
   char buffer[BUFFER_SIZE];

   va_start(arguments, format);
   if ((n = vsnprintf(buffer, max_string_length, format, arguments)) < 0) {
      throw runtime_error("Error while writing to memory buffer!");
   } else if (n > max_string_length) {
      throw runtime_error("Too small memory buffer size for writing!");
   }
   va_end(arguments);

   if (bgzf_write(f, buffer, n) < n) {
      throw runtime_error("Error while writing to output file!");
   }
}

