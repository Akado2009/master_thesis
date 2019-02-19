#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <map>
#include <typeinfo>
#include <iomanip>

struct Data {
    int patho_one;
    int patho_more;
    int benign_one;
    int benign_more;
    int benign;
    int pathogenic;
};

struct Mutation {
    int position;
    char from;
    char to;

    Mutation(int position, char from, char to): position(position), from(from), to(to){}
};

Mutation* construct_mutation(int position, char from, char to) {
    Mutation* mutation = new Mutation(position, from, to);
    return mutation;
}

std::ostream& operator<<(std::ostream& os, const Mutation& mutation) {
    return os << "Mutation: position " << mutation.position << " from " << mutation.from << " to " << mutation.to << std::endl;
}

struct Transcript {
    int exon_count;
    int cds_start;
    int cds_end;
    std::vector<int> starts;
    std::vector<int> stops;

    Transcript(int exon_count, int cds_start, int cds_end, std::vector<int> starts, std::vector<int> stops) : exon_count(exon_count), cds_start(cds_start), cds_end(cds_end), starts(starts), stops(stops){}
};


Transcript* construct_transcript(int exon_count, int cds_start, int cds_end, std::vector<int> starts, std::vector<int> stops) {
    Transcript* transcript = new Transcript(exon_count, cds_start, cds_end, starts, stops);
    return transcript;
}

std::ostream& operator<<(std::ostream& os, const Transcript& transcript) {
    return os << "Transcript: cds_start " << transcript.cds_start << ", cds_end " << transcript.cds_end << ", exonCount " << transcript.exon_count << std::endl;
}

std::vector<std::string> split_string(std::string &input_string, char delimiter) {
    std::vector<std::string> result;
    std::istringstream line_stream(input_string);
    std::string substring;
    while (std::getline(line_stream, substring, delimiter))
        result.push_back(substring);
    return result;
}

bool is_valid_mutation(std::string mutation) {

    if (mutation.find("=") == std::string::npos && mutation.find("p.") != std::string::npos && mutation.find("?") == std::string::npos && mutation.find(" ") != std::string::npos && mutation.find("del") == std::string::npos && mutation.find("*") == std::string::npos && mutation.find("))") == std::string::npos) {
        return true;
    } else {
        return false;
    }
}

std::map<std::string, std::vector<Mutation* > > parse_mutations_file(std::string filename, bool aminoacid, Data* data, char flag) {
    std::map<std::string, std::vector<Mutation* > > result;
    std::ifstream mutations(filename);
    for(std::string line; std::getline(mutations, line);) {
        if (line.rfind("Chromosome", 0) != 0) {
            std::vector <std::string> parsed_line = split_string(line, '\t');
            std::string name = parsed_line[2];
            if (result.find(name) == result.end()) {
                result[name] = std::vector<Mutation *>{};
            }
            result[name].push_back(construct_mutation(std::stoi(parsed_line[7]), parsed_line[8][0], parsed_line[9][0]));

            if (flag == 'p') data->pathogenic++;
            else data->benign++;
        }
    }
    mutations.close();
    return result;
}

std::vector<int> parse_exons(std::string exons){
    std::vector<int> exon_coords;
    std::vector<std::string> parsed_exons = split_string(exons, ',');
    for(auto coordinate: parsed_exons) {
        exon_coords.push_back(std::stoi(coordinate));
    }
    return exon_coords;
}

std::map<std::string, Transcript* > parse_support_table (std::string filename) {
    std::map<std::string, Transcript* > result;
    std::ifstream support_table(filename);
    for(std::string line; std::getline(support_table, line);) {
//#name	chrom	strand	cds_start	cds_end	exonCount	exonStarts	exonEnds
        if (line.rfind("#name", 0) != 0) {
            std::vector<std::string> parsed_line = split_string(line, '\t');
            char strand = parsed_line[2][0];
            std::vector<int> starts = parse_exons(parsed_line[6]);
            std::vector<int> stops = parse_exons(parsed_line[7]);
            int cds_start = std::stoi(parsed_line[3]);
            int cds_end = std::stoi(parsed_line[4]);
            int exon_count = std::stoi(parsed_line[5]);
            std::string name = parsed_line[0];
            if (strand == '+') {
                result[name] = construct_transcript(exon_count, cds_start, cds_end, starts, stops);
            } else {
                result[name] = construct_transcript(exon_count, cds_end, cds_start, stops, starts);
            }
        }
    }
    support_table.close();
    return result;
}

int get_position_in_alignment(Transcript* transcript, int mutation_position) {
    int exon_index = -1;
    for(int i = 0; i < transcript->starts.size(); ++i) {
        if (transcript->starts[i] <= mutation_position && transcript->stops[i] >= mutation_position) {
            exon_index = i;
            break;
        }
    }
    if (exon_index < 1) return -1;
    int start_gap = transcript->stops[0] - transcript->cds_start;
    int stop_gap = transcript->cds_end - transcript->starts.back();
    int actual_position = start_gap;
    for(int i = 1; i < exon_index; ++i) {
        actual_position += transcript->stops[i] - transcript->starts[i];
    }
    actual_position += mutation_position - transcript->starts[exon_index];
    return abs(actual_position) - 1;
}

void save_alignment(std::string name, Mutation* mutation, std::vector<std::string> sequences) {
    std::string filename = "data/examples/" + name + "_" + std::to_string(mutation->position) + "_" + mutation->from + "_" + mutation->to + ".fa";
    std::ofstream example_alignment(filename);
    int counter = 0;

    std::vector<std::string> species{"hg38", "panTro4", "panPan1", "gorGor3", "ponAbe2", "nomLeu3", "rheMac3", "macFas5", "papAnu2", "chlSab2", "nasLar1", "rhiRox1", "calJac3", "saiBol1", "tarSyr2", "micMur1", "otoGar3"};

    for(auto &seq: sequences) {
        std::transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
        char uppercase_letter = std::toupper(seq[mutation->position]);
        seq[mutation->position] = uppercase_letter;
        example_alignment << ">" << species[counter] << std::endl;
        example_alignment << seq << std::endl;
        counter++;
    }
    example_alignment.close();
}

void parse_alignment(std::string filename, std::map<std::string, std::vector<Mutation* > > &path_mutations, std::map<std::string, std::vector<Mutation* > > &benign_mutations,  std::map<std::string, Transcript* > &support, bool aminoacid, bool export_alignments, Data* info, bool get_alignments=false) {
    std::ifstream alignment(filename);
    std::string current_name;
    std::vector<std::string> sequences;

    for(std::string line; std::getline(alignment, line);) {
        if (line.rfind(">", 0) == 0) {
            // get current name
//            >NM_000299_hg38 2244 chr1:201283703-201328836+
            std::vector<std::string> parsed_line = split_string(line, ' ');
            current_name = "NM_" + split_string(parsed_line[0], '_')[1];
        } else if (line.empty() && sequences.size() > 0) {
            std::vector<Mutation *> benign_possible_mutations;
            std::vector<Mutation *> path_possible_mutations;
            if (path_mutations.find(current_name) != path_mutations.end()) {
                path_possible_mutations = path_mutations[current_name];
            }
            if (benign_mutations.find(current_name) != benign_mutations.end()) {
                benign_possible_mutations = benign_mutations[current_name];
            }
            for(auto const &mutation: path_possible_mutations) {
                int actual_position = mutation->position;
                if (!aminoacid) actual_position = get_position_in_alignment(support[current_name], mutation->position);
                if (actual_position > 0) {
                    actual_position = actual_position - 1;
                    int species_count = 0;

                    if (sequences[0][actual_position] == mutation->from) {
                        for (int i = 1; i < sequences.size(); ++i) {
                            if (sequences[i][actual_position] == mutation->to) species_count++;
                        }
                    }
                    if (species_count == 1) {
                        info->patho_one++;
                    }
                    if (species_count > 1) {
                        info->patho_more++;
//                        if (export_alignments) save_alignment(current_name, mutation, sequences);
                    }
//                    if (export_alignments) save_alignment(current_name, mutation, sequences);
                }
            }
            for(auto const &mutation: benign_possible_mutations) {
                int actual_position = mutation->position;
                if (!aminoacid) actual_position = get_position_in_alignment(support[current_name], mutation->position);
//                int actual_position = get_position_in_alignment()
                if (actual_position > 0) {
                    actual_position = actual_position - 1;
                    int species_count = 0;
                    if (sequences[0][actual_position] == mutation->from) {
                        for (int i = 1; i < sequences.size(); ++i) {

                            if (sequences[i][actual_position] == mutation->to) species_count++;
                        }
                    }
                    if (species_count == 1) info->benign_one++;
                    if (species_count > 1) {
                        info->benign_more++;
                        if (export_alignments) save_alignment(current_name, mutation, sequences);
                    }

                }
            }

            sequences.clear();
            // do logic here
        } else if (line.empty()) {
            continue;
        } else {
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            sequences.push_back(line);
        }
    }
    alignment.close();
}

void output_table(Data* info) {
    std::ofstream output_table("result_table.txt");
    output_table << "==============================================================" << std::endl;
    output_table << "ClinVar, Feb.2019, GRCh38" << std::endl;
    output_table << std::left << std::setw(21) << "Benign" << '\t' << info->benign << '\t' << double(info->benign) / (info->benign + info->pathogenic) * 100 << "%" << std::endl;
    output_table << std::left << std::setw(21) << "Pathogenic" << '\t' << info->pathogenic << '\t' << double(info->pathogenic) / (info->benign + info->pathogenic) * 100 << "%" << std::endl;
    output_table << std::left << std::setw(21) << "Total" << '\t' << info->benign + info->pathogenic << std::endl;
    output_table << std::endl;
    output_table << "Pathogenic and benign CPDs, 16 primates (19 vertebrates)" << std::endl;
    output_table << std::left << std::setw(23) << "Significance" << '\t' << "Total" << '\t' << "CPDs" << '\t' << "in 1" << '\t' << "in >1" << std::endl;
    output_table << "--------------------------------------------------------------" << std::endl;
    output_table << std::left << std::setw(23) << "Benign" << '\t' << info->benign << '\t' << info->benign_more + info->benign_one << '\t' << info->benign_one << '\t' << info->benign_more << std::endl;
    output_table << std::left << std::setw(23) << "Pathogenic" << '\t' << info->pathogenic << '\t' << info->patho_one + info->patho_more << '\t' << info->patho_one << '\t' << info->patho_more << std::endl;
    output_table << "==============================================================" << std::endl;
    output_table.close();
}

int main() {

    Data* info = new Data();
    bool aminoacid = true;
    bool export_alignments = false;
    // parse mutations into map [NAME] = {mutations}
    std::map<std::string, std::vector<Mutation* > > path_mutations = parse_mutations_file("data/results/pathogenic.tab", aminoacid, info, 'p');
    std::map<std::string, std::vector<Mutation* > > benign_mutations = parse_mutations_file("data/results/benigns.tab", aminoacid, info, 'b');
    // parse support into map [NAME] = DESCRIPTION e.g. exons
    std::map<std::string, Transcript* > support_info = parse_support_table("data/alignments/support_table.tab");
    // parse alignment
    std::cout << info->pathogenic << " vs " << info->benign << std::endl;
    parse_alignment("data/alignments/new_alignment_16", path_mutations, benign_mutations, support_info, aminoacid, export_alignments, info, true);
//    output_result();
    std::cout << "Patho_one: " << info->patho_one << ", patho_more: " << info->patho_more << std::endl;
    std::cout << "Benign_one: " << info->benign_one << ", benign_more: " << info->benign_more << std::endl;

    output_table(info);
    path_mutations.clear();
    benign_mutations.clear();
    support_info.clear();
    return 0;
}