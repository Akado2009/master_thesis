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

struct Description {
    int start;
    int stop;
    std::string chromosome;
    std::string name;

    Description(int start, int stop, std::string chromosome, std::string name): start(start), stop(stop), chromosome(chromosome), name(name){}
};

struct Mutation {
    int position;
    char from;
    char to;

    Mutation(int position, char from, char to): position(position), from(from), to(to){}
};

std::ostream& operator<<(std::ostream& os, const Mutation& mutation) {
    return os << "Mutation: position " << mutation.position << " from " << mutation.from << " to " << mutation.to << std::endl;
}

std::vector<std::string> split_string(std::string &input_string, char delimiter) {
    std::vector<std::string> result;
    std::istringstream line_stream(input_string);
    std::string substring;
    while (std::getline(line_stream, substring, delimiter))
        result.push_back(substring);
    return result;
}

Description* get_description(std::vector<std::string> raw_desc) {
    int start, stop;
    std::string chromosome = "";
    std::vector<std::string> result = split_string(raw_desc[0], '_');
    if (result[2] == "hg38") {
        std::vector <std::string> position = split_string(raw_desc[2], '-');

        stop = std::stoi(position[1]);
        std::vector <std::string> chrome_start = split_string(position[0], ':');
        start = std::stoi(chrome_start[1]);
        chromosome = chrome_start[0];

    }
    Description* desc = new Description(start, stop, chromosome, "NM_" + result[1]);
    return desc;
}

Mutation* construct_Mutation(std::vector<std::string> &input) {
//    std::map<std::string, char> aa_dict{
//            {"Ala", 'A'},
//            {"Arg", 'R'},
//            {"Asn", 'N'},
//            {"Asp", 'D'},
//            {"Cys", 'C'},
//            {"Glu", 'E'},
//            {"Gln", 'Q'},
//            {"Gly", 'G'},
//            {"His", 'H'},
//            {"Ile", 'I'},
//            {"Leu", 'L'},
//            {"Lys", 'K'},
//            {"Met", 'M'},
//            {"Phe", 'F'},
//            {"Pro", 'P'},
//            {"Ser", 'S'},
//            {"Thr", 'T'},
//            {"Trp", 'W'},
//            {"Tyr", 'Y'},
//            {"Val", 'V'}
//    };
//    char from, to;
//    std::vector<std::string> preparsed_line = split_string(input[2], ' ');
//    std::vector<std::string> temp_line = split_string(preparsed_line[1], '.');
//    int i = 0;
//    int j = temp_line[1].size() - 2; // last one is (
//    // each AA is always 3
//    std::string first, second, position;
//    for(; i < 3; ++i) first += temp_line[1][i];
//    for(; j > temp_line[1].size() - 5; --j) second += temp_line[1][j];
//    std::reverse(second.begin(), second.end());
//    for(; i <= j; ++i) position += temp_line[1][i];
//    Mutation* mutation = new Mutation(std::stoi(position), aa_dict[first], aa_dict[second]);
//    for(auto const &x: input) std::cout << x << std::endl;
    Mutation* mutation = new Mutation(std::stoi(input[6]), input[8][0], input[9][0]);
//    std::cout << *mutation << std::endl;
    return mutation;
}

std::vector<Mutation *> find_mutations(int start, int stop, std::vector<Mutation *> &mutations) {
    std::vector<Mutation *> result;
    for(auto mutation: mutations) {
        if (mutation->position >= start && mutation->position <= stop) {
//            std::cout << "Start: " << start << ", stop: " << stop << ", position: " << mutation->position << std::endl;
            // get actual index inside the sequence
            mutation->position = mutation->position - start + 1;
            result.push_back(mutation);
        }
        if (mutation->position > stop) {
            break;
        }
    }
    return result;
}

std::map<std::string, std::vector<Mutation *> > parse_mutations_file(std::string filename, int &count) {
    std::ifstream input_file(filename);
    std::map<std::string, std::vector<Mutation *> > result;

    for(std::string line; std::getline(input_file, line);) {
        if (line.rfind("Chromosome", 0) != 0) {
            std::vector<std::string> parsed_line = split_string(line, '\t');
            std::string key = "chr" + parsed_line[0];
            if (result.find(key) == result.end()) {
                std::vector<Mutation *> to_insert;
                result[key] = to_insert;
            }
//            create a function for this
            if (parsed_line[2].find("=") == std::string::npos && parsed_line[2].find("p.") != std::string::npos && parsed_line[2].find("?") == std::string::npos && parsed_line[2].find(" ") != std::string::npos && parsed_line[2].find("del") == std::string::npos && parsed_line[2].find("*") == std::string::npos && parsed_line[2].find("))") == std::string::npos) {
                result[key].push_back(construct_Mutation(parsed_line));
                count++;
            }
        }
    }
    input_file.close();
    return result;
}

void save_alignment(std::vector<Description*> &descs, std::vector<std::string>& sequences, int position, char from, char to) {
    std::string name = descs[0]->name;
    std::ofstream example_alignment("data/examples/" + name + "_" + std::to_string(position) + "_" + from + "_" + to + ".fa");
    int counter = 0;
//    std::cout << "Mutation position is: " << position << ", but seq length is: " << sequences[0].size() << std::endl;
    for(auto &seq: sequences) {
        std::transform(seq.begin(), seq.end(), seq.begin(), ::tolower);
        char uppercase_letter = std::toupper(seq[position]);
        seq[position] = uppercase_letter;
        example_alignment << ">" << counter << std::endl;
        example_alignment << seq << std::endl;
        counter++;
    }
    example_alignment.close();
}

std::map<std::string, std::pair<int, int> > parse_support_info(std::string filename) {
    std::map<std::string, std::pair<int, int> > result;
    std::ifstream support_table(filename);
    for(std::string line; std::getline(support_table, line);) {
        std::vector<std::string> parsed_line = split_string(line, '\t');
        int start;
        int stop;
        if (parsed_line[2] == "+") {
            start = std::stoi(parsed_line[3]);
            stop = std::stoi(parsed_line[4]);
        } else {
            stop = std::stoi(parsed_line[3]);
            start = std::stoi(parsed_line[4]);
        }
        result[parsed_line[0]] = std::make_pair(start, stop);
    }
    return result;
}

int main() {
    auto start = std::chrono::high_resolution_clock::now();
    int patho_one = 0;
    int patho_more = 0;
    int benign_one = 0;
    int benign_more = 0;
    int benign = 0;
    int pathogenic = 0;

//    PARSING MUTATIONS
    std::map<std::string, std::vector<Mutation *> > benign_mutations = parse_mutations_file("data/results/pathogenic.tab", pathogenic);
    std::map<std::string, std::vector<Mutation *> > pathogenic_mutations = parse_mutations_file("data/results/benigns.tab", benign);

    //sort all mutations
    for(auto const& key: pathogenic_mutations) {
        sort(pathogenic_mutations[key.first].begin(), pathogenic_mutations[key.first].end(), [ ](const auto& lhs, const auto& rhs) {
                return lhs->position < rhs->position;
        });
    }

    for(auto const& key: benign_mutations) {
        sort(benign_mutations[key.first].begin(), benign_mutations[key.first].end(), [ ](const auto& lhs, const auto& rhs) {
            return lhs->position < rhs->position;
        });
    }

    std::map<std::string, std::pair<int, int> > support_info = parse_support_info("data/alignments/support_table");

    std::ifstream alignment("data/alignments/alignments_16_NUC");

//    PARSING ALIGNMENTS + DOING ALL THE STUFF

    std::vector<Description*> descriptions;
    std::vector<std::string> sequences;

    for(std::string line; std::getline(alignment, line);) {
        if (line.rfind(">", 0) == 0) {
            descriptions.push_back(get_description(split_string(line, ' ')));
        } else if (line.empty() && descriptions.size() > 0) {
            // DESC 0 is default human since alignment starts with human seq
           // do for benigns and for pathogenic, (just copy paste?)
            if (support_info.find(descriptions[0]->name) != support_info.end()) {
                descriptions[0]->start = support_info[descriptions[0]->name].first;
                descriptions[0]->stop = support_info[descriptions[0]->name].second;
                std::vector<Mutation *> possible_mutations_pathogenic = find_mutations(descriptions[0]->start, descriptions[0]->stop, pathogenic_mutations[descriptions[0]->chromosome]);

                for(auto mutation: possible_mutations_pathogenic) {
                    int species_count = 0;
                    std::cout << "Difference is: " << descriptions[0]->stop - descriptions[0]->start + 1 << ", actual length is: " << sequences[0].size() << ", mutation position: " << mutation->position << std::endl;
                    if (sequences[0][mutation->position] == mutation->from) {
                        // norm in human

                        for(int i = 1; i < sequences.size(); ++i) {

                            if (sequences[i][mutation->position] == mutation->to)
                                species_count++;
                        }
                    }
                    if (species_count == 1) patho_one++;
                    if (species_count > 1) {
                        patho_more++;
                        save_alignment(descriptions, sequences, mutation->position, mutation->from, mutation->to);
                    }
                }

                std::vector<Mutation *> possible_mutations_benign = find_mutations(descriptions[0]->start, descriptions[0]->stop, benign_mutations[descriptions[0]->chromosome]);
                for(auto mutation: possible_mutations_benign) {
                    int species_count = 0;
                    // sequences[0] is a human one
                    if (sequences[0][mutation->position] == mutation->from) {
                        // norm in human
                        for(int i = 1; i < sequences.size(); ++i) {
                            if (sequences[i][mutation->position] == mutation->to)
                                species_count++;
                        }
                    }
                    if (species_count == 1) benign_one++;
                    if (species_count > 1) benign_more++;
                }
            }

            descriptions.clear();
            sequences.clear();
        } else if (line.empty()) {
            continue;
        } else {
            std::transform(line.begin(), line.end(), line.begin(), ::toupper);
            sequences.push_back(line);
        }
    }
    benign_mutations.clear();
    pathogenic_mutations.clear();
    descriptions.clear();
    alignment.close();

    auto finish = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed = finish - start;
    std::cout << "Type\tOne\tMore\n";
    std::cout << "Pathologic" << '\t' << patho_one << '\t' << patho_more << std::endl;
    std::cout << "Bening" << '\t' << benign_one << '\t' << benign_more << std::endl;
    std::cout << "Execution time: " << elapsed.count() << std::endl;


    std::ofstream output_table("result_table.txt");
    output_table << "==============================================================" << std::endl;
    output_table << "ClinVar, Feb.2019, GRCh38" << std::endl;
    output_table << std::left << std::setw(21) << "Benign" << '\t' << benign << '\t' << double(benign) / (benign + pathogenic) * 100 << "%" << std::endl;
    output_table << std::left << std::setw(21) << "Pathogenic" << '\t' << pathogenic << '\t' << double(pathogenic) / (benign + pathogenic) * 100 << "%" << std::endl;
    output_table << std::left << std::setw(21) << "Total" << '\t' << benign + pathogenic << std::endl;
    output_table << std::endl;
    output_table << "Pathogenic and benign CPDs, 16 primates (19 vertebrates)" << std::endl;
    output_table << std::left << std::setw(23) << "Significance" << '\t' << "Total" << '\t' << "CPDs" << '\t' << "in 1" << '\t' << "in >1" << std::endl;
    output_table << "--------------------------------------------------------------" << std::endl;
    output_table << std::left << std::setw(23) << "Benign" << '\t' << benign << '\t' << benign_more + benign_one << '\t' << benign_one << '\t' << benign_more << std::endl;
    output_table << std::left << std::setw(23) << "Pathogenic" << '\t' << pathogenic << '\t' << patho_one + patho_more << '\t' << patho_one << '\t' << patho_more << std::endl;
    output_table << "==============================================================" << std::endl;
    output_table.close();
    return 0;
}