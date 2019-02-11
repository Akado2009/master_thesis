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

std::vector<std::string> split_string(std::string input_string, char delimiter) {
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
    if (result[1] == "hg38") {
        std::vector <std::string> position = split_string(raw_desc[4], '-');

        stop = std::stoi(position[1]);
        std::vector <std::string> chrome_start = split_string(position[0], ':');
        start = std::stoi(chrome_start[1]);
        chromosome = chrome_start[0];
    }
    Description* desc = new Description(start, stop, chromosome, result[1]);
    return desc;
}

Mutation* construct_Mutation(std::vector<std::string> &input) {
    Mutation* mutation = new Mutation(std::stoi(input[1]), input[3][0], input[4][0]);
    return mutation;
}

std::vector<Mutation *> find_mutations(int start, int stop, std::vector<Mutation *> &mutations) {
    std::vector<Mutation *> result;
//    std::cout << "Mutation: " << start << " " << stop << std::endl;
    for(auto mutation: mutations) {
        if (mutation->position >= start && mutation->position <= stop) {
            // get actual index inside the sequence
            mutation->position = mutation->position - start;
            result.push_back(mutation);
        }
        if (mutation->position > stop) {
            break;
        }
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
    std::map<std::string, std::vector<Mutation *> > benign_mutations;
    std::map<std::string, std::vector<Mutation *> > pathogenic_mutations;

    std::ifstream mutations("data/vcf/clinvar_20190204.vcf");
    for(std::string line; std::getline(mutations, line); ) {
        if (line.rfind('#', 0) != 0) {
            std::vector<std::string> result = split_string(line, '\t');
            if (result[3].size() == 1 && result[4].size() == 1) {
                std::string key = "chr" + result[0];
                if (result[7].find("Pathogenic") != std::string::npos) {
                    if (pathogenic_mutations.find(key) == pathogenic_mutations.end()) {
                        std::vector<Mutation *> to_insert_pathogenic;
                        pathogenic_mutations[key] = to_insert_pathogenic;
                    }
                    pathogenic_mutations[key].push_back(construct_Mutation(result));
                    pathogenic++;
                } else {
                    if (benign_mutations.find(key) == benign_mutations.end()) {
                        std::vector<Mutation *> to_insert_benign;
                        benign_mutations[key] = to_insert_benign;
                    }
                    benign_mutations[key].push_back(construct_Mutation(result));
                    benign++;
                }
            }
        }
    }
    mutations.close();

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

    std::ifstream alignment("data/alignments/knownCanonical.exonNuc.fa");

//    PARSING ALIGNMENTS + DOING ALL THE STUFF

    std::vector<Description*> descriptions;
    std::vector<std::string> sequences;

    for(std::string line; std::getline(alignment, line);) {
        if (line.rfind(">", 0) == 0) {
            descriptions.push_back(get_description(split_string(line, ' ')));
        } else if (line.empty()) {
            // DESC 0 is default human since alignment starts with human seq
           // do for benigns and for pathogenic, (just copy paste?)
            std::vector<Mutation *> possible_mutations_pathogenic = find_mutations(descriptions[0]->start, descriptions[0]->stop, pathogenic_mutations[descriptions[0]->chromosome]);

            for(auto mutation: possible_mutations_pathogenic) {
                int species_count = 0;
//                std::cout << sequences[0][mutation->position] << " " << mutation->from << std::endl;
                // sequences[0] is a human one
                if (sequences[0][mutation->position] == mutation->from) {
                    // norm in human
                    for(int i = 1; i < sequences.size(); ++i) {
                        if (sequences[i][mutation->position] == mutation->to)
                            species_count++;
                    }
                }

//                std::cout << typeid(sequences[0][0]).name() << std::endl;
                if (species_count == 1) patho_one++;
                if (species_count > 1) patho_more++;
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
//            std::cout << "Patho: " << patho_one << " " << patho_more << std::endl;
//            std::cout << "Bening: " << benign_one << " " << benign_more << std::endl;
//            std::cout << "--------------------------" << std::endl;
//             perform all the stuff
            descriptions.clear();
            sequences.clear();
        } else {
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
    output_table << std::left << std::setw(23) << "Benign" << '\t' << benign << '\t' << benign_more + benign_one << '\t' << benign_one << '\t' << benign_one << std::endl;
    output_table << std::left << std::setw(23) << "Pathogenic" << '\t' << pathogenic << '\t' << patho_one + patho_more << '\t' << patho_one << '\t' << patho_more << std::endl;
    output_table << "==============================================================" << std::endl;
    output_table.close();
    return 0;
}