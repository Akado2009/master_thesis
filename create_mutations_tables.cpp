#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>
#include <chrono>
#include <map>
#include <iterator>

//riginSimple	ClinicalSignificance	ReviewStatus	NumberSubmitters	PhenotypeIDS	PhenotypeList	OtherIDs	IsCancer

//std::vector<int> data_to_export{18, 1, 2, 3, 4, 10, 19, 20, 21, 22, 15, 6, 24, 25, 12, 13, 28};
//0 #AlleleID
//1 Type
//2 Name
//3 GeneID
//4 GeneSymbol
//5 HGNC_ID
//6 ClinicalSignificance
//7 ClinSigSimple
//8 LastEvaluated
//9 RS# (dbSNP)
//10 nsv/esv (dbVar)
//11 RCVaccession
//12 PhenotypeIDS
//13 PhenotypeList
//14 Origin
//15 OriginSimple
//16 Assembly
//17 ChromosomeAccession
//18 Chromosome
//19 Start
//20 Stop
//21 ReferenceAllele
//22 AlternateAllele
//23 Cytogenetic
//24 ReviewStatus
//25 NumberSubmitters
//26 Guidelines
//27 TestedInGTR
//28 OtherIDs
//29 SubmitterCategories
//30 VariationID

std::vector<std::string> split_string(std::string input_string, char delimiter) {
    std::vector<std::string> result;
    std::istringstream line_stream(input_string);
    std::string substring;
    while (std::getline(line_stream, substring, delimiter))
        result.push_back(substring);
    return result;
}

struct Mutation {
    std::string position;
    char from;
    char to;

    Mutation(std::string position, char from, char to): position(position), from(from), to(to){}
};

Mutation* construct_mutation(std::string input_mutation) {
    std::map<std::string, char> aa_dict{
            {"Ala", 'A'},
            {"Arg", 'R'},
            {"Asn", 'N'},
            {"Asp", 'D'},
            {"Cys", 'C'},
            {"Glu", 'E'},
            {"Gln", 'Q'},
            {"Gly", 'G'},
            {"His", 'H'},
            {"Ile", 'I'},
            {"Leu", 'L'},
            {"Lys", 'K'},
            {"Met", 'M'},
            {"Phe", 'F'},
            {"Pro", 'P'},
            {"Ser", 'S'},
            {"Thr", 'T'},
            {"Trp", 'W'},
            {"Tyr", 'Y'},
            {"Val", 'V'},
            {"Ter", '-'}
    };
    char from, to;
    std::string first, second, position;
    int i = 0;
    std::vector<std::string> preparsed_line = split_string(input_mutation, ' ');
    std::vector<std::string> temp_line = split_string(preparsed_line[1], '.');
    int j = temp_line[1].size() - 2;
    for(; i < 3; ++i) first += temp_line[1][i];
    for(; j > temp_line[1].size() - 5; --j) second += temp_line[1][j];
    std::reverse(second.begin(), second.end());
    for(; i <= j; ++i) position += temp_line[1][i];
    Mutation* mutation = new Mutation(position, aa_dict[first], aa_dict[second]);
    return mutation;
}

std::ostream& operator<<(std::ostream& os, const Mutation& mutation) {
    return os << "Mutation: position " << mutation.position << " from " << mutation.from << " to " << mutation.to << std::endl;
}

void write_string_to_file(std::ofstream &output_stream, std::vector<int> &indice, std::vector<std::string> &data, char delimiter) {
    for(std::size_t i = 0; i < indice.size(); ++i) {
        if (i != 0) output_stream << delimiter;
        output_stream << data[indice[i]];
    }
    output_stream << std::endl;
}

std::string parse_id(std::string &input) {
    return split_string(input, '.')[0];
}

bool is_valid_mutation(std::string mutation) {
    bool valid = true;
    std::vector<std::string> no_to_find{"=", "?", "del", "*", "Ter", "))"};
    std::vector<std::string> to_find{"p.", " "};
    for(auto const &keyword: to_find) valid &= (mutation.find(keyword) != std::string::npos);
    for(auto const &keyword: no_to_find) valid &= (mutation.find(keyword) == std::string::npos);
    return valid;
}

//dobavit isklucheniya
bool is_valid_benign(std::string mutation) {
    bool valid = false;
    std::vector<std::string> keywords{"Benign", "benign"};
    for(auto const &keyword: keywords) valid |= (mutation.find(keyword) != std::string::npos);
//    valid &= (mutation.find("Conflicting") == std::string::npos);
    return valid;
}

bool is_valid_pathogenic(std::string mutation) {
    bool valid = false;
    std::vector<std::string> keywords{"Pathogenic", "pathogenic"};
    for(auto const &keyword: keywords) valid |= (mutation.find(keyword) != std::string::npos);
//    valid &= (mutation.find("Conflicting") == std::string::npos);
    return valid;
}

int main() {
    int benigns = 0;
    int pathogenic = 0;
    int overall = 0;
    std::map<std::string, bool> duplicate_ids;

    std::vector<int> data_to_export{18, 1, 2, 3, 4, 10, 19, 20, 21, 22, 15, 6, 24, 25, 12, 13, 28};
    std::ifstream read("data/vcf/variant_summary.txt");
    std::ofstream benigns_write("data/results/benigns.tab");
    std::ofstream pathogenic_write("data/results/pathogenic.tab");
    std::ofstream ids("data/results/ids.tab");

    for(std::string line; std::getline(read, line);) {
        std::vector<std::string> parsed_line = split_string(line, '\t');
        if (line.rfind("#AlleleID", 0) == 0) {
            parsed_line[19] = "GenomicCoordinates";
            parsed_line[20] = "CDSCoordinates";
            write_string_to_file(pathogenic_write, data_to_export, parsed_line, '\t');
            write_string_to_file(benigns_write, data_to_export, parsed_line, '\t');
        }
        if (parsed_line[1] == "single nucleotide variant" && parsed_line[16] == "GRCh38"){
            if (is_valid_mutation(parsed_line[2])) {

                Mutation *mut = construct_mutation(parsed_line[2]);
                // modify teh struct of the table
                parsed_line[2] = split_string(parsed_line[2], '.')[0];
                parsed_line[20] = mut->position;
                parsed_line[21] = mut->from;
                parsed_line[22] = mut->to;
                if (is_valid_pathogenic(parsed_line[6]) && !is_valid_benign(parsed_line[6])) {
                    write_string_to_file(pathogenic_write, data_to_export, parsed_line, '\t');
                    pathogenic++;
                } else if (is_valid_benign(parsed_line[6]) && !is_valid_pathogenic(parsed_line[6])) {
                    write_string_to_file(benigns_write, data_to_export, parsed_line, '\t');
                    benigns++;
                }
                delete mut;
                std::string id = parse_id(parsed_line[2]);
                if (duplicate_ids.find(id) == duplicate_ids.end()) {
                    ids << id << std::endl;
                    duplicate_ids[id] = true;
                }
                overall++;
            }
        }

    }

    read.close();

    std::ofstream write("vcf_info.txt");
    write << "Overall number: " << overall << std::endl;
    write << "Benings number: " << benigns << std::endl;
    write << "Pathogenic number: " << pathogenic << std::endl;
    write.close();
    benigns_write.close();
    pathogenic_write.close();
    ids.close();
    return 0;
}
