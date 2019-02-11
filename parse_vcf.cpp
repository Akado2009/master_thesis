#include <iostream>
#include <fstream>
#include <vector>
#include <string>
#include <sstream>
#include <algorithm>

//HEADER
//#CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO
// 1	182382427	293897	C	T	.	.	ALLELEID=278615;CLNDISDB=MedGen:C1864910,OMIM:610015,Orphanet:ORPHA71278;CLNDN=Glutamine_deficiency,_congenital;CLNHGVS=NC_000001.11:g.182382427C>T;CLNREVSTAT=criteria_provided,_single_submitter;CLNSIG=Uncertain_significance;CLNVC=single_nucleotide_variant;CLNVCSO=SO:0001483;CLNVI=Illumina_Clinical_Services_Laboratory,Illumina:479995;GENEINFO=GLUL:2752;MC=SO:0001624|3_prime_UTR_variant;ORIGIN=1;RS=868130362

int main() {
    int benings = 0;
    int pathogenic = 0;
    int overall = 0;

    std::ifstream read("data/vcf/clinvar_20190204.vcf");
    std::ofstream benings_write("data/results/benings.tab");
    std::ofstream pathogenic_write("data/results/pathogenic.tab");

    benings_write << "Chromosome\tPosition\tID\tFrom\tTo" << std::endl;
    pathogenic_write << "Chromosome\tPosition\tID\tFrom\tTo" << std::endl;

    for(std::string line; std::getline(read, line); ) {
        if (line.rfind('#', 0) != 0) {
            std::vector<std::string> result;
            std::istringstream line_stream(line);
            std::string substring;
            while (std::getline(line_stream, substring, '\t'))
                result.push_back(substring);
            overall++;
            if (result[3].size() == 1 && result[4].size() == 1) {
                if (result[7].find("Pathogenic") != std::string::npos) {
                    pathogenic++;
                    pathogenic_write << result[0] << "\t" << result[1] << "\t" << result[2] << "\t" << result[3] << "\t"
                                     << result[4] << std::endl;
                } else {
                    benings++;
                    benings_write << result[0] << "\t" << result[1] << "\t" << result[2] << "\t" << result[3] << "\t"
                                  << result[4] << std::endl;
                }
            }
        }
    }
    read.close();

    std::ofstream write("vcf_info.txt");
    write << "Overall number: " << overall << std::endl;
    write << "Benings number: " << benings << std::endl;
    write << "Pathogenic number: " << pathogenic << std::endl;
    write.close();
    benings_write.close();
    pathogenic_write.close();
    return 0;
}
