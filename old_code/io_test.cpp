#include <iostream>
#include <iomanip>


int main() {
    int ben = 89083;
    int path = 33356;
    int ben_one = 4738;
    int path_one = 650;
    int ben_more = 5019;
    int path_more = 560;

    std::cout << "ASDasdasdasdsaasd" << std::endl;
    std::cout << std::left << std::setw(21) << "Significance" << '\t' << "Total" << '\t' << "CPDs" << '\t' << "in 1" << '\t' << "in >1" << std::endl;
    std::cout << std::left << std::setw(21) << "Benign" << '\t' << ben << '\t' << ben_one + ben_more << '\t' << ben_one << '\t' << ben_more << std::endl;
    return 0;
}