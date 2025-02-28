#include <iostream>
#include <chm/sqlite3.h>
#include <chm/Core.h>
using namespace chm;


int main() {

    if(!chm::Init("C:/Users/olami/Documents/Programming/VCS/Chemilator/assets/chemilator.db"))
        std::cerr << "Unable to load chm";

    // for(auto& atom: Atom::Get())
    // {
    //     std::cout << atom.first << " - " << atom.second.name << std::endl;
    // }

    std::cout << std::endl;
    return 0;
}