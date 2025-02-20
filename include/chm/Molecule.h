#ifndef __BYTENOL_CHM_MOLECULE_H__
#define __BYTENOL_CHM_MOLECULE_H__

#include <string>
#include <map>
#include <stack>
#include <vector>

namespace chm
{

    // struct Atom
    // {
    //     std::string symbol;
    //     std::string name;
    //     std::string desc;
    // };


    class Molecule
    {
        std::string formula;
        std::map<std::string, int> atomDict;
        std::string charge; // not fully implemented

        public:
            Molecule(const std::string& formula);

            decltype(atomDict) GetAtoms();

            int GetCount(const std::string& symbol);

        private:

            void GetCharges(const std::string& formula);

            void SetupAtom();

            std::string GetAtomString();

            std::vector<std::pair<std::string, std::string>> GetAtomAsString(const std::string& text);

            /**
             * This method accepts a string of all [A-Z|a-z]+ and [0-9] only
             * @param text is the text to count
             * @param scale is how much to scale the count
             * 
             * (S2O4, 2) => S4O8
             * (He3FS2) => He3F1S2
             */
            std::string CountAtomAsText(const std::string& text, int scale = 1);

    };

}

#endif 