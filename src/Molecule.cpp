#include <chm/Molecule.h>


using namespace chm;

Molecule::Molecule(const std::string& formula)
{
    this->formula = "(("+formula+"))";  // the double parenthesis here is quite important
    SetupAtom();
}

decltype(Molecule::atomDict) Molecule::GetAtoms()
{
    return atomDict;
}

int Molecule::GetCount(const std::string& symbol)
{
    if(atomDict.find(symbol) == atomDict.end()) return 0;
    return atomDict[symbol];
}


void Molecule::GetCharges(const std::string& formula)
{
    auto delimeterStart = formula.find("{");
    for(size_t i = delimeterStart; i < formula.size(); i++) 
    {
        if(formula[i] == '}') break;
        charge += formula[i];
    }
}

void Molecule::SetupAtom()
{
    const auto atomsFlattened = GetAtomString();
    auto h = GetAtomAsString(atomsFlattened);
    for(const auto& data: h)
    {
        auto& symb = data.first;
        if(atomDict.find(symb) == atomDict.end())
            atomDict[symb] = std::stoi(data.second);
        else
            atomDict[symb] += std::stoi(data.second);
    }
}

std::string Molecule::GetAtomString()
{
    std::string atomFlatText;
    std::stack<std::string> currentNodes;
    std::stack<std::string> lastNodeData;
    // std::stac

    for(size_t i = 0; i < formula.size(); i++)
    {
        auto ch = formula[i];
        if(ch == '(')
        {
            currentNodes.push("");
            continue;
        } else if(ch == ')') {
            auto currIndex = i;
            std::string s_count = "";
            while(++currIndex < formula.size())
            {
                auto nextDigit = formula[currIndex];
                if(!std::isdigit(nextDigit)) break;
                s_count += nextDigit;
            }
            int i_count = std::stoi(s_count.empty() ? "1": s_count);
            i = currIndex - 1;

            auto lastNode = currentNodes.top();
            auto text = CountAtomAsText(lastNode, i_count);
            currentNodes.pop();

            if(currentNodes.size())
            {
                currentNodes.top() += text;
                atomFlatText = currentNodes.top();
            }
            continue;
        }

        if(std::isalnum(ch)) currentNodes.top() += ch;
    }

    return atomFlatText;
}

std::vector<std::pair<std::string, std::string>> Molecule::GetAtomAsString(const std::string& text)
{
    std::vector<std::pair<std::string, std::string>> atoms;
    for(size_t i = 0; i < text.size(); i++)
    {
        auto ch = text[i];
        if(std::isupper(ch))
            atoms.push_back({"", ""});
        
        if(atoms.empty()) continue;
        if(std::isalpha(ch)) atoms.back().first += ch;
        if(std::isdigit(ch)) atoms.back().second += ch;
    }
    return atoms;
}

/**
 * This method accepts a string of all [A-Z|a-z]+ and [0-9] only
 * @param text is the text to count
 * @param scale is how much to scale the count
 * 
 * (S2O4, 2) => S4O8
 * (He3FS2) => He3F1S2
 */
std::string Molecule::CountAtomAsText(const std::string& text, int scale)
{
    std::string res;
    auto atoms = GetAtomAsString(text);
    for(auto& atom: atoms)
        res += atom.first + std::to_string(std::stoi(atom.second.empty() ? "1" : atom.second) * scale);
    
    return res;
}