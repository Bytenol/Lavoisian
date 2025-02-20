#include <iostream>
#include <vector>
#include <string>
#include <random>
#include <emscripten/emscripten.h>

// #include <chm/Molecule.h>
#include <chm/Equation.h>

using namespace std;
using namespace chm;

extern "C"
{

    EMSCRIPTEN_KEEPALIVE
    const char* GetRandomEqn()
    {
        vector<string> equations {
            "H2 + O2 = H2O",
            "H2SO4 + NaOH = Na2SO4 + H2O",
            "Mg + N2 = Mg2N3",
            "H2O = H2 + O2",
            "HCl + Zn = ZnCl2 + H2",
            "Na2SO4 + H2SO4 = NaHSO4",
            "H2O2 = H2O + O2",
            "CH4 + O2 = CO2 + H2O",
            "C3H8 + O2 = CO2 + H2O",
            "H3PO3 = H3PO4 + PH3",
            "Fe2(SO4)3 + NH3 + H2O = Fe(OH)3 + (NH4)2SO4",
            "Ca + H2O = Ca(OH)2 + H2",
            "Na + Cl2 = NaCl",
            "H2 + O2 = H2O",
            "NaNO3 = NaNO2 + O2",
            "Cu + HNO3 = Cu(NO3)2 + H2O + NO",
            "Fe2SiO4 + Mg2SiO4 + H2O + CO2 = Mg6(Si4O10)(OH)8 + Fe2O3 + CH4"
        };

        int index = 1;
        static const char* eqn = equations[index].c_str();
        return eqn;
    }
}

int main()
{
    

    Equation eqn{ "H2 + O2 = H2O" };
    auto reactants = eqn.GetReactants();
    auto products = eqn.GetProducts();
    auto atoms = eqn.GetReactingAtoms();

    for(auto& atom: atoms) cout << atom << " = ";
    

    cout << endl;
    return 0;
}