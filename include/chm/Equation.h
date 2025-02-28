#ifndef __BYTENOL_CHM_EQUATION_H__
#define __BYTENOL_CHM_EQUATION_H__

#include <set>

#include "Core.h"

namespace chm
{

    enum class EquationType
    {
        MOLECULAR = 0x00001,
        IONIC = 0x00010,
        THERMOCHEMICAL = 0x00100, 
    };


    class Equation
    {
        std::vector<int> coeff;
        std::string text, error;
        std::vector<chm::ChemicalEntity> m_reactants, m_products;
        std::set<std::string> lhs_atoms, rhs_atoms;
        std::vector<std::vector<int>> matrices;

        public:
            explicit Equation(const std::string& text);

            const decltype(m_reactants)& GetReactants() const;
            const decltype(m_products)& GetProducts() const;
            const decltype(lhs_atoms)& GetReactingAtoms() const;

        private:

            /**
             * This method setup basic functionalities including the information
             * about the reactants, products and the atoms in the equation as a whole
             */
            void SetReactingSpecies();

            void SetMatrix();

            void BalanceMolecularEqn();
    };

}

#endif 