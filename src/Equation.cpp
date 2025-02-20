#include <chm/Equation.h>

using namespace chm;

Equation::Equation(const std::string& text)
{
    this->text = text;
    SetReactingSpecies();
    SetMatrix();
    BalanceMolecularEqn();
}


const decltype(Equation::lhs_atoms)& Equation::GetReactingAtoms() const
{
    return lhs_atoms;
}


void Equation::SetReactingSpecies()
{
    size_t iter = text.size();

    // remove all whitespace from the text
    while((iter = text.find(' ')) < text.size())
        text.erase(text.begin() + iter);

    auto eqn_side = Split(text, '=');
    if(eqn_side.size() != 2) {
        error = "Invalid Equation: Both reactants and products must be seperated by comma";
        return;
    }

    for(iter = 0; iter < 2; iter++)
    {
        auto& molecules = iter > 0 ? m_products : m_reactants;
        auto& atoms = iter > 0 ? rhs_atoms : lhs_atoms;
        auto molecularText = Split(eqn_side[iter], '+');
        for(size_t j = 0; j < molecularText.size(); j++) 
            molecules.push_back({ molecularText[j] });

        for(auto& molecule: molecules)
        {
            for(auto& h: molecule.GetAtoms())
                atoms.emplace(h.first);
        }
    }

    if(!(lhs_atoms == rhs_atoms)) {
        error = "Equation does not Obey the law of conservation of matter";
        return;
    }
}


const decltype(Equation::m_reactants)& Equation::GetReactants() const
{
    return m_reactants;
}


const decltype(Equation::m_products)& Equation::GetProducts() const
{
    return m_products;
}


void Equation::SetMatrix()
{
    auto rSum = m_reactants;
    rSum.insert(rSum.end(), m_products.begin(), m_products.end());

    matrices.clear();
    for(auto& atom: lhs_atoms)
    {
        matrices.push_back({});
        for(size_t i = 0; i < rSum.size(); i++)
            matrices.back().push_back( rSum[i].GetCount(atom) );
    }
}



void Equation::BalanceMolecularEqn()
{
    tensor_t<Fraction> reducedMatrix;
    vector_t<Fraction> coeff;
    ReduceRow(matrices, reducedMatrix, coeff);
}