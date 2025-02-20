#include <chm/Equation.h>

using namespace chm;

Equation::Equation(const std::string& text)
{
    this->text = text;
    SetReactingSpecies();
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



void Equation::BalanceMolecularEqn()
{
    auto rSum = m_reactants;
    rSum.insert(rSum.end(), m_products.begin(), m_products.end());
    
    std::vector<std::vector<int>> matrices;

    for(auto& atom: lhs_atoms)
    {
        matrices.push_back({});
        for(int i = 0; i < rSum.size(); i++)
            matrices.back().push_back( rSum[i].GetCount(atom) );
    }

    // swap rows
    // @todo check if some matrices are impossible to be diagonal
    auto swap_row = [&]() {
        decltype(matrices) m = matrices;
        bool isDiagonal = false;
        while(std::next_permutation(m.begin(), m.end()) && !isDiagonal)
        {
            isDiagonal = true;
            for(int i = 0; i < matrices[0].size() - 1; i++)
            {
                if(m[i][i] == 0) 
                {
                    isDiagonal = false;
                    break;
                }
            }

            if(isDiagonal) matrices = m;
        }
    };

    std::vector<std::pair<int, int>> lowerTriangle, upperTriangle;

    for(int i = 0; i < matrices[0].size() - 1; i++)
    {
        for(int j = 0; j < 2; j++)
        {
            int start = (j == 0) ? i + 1: 0;
            int end = (j == 0) ? matrices.size(): i;
            auto& triangle_ref = (j == 0) ? lowerTriangle: upperTriangle;
            for(int k = start; k < end; k++)
                triangle_ref.push_back({ k, i });
        }
    }

    // check if any triangles are zero
    decltype(lowerTriangle) posToReduce;

    auto getAllPosToReduce = [&]() {
        posToReduce.clear();
        bool isDone = true;
        for(int i = 0; i < 2; i++)
        {
            auto& triangle_ref = (i == 0) ? lowerTriangle : upperTriangle;
            for(const auto& pos: triangle_ref)
                if(matrices[pos.first][pos.second] != 0) {
                    isDone = false;
                    posToReduce.push_back(pos);
                }
            if(!isDone) break;
        }
    };

    // reduce row echelon
    bool isReduced = false;
    while(!isReduced)
    {
        swap_row(); //@todo: go to swap_row and prevent infinite loop if any diagonal is zero
        getAllPosToReduce();
        if(posToReduce.empty()) isReduced = true;
        for(const auto& pos: posToReduce)
        {
            auto& pivot = matrices[pos.second];
            const int a = pivot[pos.second];
            const int b = matrices[pos.first][pos.second];
            for(int i = 0; i < pivot.size(); i++)
            {
                const int x = matrices[pos.first][i];
                matrices[pos.first][i] = a * x - b * pivot[i];
            }
        }
        // break;
    }

    // std::vector<chm::Fraction> solutions;
    // int maxDenom = -1;
    // for(int i = 0; i < matrices[0].size() - 1; i++) {
    //     solutions.push_back({ matrices[i][matrices[0].size() - 1], matrices[i][i] });
    //     maxDenom = std::max(maxDenom, solutions.back().GetDenom());
    // }

    // std::vector<int> coeff;
    // for(auto& solution: solutions) 
    // {
    //     solution *= maxDenom;
    //     coeff.push_back(std::abs(solution.GetNum()));
    // }
    // coeff.push_back(maxDenom);
        
    // for(const auto& c: coeff) std::cout << c << " = ";

    // print out the matrices
    std::cout << std::endl;
    // for(auto & i: matrices)
    // {
    //     for(auto& j: i)
    //         std::cout << j << " ";
    //     std::cout << std::endl;
    // }

}