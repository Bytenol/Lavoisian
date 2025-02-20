#include <chm/Math.h>
// #include "Math.h"


namespace chm
{

    Fraction::Fraction(int n, int d)
    {
        num = n;
        denom = d;
        Simplify();
    }


    Fraction& Fraction::operator*=(const int& scale) 
    {
        num *= scale;
        Simplify();
        return *this;
    }


    void Fraction::Simplify()
    {
        int gcd = GetGCD(num, denom);
        num /= gcd;
        denom /= gcd;
    }

    
    std::string Fraction::ToString()
    {
        if(std::abs(denom) == 1) return std::to_string(num);
        return std::to_string(num) + "/" + std::to_string(denom);
    }

    int GetGCD(int a, int b)
    {
        while(b != 0) 
        {
            int temp = b;
            b = a % b;
            a = temp;
        }
        return a;
    }


    void ReduceRow(const tensor_t<int>& m_in, tensor_t<Fraction>& m_out, vector_t<Fraction>& solution)
    {
        auto matrices = m_in;

        const int COL_SIZE = matrices[0].size() - 1;
        const int ROW_SIZE = matrices.size();

        // check if a given 2d matrix is diagonal
        auto is_diagonal = [=](decltype(matrices)& m) -> bool {
            bool isDiagonal = true;
            for(int i = 0; i < COL_SIZE; i++)
            {
                if(m[i][i] == 0) 
                {
                    isDiagonal = false;
                    break;
                }
            }
            return isDiagonal;
        };

        // swap rows of gausian procedures
        auto swap_row = [&]() {
            decltype(matrices) m = matrices;
            bool isDiagonal = is_diagonal(m);
            while(std::next_permutation(m.begin(), m.end()) && !isDiagonal)
            {
                isDiagonal = is_diagonal(m);
                if(isDiagonal) matrices = m;
            }

            return isDiagonal;
        };

        // position to reduce
        std::vector<std::pair<int, int>> lowerTriangle, upperTriangle;

        for(size_t i = 0; i < matrices[0].size() - 1; i++)
        {
            for(size_t j = 0; j < 2; j++)
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
            for(size_t i = 0; i < 2; i++)
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
            if(!swap_row()) {
                std::cerr << "No solution" << std::endl;
                return;
            }
            getAllPosToReduce();
            if(posToReduce.empty()) isReduced = true;
            for(const auto& pos: posToReduce)
            {
                auto& pivot = matrices[pos.second];
                const int a = pivot[pos.second];
                const int b = matrices[pos.first][pos.second];
                for(size_t i = 0; i < pivot.size(); i++)
                {
                    const int x = matrices[pos.first][i];
                    matrices[pos.first][i] = a * x - b * pivot[i];
                }
            }
        }


        int maxDenom = -1;
        for(int i = 0; i < COL_SIZE; i++)
            maxDenom = std::max(maxDenom, matrices[i][i]);

        for(int i = 0; i < ROW_SIZE; i++)
        {
            m_out.push_back({});
            for(int j = 0; j < (COL_SIZE + 1); j++) {
                m_out[i].push_back(Fraction(matrices[i][j], maxDenom ));
                if(j == COL_SIZE) 
                {
                    m_out[i][j] *= maxDenom;
                    solution.push_back(m_out[i][j]);
                }
            }
        } // for ROW_SIZE_ENDS

        solution.push_back(Fraction(maxDenom));

    }

}