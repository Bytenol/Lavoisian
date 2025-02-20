/**
 * Make `ReduceRow` accepts fraction as default
 */
#ifndef __BYTENOL_CHM_MATH_H__
#define __BYTENOL_CHM_MATH_H__

#include <iostream>
#include <vector>
#include <string>
#include <algorithm>

namespace chm
{

    template<typename T>
    using tensor_t = std::vector<std::vector<T>>;

    template<typename T>
    using vector_t = std::vector<T>;

    class Fraction
    {
        int num;
        int denom;

        public:
            explicit Fraction(int n = 0, int d = 1);

            bool isInt(){ return (GetDenom() == 1); }

            inline int GetDenom() { return std::abs(denom); }

            inline int GetNum() { return num * (denom < 0 ? -1: 1); }

            inline int GetValue() { return (GetDenom() == 1 ? GetNum(): num/denom); }

            Fraction& operator*=(const int& scale);

            friend std::ostream& operator<<(std::ostream& os, const Fraction& f)
            {
                if(std::abs(f.denom) == 1) os << f.num;
                else os << f.num << "/" << f.denom;
                return os;
            }

        private:
            void Simplify();
            std::string ToString();
    };

    int GetGCD(int a, int b);

    /**
     * If no solution returns an empty matrix
     */
    void ReduceRow(const tensor_t<int>& m_in, tensor_t<Fraction>& m_out, vector_t<Fraction>& solution);
}

#endif 