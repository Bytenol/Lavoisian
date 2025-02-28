/**
 * Add getIonic magnitude to atom
 */
#ifndef __BYTENOL_CHM_CORE_H__
#define __BYTENOL_CHM_CORE_H__

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <cstring>
#include <stack>
#include <algorithm>

#include "sqlite3.h"
#include "Error.h"

namespace chm
{

    template<typename T>
    using tensor_t = std::vector<std::vector<T>>;

    template<typename T>
    using vector_t = std::vector<T>;

    extern std::string DB_PATH;

    bool Init(const std::string& dbPath);


    /**
     * An atom is the building block of other neccessary materials and entities
     */
    class Atom
    {
        friend bool Init(const std::string& dbPath);
        static std::vector<Atom> atoms;

        public:

            int id;
            
            std::string name,
                symbol,
                appearance,
                discovered_by,
                named_by,
                phase,
                source,
                summary,
                block,
                cpk_hex,
                image_url,
                image_attribution,
                category,
                image_title,
                electron_configuration,
                electron_configuration_semantic;

            double atomic_mass,
                boil,
                density,
                melt,
                molar_heat,
                electron_affinity,
                electronegativity_pauling,
                number,
                period,
                xpos,
                ypos,
                wxpos,
                wypos,
                group;  //@todo not implemented in the database

            /***
             * @todo Implement blob for image and .glb files for the fields
             * bohr_model
             * bohr_model_2d
             */

            int count = 1;
            
            static decltype(atoms)& Get();
            static Atom* Get(const int& id);
            static Atom* Get(const std::string& symbol);
            // static Atom* Get(const std::string& symbol);

        private:
            // * Run inside the Init function to setup basic atom's functionalities
            static bool InitStatic();
    };


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //////////////////  CHEMICAL ENTITY       //////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////

    class ChemicalEntity
    {
        std::string formula, parsedFormula;
        std::map<std::string, int> atomDict;
        
        int atomicity = 0;
        int no_of_moles = 1;

        public:
            int charge;
            /**
             * H2, H2SO4, CO3{-<charge_number>}, Na, Na{+<charge_no>}
             */
            explicit ChemicalEntity(const std::string& formula);
            bool IsIon() const;
            decltype(atomDict) GetAtom();
            Atom GetAtom(const std::string& symbol);
            double GetMolarMass() const;

            friend std::ostream& operator<<(std::ostream& os, const ChemicalEntity& molecule)
            {
                os << molecule.parsedFormula;
                return os;
            }

        private:
            int ExtractCharge();
            std::vector<std::pair<std::string, std::string>> GetAtomAsString(const std::string& text);
            void SetupAtom();
            std::string GetParsedFormula();
    };


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //////////////////             MATH       //////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////


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


    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    //////////////////             Utils       //////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    ////////////////////////////////////////////////////////////////////////
    vector_t<std::string> Split(const std::string& str, char delimeter);
}

#endif 