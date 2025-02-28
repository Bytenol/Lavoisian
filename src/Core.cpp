#include <chm/Core.h>


namespace chm
{
    std::string DB_PATH = "";


    bool Init(const std::string& dbPath)
    {
        DB_PATH = dbPath;
        Atom::InitStatic();
        return true;
    }


    std::vector<Atom> Atom::atoms = { };

    decltype(Atom::atoms)& Atom::Get()
    {
        return atoms;
    }

    Atom* Atom::Get(const int& id)
    {
        if(id < 0 || id >= atoms.size())
            return nullptr;
        return &(atoms[id]);
    }


    Atom* Atom::Get(const std::string& symbol)
    {
        auto f = std::find_if(atoms.begin(), atoms.end(), [&symbol](const Atom& a) {
            return a.symbol == symbol;
        });

        if(f == atoms.end()) return nullptr;
        return &(*f);
    }


    bool Atom::InitStatic()
    {
        atoms.clear();

        const char* dbPath = DB_PATH.c_str();
        sqlite3* db = nullptr;
        char* errMsg = nullptr;
        if(sqlite3_open(dbPath, &db) != SQLITE_OK)
            throw DBFailedException("Failed to open Sqlite3 database. Please Make sure the DB_PATH is set correctly");

        // fetch all elements
        const char* sql = "SELECT * FROM elements;";
        sqlite3_stmt* stmt;

        if(sqlite3_prepare_v2(db, sql, -1, &stmt, nullptr) != SQLITE_OK) {
            throw DBFailedException("Error Preparing Database Statement");
        }

        auto getText = [&stmt](int index) {
            const auto text = reinterpret_cast<const char*>(sqlite3_column_text(stmt, index));
            return (text ? text : "NULL");
        };

        auto getInt = [&stmt](int index) { return sqlite3_column_int(stmt, index); };

        auto getDouble = [&stmt](int index) { return sqlite3_column_double(stmt, index); };

        while(sqlite3_step(stmt) == SQLITE_ROW)
        {
            Atom atom;
            atom.id = getInt(0);
            atom.name = getText(1);
            atom.symbol = getText(2);
            atom.appearance = getText(3);
            atom.discovered_by = getText(4);
            atom.named_by = getText(5);
            atom.phase = getText(6);
            atom.source = getText(7);
            atom.summary = getText(8);
            atom.block = getText(9);
            atom.cpk_hex = getText(10);
            atom.image_url = getText(11);
            atom.image_attribution = getText(12);
            atom.category = getText(13);
            atom.image_title = getText(14);
            atom.electron_configuration = getText(15);
            atom.electron_configuration_semantic = getText(16);
            atom.atomic_mass = getDouble(17);
            atom.boil = getDouble(18);
            atom.density = getDouble(19);
            atom.melt = getDouble(20);
            atom.molar_heat = getDouble(21);
            atom.electron_affinity = getDouble(22);
            atom.electronegativity_pauling = getDouble(23);
            atom.number = getInt(24);
            atom.period = getInt(25);
            atom.xpos = getInt(26);
            atom.ypos = getInt(27);
            atom.wxpos = getInt(28); 
            atom.wypos = getInt(29); 

            atoms.push_back(atom);
        }

        sqlite3_close(db);
        return true;
    }


    ChemicalEntity::ChemicalEntity(const std::string& formula)
    {
        this->formula = formula;
        charge = ExtractCharge();
        parsedFormula = GetParsedFormula();
        SetupAtom();
    }


    bool ChemicalEntity::IsIon() const
    {
        return (charge != 0);
    }


    decltype(ChemicalEntity::atomDict) ChemicalEntity::GetAtom()
    {
        return atomDict;
    }


    Atom ChemicalEntity::GetAtom(const std::string& symbol)
    {
        auto ptr = Atom::Get(symbol);
        if(!ptr) return {};
        if(atomDict.find(symbol) == atomDict.end()) return {};
        Atom atom = *ptr;
        atom.count = atomDict[symbol];
        return atom;
    }


    double ChemicalEntity::GetMolarMass() const
    {
        double res = 0;
        for(auto [symb, count]: atomDict) {
            auto h = Atom::Get(symb);
            if(h) {
                res += h->atomic_mass * count;
            }
        }
        return res;
    }


    int ChemicalEntity::ExtractCharge()
    {
        auto index = formula.find('{');
        if(!(index < formula.size()))
            return 0;

        std::string s_charge = "";
        for(size_t i = index + 1; i < formula.size(); i++)
        {
            auto ch = formula[i];
            if(ch != '}')
                s_charge += ch;
            else 
                break;
        }

        int res = 1;
        try {
            res = std::stoi(s_charge.substr(1)) * (s_charge[0] == '-' ? -1: 1);
        } catch(...) {
            res = 0;
        }

        this->formula = this->formula.substr(0, index);

        return res;
    }


    std::vector<std::pair<std::string, std::string>> ChemicalEntity::GetAtomAsString(const std::string& text)
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

    void ChemicalEntity::SetupAtom()
    {
        auto h = GetAtomAsString(parsedFormula);
        for(const auto& data: h)
        {
            auto& symb = data.first;
            if(atomDict.find(symb) == atomDict.end())
                atomDict[symb] = std::stoi(data.second);
            else
                atomDict[symb] += std::stoi(data.second);
        }
    }

    std::string ChemicalEntity::GetParsedFormula()
    {

        auto CountAtomAsText = [this](const std::string& text, int scale) -> std::string {
            std::string res;
            auto atoms = GetAtomAsString(text);
            for(auto& atom: atoms)
                res += atom.first + std::to_string(std::stoi(atom.second.empty() ? "1" : atom.second) * scale);
            
            return res;
        };

        std::string atomFlatText;
        std::stack<std::string> currentNodes;
        std::stack<std::string> lastNodeData;

        auto formula = "((" + this->formula + "))";

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



    vector_t<std::string> Split(const std::string& str, char delimeter)
    {
        std::vector<std::string> res;
        std::string text = str;
        size_t n;
        while((n = text.find(delimeter)) < text.size()) {
            res.push_back(text.substr(0, n));
            text = text.substr(n + 1, text.size());
        }
        res.push_back(text);
        return res;
    }

}
