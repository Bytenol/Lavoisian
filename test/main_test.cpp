#include <gtest/gtest.h>
#include <iostream>
#include <vector>
#include <string>
#include <chm/Core.h>

using namespace chm;

TEST(AtomSetup, CoreTest)
{
    bool dbOpen = true;
    if(!Init("C:/Users/olami/Documents/Programming/VCS/Chemilator/test/test.db"))
        dbOpen = false;

    ASSERT_EQ(dbOpen, true) << " SQLite3 Database could not be opened ";
    EXPECT_EQ(Atom::Get().size(), 119u) << " Atom size is not up to 119 ";
}


TEST(AtomSearch, CoreTest)
{
    EXPECT_EQ(Atom::Get(0)->symbol, "H") << " Atom index does not match the symbol ";
    EXPECT_EQ(Atom::Get(1)->symbol, "He") << " Atom index does not match the symbol ";
    EXPECT_EQ(Atom::Get(0)->period, 1) << " Atom index does not match the period ";
    EXPECT_EQ(Atom::Get(1)->named_by, "NULL") << " Atom index does not match the named_by ";
    EXPECT_EQ(Atom::Get(0)->atomic_mass, 1.008) << " Atom index does not match the atomic_mass ";
    EXPECT_EQ(Atom::Get(-1), nullptr);
    EXPECT_EQ(Atom::Get(199), nullptr);

    EXPECT_NE(Atom::Get("Na"), nullptr) << " Atom by name does not match ";
    EXPECT_EQ(Atom::Get("Xyl"), nullptr) << " Atom by name does not match ";
}


TEST(ChemicalEntitySetup, ChemicalEntityTest)
{
    ChemicalEntity h2{ "H2" };
    EXPECT_EQ(h2.charge, 0) << " Atom does not have a charge ";
}


TEST(IonChargeMatches, ChemicalEntityTest)
{
    std::vector<std::string> freeAtoms{ "Na", "Mg+", "H2", "HSO4{}", "CO3{+}", "K2CR2O7{-}" };
    for(auto& cases: freeAtoms)
    {
        ChemicalEntity chm_test(cases);
        EXPECT_EQ(chm_test.charge, 0) << " Ion Must have a Charge ";
        EXPECT_NE(chm_test.IsIon(), true);
    }
        

    EXPECT_EQ((ChemicalEntity("Na{+1}")).charge, 1) << " Ion charge mismatch ";
    EXPECT_EQ((ChemicalEntity("Mg{+2}")).charge, 2) << " Ion charge mismatch ";
    EXPECT_EQ((ChemicalEntity("Al{+3}")).charge, 3) << " Ion charge mismatch ";
    EXPECT_EQ((ChemicalEntity("O{-2}")).charge, -2) << " Ion charge mismatch ";
    EXPECT_EQ((ChemicalEntity("F{-1}")).charge, -1) << " Ion charge mismatch ";
    EXPECT_EQ((ChemicalEntity("Xxx{-1234}")).charge, -1234) << " Ion charge mismatch ";

    ChemicalEntity xxx("Xxx{+1234}");
    EXPECT_EQ(xxx.charge, 1234) << " Ion charge mismatch ";
    EXPECT_EQ(xxx.IsIon(), true);
}



TEST(MolecularFormula, ChemicalEntityTest)
{
    ChemicalEntity ethanol{"CH2COOH"};
    EXPECT_EQ(ethanol.GetAtom("H").count, 3);

    ChemicalEntity naoh{"Pb(NO3)2"};
    EXPECT_EQ((int)naoh.GetMolarMass(), 331);
}