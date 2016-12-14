#include "datamodel.hpp"
#include "vcf.hpp"
#include "CppUTest/TestHarness.h"

#include <iostream>
TEST_GROUP(DataModel) {};

TEST(DataModel, GenotypeMethods) {
    Dataset d = read_vcf("unittests/test2.vcf", "AF");

    std::vector<int> expected_mac{0,1,0,2,1};
    auto observed_mac = d.individuals["A"]->chromosomes[0]->dosages();

    CHECK(expected_mac == observed_mac);
    
    expected_mac = {0,2,0,2,2};
    observed_mac = d.individuals["B"]->chromosomes[0]->dosages();
    CHECK(expected_mac == observed_mac);

    expected_mac = {0,1,1,2,0};
    observed_mac = d.individuals["C"]->chromosomes[0]->dosages();
    CHECK(expected_mac == observed_mac);    

    std::vector<int> exp_hma_a = {1,3,4};
    auto observed_hma_a = d.individuals["A"]->chromosomes[0]->has_minor_allele();
    CHECK(exp_hma_a == observed_hma_a); 

    std::vector<int> exp_hma_b = {1,3,4};
    auto observed_hma_b = d.individuals["B"]->chromosomes[0]->has_minor_allele();
    CHECK(exp_hma_b == observed_hma_b); 

    std::vector<int> exp_hma_c = {1,2,3};
    auto observed_hma_c = d.individuals["C"]->chromosomes[0]->has_minor_allele();
    CHECK(exp_hma_c == observed_hma_c); 

    std::vector<int> exp_het_a = {1,4};
    auto observed_het_a = d.individuals["A"]->chromosomes[0]->heterozygous();
    CHECK(exp_het_a == observed_het_a); 

    std::vector<int> exp_het_b = {};
    auto observed_het_b = d.individuals["B"]->chromosomes[0]->heterozygous();
    CHECK(exp_het_b == observed_het_b); 

    std::vector<int> exp_het_c = {1,2};
    auto observed_het_c = d.individuals["C"]->chromosomes[0]->heterozygous();
    CHECK(exp_het_c == observed_het_c); 

    std::vector<int> exp_hzm_a = {3};
    auto observed_hzm_a = d.individuals["A"]->chromosomes[0]->homozygous_minor();
    CHECK(exp_hzm_a == observed_hzm_a); 


    std::vector<int> exp_hzm_b = {1,3,4};
    auto observed_hzm_b = d.individuals["B"]->chromosomes[0]->homozygous_minor();
    CHECK(exp_hzm_b == observed_hzm_b); 

    std::vector<int> exp_hzm_c = {3};
    auto observed_hzm_c = d.individuals["C"]->chromosomes[0]->homozygous_minor();
    CHECK(exp_hzm_c == observed_hzm_c); 

};