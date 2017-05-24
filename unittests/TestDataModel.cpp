#include "datamodel.hpp"
#include "vcf.hpp"
#include "CppUTest/TestHarness.h"

#include <iostream>
TEST_GROUP(DataModel) {};

TEST(DataModel, GenotypeMethods) {
    VCFParams vcfp = {false, false, false, "AF"};

    Dataset d = read_vcf("unittests/data/vcf/test2.vcf", vcfp);

    std::map<std::string, Individual> dinds;
    for (auto& ind : d.individuals) { dinds[ind.label] = ind; }

    std::vector<int> expected_mac{0,1,0,0,1};
    auto observed_mac = dinds["A"].chromosomes[0].dosages();

    CHECK(expected_mac == observed_mac);
    
    expected_mac = {0,2,0,0,2};
    observed_mac = dinds["B"].chromosomes[0].dosages();
    CHECK(expected_mac == observed_mac);

    expected_mac = {0,1,1,0,0};
    observed_mac = dinds["C"].chromosomes[0].dosages();
    CHECK(expected_mac == observed_mac);    

    std::vector<int> exp_hma_a = {1,4};
    auto observed_hma_a = dinds["A"].chromosomes[0].has_minor_allele();
    CHECK(exp_hma_a == observed_hma_a); 

    std::vector<int> exp_hma_b = {1,4};
    auto observed_hma_b = dinds["B"].chromosomes[0].has_minor_allele();
    CHECK(exp_hma_b == observed_hma_b); 

    std::vector<int> exp_hma_c = {1,2};
    auto observed_hma_c = dinds["C"].chromosomes[0].has_minor_allele();
    CHECK(exp_hma_c == observed_hma_c); 

    std::vector<int> exp_het_a = {1,4};
    auto observed_het_a = dinds["A"].chromosomes[0].heterozygous();
    CHECK(exp_het_a == observed_het_a); 

    std::vector<int> exp_het_b = {};
    auto observed_het_b = dinds["B"].chromosomes[0].heterozygous();
    CHECK(exp_het_b == observed_het_b); 

    std::vector<int> exp_het_c = {1,2};
    auto observed_het_c = dinds["C"].chromosomes[0].heterozygous();
    CHECK(exp_het_c == observed_het_c); 

    std::vector<int> exp_hzm_a = {};
    auto observed_hzm_a = dinds["A"].chromosomes[0].homozygous_minor();
    CHECK(exp_hzm_a == observed_hzm_a); 


    std::vector<int> exp_hzm_b = {1,4};
    auto observed_hzm_b = dinds["B"].chromosomes[0].homozygous_minor();
    CHECK(exp_hzm_b == observed_hzm_b); 

    std::vector<int> exp_hzm_c = {};
    auto observed_hzm_c = dinds["C"].chromosomes[0].homozygous_minor();
    CHECK(exp_hzm_c == observed_hzm_c); 


    auto gta = dinds["A"].chromosomes[0];
    std::vector<int> exp_hap;
    std::vector<int> obs;
    exp_hap = {};
    obs = dinds["A"].chromosomes[0].hapa;
    CHECK(obs == exp_hap);
    
    exp_hap = {1,4};
    obs = dinds["A"].chromosomes[0].hapb;
    CHECK(obs == exp_hap);

    exp_hap = {1,4};
    auto gtb = dinds["B"].chromosomes[0];
    CHECK(gtb.hapa == exp_hap);
    CHECK(gtb.hapb == exp_hap);

    auto gtc = dinds["C"].chromosomes[0];
    exp_hap = {};
    CHECK(gtc.hapa == exp_hap);
    exp_hap = {1,2};
    CHECK(gtc.hapb == exp_hap);
};