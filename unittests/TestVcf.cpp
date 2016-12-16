#include "datamodel.hpp"
#include "vcf.hpp"
#include "CppUTest/TestHarness.h"

TEST_GROUP(VCF) {};

TEST(VCF, ReadVCF) {
    VCFParams vcfp = {false, false, false, "AF"};
    Dataset d = read_vcf("unittests/test.vcf", vcfp);
    CHECK_EQUAL(3, d.ninds());
    CHECK_EQUAL(2, d.nchrom());
    CHECK_EQUAL(6, d.nmark());

    CHECK_EQUAL("1", d.chromosomes[0]->label);
    CHECK_EQUAL(1, d.chromosomes[0]->nmark());
    std::vector<double> expected_freqs = {0.5};
    auto observed_freqs = d.chromosomes[0]->frequencies;
    CHECK(expected_freqs == observed_freqs);

    CHECK_EQUAL("20", d.chromosomes[1]->label);
    CHECK_EQUAL(5, d.chromosomes[1]->nmark());
    expected_freqs = {0.5, 0.017,0,0,0};
    CHECK(expected_freqs == d.chromosomes[1]->frequencies);

};

TEST(VCF, VCFEmpiricalFrequency) {
    VCFParams vcfp = {false, false, true, "DOESNT MATTER"};
    Dataset d = read_vcf("unittests/test.vcf", vcfp);
    CHECK_EQUAL(3, d.ninds());
    CHECK_EQUAL(2, d.nchrom());
    CHECK_EQUAL(6, d.nmark());

    CHECK_EQUAL("1", d.chromosomes[0]->label);
    CHECK_EQUAL(1, d.chromosomes[0]->nmark());
    std::vector<double> expected_freqs = {0.5};
    auto observed_freqs = d.chromosomes[0]->frequencies;
    CHECK(expected_freqs == observed_freqs);

    CHECK_EQUAL("20", d.chromosomes[1]->label);
    CHECK_EQUAL(5, d.chromosomes[1]->nmark());
    expected_freqs = {3/6.0, 1/6.0, 0, 1-(4/6.0), 0};
    observed_freqs = d.chromosomes[1]->frequencies;
    CHECK(expected_freqs == observed_freqs);

}

TEST(VCF, MinorAlleleInvert) {
    VCFRecordGenotypeContainer con(3);
    
    // std::string line = "1   10  rs100   A   C   100 PASS    AF=.5   GT  0/0 0/0 0/0";
    // VCFRecord rec(line);

    con.alts = {0,1,4,5};
    con.invert();

    std::vector<size_t> expected = {2,3};
    CHECK(con.alts == expected);


}