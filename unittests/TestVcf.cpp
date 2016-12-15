#include "datamodel.hpp"
#include "vcf.hpp"
#include "CppUTest/TestHarness.h"

TEST_GROUP(VCF) {};

TEST(VCF, ReadVCF) {
    VCFParams vcfp = {false, false, false, "AF"};
    Dataset d = read_vcf("test.vcf", vcfp);
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
    Dataset d = read_vcf("test.vcf", vcfp);
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
    expected_freqs = {3/6.0, 1/6.0, 0, 4/6.0, 1.0};
    observed_freqs = d.chromosomes[1]->frequencies;
    CHECK(expected_freqs == observed_freqs);

}