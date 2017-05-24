#include "datamodel.hpp"
#include "adios.hpp"
#include "vcf.hpp"
#include "Linalg.hpp"
#include "math.h"
#include "CppUTest/TestHarness.h"

#include <iostream>
TEST_GROUP(adios) {};

TEST(adios, InformativeSites) {
    VCFParams vcfp = {false, false, false, "AF"};

    Dataset d = read_vcf("unittests/data/vcf/test_informative_sites.vcf", vcfp);
    auto ind1 = d.individuals[0];
    auto ind2 = d.individuals[0];

    std::vector<int> rares;

    for (int i = 0; i < d.chromosomes[0]->nmark(); ++i) {
        if (d.chromosomes[0]->frequencies[i] < 0.05) { rares.push_back(i); }
    }

    auto p = adios::find_informative_sites_unphased(ind1, ind2, 0, rares);
    // std::vector<int> expected_sites = {2, 6, 11, 13, 14, 15, 16, 17};
    std::vector<int> expected_sites = {2,6,10, 11, 12, 13, 14, 15, 16, 17};
    auto observed_sites = p.second;
    CHECK(expected_sites == observed_sites);
    // std::vector<int> expected_states = {2, 6, 2, 4, 5, 6, 7, 8};
    // auto observed_states = p.first;
    // CHECK(observed_states == expected_states);
};

// TEST(adios, TransitionMatrix) {
//     double gamma = pow(10,-4);
//     double rho = pow(10,-2);
//     Matrix expected = {
//         {9.9999e-01, 1e-5,       1e-10},
//         {1e-2,       9.8999e-01, 1e-05},
//         {1e-4,       1e-2,       9.899e-01}
//     };
//     Matrix observed = adios::unphased_transition_matrix(10);
//     CHECK((expected - observed).sum() < 1e-6);
// }

TEST(adios, UnphasedEmissionMatrix) {
    Matrix expected = Matrix( {
        {  9.96005996e-01,   9.97002999e-01,   9.98001000e-01},
        {  9.97002999e-04,   9.98001000e-04,   0.00000000e+00},
        {  9.98001000e-07,   0.00000000e+00,   0.00000000e+00},
        {  9.97002999e-04,   9.98001000e-04,   0.00000000e+00},
        {  3.99200400e-06,   9.99000000e-04,   1.99800000e-03},
        {  9.99000000e-10,   9.99000000e-07,   0.00000000e+00},
        {  9.98001000e-07,   0.00000000e+00,   0.00000000e+00},
        {  9.99000000e-10,   9.99000000e-07,   0.00000000e+00},
        {  1.00000000e-12,   1.00000000e-09,   1.00000000e-06}
    });
    Matrix observed = adios::unphased_emission_matrix(0.001);
    CHECK((expected - observed).sum() < 1e-6);

}

TEST(adios, GenotypeErrorMatrix) {
    Matrix expected = {
        {  9.98001000e-01,   1.99800000e-03,   1.00000000e-06},
        {  9.99000000e-04,   9.98002000e-01,   9.99000000e-04},
        {  1.00000000e-06,   1.99800000e-03,   9.98001000e-01}
    };
    Matrix observed = adios::unphased_genotype_error_matrix(0.001);
    CHECK((expected - observed).sum() < 1e-6);
}