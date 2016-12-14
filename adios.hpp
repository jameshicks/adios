#ifndef ADIOS_HPP
#define ADIOS_HPP

#include <map>
#include <set>
#include <utility> // for std::pair 
#include <cmath> // For log
#include <memory>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>



#include "combinatorics.hpp"
#include "Linalg.hpp"
#include "HiddenMarkov.hpp"
#include "datamodel.hpp"
#include "common.hpp"
// using AlleleSites;

namespace adios {


typedef std::pair<std::vector<int>, std::vector<int>> adios_sites;
typedef std::shared_ptr<Individual> Indptr;
typedef std::pair<Indptr, Indptr> Indptr_pair;
typedef std::shared_ptr<ChromInfo> Chromptr;
using std::pair;

// Parameters for ADIOS.
struct adios_parameters {
    double rare_thresh;                             // Rare variant frequency threshold 
    double err_rate;                                // Allele error rate
    Matrix allele_error_mat;                        // Matrix of genotyping error probabilitiez
    Matrix transition_mat;                          // HMM transition matrix
    std::vector<std::vector<int>> rare_sites;       // The set of sites with rare variation
    int gamma_;                                     // Probability to enter IBD (10^(-gamma))
    int rho;                                        // Probability of exitiing IBD (10^(-rho))
    int min_length;                                 // Minimum segment length of to consider
    size_t min_mark;                                // Minimum number of markers allowed in a segment
    double min_lod;                                 // Minimum allowed quality score
    void get_rare_sites(Dataset& data);             // Get the rare sites
    void calculate_emission_mats(const Dataset& d); // Precompute emission matrices
    std::map<double, Matrix> emission_mats;         // Precomputer emission mats
};


class Segment {
private:
    Indptr ind1;  // The first individual
    Indptr ind2;  // The second individual
    Chromptr chrom;  // Chromosome info

    size_t start;
    size_t stop;
public:

    size_t full_start;
    size_t full_stop;
    int state;

    size_t nmark;
    double lod;


    inline int length(void) const {
        return chrom->positions[full_stop] - chrom->positions[full_start];
    }

    Segment(Indptr a, Indptr b, ValueRun& run, Chromptr c,
            std::vector<int>& obs, std::vector<Matrix>& emissions,
            std::vector<int>& adiossites, const adios_parameters& params);

    void trim(std::vector<int>& observations);
    double calculate_lod(std::vector<int>& observations,
                         std::vector<Matrix>& emissions,
                         const adios::adios_parameters& params) const;
    bool passes_filters(const adios::adios_parameters& params) const;
    std::string record_string(void) const;

};

// Make the model parameters from the command line args 
adios_parameters params_from_args(std::map<std::string, std::vector<std::string>> args);

// Returns true if the observed genotype configuration includes a
// shared rare variant.
inline bool is_shared_rv(int obs) { return ((obs >= 4) && (obs != 6)); }

// Creates the transition matrix for the HMM 
Matrix transition_matrix(unsigned int l10gamma, unsigned int l10rho);

// A matrix describing the probability of a having a genotype
// given the observed genotype using an allele miscall error rate (eps)
Matrix genotype_error_matrix(double eps);

// Creates the HMM emission matrix for a genotype with minor allele frequency q
Matrix emission_matrix(double q);

// ADIOS works on a subset of markers. This function translates
// the indices of the markers used back to the full indices.
AlleleSites update_indices(const AlleleSites& inp,
                           const std::map<int, int>& translator);

// For two individuals, find a set of informative sites
adios_sites find_informative_sites(const Indptr& ind1,
                                   const Indptr& ind2,
                                   int chromidx,
                                   const AlleleSites& rares);

 
// Perform adios on the entire dataset d using parameters `params`
void adios(Dataset& d, const adios_parameters& params);

// Perform adios on a pair of individuals on one chromosome
void adios_pair(const Indptr_pair& pair,
                int chromidx,
                const adios_parameters& params);
};

#endif