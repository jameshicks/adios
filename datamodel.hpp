#ifndef DATAMODEL_HPP
#define DATAMODEL_HPP

#include <vector>
#include <map>
#include <set>
#include <unordered_set>
#include <string>
#include <memory>
#include <algorithm>
#include <iterator>
#include <iostream>

#include "setops.hpp"

using std::shared_ptr;
// #define protected public
// #define private public

// typedef std::vector<int> AlleleSites;
typedef std::vector<int> AlleleSites;

class Variant
{
public:
    std::string label;
    int position;
    double frequency;
    Variant(const std::string& lab, int bp, double maf);
};

class ChromInfo
{
public:
    std::vector<Variant> variants;
    std::string label;
    std::map<std::string, size_t> exclusions;
    void add_variant(const std::string& lab, int bp, double maf);
    size_t nmark(void) const;
    int size(void) const;
    std::vector<std::string> variant_labels;
    std::vector<int> positions;
    std::vector<double> frequencies;


    ChromInfo(void);
    ChromInfo(const std::string& lab);
    // ~ChromInfo(void);
};

class Genotypes
{
public:
    AlleleSites hapa;
    AlleleSites hapb;
    AlleleSites missing;
    shared_ptr<ChromInfo> info;

    std::vector<int> todense(int haplotype);
    void set_allele(int markidx, int hapidx, int allele);
    int get_minor_allele_count(int markidx);
    AlleleSites heterozygous(void) const;
    AlleleSites homozygous_minor(void) const;
    AlleleSites has_minor_allele(void) const;
    std::vector<int> dosages(void);
    Genotypes(shared_ptr<ChromInfo> c);

};

class Individual
{
public:
    std::string label;
    std::vector<shared_ptr<Genotypes>> chromosomes;

    void add_empty_chromosome(shared_ptr<ChromInfo> inf);
    void set_allele(int chromidx, int markidx, int hapidx, int allele);
    int get_minor_allele_count(int chromidx, int markidx);
    Individual(void);
    Individual(const std::string& lab);
    // ~Individual(void);
};

class Dataset
{
public:
    std::map<std::string, shared_ptr<Individual>> individuals;
    std::vector<shared_ptr<ChromInfo>> chromosomes;
    size_t ninds(void) const;
    size_t nchrom(void) const;
    size_t nmark(void) const;
    // int ngeno(void);
    shared_ptr<Individual> add_individual(const std::string& label);
    shared_ptr<ChromInfo> add_chromosome(const std::string& label);

};

// Inline funcs

inline void Individual::set_allele(int chromidx, int markidx, int hapidx, int allele)
{
    chromosomes[chromidx]->set_allele(markidx, hapidx, allele);
}

inline void Genotypes::set_allele(int markidx, int hapidx, int allele)
{
    if (allele < 0) {
        missing.push_back(markidx);
    } else {
        AlleleSites& chromatid = hapidx ? hapa : hapb;
        chromatid.push_back(markidx);
    }
}

#endif
