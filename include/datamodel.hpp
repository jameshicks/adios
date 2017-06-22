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

#include "utility.hpp"
#include "setops.hpp"

using std::shared_ptr;
// #define protected public
// #define private public

typedef std::vector<int> AlleleSites;


struct chromspan {
    int start;
    int stop;
    inline int length(void) const { return stop - start; }
};


class Dataset;



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

    AlleleSites het; 
    AlleleSites hzm;
    AlleleSites hma;

    shared_ptr<ChromInfo> info;

    std::vector<int> todense(int haplotype);
    void set_allele(int markidx, int hapidx, int allele);
    int get_minor_allele_count(int markidx);
    AlleleSites heterozygous(void) const;
    AlleleSites homozygous_minor(void) const;
    AlleleSites has_minor_allele(void) const;
    void finalize(void);

    std::vector<int> dosages(void);
    
    Genotypes(shared_ptr<ChromInfo> c);

};

class Individual
{
public:
    std::string label;
    std::vector<Genotypes> chromosomes;

    void add_empty_chromosome(shared_ptr<ChromInfo> inf);
    void get_empty_chromosomes(const Dataset& d);
    void set_allele(int chromidx, int markidx, int hapidx, int allele);
    int get_minor_allele_count(int chromidx, int markidx);
    void finalize(void);

    Individual(void);
    Individual(const std::string& lab);
    // ~Individual(void);
    Individual(const Individual& ind);
};

class Dataset
{
public:
    std::vector<Individual> individuals;
    std::vector<shared_ptr<ChromInfo>> chromosomes;
    
    // Number of individuals
    size_t ninds(void) const;

    // Number of chromosomes
    size_t nchrom(void) const;

    // Number of retained markers 
    size_t nmark(void) const;

    // Number of excluded markers
    size_t nexcluded(void) const;
    void add_individual(const std::string& label);
    void add_chromosome(const std::string& label);
    void round_frequencies(unsigned int places);
    void floor_frequencies(double floor);
    void subset(std::set<std::string> indlabs);
    void finalize(void);

};

void copy_genospan(const Individual& from, int hapfrom, Individual& to, int hapto, 
    int chromidx, const chromspan& cs);

// Inline funcs

inline void Individual::set_allele(int chromidx, int markidx, int hapidx, int allele)
{
    chromosomes[chromidx].set_allele(markidx, hapidx, allele);
}

inline void Genotypes::set_allele(int markidx, int hapidx, int allele)
{
    if (allele < 0) {
        missing.push_back(markidx);
    } else {
        AlleleSites& chromatid = hapidx ? hapb : hapa;
        chromatid.push_back(markidx);
    }
}



void add_error(Individual& ind, double error_rate);

#endif
