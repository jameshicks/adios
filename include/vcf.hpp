#ifndef VCF_HPP
#define VCF_HPP
#include <iostream>
#include <fstream>


#include <string>
#include <vector>
#include <map>
#include <utility>
#include <memory>
#include <numeric>
#include <string.h>
#include <stdlib.h>
#include <iterator>
#include <algorithm>


#include "common.hpp"
#include "stringops.hpp"
#include "datamodel.hpp"
#include "FileIOManager.hpp"


// Receives variant calls from VCFRecord
struct VCFRecordGenotypeContainer {
    size_t ninds; 
    std::vector<size_t> missing;
    std::vector<size_t> alts;
    VCFRecordGenotypeContainer(size_t n);

    // Set minor alleles to major and vice-versa. Ignores missing sites.
    void invert(void);

    // Calculates frequency of minor allele (missings excluded)
    double allele_frequency(void) const;
    inline bool monomorphic(void) const { return alts.size() == 0; }
    inline bool singleton(void) const { return alts.size() == 1; }
    inline void clear(void) { missing.clear(); alts.clear(); }


};

struct VCFParams {
    bool drop_singletons;
    bool drop_monomorphs;

    // Calculate allele frequencies from data?
    bool empirical_freqs;

    // VCF Info field to retrieve frequency from if not 
    // calculating empirical frequencies.
    std::string freq_field;
};

class VCFRecord {
public:
    std::string chrom;
    int pos;
    double freq;
    std::string label;
    std::string format;
    std::vector<std::string> alleles;
    std::string data; 
    std::string infostr;

    VCFRecord(const std::string& line);
    std::map<std::string, std::string> infomap(void) const;
    std::string get_info_by_key(const char* key);
    void get_minor_alleles(VCFRecordGenotypeContainer& container) const;
    
    // Returns the frequency specified from an INFO field, otherwise 0.0;
    double get_info_freq(const std::string& info_field);
    inline int nalleles(void) const;
    inline bool is_snv(void) const;
};



Dataset read_vcf(const std::string& filename, const VCFParams& fileparams);
#endif
