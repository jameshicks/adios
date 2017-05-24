#include "datamodel.hpp"

//Variant

Variant::Variant(const std::string& lab, int bp, double maf) {
    label.assign(lab);
    position = bp;
    frequency = maf;
}

//ChromInfo
void ChromInfo::add_variant(const std::string& lab, int bp, double maf) {
    auto v = Variant(lab, bp, maf);
    positions.push_back(bp);
    frequencies.push_back(maf);
    variants.push_back(v);
}

size_t ChromInfo::nmark(void) const {
    return positions.size();
}

int ChromInfo::size(void) const {
    if (!nmark()) return 0;
    return positions.back() - positions.front();
}


ChromInfo::ChromInfo(void) {
    label.assign("");
}

ChromInfo::ChromInfo(const std::string& lab) {
    label.assign(lab);
}


// Genotypes
std::vector<int> Genotypes::todense(int haplotype) {
    int n = info->nmark();
    std::vector<int> outp(n);
    auto selected = haplotype ? hapb : hapa;
    for (auto it = selected.begin(); it != selected.end(); ++it) {
        outp[*it] = 1;
    }
    return outp;
}


int Genotypes::get_minor_allele_count(int markidx) {
    int a = std::find(hapa.begin(), hapa.end(), markidx) != hapa.end();
    int b = std::find(hapb.begin(), hapb.end(), markidx) != hapb.end();

    return a + b;
}

AlleleSites Genotypes::heterozygous(void) const {
    return setops::symmetric_difference(hapa, hapb);
}

AlleleSites Genotypes::homozygous_minor(void) const {
    return setops::intersection(hapa, hapb);
}

AlleleSites Genotypes::has_minor_allele(void) const {
    return setops::union_(hapa, hapb);
}

std::vector<int> Genotypes::dosages(void) {
    std::vector<int> outp(info->nmark());

    for (auto i : hapa) {outp[i] += 1;}
    for (auto i : hapb) {outp[i] += 1;}
 

    return outp;
}
Genotypes::Genotypes(shared_ptr<ChromInfo> c) {
    info = c;
}

// Individual

Individual::Individual(void) {
    label.assign("");
}

Individual::Individual(const std::string& lab) {
    label.assign(lab);
}

Individual::Individual(const Individual& ind) {
    label = ind.label;
    chromosomes = ind.chromosomes;
}

void Individual::add_empty_chromosome(shared_ptr<ChromInfo> inf) {
    Genotypes g(inf);
    chromosomes.push_back(g);
}

void Individual::get_empty_chromosomes(const Dataset& d) {
    for (shared_ptr<ChromInfo> inf : d.chromosomes) add_empty_chromosome(inf);
}

int Individual::get_minor_allele_count(int chromidx, int markidx) {
    return chromosomes[chromidx].get_minor_allele_count(markidx);
}

// Dataset

size_t Dataset::ninds(void) const {
    return individuals.size();
}

size_t Dataset::nchrom(void) const {
    return chromosomes.size();
}

size_t Dataset::nmark(void) const {
    size_t tot = 0;
    for (size_t i = 0; i < nchrom(); ++i) {
        tot += chromosomes[i]->nmark();
    }
    return tot;
}

size_t Dataset::nexcluded(void) const {
    size_t tot = 0;
    for (size_t i = 0; i < nchrom(); ++i) {
        auto c = chromosomes[i];
        for (auto kv : c->exclusions) {
            tot += kv.second;
        }
    }
    return tot;
}

void Dataset::add_individual(const std::string& lab) {
    individuals.push_back(Individual(lab));

}

void Dataset::add_chromosome(const std::string& lab) {
    shared_ptr<ChromInfo> c(new ChromInfo(lab));
    chromosomes.push_back(c);

    for (Individual& ind : individuals) {
        ind.add_empty_chromosome(c);
    }
 
}

void Dataset::round_frequencies(unsigned int places) {
    for (auto c : chromosomes) {
        for (size_t i=0; i<c->frequencies.size(); ++i) {
            c->frequencies[i] = round(c->frequencies[i], places);
        }
    }
}

void Dataset::floor_frequencies(double floor) {
    for (auto c : chromosomes) {
        for (size_t i=0; i<c->frequencies.size(); ++i) {
            double fq = c->frequencies[i];
            c->frequencies[i] = fq < floor ? floor : fq;
        }
    }
}



void copy_genospan(const Individual& from, int hapfrom, 
                          Individual& to, int hapto,
                          int chromidx, const chromspan& cs) {
    using std::distance;
    using std::lower_bound;
    using std::upper_bound;

    ChromInfo& info = *(from.chromosomes[chromidx].info);


    int startidx = distance(info.positions.begin(),
                            lower_bound(info.positions.begin(), 
                                        info.positions.end(), 
                                        cs.start)); 
    
    int stopidx = distance(info.positions.begin(),
                           lower_bound(info.positions.begin(), 
                                       info.positions.end(), 
                                       cs.stop));
    
    AlleleSites templatechrom = hapfrom ? 
                                from.chromosomes[chromidx].hapb : 
                                from.chromosomes[chromidx].hapa;

    AlleleSites newchrom;
    AlleleSites oldchrom = hapto ? 
                           to.chromosomes[chromidx].hapb : 
                           to.chromosomes[chromidx].hapa;
    int maxold = oldchrom.back();
    int i = 0;
    while (oldchrom[i] < startidx) {
        newchrom.push_back(oldchrom[i]);
        i++;
    }

    auto copystartit = std::lower_bound(templatechrom.begin(),
                                        templatechrom.end(),
                                        startidx);
    auto copystopit = std::upper_bound(templatechrom.begin(),
                                       templatechrom.end(),
                                       stopidx);

    for (auto shared_it = copystartit; shared_it != copystopit; shared_it++) 
        newchrom.push_back(*shared_it);

    for (auto it = std::lower_bound(oldchrom.begin(), oldchrom.end(), stopidx);
            it != oldchrom.end();
            ++it) {
        newchrom.push_back(*it);
    }

    if (hapto) {
        to.chromosomes[chromidx].hapb = newchrom;
    } else {
        to.chromosomes[chromidx].hapa = newchrom;
    }
}

