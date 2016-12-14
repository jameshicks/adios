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
    variant_labels.push_back(lab);
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


void Individual::add_empty_chromosome(shared_ptr<ChromInfo> inf) {
    shared_ptr<Genotypes> g(new Genotypes(inf));
    chromosomes.push_back(g);
}

int Individual::get_minor_allele_count(int chromidx, int markidx) {
    return chromosomes[chromidx]->get_minor_allele_count(markidx);
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

shared_ptr<Individual> Dataset::add_individual(const std::string& lab) {
    shared_ptr<Individual> ind(new Individual(lab));
    individuals[lab] = ind;
    return ind;
}

shared_ptr<ChromInfo> Dataset::add_chromosome(const std::string& lab) {
    shared_ptr<ChromInfo> c(new ChromInfo(lab));
    chromosomes.push_back(c);

    for (auto it = individuals.begin(); it != individuals.end(); ++it) {
        it->second->add_empty_chromosome(c);
    }
    return c;
}
