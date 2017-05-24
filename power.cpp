#include "power.hpp"



bool overlaps(const chromspan& s, const chromspan& cs) {
    int a = s.start;
    int b = s.stop;

    int c = cs.start;
    int d = cs.stop;

    return ((a <= c && c <= b) || (c <= a && a <= d));
}

int overlap_amount(const chromspan& a, const chromspan& b) {
    return std::min(a.stop, b.stop) - std::max(a.start, b.start);
}

namespace adios {

Indpair dummy_indpair(const Dataset& d, int chromidx, const chromspan& cs) {
    int ninds = d.ninds();

    // Step 1: select 3 random individuals
    int aidx, bidx, tidx;
    aidx = randint(0, ninds-1);

    do {
        bidx = randint(0, ninds-1);
    } while (bidx == aidx);
    
    do {
        tidx = randint(0, ninds-1);
    } while (tidx == aidx || tidx == bidx);

    // Step 2: make copies of two of them
    Individual ind_a = d.individuals[aidx];
    Individual ind_b = d.individuals[bidx];
    const Individual& ind_template = d.individuals[tidx];



    // Step 3: from the third individual copy a genome region to one chromosome of the
    //         two others
    int shared_hapidx = randint(0, 1);
    copy_genospan(ind_template, shared_hapidx, ind_a, randint(0, 1), chromidx, cs);
    copy_genospan(ind_template, shared_hapidx, ind_b, randint(0, 1), chromidx, cs);

    // Step 4: return the pair

    Indpair pr = {ind_a, ind_b};
    return pr;
}


chromspan random_range(const std::shared_ptr<ChromInfo>& c, int size) {
    // Step 1: randomly select a start point
    int chromsize = c->size();
    int start = randint(0, chromsize - size);
    // Step 2: return start + size
    chromspan cs = {start, start + size};
    return cs;
}


PowerReplicateResult power_replicate(const Dataset& d,
                                     const adios_parameters& params,
                                     int chromidx,
                                     unsigned int segsize) {
    PowerReplicateResult prr;

    // Step 1: Pick a span;
    chromspan synthetic = random_range(d.chromosomes[chromidx], segsize);
    // Step 2: create the dummy indpair
    Indpair dummy_pair = dummy_indpair(d, chromidx, synthetic);

    // Step 3: Run adios_pair_unphased on the pair
    adios_result res = adios_pair_unphased(dummy_pair.ind1, 
                                           dummy_pair.ind2, 
                                           chromidx, 
                                           params);

    auto& chrominf = d.chromosomes[chromidx];
    std::vector<chromspan> result_spans;
    for (auto& s : res.segments) {
        chromspan rs = { 
                         chrominf->positions[s.full_start], 
                         chrominf->positions[s.full_stop]
                       };
        
        result_spans.push_back(rs);
    }

    // Step 4: Find segments overlapping the true segment
    int oa = 0;

    for (int i = 0; i < res.segments.size(); ++i) {
        auto& rs = result_spans[i];
        if (overlaps(rs, synthetic)) {
            prr.segments.push_back(res.segments[i]);
            oa += overlap_amount(rs, synthetic);
        }
    }

    prr.success = prr.segments.size() > 0;
    prr.nseg = prr.segments.size();
    prr.overlap_amount = oa;

    return prr;
}

PowerResult calc_power(Dataset& d, const adios_parameters& params,
                int chromidx,
                unsigned int segsize,
                unsigned int nrep) {

    PowerResult results;
    results.segsize = segsize;

    #pragma omp parallel for 
    for (unsigned long int i = 0; i < nrep; ++i) {
        PowerReplicateResult res = power_replicate(d, params, chromidx, segsize);
        results.replicates.push_back(res);
    }

    return results;
}

double PowerResult::mean_num_segments(void) const {
    int i = 0;
    int n = 0;
    for (auto& rep : replicates) {
        n += rep.success;
        i += rep.success ? rep.segments.size() : 0; 
    }

    return (double)i / n;
}

std::vector<int> PowerResult::length_diffs(void) const {
    std::vector<int> diffs;

    for (auto& rep : replicates) {
        if (rep.success) {
            diffs.push_back(rep.segments[0].length() - segsize);
        }
    }

    return diffs;
}

double PowerResult::prop_detected(void) const {
    long long int detected  = 0;
    long long int n = segsize * replicates.size();

    for (auto& rep : replicates) {
        if (rep.success) {
            detected += rep.overlap_amount;
        }
    }

    return detected / (double)n;
}

}