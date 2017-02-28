#include "power.hpp"

int randint(int lo, int hi) {
    // Produce a random interval in the closed interval [lo, hi]
    return lrand48() % (hi + 1 - lo) + lo;
}

bool overlaps(const adios::chromspan& s, const adios::chromspan& cs) {
    int a = s.start;
    int b = s.stop;

    int c = cs.start;
    int d = cs.stop;

    return ((a <= c && c <= b) || (c <= a && a <= d));
}

int overlap_amount(const adios::chromspan& a, const adios::chromspan& b) {
    return std::min(a.stop, b.stop) - std::max(a.start, b.start);
}

namespace adios {

Indpair dummy_indpair(const Dataset& d, int chromidx, const chromspan& cs) {
    int ninds = d.ninds();

    std::vector<const Individual*> inds;
    for (auto it = d.individuals.begin(); it != d.individuals.end(); ++it) {
        inds.push_back(&(it->second));
    }
    // Step 1: select 3 random individuals
    int aidx = randint(0, ninds-1);
    int bidx = randint(0, ninds-1);
    int tidx = randint(0, ninds-1);

    Individual ind_a = *(inds[aidx]);
    Individual ind_b = *(inds[bidx]);
    const Individual& ind_template = *(inds[tidx]);
    // Step 2: make copies of two of them


    // Step 3: from the third individual copy a genome region to one chromosome of the
    //         two others
    int shared_hapidx = randint(0, 1);
    copy_genospan(ind_template, shared_hapidx, ind_a, randint(0, 1), chromidx, cs);
    copy_genospan(ind_template, shared_hapidx, ind_b, randint(0, 1), chromidx, cs);

    // Step 4: return the pair

    Indpair pr = {ind_a, ind_b};
    return pr;
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
        chromspan rs = { chrominf->positions[s.full_start], 
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