#ifndef POWER_HPP
#define POWER_HPP

#include <stdlib.h>
#include <initializer_list>


#include "datamodel.hpp"
#include "adios.hpp"
namespace adios {

struct chromspan {
    int start;
    int stop;
    inline int length(void) const { return stop - start; }
};


struct PowerReplicateResult {
    bool success;
    int nseg;
    unsigned long overlap_amount; 
    std::vector<Segment> segments;

};

struct PowerResult {
    int segsize;
    std::vector<PowerReplicateResult> replicates; 
    
    inline int nsuccess(void) const {
        int i = 0;
        for (auto& rep : replicates) i += rep.success;
        return i; 
    }

    inline double power(void) const { 
        return (double)nsuccess() / replicates.size(); 
    }

    double mean_num_segments(void) const;
    std::vector<int> length_diffs(void) const;
    double prop_detected(void) const; 

};

void copy_genospan(const Individual& from, int hapfrom, Individual& to, int hapto, 
    int chromidx, const chromspan& cs);


Indpair dummy_indpair(const Dataset& d, int chromidx, const chromspan& cs);
chromspan random_range(std::shared_ptr<ChromInfo>& c,  int size);
PowerReplicateResult power_replicate(const Dataset& d,
                                     const adios_parameters& params,
                                     int chromidx,
                                     unsigned int segsize);
PowerResult calc_power(Dataset& d, const adios_parameters& params,
                int chromidx,
                unsigned int segsize,
                unsigned int nrep);


}
#endif