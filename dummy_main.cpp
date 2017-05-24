#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>

#include "config.h"
#include "adios.hpp"
#include "datamodel.hpp"
#include "vcf.hpp"
#include "ArgumentParser.hpp"
#include "utility.hpp"

struct synthseg {
    std::string ind_a;
    std::string ind_b;
    std::string chromlab;
    int start;
    int stop;

    bool overlaps(const synthseg& other) {
        if (chromlab.compare(other.chromlab) != 0) return false;
        return  ((start <= other.start && other.start <= stop) ||
                 (other.start <= start && start <= stop));
    }

    std::string to_string(void) {
        std::stringstream ss;
        ss << ind_a << '\t' << ind_b << '\t' << chromlab;
        ss << '\t' << start << '\t' << stop;
        return ss.str();
    }
};

AlleleSites synthetic_chromosome(Dataset& d, int chromidx, int chunksize) {

    AlleleSites newchrom;
    int cur_pos = 0;
    int chromsize = d.chromosomes[chromidx]->size();  
    while (cur_pos < chromsize) { 
        int indidx = randint(0, d.ninds()-1);
        Individual& template_ind = d.individuals[indidx];

        Genotypes& template_pair = template_ind.chromosomes[chromidx];
        AlleleSites& template_chrom = randint(0,1) ? template_pair.hapa : template_pair.hapb;
        auto copystartit = std::upper_bound(template_chrom.begin(), template_chrom.end(), cur_pos);
        auto copystopit = std::upper_bound(template_chrom.begin(), template_chrom.end(), cur_pos+chunksize);
        
        for (auto it = copystartit; it != copystopit; it++) newchrom.push_back(*it);

        cur_pos += chunksize;
    } 

    return newchrom; 
}

Individual dummy_ind(Dataset& d,  int chunksize) {

    Individual dummy;
    dummy.get_empty_chromosomes(d);

    int ninds = d.ninds();

    for (int chromidx = 0; chromidx < d.nchrom(); chromidx++) {
        dummy.chromosomes[chromidx].hapa = synthetic_chromosome(d, chromidx, chunksize);
        dummy.chromosomes[chromidx].hapb = synthetic_chromosome(d, chromidx, chunksize);
    }

    return dummy; 
}



int main(int argc, char** argv) {
    using std::string;
    using std::vector;

    ArgumentParser parser;
    vector<string> rawargs;
    for (int argidx = 0; argidx < argc; argidx++) {
        rawargs.push_back(string(argv[argidx]));
    }

    vector<CommandLineArgument> arginfo = {
        //                  Argument      action       default        nargs   help
        CommandLineArgument{"vcf",        "store",     {""},          1,     "Input VCF"},
        CommandLineArgument{"out",        "store",     {""},          1,     "Output VCF"},
        CommandLineArgument{"synthetics", "store",     {""},          1,     "Number of synthetic individuals"},
        CommandLineArgument{"prefix",     "store",     {"SYNTH"},     1,     "Prefix for synthetics"},
        CommandLineArgument{"chunksize",  "store",     {"20000"},     1,     "Null chunk size"},
        CommandLineArgument{"seed",       "store",     {"TIME"},      1,     "RNG seed"},
        CommandLineArgument{"nseg",       "store",     {"0"},         1,     "Number of synthetic segments"},
        CommandLineArgument{"seglen",     "store",     {"0"},         1,     "Length of synthetic segment"},
        CommandLineArgument{"help",       "store_yes", {"NO"},        0,     "Print this help message" }
    };

    for (auto& argi : arginfo) { parser.add_argument(argi); }

    try {
        parser.update_args(rawargs);
        
    } catch (std::out_of_range& e) {
        std::cerr << e.what() << '\n';
    }
    auto args = parser.args;

    if (parser.has_arg("help")) {
        parser.print_help(); 
        return 0;
    }

    int rseed;
    if (args["seed"][0] != "TIME") {
        rseed = std::stoi(args["seed"][0]);
        srand48(rseed);
        std::cout  << "Random Seed: " << rseed << '\n';
    } 

    std::cout << "Reading data\n";
    VCFParams vcfp = {false, false, true, "AF"};

    Dataset data;
    try {
        data = read_vcf(args["vcf"][0], vcfp);
    } catch (const std::exception& e) {
        std::cout << "Could not process file: " << args["vcf"][0] << ": ";
        std::cout << e.what() << '\n';
        return 1;
    }

    std::cout << data.ninds() << " individuals\n";

    for (auto& c : data.chromosomes) {
        for (auto& kv : c->exclusions) {
            std::cout << "Excluded: " << kv.second << " " << kv.first << "\n";
        }
    }

    int chunksize = std::stod(args["chunksize"][0]);

    int nsynth = std::stol(args["synthetics"][0]);
    
    std::vector<Individual> dummies;
    for (int rep = 0; rep < nsynth; rep++) {
        std::stringstream ss; 
        ss << args["prefix"][0] << "_" << rep; 

        Individual dummy = dummy_ind(data, chunksize);
        dummy.label = ss.str();

        dummies.push_back(dummy);
    }

    int nseg = std::stol(args["nseg"][0]);
    double segsize = std::stod(args["seglen"][0]);

    std::vector<synthseg> synthsegs;


    for (int i = 0; i < nseg; i++) {
        int aidx, bidx, chridx; 
        aidx = randint(0, nsynth-1);

        do {
            bidx = randint(0, nsynth-1);
        } while (aidx != bidx);

        Individual& inda = dummies[aidx];
        Individual& indb = dummies[bidx];

        chridx = 0;

        if (data.chromosomes[chridx]->size() < segsize) throw std::out_of_range("Chromosome too small");

        bool nooverlaps = true;
        synthseg s;
        do {
            int start = randint(1, data.chromosomes[chridx]->size());
            s = {inda.label, indb.label, data.chromosomes[chridx]->label, start, (int)(start + segsize)};

            for (auto& existing : synthsegs) nooverlaps &= !(s.overlaps(existing));

        } while (!nooverlaps);

        chromspan cs = {s.start, s.stop};
        copy_genospan(inda, 0, indb, 0, chridx, cs);
        synthsegs.push_back(s);
    }

    if (nseg > 0) {
        std::stringstream ss;
        ss << args["out"][0] << ".trueseg";
        std::string ofn = ss.str();
        DelimitedFileWriter ts(ofn, '\t');
        for (auto& s : synthsegs) {
            ts.writeline(s.to_string());
        }
    } 

    // Add the synthetic individuals to the dataset
    for (auto& dummy : dummies) {
        data.individuals.push_back(dummy);
    }

    std::cout << "Writing output VCF\n";
    write_vcf(data, args["out"][0]);

    return 0;
}
