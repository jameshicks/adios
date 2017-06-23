#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>
#include <set>
#include <fstream>


#include "config.h"
#include "adios.hpp"
#include "datamodel.hpp"
#include "vcf.hpp"
#include "ArgumentParser.hpp"
#include "utility.hpp"
#include "setops.hpp"

struct synthseg {
    std::string ind_a;
    std::string ind_b;
    std::string chromlab;
    int start;
    int stop;

    bool overlaps(const synthseg& other) {
        if (chromlab.compare(other.chromlab) != 0) return false;
        return  ((start <= other.start && other.start <= stop) ||
                 (other.start <= start && start <= other.stop));
    }

    std::string to_string(void) {
        std::stringstream ss;
        ss << ind_a << '\t' << ind_b << '\t' << chromlab;
        ss << '\t' << start << '\t' << stop;
        return ss.str();
    }
};

AlleleSites errored_chromosome(AlleleSites& c, int nmark, double error_rate) {
    std::vector<int> errsites;
    for (int i = 0; i < nmark; i++) {
        if (drand48() < error_rate) errsites.push_back(i);
    }


    // Use symmetric difference: when a site is in both the chromosome and the
    // error sites, it is changing to the major allele and must be removed from
    // the chromosome. Otherwise it must be added to the chromosome. Symmetric
    // difference will do both easily.
    AlleleSites nc = setops::symmetric_difference(c, errsites);
    return nc;
}

void add_error(Individual& ind, double error_rate) {
    for (int cidx = 0; cidx < ind.chromosomes.size(); cidx++) {
        auto& gt = ind.chromosomes[cidx];
        int nmark = gt.info->nmark();
        gt.hapa = errored_chromosome(gt.hapa, nmark, error_rate);
        gt.hapb = errored_chromosome(gt.hapb, nmark, error_rate);

    }
}

AlleleSites synthetic_chromosome(Dataset& d, int chromidx, int chunksize, std::string fn) {

    AlleleSites newchrom;
    int cur_pos = 0;
    int chromsize = d.chromosomes[chromidx]->size();

    std::ofstream chunkfile; chunkfile.open(fn);

    while (cur_pos < chromsize) {
        int indidx = randint(0, d.ninds() - 1);
        Individual& template_ind = d.individuals[indidx];

        Genotypes& template_pair = template_ind.chromosomes[chromidx];
        int whichhap = randint(0,1);
        AlleleSites& template_chrom = whichhap ? template_pair.hapa : template_pair.hapb;

        chunkfile << cur_pos << '\t' << cur_pos + chunksize << '\t' << template_ind.label << '\t' << whichhap << '\n'; 
        auto copystartit = std::lower_bound(template_chrom.begin(), template_chrom.end(), cur_pos);
        auto copystopit = std::upper_bound(template_chrom.begin(), template_chrom.end(), cur_pos + chunksize);

        for (auto it = copystartit; it != copystopit; it++) newchrom.push_back(*it);

        cur_pos += chunksize;
    }

    chunkfile.close();
    return newchrom;
}

Individual dummy_ind(Dataset& d,  int chunksize, std::string label) {

    Individual dummy;
    dummy.get_empty_chromosomes(d);
    dummy.label = label;

    int ninds = d.ninds();

    for (int chromidx = 0; chromidx < d.nchrom(); chromidx++) {
        std::stringstream fn1;
        fn1 << label << "_chr" << d.chromosomes[chromidx]->label << "_0"  << ".chunks"; 
        dummy.chromosomes[chromidx].hapa = synthetic_chromosome(d, chromidx, chunksize, fn1.str());

        std::stringstream fn2;
        fn2 << label << "_chr" << d.chromosomes[chromidx]->label << "_1"  << ".chunks";
        dummy.chromosomes[chromidx].hapb = synthetic_chromosome(d, chromidx, chunksize, fn2.str());
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
        CommandLineArgument{"error",      "store",     {"0"},         1,     "Genotype error rate"},
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


    std::cout << "Input file: " << args["vcf"][0] << '\n';
    std::cout << "Output prefix: " << args["out"][0] << '\n';
    std::cout << "Synthetic individuals: " << args["synthetics"][0] << '\n';
    std::cout << "Synthetic segments: " << args["nseg"][0] << '\n';

    if (std::stol(args["nseg"][0]) > 0) {
        std::cout << "Synthetic segment length: " << args["seglen"][0] << '\n';
    }

    std::cout << "Random seed: " << args["seed"][0] << "\n\n";

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
        std::cout  << "Chromosome " <<  c->label << ": " << c->size() << '\n';
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

        Individual dummy = dummy_ind(data, chunksize, ss.str());

        dummies.push_back(dummy);
    }

    int nseg = std::stol(args["nseg"][0]);
    double segsize = std::stod(args["seglen"][0]);

    std::vector<synthseg> synthsegs;


    for (int i = 0; i < nseg; i++) {
        int aidx, bidx, chridx;
        aidx = randint(0, nsynth - 1);

        do {
            bidx = randint(0, nsynth - 1);
        } while (aidx == bidx);

        Individual& inda = dummies[aidx];
        Individual& indb = dummies[bidx];

        chridx = 0;

        if (data.chromosomes[chridx]->size() < segsize) throw std::out_of_range("Chromosome too small");

        bool any_overlaps = false;
        synthseg s;
        do {
            int start = randint(1, data.chromosomes[chridx]->size() - segsize);
            s = {inda.label, indb.label, data.chromosomes[chridx]->label, start, (int)(start + segsize)};

            // Check if this segment overlaps anything we've made previously
            bool any_overlaps = false;

            for (auto& existing : synthsegs) {
                // To overlap, one of the individuals must be in the existing segment
                bool inds_in_common ((s.ind_a.compare(existing.ind_a) == 0) ||
                                     (s.ind_a.compare(existing.ind_b) == 0) ||
                                     (s.ind_b.compare(existing.ind_a) == 0) ||
                                     (s.ind_b.compare(existing.ind_b) == 0));
                any_overlaps = any_overlaps || (inds_in_common && s.overlaps(existing));

            }


        } while (any_overlaps);

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

    double error_rate = std::stod(args["error"][0]);
    if (error_rate > 0) {
        std::cout << "Adding error\n";
        for (auto& ind : dummies) add_error(ind, error_rate);
    }
    // Add the synthetic individuals to the dataset
    for (auto& dummy : dummies) {
        data.individuals.push_back(dummy);
    }

    std::cout << "Writing output VCF\n";

    std::stringstream ofnss;
    ofnss << args["out"][0] << ".vcf";
    write_vcf(data, ofnss.str());

    return 0;
}
