#include <string>
#include <vector>
#include <stdlib.h>
#include <algorithm>

#include "config.h"
#include "datamodel.hpp"
#include "vcf.hpp"
#include "ArgumentParser.hpp"
#include "utility.hpp"


AlleleSites synthetic_chromosome(Dataset& d, int chromidx, int chunksize) {
    std::vector<std::string> indlabs;
    for (auto& kv : d.individuals) { indlabs.push_back(kv.first); } 

    AlleleSites newchrom;
    int cur_pos = 0;
    int chromsize = d.chromosomes[chromidx]->size();  
    while (cur_pos < chromsize) { 
        int indidx = randint(0, d.ninds()-1);
        Individual& template_ind = d.individuals[indlabs[indidx]];

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
        //                  Argument      action     default nargs   help
        CommandLineArgument{"vcf",        "store",   {""},          1,     "Input VCF"},
        CommandLineArgument{"out",        "store",   {""},          1,     "Output VCF"},
        CommandLineArgument{"synthetics", "store",   {""},          1,     "Number of synthetic individuals"},
        CommandLineArgument{"prefix",     "store",   {"SYNTH"},     1,     "Prefix for synthetics"},
        CommandLineArgument{"chunksize",  "store",   {"20000"},     1,     "Null chunk size"},
        CommandLineArgument{"seed",       "store",   {"TIME"},      1,     "RNG seed"}
    };
    for (auto& argi : arginfo) { parser.add_argument(argi); }

    try {
        parser.update_args(rawargs);
        
    } catch (std::out_of_range& e) {
        std::cerr << e.what() << '\n';
    }
    auto args = parser.args;

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

    int nsynth = std::stod(args["synthetics"][0]);
    for (int rep = 0; rep < nsynth; rep++) {
        std::stringstream ss; 
        ss << args["prefix"][0] << "_" << rep; 

        Individual dummy = dummy_ind(data, chunksize);
        dummy.label = ss.str();

        data.individuals[dummy.label] = dummy;
    }

    std::cout << "Writing output VCF\n";
    write_vcf(data, args["out"][0]);

    return 0;
}
