#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <stdlib.h>

#include "config.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#endif

#include "vcf.hpp"
#include "adios.hpp"
#include "ArgumentParser.hpp"
#include "datamodel.hpp"
#include "FileIOManager.hpp"



int main(int argc, char** argv) {
    using std::string;
    using std::cout;
    using std::endl;


    std::vector<string> rawargs;
    for (int argidx = 0; argidx < argc; ++argidx) {
        rawargs.push_back(string(argv[argidx]));
    }

    ArgumentParser parser;
    std::vector<CommandLineArgument> arginfo = {
        //                  Argument           Action       Default               narg  help string
        CommandLineArgument{"vcf",               "store",     {""},               1,    "VCF input file"},
        CommandLineArgument{"vcf_freq",          "store",     {"AF"},             1,    "VCF INFO field containing allele frequency"},
        CommandLineArgument{"out",               "store",     {"-"},              1,    "Write output to file"},
        CommandLineArgument{"empirical_freqs",   "store_yes", {"NO"},             0,    "Calculate allele frequencies from data"},
        CommandLineArgument{"keep_singletons",   "store_yes", {"NO"},             0,    "Include singleton variants from dataset"},
        CommandLineArgument{"keep_monomorphic",  "store_yes", {"NO"},             0,    "Include monomorphic positions in dataset"},
        CommandLineArgument{"rare",              "store",     {"0.05"},           1,    "Rare frequency threshold"},
        CommandLineArgument{"freq_floor",        "store",     {"0.001"},          1,    "Variants with frequencies below this are set to this"},
        CommandLineArgument{"minlod",            "store",     {"3.0"},            1,    "Minimum LOD to report a segment"},
        CommandLineArgument{"minlength",         "store",     {"1.0"},            1,    "Miniumum length to report a segment (Mb)"},
        CommandLineArgument{"minmark",           "store",     {"16"},             1,    "Minimum number of shared rare variants to report an IBD segment"},
        CommandLineArgument{"err",               "store",     {"0.001", "0.005"}, 2,    "Allele error rates (rare, common)"},
        CommandLineArgument{"transition",        "store",     {"6", "3"},         2,    "IBD entrance/exit penalty: P(Transition) = 10^(-x))"},
        CommandLineArgument{"threads",           "store",     {"1"},              1,    "Number of threads"},
        CommandLineArgument{"help",              "store_yes", {"NO"},             0,    "Display this help message"   },
        CommandLineArgument{"version",           "store_yes", {"NO"},             0,    "Print version information"   },
        CommandLineArgument{"viterbi",           "store_yes", {"NO"},             0,    "Use maximum a posteriori decoding"}
    };
    for (auto argi : arginfo) { parser.add_argument(argi); }

    parser.update_args(rawargs);


    if (parser.has_arg("version")) {
        std::cout << "adios v0.8" << std::endl << std::endl;
        return 0;
    }

    if (parser.has_arg("help")) {
        parser.print_help();
        return 0;
    }

    std::vector<std::string> errors = parser.validate_args();

    if (!errors.empty()) {
        for (auto e : errors) {
            std::cerr << e << '\n';
        }
        return 64;
    }
    auto args = parser.args;

    // int nthreads = 1;
    int nthreads = atoi(args["threads"][0].c_str());


    bool empirical_freqs = (args["empirical_freqs"][0].compare("YES") == 0);

    cout << "adios v0.8" << std::endl << std::endl;


#ifdef HAVE_OPENMP
    omp_set_num_threads(nthreads);
#else
    if (nthreads != 1) {
        cout << '\n';
        cout << "adios was not compiled with multithreading support\n";
        cout << "Defaulting to single threading\n";
        cout << '\n';
    }
#endif


    adios::adios_parameters params = adios::params_from_args(args);

    cout << "VCF file: " << args["vcf"][0] << '\n';
    cout << "Frequencies: " << (empirical_freqs ? std::string("Calculated from dataset") : args["vcf_freq"][0]) << '\n';
    cout << "Rare frequency threshold: " << params.rare_thresh << '\n';
    cout << "Transition costs: " << params.gamma_ << " (IBD entry), " << params.rho << " (IBD exit)\n";
    cout << "Minimum segment LOD: " << params.min_lod << '\n';
    cout << "Minimum segment length: " << bp_formatter(params.min_length) << '\n';
    cout << "Minimum markers to declare IBD: " << params.min_mark << '\n';
    cout << "Genotype error rate: " << params.err_rate_common << " (common variants), " << params.err_rate_rare << " (rare variants)" << '\n';
    cout << "Decoding: " << (params.viterbi ? "MAP" : "ML") << '\n';
#ifdef HAVE_OPENMP
    cout << "Threads: " << nthreads << '\n';
#endif

    cout << endl;

    VCFParams vcfp = {!(args["keep_singletons"][0].compare("YES") == 0),
                      !(args["keep_monomorphic"][0].compare("YES") == 0),
                      empirical_freqs,
                      args["vcf_freq"][0]
                     };


    auto start = std::chrono::steady_clock::now();
    Dataset data;
    try {
        data = read_vcf(args["vcf"][0], vcfp);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Could not open file: " << args["vcf"][0] << ": ";
        std::cerr << e.what() << std::endl;
        return 1;
    }

    params.get_rare_sites(data);

    // If we calculated the data we can round them to a sensible place too.
    if (empirical_freqs) {
        data.round_frequencies(4);
    }

    data.floor_frequencies(std::stod(args["freq_floor"][0]));
    auto end = std::chrono::steady_clock::now();
    double elapsedSeconds = ((end - start).count()) * std::chrono::steady_clock::period::num / static_cast<double>(std::chrono::steady_clock::period::den);

    size_t nmark_total = data.nmark() + data.nexcluded();
    std::cout << "Processed " << nmark_total << " variants ";
    std::cout << "in " << elapsedSeconds << "s ";
    std::cout << "(" << (nmark_total / elapsedSeconds) << " variants/sec)\n\n";

    std::cout << data.ninds() << " individuals" << std::endl;
    for (size_t chridx = 0; chridx < data.nchrom(); ++chridx) {
        auto c = data.chromosomes[chridx];
        std::cout << "Chromosome " << c->label << " (" << c->size() / 1000000 << "Mb)";
        std::cout << ": " << c->nmark() << " variants, ";
        std::cout << params.rare_sites[chridx].size() << " rare. ";

        if (!(c->exclusions.empty())) {
            std::cout << "Excluded";
            for (auto kv : (c->exclusions)) {
                std::cout << " " << kv.second << ' ' << kv.first << ',';
            }
            std::cout << '\n';
        }

        std::cout << std::endl;
    }

    // Precompute the emission matrices.
    params.calculate_emission_mats(data);

    DelimitedFileWriter output(args["out"][0], '\t');
    adios::adios(data, params, output);
}







