#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <sstream>
#include <memory>
#include <stdexcept>
#include <chrono>
#include <stdlib.h>
#include <time.h>
#include <sys/resource.h>

#include "config.h"

#ifdef HAVE_OPENMP
#include <omp.h>
#define OMP_AVAILABLE 1 
#else 
#define OMP_AVAILABLE 0
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

    struct timespec prog_start;
    struct timespec prog_stop;

    if (clock_gettime(CLOCK_MONOTONIC, &prog_start)) {
        std::cerr << "Error getting clock time\n";
    }

    string sprog_start = current_time_string();


    std::vector<string> rawargs;
    for (int argidx = 0; argidx < argc; ++argidx) {
        rawargs.push_back(string(argv[argidx]));
    }

    ArgumentParser parser;
    std::vector<CommandLineArgument> arginfo = {
        //                  Argument           Action       Default               narg  help string
        CommandLineArgument{"vcf",               "store",     {""},               1,    "VCF input file"},
        CommandLineArgument{"vcf_freq",          "store",     {"-"},              1,    "VCF INFO field containing allele frequency"},
        CommandLineArgument{"out",               "store",     {"-"},              1,    "Output file prefix"},
        CommandLineArgument{"keep_singletons",   "store_yes", {"NO"},             0,    "Include singleton variants from dataset"},
        CommandLineArgument{"keep_monomorphic",  "store_yes", {"NO"},             0,    "Include monomorphic positions in dataset"},
        CommandLineArgument{"rare",              "store",     {"0.05"},           1,    "Rare frequency threshold"},
        CommandLineArgument{"freq_floor",        "store",     {"0.001"},          1,    "Variants with frequencies below this are set to this"},
        CommandLineArgument{"minlod",            "store",     {"3.0"},            1,    "Minimum LOD to report a segment"},
        CommandLineArgument{"minlength",         "store",     {"1.0"},            1,    "Miniumum length to report a segment (Mb)"},
        CommandLineArgument{"minmark",           "store",     {"16"},             1,    "Minimum number of shared rare variants to report an IBD segment"},
        CommandLineArgument{"err",               "store",     {"0.001", "0.005"}, 2,    "Allele error rates (rare, common)"},
        CommandLineArgument{"transition",        "store",     {"6", "3"},         2,    "IBD entrance/exit penalty: P(Transition) = 10^(-x))"},
        CommandLineArgument{"threads",           "store",     {"1"},              1,    OMP_AVAILABLE ? "Number of threads" : "SUPPRESS"},
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

    std::string logfilename;
    if (!args["out"][0].compare("-")) {
        logfilename = "adios.log";
    } else {
        std::stringstream ss; 
        ss << args["out"][0] << ".log";
        logfilename = ss.str();
    }
     
    

    Logstream log(logfilename);

    int nthreads = atoi(args["threads"][0].c_str());


    bool empirical_freqs = !(args["vcf_freq"][0].compare("-"));

    log << "adios v0.8\n";
    log << "Started at " << sprog_start << "\n\n";


#ifdef HAVE_OPENMP
    omp_set_num_threads(nthreads);
#else
    if (nthreads != 1) {
        log << '\n';
        log << "adios was not compiled with multithreading support\n";
        log << "Defaulting to single threading\n";
        log << '\n';
    }
#endif


    adios::adios_parameters params = adios::params_from_args(args);

    log << "VCF file: " << args["vcf"][0] << '\n';
    log << "Frequencies: " << (empirical_freqs ? std::string("Calculated from dataset") : args["vcf_freq"][0]) << '\n';
    log << "Rare frequency threshold: " << params.rare_thresh << '\n';
    log << "Transition costs: " << params.gamma_ << " (IBD entry), " << params.rho << " (IBD exit)\n";
    log << "Minimum segment LOD: " << params.min_lod << '\n';
    log << "Minimum segment length: " << bp_formatter(params.min_length) << '\n';
    log << "Minimum markers to declare IBD: " << params.min_mark << '\n';
    log << "Genotype error rate: " << params.err_rate_common << " (common variants), " << params.err_rate_rare << " (rare variants)" << '\n';
    log << "Decoding: " << (params.viterbi ? "MAP" : "ML") << '\n';
#ifdef HAVE_OPENMP
    log << "Threads: " << nthreads << '\n';
#endif

    log << '\n';

    VCFParams vcfp = {!(args["keep_singletons"][0].compare("YES") == 0),
                      !(args["keep_monomorphic"][0].compare("YES") == 0),
                      empirical_freqs,
                      args["vcf_freq"][0]
                     };


    auto start = std::chrono::steady_clock::now();
    
    Dataset data;
    try {
        data = read_vcf(args["vcf"][0], vcfp);
    } catch (const std::exception& e) {
        log << "Could not process file: " << args["vcf"][0] << ": ";
        log << e.what() << '\n';
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
    log << "Processed " << nmark_total << " variants ";
    log << "in " << elapsedSeconds << "s ";
    log << "(" << (nmark_total / elapsedSeconds) << " variants/sec)\n\n";

    log << data.ninds() << " individuals\n";
    for (size_t chridx = 0; chridx < data.nchrom(); ++chridx) {
        auto c = data.chromosomes[chridx];
        log << "Chromosome " << c->label << " (" << c->size() / 1000000 << "Mb)";
        log << ": " << c->nmark() << " variants, ";
        log << params.rare_sites[chridx].size() << " rare. ";

        if (!(c->exclusions.empty())) {
            log << "Excluded";
            for (auto kv : (c->exclusions)) {
                log << " " << kv.second << ' ' << kv.first << ',';
            }
            log << '\n';
        }

        log << '\n';
    }

    // Precompute the emission matrices.
    params.calculate_emission_mats(data);

    std::string output_filename = !(args["out"][0].compare("-")) ? 
                                   "-" : (args["out"][0] + ".ibd");  
    DelimitedFileWriter output(output_filename, '\t');
    adios::adios(data, params, output);

    log << "\n\nCompleted at " << current_time_string() << '\n';
    
    if (clock_gettime(CLOCK_MONOTONIC, &prog_stop)) {
        std::cerr << "Error getting clock time\n";
    }

    timeval real = {prog_stop.tv_sec - prog_start.tv_sec, 0};

    struct rusage usage; 
    getrusage(RUSAGE_SELF, &usage);
    log << print_elapsed(real) << "\tReal time\n";
    log << print_elapsed(usage.ru_utime) << "\tCPU time\n";
    log << print_elapsed(usage.ru_stime) << "\tSystem time\n";
}







