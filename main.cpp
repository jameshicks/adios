#include <map>
#include <vector>
#include <iostream>
#include <string>
#include <memory>
#include <stdexcept>
#include "vcf.hpp"
#include "adios.hpp"
#include "ArgumentParser.hpp"
#include "datamodel.hpp"



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
        //                  Argument        Action       Default     narg  help string
        CommandLineArgument{"vcf",          "store",     {""},       1,    "VCF input file"},
        CommandLineArgument{"rare",         "store",     {"0.05"},   1,    "Rare frequency threshold"},
        CommandLineArgument{"minlod",       "store",     {"3.0"},    1,    "Minimum IBDLOD"},
        CommandLineArgument{"minlength",    "store",     {"1.0"},    1,    "Miniumum segment length (Mb)"},
        CommandLineArgument{"minmark",      "store",     {"16"},     1,    "Minimum number of markers to establish IBD"},
        CommandLineArgument{"err",          "store",     {"0.001"},  1,    "Allele error rate"},
        CommandLineArgument{"transition",   "store",     {"5", "3"}, 2,    "Transition IBD penalty (10^(-x))"},
        CommandLineArgument{"vcf_freq",     "store",     {"AF"},     1,    "VCF INFO field containing allele frequency"},
        CommandLineArgument{"help",         "store_yes", {"NO"},     0,    "Display this help message"   },
        CommandLineArgument{"version",      "store_yes", {"NO"},     0,    "Print version information"   }
    };
    for (auto argi : arginfo) { parser.add_argument(argi); }

    parser.update_args(rawargs);

    std::vector<std::string> errors = parser.validate_args();

    if (!errors.empty()) {
        for (auto e : errors) {
            std::cerr << e << '\n';
        }
        return 64;
    }

    if (parser.has_arg("version")) {
        std::cout << "adios v0.8" << std::endl << std::endl;
        return 0;
    }

    if (parser.has_arg("help")) {
        parser.print_help();
        return 0;
    }

    auto args = parser.args;



    std::cout << "adios v0.8" << std::endl << std::endl;

    adios::adios_parameters params = adios::params_from_args(args);

    cout << "VCF file: " << args["vcf"][0] << '\n';
    cout << "Frequencies: " << args["vcf_freq"][0] << '\n';
    cout << "Rare frequency threshold: " << params.rare_thresh << '\n';
    cout << "Transition costs: " << params.gamma_ << " (IBD entry), " << params.rho << " (IBD exit)\n";
    cout << "Minimum segment LOD: " << params.min_lod << '\n';
    cout << "Minimum segment length: " << bp_formatter(params.min_length) << '\n';
    cout << "Minimum markers to declare IBD: " << params.min_mark << '\n';
    cout << "Genotype error rate: " << params.err_rate << '\n';
    cout << endl;

    Dataset data;
    try {
        data = read_vcf(args["vcf"][0], args["vcf_freq"][0]);
    } catch (const std::invalid_argument& e) {
        std::cerr << "Could not open file: " << args["vcf"][0] << std::endl;
        return 1;
    }

    params.get_rare_sites(data);

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


    adios::adios(data, params);
}







