#include "ArgumentParser.hpp"

void ArgumentParser::add_argument(CommandLineArgument arg) {
    arg_info[arg.label] = arg;
    args[arg.label] = arg.default_;
}

void ArgumentParser::update_args(std::vector<std::string> newargs) {
    using std::string;
    std::map<string, string> parsed_args;
    string cur_flag = "";
    size_t cur_narg = 0;
    for (auto it = newargs.begin() + 1; it != newargs.end(); ++it) {
        if (stringops::startswith(*it, "--")) {
            cur_flag = it->substr(2);
            
            if (args.find(cur_flag) == args.end()) { 
                throw std::out_of_range("Unknown argument: " + cur_flag); 
            }
            CommandLineArgument a = arg_info.at(cur_flag);
            
            if (a.action.compare("store") == 0) {
                args[cur_flag] = {};
            } else if (a.action.compare("store_yes") == 0) {
                args[cur_flag] = {"YES"};
            }

        } else {
            if (cur_narg < arg_info[cur_flag].nargs) {
                args[cur_flag].push_back(*it);
            } else {
                extras.push_back(*it);
            }
        }

    }
}

std::vector<std::string> ArgumentParser::validate_args(void) const {
    std::vector<std::string> errors;
    for (auto it = args.begin(); it != args.end(); ++it) {
        CommandLineArgument ai = arg_info.at(it->first);
        if (it->second.front().empty()) {
            errors.push_back("Argument required: " + it->first);
        }
        if (ai.nargs && it->second.size() != ai.nargs) {
            errors.push_back("Not enough arguments: " + it->first + " (" + std::to_string(ai.nargs) + " required)");
        }
    }
    return errors;
}

bool ArgumentParser::has_arg(const std::string& arg) const {
    std::vector<std::string> val = args.at(arg);
    return !val.empty() && !(val[0].compare("NO") == 0);
}

void ArgumentParser::print_help(void) const {
    using std::cout;
    using std::endl;
    using std::left;

    std::cout.width(25); cout << std::left << "Flag";
    std::cout.width(30); cout << "Description" << endl;
    cout << endl;
    for (auto kv : arg_info) {
        CommandLineArgument cla = kv.second;
        if (!cla.help.compare("SUPPRESS")) continue;
        std::cout.width(25); std::cout << std::left << "--" + cla.label;
        std::cout.width(30); std::cout << cla.help << std::endl;
    }
    cout << endl;
}