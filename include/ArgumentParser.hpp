#include <map>
#include <vector>
#include <string>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

#include "stringops.hpp"

struct CommandLineArgument {
    std::string label;
    std::string action;
    std::vector<std::string> default_;
    size_t nargs;
    std::string help;
};

struct ArgumentParser {
    std::map<std::string, CommandLineArgument> arg_info;
    std::map<std::string, std::vector<std::string>> args;
    std::vector<std::string> extras;

    void add_argument(CommandLineArgument arg);
    bool has_arg(const std::string& arg) const;
    void update_args(std::vector<std::string> args);
    std::vector<std::string> validate_args(void) const;
    void print_help(void) const;
};