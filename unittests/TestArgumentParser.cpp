#include <vector>
#include <map>
#include <string>
#include "ArgumentParser.hpp"
#include "CppUTest/TestHarness.h"

TEST_GROUP(ArgumentParser) {};

TEST(ArgumentParser, store_yes) {
    ArgumentParser parser;
    parser.add_argument(CommandLineArgument{"test", "store_yes", {"NO"}, 0, "testyes"  });
    CHECK(parser.args.at("test")[0].compare("NO") == 0);

    std::vector<std::string> argline = {"progname", "--test"}; 
    parser.update_args(argline);
    CHECK(parser.args.at("test")[0].compare("YES") == 0);
}