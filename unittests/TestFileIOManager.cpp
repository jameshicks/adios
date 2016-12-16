#include <iostream>
#include "FileIOManager.hpp"
#include "CppUTest/TestHarness.h"

TEST_GROUP(FileIO) {};

TEST(FileIO, readline) {
    FileReader f("unittests/data/abc.txt");
    std::string l;
    CHECK(f.getline().compare("a b c") == 0);
    CHECK(f.getline().compare("d e") == 0);
    CHECK(f.getline().compare("f") == 0);
    CHECK(f.eof());

    FileReader lf("unittests/data/longline.txt");
    std::string expected;
    for (size_t i = 0; i < 100000; ++i) { expected.push_back('a'); }
    std::string observed = lf.getline();
    CHECK(expected.compare(observed) == 0);
    CHECK(lf.eof());
}