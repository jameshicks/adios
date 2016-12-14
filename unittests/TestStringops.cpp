#include "stringops.hpp"
#include "CppUTest/TestHarness.h"


TEST_GROUP(Stringops) {
};


TEST(Stringops, StringSplitter) {
    using stringops::split;
    std::string data("abc def ghi");

    auto expected = std::vector<std::string>{"abc", "def", "ghi"};
    CHECK(split(data, " ") == expected);

    expected = {"abc", "def ghi"};
    auto observed = split(data, " ", 1);
    CHECK(observed == expected);

    data = std::string("abcdef");
    expected = std::vector<std::string>{"abcdef"};
    CHECK(split(data, " ") == expected);

    data = std::string("");
    // expected = std::vector<std::string>;
    CHECK(split(data, " ").empty());

};

TEST(Stringops, Startswith) {
    using stringops::startswith;
    CHECK(startswith("ADIOS", "AD"));
    CHECK(!startswith("ADIOS", "SO"));
    CHECK(startswith("ADIOS", "ADIOS"));
    CHECK(!startswith("AD", "ADIOS"));
    CHECK(!startswith("", ""));
}


TEST(Stringops, StringJoin) {
    using stringops::join;
    using std::string;
    std::vector<std::string> test1 = {"hello", "yes", "this", "is", "dog"};
    string observed = join(test1, " ");
    string expected = "hello yes this is dog";
    CHECK(observed.compare(expected) == 0); 

}