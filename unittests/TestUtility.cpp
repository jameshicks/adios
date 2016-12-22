#include "utility.hpp"
#include "CppUTest/TestHarness.h"

#include <iostream>

TEST_GROUP(UtilityFunctions) {
};


TEST(UtilityFunctions, ValueRuns) {
    ValueRun v = {1, 2, 3};
    CHECK(v == v);

    std::vector<int> exp_none = {0,0,0,0,0};
    auto observed = runs_gte(exp_none, 1);
    CHECK(observed.empty());

    std::vector<int> exp_one = {0,2,2,2,2,2,2,2,0,0,0};
    observed = runs_gte(exp_one,1);
    CHECK_EQUAL(1, observed.size());
    v = {1, 8, 2};
    CHECK(v == observed[0]);


    std::vector<int> exp_runtoend = {0,0,0,2,2,2,2,2};
    observed = runs_gte(exp_runtoend, 1);
    v ={3,8,2};
}

TEST(UtilityFunctions, arange) {
    std::vector<double> observed = arange(0,1,.25);
    std::vector<double> expected = {0, 0.25, 0.5, 0.75, 1.0};
    CHECK(observed==expected);

}