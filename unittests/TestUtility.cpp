#include "common.hpp"
#include "CppUTest/TestHarness.h"


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

}