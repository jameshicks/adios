#include <vector>
#include "combinatorics.hpp"

#include "CppUTest/TestHarness.h"

TEST_GROUP(Combinatorics)
{
};

TEST(Combinatorics, PairCombinations) {
    using std::make_pair;
    std::vector<int> testvals = {1,2,3};
    auto pairs = combinatorics::pair_combinations(testvals);

    CHECK_EQUAL(3, testvals.size());
    
    auto expected = make_pair(1,2);
    CHECK(expected == pairs[0]);
    
    expected = make_pair(1,3);
    CHECK(expected == pairs[1]);
    
    expected = make_pair(2,3);
    CHECK(expected == pairs[2]);
}
