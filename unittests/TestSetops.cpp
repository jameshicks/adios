#include <iostream>
#include <set>
#include <vector>
#include "setops.hpp"

#include "CppUTest/TestHarness.h"
using std::set;
using std::vector;
TEST_GROUP(Setops)
{
};

TEST(Setops, BasicUnion) {
    set<int> a = {1,2,3,4,5};
    set<int> b = {3,4,5,6,7};
    set<int> emptyset;

    set<int> expected_union = {1,2,3,4,5,6,7};

    auto observed_union = setops::union_(a,b);
    CHECK(expected_union == observed_union);

    observed_union = setops::union_(a,emptyset);
    CHECK(observed_union == a);

    observed_union = setops::union_(emptyset, emptyset);
    CHECK(observed_union == emptyset);
};

TEST(Setops, BasicIntersection) {
    set<int> a = {1,2,3,4,5};
    set<int> b = {3,4,5,6,7};
    set<int> emptyset;

    set<int> expected_intersection = {3,4,5};
    auto observed_intersection = setops::intersection(a,b);
    CHECK(expected_intersection == observed_intersection);

    CHECK(setops::intersection(a, emptyset) == emptyset);
    CHECK(setops::intersection(emptyset, emptyset) == emptyset);

};

TEST(Setops, BasicSymmetricDiff) {

    set<int> a = {1,2,3,4,5};
    set<int> b = {3,4,5,6,7};
    set<int> emptyset;

    set<int> expected_symdiff = {1,2,6,7};
    auto observed_symdiff = setops::symmetric_difference(a,b);
    CHECK(expected_symdiff == observed_symdiff);

    CHECK(setops::symmetric_difference(a, emptyset) == a);
    CHECK(setops::symmetric_difference(b, emptyset) == b);
    CHECK(setops::symmetric_difference(emptyset, emptyset) == emptyset);
    CHECK(setops::symmetric_difference(a,a) == emptyset);
    CHECK(setops::symmetric_difference(b,b) == emptyset);
};

TEST(Setops, BasicDiffSet) {

    set<int> a = {1,2,3,4,5};
    set<int> b = {3,4,5,6,7};
    set<int> emptyset;

    set<int> expected_diff = {1,2};
    CHECK(setops::difference(a,b) == expected_diff);

    expected_diff = {6,7};
    CHECK(setops::difference(b,a) == expected_diff);

    CHECK(setops::difference(a,emptyset) == a);
    CHECK(setops::difference(b,emptyset) == b);
    CHECK(setops::difference(emptyset, a) == emptyset);
    CHECK(setops::difference(emptyset, b) == emptyset);
    CHECK(setops::difference(emptyset, emptyset) == emptyset);
    

};

TEST(Setops, BasicDiffVec) {

    vector<int> a = {1,2,3,4,5};
    vector<int> b = {3,4,5,6,7};
    vector<int> emptyset;

    vector<int> expected_diff = {1,2};
    CHECK(setops::difference(a,b) == expected_diff);

    expected_diff = {6,7};
    CHECK(setops::difference(b,a) == expected_diff);

    CHECK(setops::difference(a,emptyset) == a);
    CHECK(setops::difference(b,emptyset) == b);
    CHECK(setops::difference(emptyset, a) == emptyset);
    CHECK(setops::difference(emptyset, b) == emptyset);
    CHECK(setops::difference(emptyset, emptyset) == emptyset);
    

};


TEST(Setops, MultiUnion) {

    set<int> a = {1,2,3,4,5};
    set<int> b = {3,4,5,6,7};
    set<int> c = {6,7,8,9,10};

    set<int> expected_union = {1,2,3,4,5,6,7,8,9,10};
    CHECK(setops::multi_union({a,b,c}) == expected_union);
    

};


TEST(Setops, MultiIntersect) {

    set<int> a = {1,2,3,4};
    set<int> b = {3,4,5};
    set<int> c = {3,4,5,6};

    set<int> expected_intersection = {3,4};
    CHECK(setops::multi_intersection({a,b,c}) == expected_intersection);
    

};


