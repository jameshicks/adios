#include <iostream>

#include <vector>
#include <stdexcept>
#include <iostream>
#include <cmath>
#include "Linalg.hpp"

#include "CppUTest/TestHarness.h"
TEST_GROUP(Vector)
{
};

TEST(Vector, VectorBasic)
{
    std::vector<double> d(3);
    // CHECK((d == Linalg::Vector{0, 0, 0}));
    Linalg::Vector v(d);
    Linalg::Vector v2(d);
    CHECK(v == v2);
    CHECK(v.sum() == 0);

    v2 = Linalg::Vector{1, 2, 3};

    Linalg::Vector q = {1,2,3};
    CHECK(q.sum() == 6);
    CHECK(!(q == v));
    CHECK(q != v);
    CHECK(q.get(1) == 2);
    CHECK_EQUAL(1, q.get(0));
    CHECK(q.get(0) != 2);
    CHECK(q.size == 3);
    CHECK_THROWS(std::invalid_argument, q.get(100));

    CHECK_EQUAL(0, Linalg::dot_product(q,v));
    DOUBLES_EQUAL(14.0, Linalg::dot_product(q, q), 0.001);
    CHECK((q == Linalg::Vector{1, 2, 3}));
    CHECK((q != Linalg::Vector{3, 2, 1}));

    v.set_all(2);
    CHECK(v == (Linalg::Vector{2,2,2}));
};

TEST(Vector, VectorMinMax) {
    Linalg::Vector v{1, 100, 3};
    CHECK_EQUAL(1, v.argmax());
    CHECK_EQUAL(100, v.max());
    CHECK_EQUAL(0, v.argmin());
    CHECK_EQUAL(1, v.min())

    v = Linalg::Vector{100, 1, 100};
    CHECK_EQUAL(100, v.max());
    CHECK_EQUAL(0, v.argmax());
    CHECK_EQUAL(1, v.argmin());
    CHECK_EQUAL(1, v.min());

    v = Linalg::Vector{1, 100, 1};
    CHECK_EQUAL(1, v.argmax());
    CHECK_EQUAL(100, v.max());
    CHECK_EQUAL(0, v.argmin());
    CHECK_EQUAL(1, v.min());

};

TEST(Vector, VectorApply) {
    Linalg::Vector v{10, 100, 1000};
    Linalg::Vector expected{1, 2, 3};
    Linalg::Vector observed = v.apply(&log10);
    CHECK(observed == expected);

    v = Linalg::Vector{2.5, 2.5, 2.5};
    expected = Linalg::Vector{3, 3, 3};
    observed = v.apply(&ceil);
    CHECK(observed == expected);

};


TEST(Vector, VectorOperator) {
    auto v1 = Linalg::Vector{0, 0, 0};
    auto v2 = Linalg::Vector{3, 3, 3};
    auto v3 = Linalg::Vector{6, 6, 6};
    CHECK((v1 + v2) ==  v2);
    CHECK((v2 - v1) ==  v2);

    CHECK((v1 + 3) == v2);
    CHECK((v2 - 3) == v1);
    CHECK((v2 * 2) == v3)
    CHECK((v3 / 2) == v2);

};

TEST(Vector, VectorMisc) {
    using Linalg::Matrix;
    Linalg::Vector v1 = {1,2,3};
    Matrix m = Linalg::diag(v1);

    Matrix expected = {{1,0,0},{0,2,0}, {0,0,3}};
    CHECK(m==expected);
}

