#include <iostream>

#include <vector>
#include <stdexcept>
#include <cmath>
#include "Linalg.hpp"

#include "CppUTest/TestHarness.h"
TEST_GROUP(Matrix)
{
};

TEST(Matrix, MatrixBasic)
{
    Linalg::Matrix m1{{1, 2, 3}, {3, 2, 1}};

    CHECK((m1 == Linalg::Matrix{{1, 2, 3}, {3, 2, 1}}));

    // Cell getters and setters
    CHECK_EQUAL(1, m1.get(0, 0));
    CHECK_EQUAL(2, m1.get(0, 1));
    CHECK_EQUAL(3, m1.get(0, 2));
    CHECK_EQUAL(3, m1.get(1, 0));
    CHECK_EQUAL(2, m1.get(1, 1));
    CHECK_EQUAL(1, m1.get(1, 2));

    m1.set(0, 0, 100.0);
    CHECK_EQUAL(100.0, m1.get(0, 0));
    m1.set(0, 0, 0.0);
    CHECK_EQUAL(0, m1.get(0, 0));
    m1.set(0, 0, 1.0);
    CHECK_EQUAL(1, m1.get(0, 0));

    // Vector getters and setters
    auto o = m1.get_row(0);
    auto r = Linalg::Vector{1, 2, 3};

    CHECK_EQUAL(o.size, r.size);
    CHECK((m1.get_row(0) == r));
    CHECK((m1.get_row(1) == Linalg::Vector{3, 2, 1}));

    o = m1.get_column(0);
    r = Linalg::Vector{1, 3};
    CHECK_EQUAL(o.size, r.size);
    CHECK(m1.get_column(0) == (Linalg::Vector{1, 3}));
    CHECK(m1.get_column(1) == (Linalg::Vector{2, 2}));
    CHECK(m1.get_column(2) == (Linalg::Vector{3, 1}));

    m1.set_row(0, Linalg::Vector{10, 10, 10});
    CHECK((m1.get_row(0) == Linalg::Vector{10, 10, 10}));
    m1.set_column(1, Linalg::Vector{20, 20});
    CHECK((m1.get_column(1) == Linalg::Vector{20, 20}))

    CHECK_THROWS(std::invalid_argument,
                 m1.set_column(1, Linalg::Vector{1, 2, 3, 4, 5}));
    CHECK_THROWS(std::invalid_argument,
                 m1.set_row(1, Linalg::Vector{1, 2, 3, 4, 5}));

    Linalg::Matrix m2{{1, 2}, {3, 4}};
    CHECK(!m1.is_square());
    CHECK(m2.is_square());

    auto m3 = Linalg::Matrix(4, 4);
    m3 = 0;
    Linalg::Matrix m4 = {{1, 2}, {3, 4}};
    m3.set_submatrix(1, 1, m4);

    Linalg::Matrix expected = {
        {0, 0, 0, 0},
        {0, 1, 2, 0},
        {0, 3, 4, 0},
        {0, 0, 0, 0}
    };
    CHECK(m3 == expected);


};

TEST(Matrix, MatrixMinMax)
{
    auto m1 = Linalg::Matrix{{1, 10, 4}, {50, 2, 20}};
    Linalg::Vector expected{1, 0};
    CHECK(expected == m1.argmax_row());

    expected = Linalg::Vector{10, 50};
    CHECK(expected == m1.max_row());

    expected = Linalg::Vector{0, 1};
    CHECK(expected == m1.argmin_row());

    expected = Linalg::Vector{1, 2};
    CHECK(expected == m1.min_row());

    expected = Linalg::Vector{1, 0, 1};
    CHECK(expected == m1.argmax_col());

    expected = Linalg::Vector{50, 10, 20};
    CHECK(expected == m1.max_col());

    expected = Linalg::Vector{0, 1, 0};
    CHECK(expected == m1.argmin_col());

    expected = Linalg::Vector{1, 2, 4};
    CHECK(expected == m1.min_col());

}

// TEST(Matrix, MatrixOperations) {
//     auto m = Linalg::Matrix{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

//     auto expected = Linalg::Matrix{
//         {1, 4, 7},
//         {2, 5, 8},
//         {3, 6, 9}};
//     CHECK(expected == m.transpose());
// }

TEST(Matrix, MatrixProduct)
{
    auto m1 = Linalg::Matrix{{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};
    auto m2 = Linalg::Matrix{{3, 2, 1}, {6, 5, 4}, {9, 8, 7}};

    auto expected = Linalg::Matrix{
        {42, 36, 30},
        {96, 81, 66},
        {150, 126, 102}};
    CHECK(Linalg::matrix_product(m1, m2) == expected);

    expected = Linalg::Matrix{{3, 4, 3}, {24, 25, 24}, {63, 64, 63}};
    CHECK(Linalg::direct_product(m1, m2) == expected);

    auto m3 = Linalg::Matrix{{1, 2}, {3, 4}};
    CHECK_THROWS(std::invalid_argument, Linalg::matrix_product(m1, m3));
    CHECK_THROWS(std::invalid_argument, Linalg::direct_product(m1, m3));

    auto a = Linalg::Matrix{{1, 2}, {3, 4}};
    auto b = Linalg::Matrix{{0, 5}, {6, 7}};

    expected = Linalg::Matrix{
        {0, 5, 0, 10},
        {6, 7, 12, 14},
        {0, 15, 0, 20},
        {18, 21, 24, 28}};

    auto observed = kronecker_product(a, b);
    CHECK(expected == observed);
};

TEST(Matrix, MatrixVectorProducts) {
    using namespace Linalg;
    Matrix A = {{1,2,3}, {4,5,6}, {7,8,9}};
    Vector y = {2,1,3};

    Vector expected = {13,31,49};
    Vector observed = matrix_vector_product(A,y);
    CHECK(observed == expected);
    CHECK_EQUAL(observed.size, A.ncol);

    Matrix B = {{ 4,  5,  6}, { 7,  8,  9}, {10, 11, 12}};
    Vector v = {1,2,3};

    expected = {48,54,60};
    CHECK(vector_matrix_product(v,B) ==  expected);
}

TEST(Matrix, ViewObjects) {
    using namespace Linalg;
    Matrix A = {{1,2,3}, {4,5,6}, {7,8,9}};
    auto a_rv0 = A.row_view(0);
    CHECK(a_rv0 == (Vector{1,2,3}));

    auto a_rv1 = A.row_view(1);
    CHECK(a_rv1 == (Vector{4,5,6}));

    auto a_rv2 = A.row_view(2);
    CHECK(a_rv2 == (Vector{7,8,9}));

    auto a_cv0 = A.col_view(0);
    CHECK(a_cv0 == (Vector{1,4,7}));

    auto a_cv1 = A.col_view(1);
    CHECK(a_cv1 == (Vector{2,5,8}));

    auto a_cv2 = A.col_view(2);
    CHECK(a_cv2 == (Vector{3,6,9}));


}
// TEST(Matrix, MatrixSlicing) {
//     Linalg::Matrix m1 = {{1, 2, 3, 4}, {1, 2, 3, 4}};

//     Linalg::Matrix expected = {{2, 3}, {2, 3}};
//     auto observed = m1.slice_columns(1, 3);
//     CHECK(observed == expected);

//     Linalg::Matrix m2 = {{1, 1}, {2, 2}, {3, 3}, {4, 4}};
//     observed = m2.slice_rows(1, 3);
//     expected = {{2, 2}, {3, 3}};
//     CHECK(observed == expected);

//     CHECK_THROWS(std::invalid_argument, m1.slice_columns(3, 2));
//     CHECK_THROWS(std::invalid_argument, m2.slice_rows(3, 2));
//     CHECK_THROWS(std::out_of_range, m1.slice_columns(1, 100));
//     CHECK_THROWS(std::out_of_range, m1.slice_rows(1, 100));
// }


