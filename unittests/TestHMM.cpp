#include <vector>
#include "Linalg.hpp"
#include "HiddenMarkov.hpp"
#include "CppUTest/TestHarness.h"

#include <iostream>
TEST_GROUP(HiddenMarkov) {};

TEST(HiddenMarkov, ForwardsBackwards)
{
    // A quick sanity check. The hidden states should just be the observations
    using Linalg::Matrix;
    std::vector<int> o = {0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
    
    // Transition matrix
    Matrix t = {
        {.75, .25},
        {.25, .75}
    };

    Matrix e = {{.99, .01}, {0.01, .99}};
    // GenotypeHMM uses different emission tables for each observation
    std::vector<Matrix> ev = {};
    for (size_t i=0; i<o.size(); ++i) { ev.push_back(e); } 

    GenotypeHMM hmmfwbw = GenotypeHMM(o, ev, t);
    auto pred_states = hmmfwbw.decode(false);

    CHECK(o == pred_states);

}

// TEST(HiddenMarkov, Viterbi)
// {
//     // A quick sanity check. The hidden states should just be the observations
//     using Linalg::Matrix;
//     std::vector<int> o = {0, 0, 0, 1, 1, 1, 1, 0, 0, 0};
    
//     // Transition matrix
//     Matrix t = {
//         {.75, .25},
//         {.25, .75}
//     };

//     Matrix e = {{.99, .01}, {0.01, .99}};
//     // GenotypeHMM uses different emission tables for each observation
//     std::vector<Matrix> ev = {};
//     for (size_t i=0; i<o.size(); ++i) { ev.push_back(e); } 

//     GenotypeHMM hmmfwbw = GenotypeHMM(o, ev, t);
//     auto pred_states = hmmfwbw.decode(true);

//     CHECK(o == pred_states);

// }
