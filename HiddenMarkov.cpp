#include "HiddenMarkov.hpp"


GenotypeHMM::GenotypeHMM(const std::vector<int>& obs,
                         const std::vector<Matrix>& emissions,
                         const Linalg::Matrix& transition)
{
    observations = obs;
    nstates = transition.nrow;
    emission_matrices = emissions;
    // Matrix transition_matrix(nstates, nstates);
    transition_matrix = transition;

}



std::vector<int> GenotypeHMM::forwards_backwards(void) const
{
    size_t nobs = observations.size();
    size_t nemiss = emission_matrices[0].nrow;
    size_t nstate = emission_matrices[0].ncol;

    std::vector<int> outp(nobs);

    Matrix fwmat = Linalg::Matrix(nstates, nobs + 1, 0.0);
    Matrix bwmat = Linalg::Matrix(nstates, nobs + 1, 0.0);

    auto first_probs = fwmat.col_view(0);
    first_probs.set_all(1.0 / nstate);

    Matrix first_emiss(nemiss, nstate, 1.0 / nemiss);

    for (size_t obsidx = 0; obsidx < nobs; ++obsidx) {
        int obs = observations[obsidx];
        const Matrix& cur_obs_probs = (obsidx == 0) ? first_emiss : emission_matrices[obsidx-1];
        auto fwcv = fwmat.col_view(obsidx);


        Matrix b = Linalg::matrix_product(transition_matrix,
                                          Linalg::diag(cur_obs_probs.row_view(obs)));

        Vector col = Linalg::vector_matrix_product(fwcv, b);
        col = col / col.sum(); // Normalize column
        fwmat.set_column(obsidx + 1, col);

    }

    auto last_probs = bwmat.col_view(bwmat.ncol - 1);
    last_probs.set_all(1.0);

    for (int obsidx=nobs; obsidx>0; obsidx--) 
    {

        const Matrix& cur_obs_probs = (obsidx == 0) ? first_emiss : emission_matrices.at(obsidx-1);
        auto bwrv = bwmat.get_column(obsidx);

        ///                                                      +++++++++++++++++++++
        auto emission_diag = Linalg::diag(cur_obs_probs.get_row(observations[obsidx - 1]));
        Matrix a = Linalg::matrix_product(transition_matrix, emission_diag);

        auto col = Linalg::matrix_vector_product(a, bwrv);
        col = col / col.sum(); // Normalize again;

        bwmat.set_column(obsidx - 1, col);
        // obsidx--;
    }


    Matrix prob_mat = Linalg::direct_product(fwmat, bwmat);

// Normalize columns once more
    for (size_t i = 0; i < prob_mat.ncol; ++i) {
        auto colv = prob_mat.col_view(i);
        double s = colv.sum();
        for (size_t j = 0; j < colv.size; ++j) {
            colv.set(j, colv.get(j) / s);
        }
    }

    Vector states = prob_mat.argmax_col();
    for (size_t i = 1; i < states.size; ++i) {outp[i - 1] = states.get(i); }

    return outp;
}




std::vector<int> GenotypeHMM::viterbi(void) const
{
    using std::log;
    int nobs = observations.size();

    Linalg::Vector log_probs(nstates);

    Linalg::Vector starts(nstates);
    for (size_t i = 0; i < starts.size; ++i) { starts.set(i, i); }

    std::vector<int> outp(nobs);

    Matrix ln_transition_matrix = transition_matrix.apply(&log);

    Matrix paths = Linalg::Matrix(nstates, nobs + 1);
    paths = 0; // set all to zero;
    paths.set_column(0, starts);

    for (int obsidx = 0; obsidx < nobs; ++obsidx) {
        auto obs = observations[obsidx];

        Linalg::Vector new_log_probs(nstates);
        for (int state = 0; state < nstates; ++state) {
            Linalg::Vector temp_prob(nstates);

            temp_prob = (log_probs +
                         log(emission_matrices[obsidx].get(obs, state)) +
                         ln_transition_matrix.col_view(state));

            auto best = temp_prob.argmax();
            auto best_path = paths.row_view(best);
            paths.set_row(state, best_path);
            paths.set(state, obsidx + 1, state);
            new_log_probs.data[state] = temp_prob.get(best);
        }
        log_probs = new_log_probs;

    }

    int best_final_state = log_probs.argmax();
    Linalg::Vector hidden_states = paths.get_row(best_final_state);
    for (size_t i = 0; i < hidden_states.size; ++i) {
        outp[i] = (int)hidden_states.get(i + 1);
    }
    return outp;
}