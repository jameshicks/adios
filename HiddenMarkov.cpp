#include "HiddenMarkov.hpp"


GenotypeHMM::GenotypeHMM(const std::vector<int>& obs,
                         const std::vector<Matrix*>& emissions,
                         const Linalg::Matrix& transition)
{
    observations = obs;
    nstates = transition.nrow;
    emission_matrices = emissions;
    transition_matrix = transition;

}



std::vector<int> GenotypeHMM::forwards_backwards(void) const
{
    size_t nobs = observations.size();
    size_t nemiss = emission_matrices[0]->nrow;
    size_t nstate = emission_matrices[0]->ncol;

    std::vector<int> outp(nobs);

    Matrix fwmat = Linalg::Matrix(nstates, nobs + 1, 0.0);
    Matrix bwmat = Linalg::Matrix(nstates, nobs + 1, 0.0);


    // Temporary variables that we're gonna keep using to avoid constantly malloc/freeing new ones
    Vector col(nstate); 
    Vector v(nstate);

    Vector fw(nstates, 1.0 / nstates);
    fwmat.set_column(0, fw);
    for (size_t obsidx = 1; obsidx < (nobs+1); ++obsidx) {
        int obs = observations[obsidx-1];
        const Matrix& cur_obs_probs = *(emission_matrices[obsidx-1]);

        auto d = cur_obs_probs.row_view(obs);

        vector_matrix_product(fw, transition_matrix, &v);
        dmatrix_vector_product(d, v, &col);

        col /= col.sum(); // Normalize column
        fwmat.set_column(obsidx, col);

        fw.swap(col);
    }

    auto last_probs = bwmat.col_view(bwmat.ncol - 1);
    last_probs.set_all(1.0);

    Vector bw(nstates, 1.0);
    bwmat.set_column(bwmat.ncol - 1, bw);
    for (int obsidx=nobs; obsidx>0; obsidx--) 
    {

        const Matrix& cur_obs_probs = *(emission_matrices[obsidx-1]);

        auto d = cur_obs_probs.row_view(observations[obsidx - 1]);
        
        Linalg::dmatrix_vector_product(d, bw, &v);
        Linalg::matrix_vector_product(transition_matrix, v, &col);
        
        col /= col.sum(); // Normalize again;

        bwmat.set_column(obsidx - 1, col);
        bw.swap(col);
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
                         log(emission_matrices[obsidx]->get(obs, state)) +
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
    for (size_t i = 1; i < hidden_states.size; ++i) {
        outp[i-1] = (int)hidden_states.get(i);
    }
    return outp;
}
