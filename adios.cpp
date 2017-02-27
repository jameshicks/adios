#include "adios.hpp"

namespace adios
{
Matrix unphased_transition_matrix(int l10gamma, int l10rho)
{
    const double gamma = pow(10, -l10gamma);
    const double rho   = pow(10, -l10rho);
    const double gamma_2 = pow(gamma, 2);
    const double rho_2 = pow(rho, 2);

    Matrix t = {
        { 1 - gamma - gamma_2,   gamma,             gamma_2 },
        { rho,                   1 - gamma - rho,   gamma },
        { rho_2,                 rho,               1 - rho - rho_2 }
    };
    return t;
}

Matrix unphased_genotype_error_matrix(double eps)
{
    const double eta = 1 - eps;

    const double eta_sq = pow(eta, 2);
    const double eps_sq = pow(eps, 2);
    const double eta_eps = eta * eps;

    Matrix outp = {
        { eta_sq,   2 * eta_eps,      eps_sq },
        { eta_eps,  eta_sq + eps_sq,  eta_eps },
        { eps_sq,   2 * eta_eps,      eta_sq }
    };

    return outp;
}

Matrix unphased_emission_matrix(double q)
{
    Matrix mat(9, 3);

    const double p = 1 - q;
    const double p_2 = pow(p, 2);
    const double p_3 = pow(p, 3);
    const double q_2 = pow(q, 2);
    const double q_3 = pow(q, 3);
    const double pq = p * q;
    const double p_2q_2 = pow(pq, 2);

    mat = {
        // P(x|IBD=0)     P(x|IBD=1)      P(x|IBD=2)
        { pow(p, 4),      p_3,        p_2 },       // 0: AA,AA
        { p_3 * q,        p_2 * q,    0 },         // 1: AA,AB
        { p_2 * q_2,      0,          0 },         // 2: AA,BB
        { p_3 * q,        p_2 * q,    0 },         // 3: AB,AA
        { 4 * p_2q_2,     pq,         2 * pq },    // 4: AB,AB
        { p * q_3,        p * q_2,    0 },         // 5: AB,BB
        { p_2q_2,         0,          0 },         // 6: BB,AA
        { p * q_3,        p * q_2,    0 },         // 7: BB,AB
        { pow(q, 4),      q_3,        q_2 }        // 8: BB,BB
    };

    return mat;
}


adios_parameters params_from_args(std::map<std::string, std::vector<std::string>> args)
{
    using std::stod;
    using std::stoi;
    adios_parameters params;

    // The error matrix we store is the genotype error matrix for an individual pair
    params.err_rate = stod(args["err"][0]);
    Matrix single_error = unphased_genotype_error_matrix(params.err_rate);
    Matrix pair_error = Linalg::kronecker_product(single_error, single_error);
    params.unphased_error_mat = pair_error;

    params.gamma_ = stoul(args["transition"][0]);
    params.rho  =   stoul(args["transition"][1]);
    params.unphased_transition_mat = unphased_transition_matrix(params.gamma_, params.rho);

    params.rare_thresh = stod(args["rare"][0]);
    params.min_length = stod(args["minlength"][0]) * 1e6;
    params.min_mark = stoi(args["minmark"][0]);
    params.min_lod = stod(args["minlod"][0]);

    params.viterbi = (args["viterbi"][0].compare("YES") == 0);

    return params;

}

void adios_parameters::get_rare_sites(Dataset& data)
{
    for (size_t chridx = 0; chridx < data.nchrom(); ++chridx) {
        std::shared_ptr<ChromInfo> chrom = data.chromosomes[chridx];
        std::vector<int> rares;
        for (size_t markidx = 0; markidx < chrom->nmark(); ++markidx) {
            if (chrom->frequencies[markidx] < rare_thresh) {
                rares.push_back(markidx);
            }
        }
        rare_sites.push_back(rares);
    }
}

void adios_parameters::calculate_emission_mats(const Dataset& data)
{
    std::set<double> fqs;
    for (auto c : data.chromosomes) {
        for (auto f : c->frequencies) {
            fqs.insert(f);
        }
    }
    for (auto fq : fqs) {
        Matrix& err = unphased_error_mat;
        Matrix m = Linalg::matrix_product(err, unphased_emission_matrix(fq));
        emission_mats[fq] = m;
    }
}


adios_sites find_informative_sites_unphased(const Individual& ind1,
    const Individual& ind2,
    const int chromidx,
    const std::vector<int>& rares)
    {
        auto chromobj = ind1.chromosomes[chromidx].info;
        int max_pos = chromobj->positions.back() + 100;

        auto hapa = ind1.chromosomes[chromidx].hapa;
        auto hapb = ind1.chromosomes[chromidx].hapb;
        auto hapc = ind2.chromosomes[chromidx].hapa;
        auto hapd = ind2.chromosomes[chromidx].hapb;

        std::vector<int> cur_idx = {0, 0, 0, 0};
        std::vector<int> cur_vars = {0, 0, 0, 0};

        int rareidx = 0;
        int nrare = rares.size();

        auto miss = setops::union_(ind1.chromosomes[chromidx].missing,
                                   ind2.chromosomes[chromidx].missing);

        int missidx = 0;
        int nmiss = miss.size();

        const int nmark = chromobj->nmark();

        std::vector<int> informatives;
        std::vector<int> states;

        // You cant know ahead of time how many sites are going to be useful
        // but in my experience it's less than 5%
        const int reserve_amount = (int)(0.05 * nmark);
        informatives.reserve(reserve_amount);
        states.reserve(reserve_amount);

        int current_position = 0;
        do {
            cur_vars[0] = cur_idx[0] >= hapa.size() ? max_pos : hapa[cur_idx[0]];
            cur_vars[1] = cur_idx[1] >= hapb.size() ? max_pos : hapb[cur_idx[1]];
            cur_vars[2] = cur_idx[2] >= hapc.size() ? max_pos : hapc[cur_idx[2]];
            cur_vars[3] = cur_idx[3] >= hapd.size() ? max_pos : hapd[cur_idx[3]];

            int current_position = *(std::min_element(cur_vars.begin(),
                                     cur_vars.end()));

            if (current_position == max_pos) { break; }

            while (rareidx < nrare && rares[rareidx] < current_position) { rareidx++; }
            while (missidx < nmiss && miss[missidx]  < current_position) { missidx++; }

            bool is_rare = current_position == rares[rareidx];
            bool is_miss = (nmiss > 0 && current_position == miss[missidx]);

            int s1 = (cur_vars[0] == current_position) + (cur_vars[1] == current_position);
            int s2 = (cur_vars[2] == current_position) + (cur_vars[3] == current_position);
            int state = 3 * s1 + s2;

            if ((state == 2 || state == 6 || is_rare) && !is_miss) {
                informatives.push_back(current_position);
                states.push_back(state);
            }

            for (int i = 0; i < 4; i++) {
                if (cur_vars[i] == current_position) { cur_idx[i]++; }
            }



        } while (current_position != max_pos);

        return make_pair(states, informatives);

    }



void adios(Dataset& d, const adios_parameters& params, DelimitedFileWriter& out)
{
    using namespace combinatorics;

    std::vector<Individual> inds;
    for (auto pr = d.individuals.begin(); pr != d.individuals.end(); ++pr) { inds.push_back(pr->second);}

    std::vector<std::string> header = {
        "IND_1", "IND_2", "CHROM", "START", "END", "LENGTH",
        "STATE", "NMARK", "NRARE", "NERR", "LOD"
    };
    out.writetoks(header);

    long ninds = inds.size();
    long npairs = nCk(ninds, 2);

    for (size_t chridx = 0; chridx < d.nchrom(); chridx++) {
        double signpost = 0.0; 
        double signpost_step = npairs > 100000 ? 0.001 : 0.01;
        int completed = 0;
        unsigned long markers_used = 0;
        unsigned long total_mark = d.chromosomes[chridx]->nmark();

        // If openmp is available, this is the loop we want to parallelize.
        // This gives each thread a set of individual pairs to compute.
        #pragma omp parallel for
        for (size_t pairidx = 0; pairidx < npairs; ++pairidx) {
            std::vector<long> indices = combination_at_index(pairidx,
                                        ninds,
                                        2);

            Individual& ind1 = inds[indices[0]];
            Individual& ind2 = inds[indices[1]];

            adios_result res = adios_pair_unphased(ind1, ind2, chridx, params);

            #pragma omp critical
            {
                // File IO needs to be locked. This block is OMP critical
                // to prevent output (both to stdout and file) from being
                // garbled.
                for (Segment s : res.segments) { 
                    out.writetoks(s.record()); 
                }
                
                markers_used += res.nmark;
                completed++;
                double progress = (double)completed / (double)(npairs);
                
                if (progress > signpost) {
                    double mean_mark = markers_used / (double)completed;

                    if (!out.is_stdout()) {
                        std::cout << "\rChromosome " << d.chromosomes[chridx]->label;
                        std::cout << ": " << sfloat(progress * 100, 1) << "% complete. ";
                        std::cout << "Average markers per pair " << sfloat(mean_mark, 2);
                        std::cout << " (" << sfloat(100 * mean_mark / total_mark, 3) << "%)";
                        std::cout << std::flush;
                    }
                    while (progress > signpost) { signpost += signpost_step; }
                }

            }

        }
    }
    if (!out.is_stdout()) { std::cout << '\n' << std::flush;  }
}


adios_result adios_pair_unphased(const Individual& ind1, const Individual& ind2,
        int chromidx,
        const adios_parameters& params)
{
    using Linalg::Matrix;

    adios_result res; 

    auto chromobj = ind1.chromosomes[chromidx].info;

    auto useful = find_informative_sites_unphased(ind1,
                                                  ind2,
                                                  chromidx,
                                                  params.rare_sites[chromidx]);

    auto observations = useful.first;
    std::vector<int> informative_sites(useful.second.begin(), useful.second.end());
    int nmark = informative_sites.size();
    res.nmark = nmark;

    // Make the emission matrices
    std::vector<Matrix*> emissions;
    emissions.reserve(nmark);
    for (int i = 0; i < nmark; ++i) {
        // Save the precomputed matrices
        Matrix* mp = const_cast<Matrix*>(&(params.emission_mats.at(chromobj->frequencies[i])));
        emissions.push_back(mp);
    }

    GenotypeHMM model(observations, emissions, params.unphased_transition_mat);
    std::vector<int> hidden_states = model.decode(params.viterbi);

    std::vector<ValueRun> runs = runs_gte_classic(hidden_states, 1, 5);

    for (ValueRun r : runs) {
        Segment seg(ind1, ind2, r, chromobj,
                    observations, emissions,
                    informative_sites, params);

        if (seg.passes_filters(params)) {
            res.segments.push_back(seg);
        }
    }

    return res;
}

Segment::Segment(const Individual& a,
                 const Individual& b,
                 ValueRun& run,
                 Chromptr c,
                 std::vector<int>& obs,
                 std::vector<Matrix*>& emissions,
                 std::vector<int>& adiossites,
                 const adios_parameters& params)
{

    ind1 = a.label;
    ind2 = b.label;
    chrom = c;

    start = run.start;
    stop = run.stop - 1; // Segments are end-inclusive for convenience.
    state = run.value;

    trim(obs); // Trim back to the first shared_rvs
    nmark = stop - start;
    full_start = adiossites[start];
    full_stop = adiossites[stop];

    lod = calculate_lod(obs, emissions, params);

    nerr = 0;
    nrare = 0;
    for (size_t i = start; i <= stop; ++i) {
        nerr += (int)(obs[i] == 2 || obs[i] == 6);
        nrare += (int)(is_shared_rv(obs[i]));
    }

}


void Segment::trim(std::vector<int>& observations)
{
    // return;
    using adios::is_shared_rv;
    while (!is_shared_rv(observations[stop])) {

        if (stop == start) { stop = 0; start = 0; return; }
        stop--;
    }
    while (!is_shared_rv(observations[start])) {

        if (stop == start) { stop = 0; start = 0; return; }
        start++;
    }
}

double Segment::calculate_lod(std::vector<int>& observations,
                              std::vector<Matrix*>& emissions,
                              const adios::adios_parameters& params) const
{
    using std::log10;
    if (start >= stop) { return -1e99; }
    // Transition probs for entering and exiting a state
    double entry = params.unphased_transition_mat.get(0, state);
    double exit = params.unphased_transition_mat.get(state, 0);

    // Transition probs for remaining in a state
    double remain_state = params.unphased_transition_mat.get(state, state);
    double remain_null = params.unphased_transition_mat.get(0, 0);
    double log_ibd_prob = log10(remain_state) * (nmark - 1);
    double log_null_prob = log10(remain_null) * (nmark - 1);

    for (size_t i = start; i < stop; ++i) {
        int obs = observations[i];
        log_ibd_prob += log10(emissions[i]->get(obs, state));
        log_null_prob += log10(emissions[i]->get(obs, 0));
    }

    log_ibd_prob += log10(entry) + log10(exit);
    log_null_prob += 2 * log10(remain_state);

    return log_ibd_prob - log_null_prob;
}


bool Segment::passes_filters(const adios::adios_parameters& params) const
{
    return ((length() >= params.min_length) &&
            (nrare >= params.min_mark) &&
            (lod >= params.min_lod));

}

std::vector<std::string> Segment::record(void) const
{
    using std::to_string;

    std::vector<std::string> s = {
        ind1,
        ind2,
        chrom->label,
        to_string(chrom->positions[full_start]),
        to_string(chrom->positions[full_stop]),
        bp_formatter(length()),
        to_string(state),
        to_string(nmark),
        to_string(nrare),
        to_string(nerr),
        sfloat(lod, 2)
    };
    return s;
}

std::string Segment::record_string(void) const
{

    std::stringstream s;
    s << ind1 << '\t';
    s << ind2 << '\t';
    s << chrom->label << '\t';
    s << chrom->positions[full_start] << '\t';
    s << chrom->positions[full_stop]  << '\t';
    s << bp_formatter(length()) << '\t';
    s << state << '\t';
    s << nmark << '\t';
    s << nrare << '\t';
    s << nerr << '\t';
    s << sfloat(lod, 2);

    return s.str();
}



}

