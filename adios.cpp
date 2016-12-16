#include "adios.hpp"

namespace adios
{
Matrix transition_matrix(int l10gamma, int l10rho)
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

Matrix genotype_error_matrix(double eps)
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

Matrix emission_matrix(double q)
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
    params.err_rate_common = stod(args["err"][0]);
    params.err_rate_rare = stod(args["err"][1]);
    Matrix single_errorc = genotype_error_matrix(params.err_rate_common);
    Matrix pair_errorc = Linalg::kronecker_product(single_errorc, single_errorc);
    params.allele_error_mat_common = pair_errorc;

    Matrix single_errorr = genotype_error_matrix(params.err_rate_rare);
    Matrix pair_errorr = Linalg::kronecker_product(single_errorr, single_errorr);
    params.allele_error_mat_rare = pair_errorr;

    params.gamma_ = stoul(args["transition"][0]);
    params.rho  =   stoul(args["transition"][1]);
    params.transition_mat = transition_matrix(params.gamma_, params.rho);

    params.rare_thresh = stod(args["rare"][0]);
    params.min_length = stoi(args["minlength"][0]) * 1e6;
    params.min_mark = stoi(args["minmark"][0]);
    params.min_lod = stod(args["minlod"][0]);

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

void adios_parameters::calculate_emission_mats(const Dataset& data) {
    std::set<double> fqs;
    for (auto c : data.chromosomes) {
        for (auto f : c->frequencies) {
            fqs.insert(f);
        }
    }
    for (auto fq : fqs) {
        Matrix& err = fq < rare_thresh ? allele_error_mat_rare : allele_error_mat_rare;
        Matrix m = Linalg::matrix_product(err, emission_matrix(fq));
        emission_mats[fq] = m;
    }
}

adios_sites find_informative_sites(const Indptr& ind1,
                                   const Indptr& ind2,
                                   const int chromidx,
                                   const std::vector<int>& rares)
{

    using namespace setops;

    // rv_positions are  Q & (A | B) & (C | D)
    // std::vector<int> rv_sites = multi_intersection({
    //     rares,
    //     ind1->chromosomes[chromidx]->has_minor_allele(),
    //     ind2->chromosomes[chromidx]->has_minor_allele()
    // });

    // A more inclusive is Q & (A | B | C | D)
    std::vector<int> any_minor = union_(
                                     ind1->chromosomes[chromidx]->has_minor_allele(),
                                     ind2->chromosomes[chromidx]->has_minor_allele());

    std::vector<int> rv_sites = intersection(rares, any_minor);


    // opposite_homozygotes are ((A & B) - (C | D)) | ((C & D) - (A | B))
    AlleleSites opp1 = difference(
                           ind1->chromosomes[chromidx]->homozygous_minor(),
                           ind2->chromosomes[chromidx]->has_minor_allele());

    AlleleSites opp2 = difference(
                           ind2->chromosomes[chromidx]->homozygous_minor(),
                           ind1->chromosomes[chromidx]->has_minor_allele());

    AlleleSites opposing_homozygotes = union_(opp1, opp2);


    AlleleSites missing = union_(ind1->chromosomes[chromidx]->missing,
                                 ind2->chromosomes[chromidx]->missing);
    // Informative sites are the union of opposing homozygotes and shared rv sites
    AlleleSites informative_sites = union_(opposing_homozygotes, rv_sites);
    informative_sites = difference(informative_sites, missing);

    // Since we're using a subset of markers, we should translate the
    // informative indexes into a newer, more convenient index space.
    // We'll translate back later.
    std::map<int, int> index_translation;

    int newidx = 0;
    for (auto it = informative_sites.begin();
            it != informative_sites.end();
            ++it) {
        index_translation[*it] = newidx;
        newidx++;
    }




    // Restrict to the informative sites
    AlleleSites a = intersection(ind1->chromosomes[chromidx]->hapa,
                                 informative_sites);
    for (size_t i = 0; i < a.size(); ++i) {a[i] = index_translation.at(a[i]); }
    AlleleSites b = intersection(ind1->chromosomes[chromidx]->hapb,
                                 informative_sites);
    for (size_t i = 0; i < b.size(); ++i) {b[i] = index_translation.at(b[i]); }
    AlleleSites c = intersection(ind2->chromosomes[chromidx]->hapa,
                                 informative_sites);
    for (size_t i = 0; i < c.size(); ++i) {c[i] = index_translation.at(c[i]); }
    AlleleSites d = intersection(ind2->chromosomes[chromidx]->hapb,
                                 informative_sites);
    for (size_t i = 0; i < d.size(); ++i) {d[i] = index_translation.at(d[i]); }


    // AuB and CuD are A ∪ B and C ∪ D
    auto AuB = union_(a, b);
    auto CuD = union_(c, d);

    // AnB, CnD are A ∩ B, C ∩ D
    auto AnB = intersection(a, b);
    auto CnD = intersection(c, d);

    // AsB,CsD are A ⊖ B, C ⊖ D
    auto AsB = symmetric_difference(a, b);
    auto CsD = symmetric_difference(c, d);


    // Generate the positions for each genotype configuration
    std::vector<AlleleSites> configurations = {
        std::vector<int>(),       // 0 AA,AA (Not selected)
        difference(CsD, AuB),   // 1 AA,AB
        difference(CnD, AuB),   // 2 AA,BB
        difference(AsB, CuD),   // 3 AB,AA
        intersection(AsB, CsD), // 4 AB,AB
        intersection(AsB, CnD), // 5 AB,BB
        difference(AnB, CuD),   // 6 BB,AA
        intersection(AnB, CsD), // 7 BB,AB
        intersection(AnB, CnD)  // 8 BB,BB
    };


    // Put these into a vector called states.
    std::vector<int> states(informative_sites.size(), 0);
    for (size_t state = 0; state < configurations.size(); ++state) {
        auto sites = configurations[state];

        for (auto site : sites) {

            states[site] = state;
        }
    }


    // We need to send the states back with the sites, so that we know
    // what they correspond to.
    auto p = make_pair(states, informative_sites);

    return p;
}

AlleleSites update_indices(const AlleleSites& inp,
                           const std::map<int, int>& translator)
{
    AlleleSites output;
    for (auto it = inp.begin();
            it != inp.end();
            ++it) {
        output.push_back(translator.at(*it));
    }
    return output;
}

void adios(Dataset& d, const adios_parameters& params)
{
    std::vector<Indptr> inds;
    for (auto pr = d.individuals.begin(); pr != d.individuals.end(); ++pr) { inds.push_back(pr->second);}

    std::cout << "IND_1\tIND_2\tCHROM\tSTART\tEND\tLENGTH\tSTATE\tNMARK\tNRARE\tNERR\tLOD\n";
    std::vector<Indptr_pair> pairs = combinatorics::pair_combinations(inds);
    for (size_t chridx = 0; chridx < d.nchrom(); chridx++) {
        for (size_t pairidx = 0; pairidx < pairs.size(); ++pairidx) {
            Indptr_pair pair = pairs[pairidx];

            adios_pair(pair, chridx, params);

        }
    }
}


void adios_pair(const Indptr_pair& inds,
                int chromidx,
                const adios_parameters& params)
{
    using Linalg::Matrix;

    Indptr ind1 = inds.first;
    Indptr ind2 = inds.second;

    auto chromobj = ind1->chromosomes[chromidx]->info;

    auto useful = find_informative_sites(ind1,
                                         ind2,
                                         chromidx,
                                         params.rare_sites[chromidx]);

    auto observations = useful.first;
    std::vector<int> informative_sites(useful.second.begin(), useful.second.end());
    int nmark = informative_sites.size();


    // Get the frequencies
    std::vector<double> freqs(nmark);
    for (int i = 0; i < nmark; ++i) {
        freqs[i] = chromobj->frequencies[i];
    }

    // Make the emission matrices
    std::vector<Matrix> emissions;
    for (int i = 0; i < nmark; ++i) {
        // Lookup the precomputed matrix
        Matrix emiss = params.emission_mats.at(freqs[i]);
        emissions.push_back(emiss);
    }

    GenotypeHMM model(observations, emissions, params.transition_mat);
    std::vector<int> hidden_states = model.decode(false);

    std::vector<ValueRun> runs = runs_gte_classic(hidden_states, 1, 5);

    for (ValueRun r : runs) {
        Segment seg(ind1, ind2, r, chromobj,
                    observations, emissions,
                    informative_sites, params);

        if (seg.passes_filters(params)) {
            std::cout << seg.record_string() << '\n';
        }
    }

}

Segment::Segment(Indptr a,
                 Indptr b,
                 ValueRun& run,
                 Chromptr c,
                 std::vector<int>& obs,
                 std::vector<Matrix>& emissions,
                 std::vector<int>& adiossites,
                 const adios_parameters& params)
{

    ind1 = a;
    ind2 = b;
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
    for (size_t i = start; i<=stop; ++i) {
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
                              std::vector<Matrix>& emissions,
                              const adios::adios_parameters& params) const
{
    using std::log;
    if (start >= stop) { return -1e99; }
    // Transition probs for entering and exiting a state
    double entry = params.transition_mat.get(0, state);
    double exit = params.transition_mat.get(state, 0);

    // Transition probs for remaining in a state
    double remain_state = params.transition_mat.get(state, state);
    double remain_null = params.transition_mat.get(0, 0);
    double log_ibd_prob = log(remain_state) * (nmark - 1);
    double log_null_prob = log(remain_null) * (nmark - 1);

    for (size_t i = start; i < stop; ++i) {
        int obs = observations[i];
        log_ibd_prob += log(emissions[i].get(obs, state));
        log_null_prob += log(emissions[i].get(obs, 0));
    }

    log_ibd_prob += log(entry) + log(exit);
    log_null_prob += 2 * log(remain_state);

    return log_ibd_prob - log_null_prob;
}


bool Segment::passes_filters(const adios::adios_parameters& params) const
{
    return ((length() >= params.min_length) &&
            (nmark >= params.min_mark) &&
            (lod >= params.min_lod));

}

std::string Segment::record_string(void) const
{
    std::stringstream s;
    s << ind1->label << '\t';
    s << ind2->label << '\t';
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

