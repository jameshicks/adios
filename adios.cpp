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
    params.err_rate_common = stod(args["err"][0]);
    params.err_rate_rare = stod(args["err"][1]);
    Matrix single_errorc = unphased_genotype_error_matrix(params.err_rate_common);
    Matrix pair_errorc = Linalg::kronecker_product(single_errorc, single_errorc);
    params.unphased_error_mat_common = pair_errorc;

    Matrix single_errorr = unphased_genotype_error_matrix(params.err_rate_rare);
    Matrix pair_errorr = Linalg::kronecker_product(single_errorr, single_errorr);
    params.unphased_error_mat_rare = pair_errorr;

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

void adios_parameters::calculate_emission_mats(const Dataset& data) {
    std::set<double> fqs;
    for (auto c : data.chromosomes) {
        for (auto f : c->frequencies) {
            fqs.insert(f);
        }
    }
    for (auto fq : fqs) {
        Matrix& err = fq < rare_thresh ? unphased_error_mat_rare : unphased_error_mat_rare;
        Matrix m = Linalg::matrix_product(err, unphased_emission_matrix(fq));
        emission_mats[fq] = m;
    }
}

// We need to get the indices in hap to correspond to the values in informative_sites
AlleleSites update_indices(const AlleleSites& hap, const AlleleSites& informative_sites) {
    AlleleSites a;

    size_t hapsize = hap.size();
    size_t ninform = informative_sites.size();
    size_t j = 0; 
    for (size_t i = 0; i < ninform; ++i) {
        while (hap[j] < informative_sites[i] && j < hapsize) ++j;
        if (j >= hapsize) break;
        if (hap[j] == informative_sites[i]) a.push_back(i);
    }
    return a;
}


adios_sites find_informative_sites_unphased(const Indptr& ind1,
                                            const Indptr& ind2,
                                            const int chromidx,
                                            const std::vector<int>& rares)
{

    using namespace setops;

    auto hma1 = ind1->chromosomes[chromidx]->has_minor_allele();
    auto hma2 = ind2->chromosomes[chromidx]->has_minor_allele();

    // we're looking for Q & (A | B | C | D) which is the same as 
    // (Q&A | Q&B | Q&C | Q&D), but we'll have to do some testing to see which
    // is faster
    auto hra1 = intersection(rares, hma1);
    auto hra2 = intersection(rares, hma2);

    std::vector<int> any_rvs = union_(hra1, hra2);

    // opposite_homozygotes are ((A & B) - (C | D)) | ((C & D) - (A | B))
    AlleleSites opp1 = difference(
                           ind1->chromosomes[chromidx]->homozygous_minor(),
                           hma2);

    AlleleSites opp2 = difference(
                           ind2->chromosomes[chromidx]->homozygous_minor(),
                           hma1);

    AlleleSites opposing_homozygotes = union_(opp1, opp2);


    AlleleSites missing = union_(ind1->chromosomes[chromidx]->missing,
                                 ind2->chromosomes[chromidx]->missing);
    // Informative sites are the union of opposing homozygotes and shared rv sites
    AlleleSites informative_sites = union_(opposing_homozygotes, any_rvs);
    informative_sites = difference(informative_sites, missing);


    AlleleSites a = update_indices(ind1->chromosomes[chromidx]->hapa, informative_sites);
    AlleleSites b = update_indices(ind1->chromosomes[chromidx]->hapb, informative_sites);
    AlleleSites c = update_indices(ind2->chromosomes[chromidx]->hapa, informative_sites);
    AlleleSites d = update_indices(ind2->chromosomes[chromidx]->hapb, informative_sites);

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


void adios(Dataset& d, const adios_parameters& params, DelimitedFileWriter& out)
{
    std::vector<Indptr> inds;
    for (auto pr = d.individuals.begin(); pr != d.individuals.end(); ++pr) { inds.push_back(pr->second);}

        std::vector<std::string> header = {
            "IND_1", "IND_2", "CHROM", "START", "END", "LENGTH",
            "STATE", "NMARK", "NRARE", "NERR", "LOD"
        };
        out.writetoks(header);
        std::vector<Indptr_pair> pairs = combinatorics::pair_combinations(inds);
    
    for (size_t chridx = 0; chridx < d.nchrom(); chridx++) {
        
        double signpost = 0.0; double signpost_step = 0.01; 
        for (size_t pairidx = 0; pairidx < pairs.size(); ++pairidx) {
            Indptr_pair pair = pairs[pairidx];

            double progress = pairidx / (double)(pairs.size());
            if (progress > signpost) {
                if (!out.is_stdout()) std::cout << sfloat(progress * 100, 1) << '%' << "..." << std::endl;
                while (progress > signpost) signpost += signpost_step;
            }

            std::vector<Segment> segs = adios_pair_unphased(pair, chridx, params);  
            for (Segment s : segs) out.writetoks(s.record());
        }
    }
}


std::vector<Segment> adios_pair_unphased(const Indptr_pair& inds,
                int chromidx,
                const adios_parameters& params)
{
    using Linalg::Matrix;

    Indptr ind1 = inds.first;
    Indptr ind2 = inds.second;

    auto chromobj = ind1->chromosomes[chromidx]->info;

    auto useful = find_informative_sites_unphased(ind1,
                                                  ind2,
                                                  chromidx,
                                                  params.rare_sites[chromidx]);

    auto observations = useful.first;
    std::vector<int> informative_sites(useful.second.begin(), useful.second.end());
    int nmark = informative_sites.size();


    // Make the emission matrices
    std::vector<Matrix> emissions;
    emissions.reserve(nmark);
    for (int i = 0; i < nmark; ++i) {
        // Lookup the precomputed matrix
        Matrix emiss = params.emission_mats.at(chromobj->frequencies[i]);
        emissions.push_back(emiss);
    }

    GenotypeHMM model(observations, emissions, params.unphased_transition_mat);
    std::vector<int> hidden_states = model.decode(params.viterbi);

    std::vector<ValueRun> runs = runs_gte_classic(hidden_states, 1, 5);

    std::vector<Segment> goodsegs;
    for (ValueRun r : runs) {
        Segment seg(ind1, ind2, r, chromobj,
                    observations, emissions,
                    informative_sites, params);

        if (seg.passes_filters(params)) {
            goodsegs.push_back(seg);
        }
    }
    return goodsegs;
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
        log_ibd_prob += log10(emissions[i].get(obs, state));
        log_null_prob += log10(emissions[i].get(obs, 0));
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

std::vector<std::string> Segment::record(void) const {
        using std::to_string;

    std::vector<std::string> s = {
        ind1->label,
        ind2->label,
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

