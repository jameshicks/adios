#include "vcf.hpp"

// VCF Parsing

VCFRecordGenotypeContainer::VCFRecordGenotypeContainer(size_t n) {
    alts.reserve(2 * n);
    missing.reserve(n);
}


void VCFRecord::set_freq(const std::string& info_field) {
    std::string val = get_info_by_key(info_field.c_str());
    if (val.empty()) freq = 0.0;
    freq = atof(val.c_str());
}

std::string VCFRecord::get_info_by_key(const char* key) {
    size_t keylen = strlen(key);
    std::string retval("");
    char* orig_infodup = strdup(infostr.c_str());
    char* infodup = orig_infodup;

    while (char* token = strsep(&infodup, ";")) {

        if (token == NULL) { break; }
        if (strlen(token) < keylen) { continue; }

        char* subtok = strsep(&token, "=");
        if (strcmp(subtok, key) == 0) {
            char* next_subtok = strsep(&token, "=");
            if (next_subtok != NULL) { retval = next_subtok; }
            else { retval = "YES"; }
            break;
        }

    }
    if (retval.empty()) { retval = ""; }
    free(orig_infodup);

    return retval;
}

VCFRecord::VCFRecord(const std::string& vcfline)  {
    using stringops::split;
    using std::string;

    char* const orig_cstr = strdup(vcfline.c_str());

    char* cstr = orig_cstr;

    char* token = strsep(&cstr, "\t"); // Field 0: Chrom
    chrom = string(token);

    token = strsep(&cstr, "\t"); // Field 1: Pos
    pos = atoi(token);

    token = strsep(&cstr, "\t"); // Field 2: label
    label = string(token);

    token = strsep(&cstr, "\t"); // Field 3: Ref
    alleles.reserve(4);
    alleles.push_back(string(token));

    token = strsep(&cstr, "\t"); // Field 4: alts
    std::vector<std::string> alts = split(string(token), ",");
    for (size_t i = 0; i < alts.size(); ++i) {
        alleles.push_back(alts[i]);
    }

    token = strsep(&cstr, "\t"); // Field 5; Qual
    token = strsep(&cstr, "\t"); // Field 6: Filter

    token = strsep(&cstr, "\t"); // Field 7: INFO
    infostr = string(token);

    token = strsep(&cstr, "\t"); // Field 8: Format
    format = string(token);

    // All the rest of the data is hiding in the rest of the string, but we
    // dont want keep tokenizing, so we'll just set our token to where it's
    // already at, plus strlen(token), so we'll be on the \0 that used to be
    // the delimiter. plus one more and we're right where we want to be.
    char* dataptr = token + strlen(token) + 1;
    data = string(dataptr);

    free(orig_cstr); // don't need this guy anymore.

    freq = 0.0;
}

std::map<std::string, std::string> VCFRecord::infomap(void) const {
    // Make the info map

    using stringops::split;
    std::map<std::string, std::string> info;

    std::vector<std::string> infotoks = split(infostr, ";");
    std::string tok;
    for (size_t i = 0; i < infotoks.size(); ++i) {
        tok = infotoks[i];
        std::vector<std::string> kv = split(tok, "=");
        if (kv.size() == 1) {
            info[kv[0]] = std::string("YES");
        } else if (kv.size() == 2) {
            info[kv[0]] = kv[1];
        } else {
            std::cout << "welp,\n";
        }
    }
    return info;
}

void VCFRecord::get_minor_alleles(VCFRecordGenotypeContainer& con) const {
    int gtidx = 0;
    auto gtfpos = format.find("GT");
    for (size_t i = 0; i < gtfpos; ++i) {
        if (format[i] == ':') gtidx++;
    }

    // std::cout << "gtidx is " << gtidx << std::endl;
    char* original_data_pointer = strdup(data.c_str());
    char* cdata = original_data_pointer;
    int indidx = -1;
    while (char* token = strsep(&cdata, " \t")) {
        indidx++;

        int subtokidx = 0;
        char* gttok = strsep(&token, ":");
        while (subtokidx < gtidx) {
            gttok = strsep(&token, ":");
            subtokidx++;
        }


        int gtlen = strlen(gttok);
        if (gtlen == 3) {

            if (gttok[0] == '.' || gttok[2] == '.') {
                con.missing.push_back(indidx);
            } else {
                const int allele_a = ((int)gttok[0] - 48);
                const int allele_b = ((int)gttok[2] - 48);

                if (allele_a != 0) {
                    con.alts.push_back(2 * indidx);
                }

                if (allele_b != 0) {
                    con.alts.push_back(2 * indidx + 1);
                }
            }
        } else  {
            std::cerr << "Malformed genotype: \'" << gttok << "\' ";
            std::cerr << "at " << chrom << ':' << pos << " (" << label << ")";
            std::cerr << " for individual at index " << indidx << ". ";
            std::cerr << "Marked as missing." << std::endl;

            con.missing.push_back(indidx);

        }

    }
    free(original_data_pointer);

}

int VCFRecord::nalleles(void) const {
    return alleles.size();
}

bool VCFRecord::is_snv(void) const {
    for (int i = 0; i < nalleles(); ++i) {
        if (alleles[i].length() > 2) {
            return false;
        }
    }
    return true;
}

Dataset read_vcf(const std::string & filename, const std::string & freq_field) {
    using stringops::split;
    using std::vector;

    Dataset data;
    std::ifstream vcffile(filename);

    if (!vcffile) { throw std::invalid_argument("Unable to open file"); }

    std::vector<std::string> indlabs;
    std::string line;
    while (getline(vcffile, line)) {
        if (stringops::startswith(line, "##")) {
            continue;
        } else if (stringops::startswith(line, "#")) {
            indlabs = split(line, "\t");
            indlabs = slice(indlabs, 9, indlabs.size());
            for (size_t i = 0; i < indlabs.size(); ++i) {
                data.add_individual(indlabs[i]);
            }
            break;
        } else {
            std::cout << "No header line??" << std::endl;
        }
    }

    std::vector<std::shared_ptr<Individual>> inds;
    for (auto v : indlabs) { inds.push_back(data.individuals[v]); }
    int ninds = inds.size();

    VCFRecordGenotypeContainer con(ninds);

    std::string last_chromid("");
    int chromidx = -1;
    int markidx = 0;
    unsigned long rawidx = 0;
    while (getline(vcffile, line)) {
        // if (rawidx % 10000 == 0) { std::cout << "Line " << rawidx << std::endl;}
        if (!line.length()) continue;
        VCFRecord rec(line);
        con.clear();

        if (rec.chrom.compare(last_chromid)) {
            if (chromidx >-1 && !data.chromosomes[chromidx]->size()) {
                data.chromosomes[chromidx]->label = rec.chrom;
            } else { 
                
            data.add_chromosome(rec.chrom);
                chromidx++;
            }
            markidx = 0;
        }

        if (!rec.is_snv()) {
            data.chromosomes[chromidx]->exclusions["Non-SNV"]++;
            continue;
        }


        if (rec.nalleles() > 2) {
            data.chromosomes[chromidx]->exclusions["Non-diallelic"]++;
            continue;
        }


        rec.get_minor_alleles(con);

        // if (con.monomorphic()) {
        //     data.chromosomes[chromidx]->exclusions["monomorphic"]++;
        //     continue;
        // }

        // if (con.singleton()) {
        //     data.chromosomes[chromidx]->exclusions["Singleton"]++;
        //     continue;
        // }

        rec.set_freq(freq_field);
        data.chromosomes[chromidx]->add_variant(rec.label, rec.pos, rec.freq);


        for (size_t i = 0; i < con.missing.size(); ++i) {
            inds[i]->set_allele(chromidx, markidx, 0, -1);
        }

        // std::cout << "minor alleles at : ";
        // for (auto m : con.alts) { std::cout << m << ", "; } std::cout << '\n';

        for (size_t i = 0; i < con.alts.size(); ++i) {
            std::div_t divres = std::div(con.alts[i], 2);
            auto ind = inds[divres.quot];
            // std::cout << "setting allele: " << divres.quot << ',' << divres.rem << '\n';
            ind->set_allele(chromidx, markidx, divres.rem, 1);
        }

        markidx++;
        rawidx++;

        last_chromid = rec.chrom;
    }
    return data;
}
