// guess phase in imputed genotypes

#include "guess_phase.h"
#include <Rcpp.h>
#include "cross.h"
#include "cross_util.h"

// guess phase, F2 autosome
// [[Rcpp::export(".guess_phase_f2A")]]
IntegerVector guess_phase_f2A(const IntegerMatrix& geno) // pos x ind
{
    const int n_pos = geno.rows();
    const int n_ind = geno.cols();
    const int matsize = n_pos*2;

    IntegerVector result(n_ind*n_pos*2); // this will be 2 x n_pos x n_ind

    for(int ind=0; ind < n_ind; ind++) {
        IntegerVector g1(n_pos), g2(n_pos);
        for(int pos=0; pos < n_pos; pos++) {
            int g = geno(pos,ind);
            if(IntegerVector::is_na(g) || g==0)
                g1[pos] = g2[pos] = NA_INTEGER;
            else if(g==1)
                g1[pos] = g2[pos] = 1;
            else if(g==3)
                g1[pos] = g2[pos] = 2;
            else {
                g1[pos] = 1;
                g2[pos] = 2;
            }
        }

        IntegerVector phased_geno = phase_geno(g1, g2);
        for(int i=0; i<matsize; i++)
            result[ind*matsize+i] = phased_geno[i];

    }

    result.attr("dim") = Dimension(2, n_pos, n_ind);
    return result;
}


// guess phase, F2 X chr
// [[Rcpp::export(".guess_phase_f2X")]]
IntegerVector guess_phase_f2X(const IntegerMatrix& geno) // pos x ind
{
    const int n_pos = geno.rows();
    const int n_ind = geno.cols();

    IntegerVector result(n_ind*n_pos*2); // this will be 2 x n_pos x n_ind

    for(int ind=0, offset=0; ind < n_ind; ind++) {
        for(int pos=0; pos < n_pos; pos++, offset += 2) {
            int g = geno(pos,ind);
            if(IntegerVector::is_na(g) || g==0)
                result[offset] = result[offset+1] = NA_INTEGER;
            else if(g==1) {
                result[offset] = result[offset+1] = 1;
            }
            else if(g==2) {
                result[offset] = 2;
                result[offset+1] = 1;
            }
            else if(g==3) {
                result[offset] = 1;
                result[offset+1] = 2;
            }
            else if(g==4) {
                result[offset] = result[offset+1] = 2;
            }
            else if(g==5) {
                result[offset] = 1;
                result[offset+1] = NA_INTEGER;
            }
            else if(g==6) {
                result[offset] = 2;
                result[offset+1] = NA_INTEGER;
            }
        }
    }

    result.attr("dim") = Dimension(2, n_pos, n_ind);
    return result;
}

// guess phase, MPP autosome
// [[Rcpp::export(".guess_phase_A")]]
IntegerVector guess_phase_A(const IntegerMatrix& geno, const String& crosstype)
{
    QTLCross* cross = QTLCross::Create(crosstype);
    const int n_gen = cross->ngen(false);
    const int n_alleles = cross->nalleles();
    delete cross;

    const int n_pos = geno.rows();
    const int n_ind = geno.cols();
    const int matsize = n_pos*2;

    IntegerVector result(n_ind*n_pos*2); // this will be 2 x n_pos x n_ind

    for(int ind=0; ind < n_ind; ind++) {
        IntegerVector g1(n_pos), g2(n_pos);
        for(int pos=0; pos < n_pos; pos++) {

            IntegerVector this_g = mpp_decode_geno(geno(pos,ind), n_alleles, false);

            g1[pos] = this_g[0];
            g2[pos] = this_g[1];
        }

        IntegerVector phased_geno = phase_geno(g1, g2);
        for(int i=0; i<matsize; i++)
            result[ind*matsize+i] = phased_geno[i];
    }

    result.attr("dim") = Dimension(2, n_pos, n_ind);
    return result;
}


// guess phase, MPP X chr
// [[Rcpp::export(".guess_phase_X")]]
IntegerVector guess_phase_X(const IntegerMatrix& geno, const String& crosstype,
                            const LogicalVector& is_female) // pos x ind
{
    QTLCross* cross = QTLCross::Create(crosstype);
    const int n_gen_A = cross->ngen(false);
    const int n_alleles = cross->nalleles();
    delete cross;

    const int n_pos = geno.rows();
    const int n_ind = geno.cols();
    const int matsize = n_pos*2;

    IntegerVector result(n_ind*n_pos*2); // this will be 2 x n_pos x n_ind

    for(int ind=0; ind < n_ind; ind++) {
        if(is_female[ind]) { // female
            IntegerVector g1(n_pos), g2(n_pos);
            for(int pos=0; pos < n_pos; pos++) {

                IntegerVector this_g = mpp_decode_geno(geno(pos,ind), n_alleles, false);

                g1[pos] = this_g[0];
                g2[pos] = this_g[1];
            }

            IntegerVector phased_geno = phase_geno(g1, g2);
            for(int i=0; i<matsize; i++)
                result[ind*matsize+i] = phased_geno[i];
        }
        else { // male
            for(int pos=0, offset=ind*matsize; pos < n_pos; pos++, offset += 2) {
                int g = geno(pos,ind);
                if(IntegerVector::is_na(g)) result[offset] = NA_INTEGER;
                else result[offset] = geno(pos,ind)-n_gen_A;
                result[offset+1] = NA_INTEGER;
            }
        }
    }

    result.attr("dim") = Dimension(2, n_pos, n_ind);
    return result;
}



// guess phase of genotypes for one individual along one chromosome
IntegerVector phase_geno(IntegerVector g1, IntegerVector g2)
{
    const int n_pos = g1.size();
    if(n_pos != g2.size())
        throw std::invalid_argument("length(g1) != length(g2)");

    IntegerVector result(2*n_pos);

    int cur1 = NA_INTEGER, cur2 = NA_INTEGER;
    for(int pos=0, offset=0; pos < n_pos; pos++, offset+=2) {
        if(IntegerVector::is_na(g1[pos]) ||
           IntegerVector::is_na(g2[pos])) {
            result[offset] = result[offset+1] = NA_INTEGER;
        }
        else if(g1[pos]==g2[pos]) { // homozygous so no need to guess
            result[offset] = result[offset+1] = cur1 = cur2 = g1[pos];
        }
        else {
            if(IntegerVector::is_na(cur1) ||
               IntegerVector::is_na(cur2)) { // not yet determined, so randomize
                if(R::runif(0.0, 1.0) < 0.5) {
                    result[offset] = cur1 = g1[pos];
                    result[offset+1] = cur2 = g2[pos];
                }
                else {
                    result[offset] = cur1 = g2[pos];
                    result[offset+1] = cur2 = g1[pos];
                }
            }
            else {
                if(cur1 == g1[pos] || cur2 == g2[pos]) { // assume no recombination
                    result[offset] = cur1 = g1[pos];
                    result[offset+1] = cur2 = g2[pos];
                }
                else if(cur2 == g1[pos] || cur1 == g2[pos]) { // assume no recombination
                    result[offset] = cur1 = g2[pos];
                    result[offset+1] = cur2 = g1[pos];
                }
                else { // randomize
                    if(R::runif(0.0, 1.0) < 0.5) {
                        result[offset] = cur1 = g1[pos];
                        result[offset+1] = cur2 = g2[pos];
                    }
                    else {
                        result[offset] = cur1 = g2[pos];
                        result[offset+1] = cur2 = g1[pos];
                    }
                }
            }
        }

    }

    return result;
}
