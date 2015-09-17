// Utility functions for Diversity Outcross transition probabilities

#include <math.h>
#include <Rcpp.h>
using namespace Rcpp;
#include "cross.h"
#include "cross_do.h"
#include "cross_do_util.h"

//////////////////////////////////////////////////////////////////////
// 1. Stuff related to partially inbred 4-way RIL
//////////////////////////////////////////////////////////////////////

/**********************************************************************
 * probability of AA haplotype on autosome at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_auto_hapAA(double r, int k)
{
    double result, s, rsq;

    rsq = r*r;
    s = sqrt(4.0*rsq-12.0*r+5.0);

    if(r==0.5) {
        if(k==1) result = 1.0/8.0;
        else result = 1.0/16.0;
    }
    else {
        result = (0.25)*(1.0/(1.0+6.0*r) +
                         (6.0*rsq-7.0*r+3.0*r*s)/((1.0+6.0*r)*s)*
                         pow((1.0 - 2.0*r - s)/4.0, (double)k) -
                         (6.0*rsq-7.0*r-3.0*r*s)/((1.0+6.0*r)*s)*
                         pow((1.0 - 2.0*r + s)/4.0, (double)k));
    }

    return(result);
}


/**********************************************************************
 * probability of AA haplotype on X chr in female at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_femX_hapAA(double r, int k)
{
    double result, t, rsq, rcube;

    rsq = r*r;
    rcube = r*rsq;
    t = sqrt(rsq-10.0*r+5.0);

    result = (0.5)*(2.0/(12.0*r+3.0) + 1.0/(3.0*r+3.0)*pow(-0.5,(double)k) -
                    (4.0*rcube - t*(4.0*rsq+3.0*r)+3.0*rsq-5.0*r)/((8.0*rsq+10.0*r+2.0)*t)*
                    pow((1.0 - r + t)/4.0, (double)k) +
                    (4.0*rcube + t*(4.0*rsq+3.0*r)+3.0*rsq-5.0*r)/((8.0*rsq+10.0*r+2.0)*t)*
                    pow((1.0 - r - t)/4.0,(double)k));

    return(result);
}


/**********************************************************************
 * probability of CC haplotype on X chr in female at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_femX_hapCC(double r, int k)
{
    double result, t, rsq;

    rsq = r*r;
    t = sqrt(rsq-10.0*r+5.0);

    result = 1.0/(12.0*r+3.0) - 1.0/(3.0*r+3.0)*pow(-0.5,(double)k) +
        (9.0*rsq +5.0*r + r*t)/((8.0*rsq+10.0*r+2.0)*t)*pow((1.0 - r + t)/4.0,(double)k) -
        (9.0*rsq +5.0*r - r*t)/((8.0*rsq+10.0*r+2.0)*t)*pow((1.0 - r - t)/4.0,(double)k);

    return(result);
}


/**********************************************************************
 * probability of AA haplotype on X chr in male at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_malX_hapAA(double r, int k)
{
    double result, t, rsq, rcube, rfourth;

    rsq = r*r;
    rcube = r*rsq;
    rfourth = rsq*rsq;
    t = sqrt(rsq-10.0*r+5.0);

    result = 1.0/(12.0*r+3.0) - 1.0/(3.0*r+3.0)*pow(-0.5, (double)k) +
        (rcube - t*(8.0*rcube+rsq-3.0*r)-10.0*rsq+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5.0)/2 *
        pow((1.0 - r + t)/4.0,(double)k) +
        +(rcube + t*(8.0*rcube+rsq-3.0*r)-10.0*rsq+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5.0)/2 *
        pow((1.0 - r - t)/4.0,(double)k);

    return(result);
}


/**********************************************************************
 * probability of CC haplotype on X chr in male at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_malX_hapCC(double r, int k)
{
    double result, t, rsq, rcube, rfourth;

    rsq = r*r;
    rcube = r*rsq;
    rfourth = rsq*rsq;
    t = sqrt(rsq-10.0*r+5.0);

    result = 1.0/(12.0*r+3.0) + 2.0/(3.0*r+3.0)*pow(-0.5, (double)k) +
        (2.0*rfourth + t*(2.0*rcube-rsq+r)-19.0*rcube+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5) *
        pow((1.0 - r + t)/4,(double)k) +
        (2.0*rfourth - t*(2.0*rcube-rsq+r)-19.0*rcube+5.0*r)/(4.0*rfourth-35.0*rcube-29.0*rsq+15.0*r+5) *
        pow((1.0 - r - t)/4,(double)k);

    return(result);
}

//////////////////////////////////////////////////////////////////////
// 2. Probability of recombinant haplotypes in DO
//////////////////////////////////////////////////////////////////////

/**********************************************************************
 * probability of recombinant haplotype on autosome at generation s
 * of the diversity outcross population, where at generation 1 the
 * mice are the result of random crosses among pre-CC mice with precc_alpha[k]
 * denoting the proportion of those mice that are at generation k
 *
 * From eqn 6 in Broman (2011) Haplotype probabilities in advanced
 * intercross populations.  Technical report #223, Department of
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_auto(double r, int s, IntegerVector precc_gen, NumericVector precc_alpha)
{
    double hapAA;
    int n_precc = precc_gen.size();
    #ifndef NDEBUG
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
    #endif

    // calculate probability of AA haplotype at generation s=1
    hapAA = 0.0;
    for(int i=0; i<n_precc; i++)
        hapAA += (precc_alpha[i] * ri4way_auto_hapAA(r, precc_gen[i]+1) * (1.0-r)/2.0);

    // later generations (using eqn (6))
    if(s > 1)
        hapAA = 1.0/64.0 + pow(1.0-r, (double)(s-1)) * (hapAA - 1.0/64.0);

    // probability of recombinant, assuming random order of initial crosses
    return( 1.0 - 8.0*hapAA );
}

/**********************************************************************
 * probability of recombinant haplotype on female X chr at generation s
 * of the diversity outcross population, where at generation 1 the
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 *
 * From eqn 7 in Broman (2011) Haplotype probabilities in advanced
 * intercross populations.  Technical report #223, Department of
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_femX_s1(double r, IntegerVector precc_gen, NumericVector precc_alpha)
{
    double result;
    int n_precc = precc_gen.size();
    #ifndef NDEBUG
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
    #endif

    // calculate probability of AA haplotype at generation s=1
    result = 0.0;
    for(int i=0; i<n_precc; i++)
        result += (precc_alpha[i] * (ri4way_femX_hapAA(r, precc_gen[i]+1) * (2.0-r) +
                                     ri4way_femX_hapCC(r, precc_gen[i]+1) * (1.0-r)));

    return(result / 8.0);
}

/**********************************************************************
 * probability of recombinant haplotype on female X chr at generation s
 * of the diversity outcross population, where at generation 1 the
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 *
 * From eqn 7 in Broman (2011) Haplotype probabilities in advanced
 * intercross populations.  Technical report #223, Department of
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_malX_s1(double r, IntegerVector precc_gen, NumericVector precc_alpha)
{
    double result;
    int n_precc = precc_gen.size();
    #ifndef NDEBUG
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
    #endif

    // calculate probability of AA haplotype at generation s=1
    result = 0.0;
    for(int i=0; i<n_precc; i++)
        result += (precc_alpha[i] * (ri4way_malX_hapAA(r, precc_gen[i]+1) * (2.0-r) +
                                     ri4way_malX_hapCC(r, precc_gen[i]+1) * (1.0-r)));

    return(result / 8.0);
}



/**********************************************************************
 * probability of recombinant haplotype on female X chr at generation s
 * of the diversity outcross population, where at generation s the
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 *
 * From eqn 9 in Broman (2011) Haplotype probabilities in advanced
 * intercross populations.  Technical report #223, Department of
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_femX(double r, int s, IntegerVector precc_gen, NumericVector precc_alpha)
{
    double result;
    #ifndef NDEBUG
    int n_precc = precc_gen.size();
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
    #endif

    if(s==1)
        result = DOrec_femX_s1(r, precc_gen, precc_alpha);
    else {
        double z = sqrt((1.0-r)*(9.0-r));
        double ws = pow((1.0-r+z)/4.0, (double)(s-1));
        double ys = pow((1.0-r-z)/4.0, (double)(s-1));

        // calculate probability of AA haplotype at generation s=1
        double f1 = DOrec_femX_s1(r, precc_gen, precc_alpha);
        double m1 = DOrec_malX_s1(r, precc_gen, precc_alpha);

        result = (2.0 + (-64.0*f1*(1.0-r) - 128.0*m1 + 3.0 - r)/z * (ys - ws) -
                  (1.0 - 64*f1)*(ws+ys))/128.0;
    }

    return( 1.0 - 8.0*result );
}

/**********************************************************************
 * probability of recombinant haplotype on female X chr at generation s
 * of the diversity outcross population, where at generation s the
 * mice are the result of random crosses among pre-CC mice with p[k]
 * denoting the proportion of those mice that are at generation k
 *
 * From eqn 9 in Broman (2011) Haplotype probabilities in advanced
 * intercross populations.  Technical report #223, Department of
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_malX(double r, int s, IntegerVector precc_gen, NumericVector precc_alpha)
{
    double result;
#ifndef NDEBUG
    int n_precc = precc_gen.size();
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
#endif

    if(s==1)
        result = DOrec_malX_s1(r, precc_gen, precc_alpha);
    else {
        double z = sqrt((1.0-r)*(9.0-r));
        double ws = pow((1.0-r+z)/4.0, (double)(s-1));
        double ys = pow((1.0-r-z)/4.0, (double)(s-1));

        // calculate probability of AA haplotype at generation s=1
        double f1 = DOrec_femX_s1(r, precc_gen, precc_alpha);
        double m1 = DOrec_malX_s1(r, precc_gen, precc_alpha);

        result = (2.0 + (64.0*m1 - 256.0*f1 + 3.0)*(1.0 - r)/z * (ys - ws) -
                  (1.0 - 64*m1)*(ws+ys))/128.0;
    }

    return( 1.0 - 8.0*result );
}



//////////////////////////////////////////////////////////////////////
// 3. transition probabilties for diversity outcross
//////////////////////////////////////////////////////////////////////

/**********************************************************************
 * transition probability for DO, autosome
 *
 * left = genotype at left locus
 * right = genotype at right locus
 * r = recombination fraction
 * s = generation of DO
 *
 * precc_alpha = proportion of preCC progenitors at generation precc_gen
 *
 * This calculates log Pr(right | left) for phase-unknown case
 *
 **********************************************************************/
double DOstep_auto(int left, int right, double r, int s,
                   IntegerVector precc_gen, NumericVector precc_alpha)
{
    double recprob;
    #ifndef NDEBUG
    int n_precc = precc_gen.size();
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
    #endif

    // pull out alleles for left and right loci
    DO* do_obj = new DO(); // need to create class object to use decode_geno
    IntegerVector leftv = do_obj->decode_geno(left);
    IntegerVector rightv = do_obj->decode_geno(right);
    int left1 = leftv[0];
    int left2 = leftv[1];
    int right1 = rightv[0];
    int right2 = rightv[1];

    // probability of recombinant haplotype
    recprob = DOrec_auto(r, s, precc_gen, precc_alpha);

    if(left1 == left2) {
        if(right1 == right2) {
            if(left1 == right1) { // AA -> AA
                return( 2.0*log(1.0 - recprob) );
            }
            else { // AA -> BB
                return( 2.0*log(recprob) - log(49.0) );
            }
        }
        else {
            if(left1 == right1 || left1 == right2) { // AA -> AB
                return( log(2.0) + log(recprob) + log(1.0-recprob) - log(7.0) );
            }
            else { // AA -> BC
                return( 2.0*log(recprob) + log(2.0) - log(49.0) );
            }
        }
    }
    else { // AB
        if(right1 == right2) {
            if(left1 == right1 || left2 == right1) { // AB -> AA
                return( log(recprob) + log(1.0 - recprob) - log(7.0) );
            }
            else { // AB -> CC
                return( 2.0*log(recprob) - log(49.0) );
            }
        }
        else {
            if((left1==right1 && left2==right2) ||
               (left1==right2 && left2==right1)) { // AB -> AB
                return( log(recprob*recprob/49.0 + (1-recprob)*(1-recprob)) );
            }
            else if(left1==right1 || left1==right2 ||
                    left2==right1 || left2==right2) { // AB -> AC
                return( log(recprob*(1.0-recprob)/7.0 + recprob*recprob/49.0) );
            }
            else { // AB -> CD
                return( 2.0*log(recprob) + log(2.0) - log(49.0) );
            }
        }
    }

}

/**********************************************************************
 * transition probability for DO, female X chr
 **********************************************************************/
double DOstep_femX(int left, int right, double r, int s,
                   IntegerVector precc_gen, NumericVector precc_alpha)
{
    double recprob;
    #ifndef NDEBUG
    int n_precc = precc_gen.size();
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
    #endif

    // pull out alleles for left and right loci
    DO* do_obj = new DO(); // need to create class object to use decode_geno
    IntegerVector leftv = do_obj->decode_geno(left);
    IntegerVector rightv = do_obj->decode_geno(right);
    int left1 = leftv[0];
    int left2 = leftv[1];
    int right1 = rightv[0];
    int right2 = rightv[1];

    // probability of recombinant haplotype
    recprob = DOrec_femX(r, s, precc_gen, precc_alpha);

    if(left1 == left2) {
        if(right1 == right2) {
            if(left1 == right1) { // AA -> AA
                return( 2.0*log(1.0 - recprob) );
            }
            else { // AA -> BB
                return( 2.0*log(recprob) - log(49.0) );
            }
        }
        else {
            if(left1 == right1 || left1 == right2) { // AA -> AB
                return( log(2.0) + log(recprob) + log(1.0-recprob) - log(7.0) );
            }
            else { // AA -> BC
                return( 2.0*log(recprob) + log(2.0) - log(49.0) );
            }
        }
    }
    else { // AB
        if(right1 == right2) {
            if(left1 == right1 || left2 == right1) { // AB -> AA
                return( log(recprob) + log(1.0 - recprob) - log(7.0) );
            }
            else { // AB -> CC
                return( 2.0*log(recprob) - log(49.0) );
            }
        }
        else {
            if((left1==right1 && left2==right2) ||
               (left1==right2 && left2==right1)) { // AB -> AB
                return( log(recprob*recprob/49.0 + (1-recprob)*(1-recprob)) );
            }
            else if(left1==right1 || left1==right2 ||
                    left2==right1 || left2==right2) { // AB -> AC
                return( log(recprob*(1.0-recprob)/7.0 + recprob*recprob/49.0) );
            }
            else { // AB -> CD
                return( 2.0*log(recprob) + log(2.0) - log(49.0) );
            }
        }
    }
}

/**********************************************************************
 * transition probability for DO, male X chr
 **********************************************************************/
double DOstep_malX(int left, int right, double r, int s,
                   IntegerVector precc_gen, NumericVector precc_alpha)
{
    double recprob;
    #ifndef NDEBUG
    int n_precc = precc_gen.size();
    if(n_precc != precc_alpha.size())
        throw std::invalid_argument("precc_gen and precc_alpha should be the same length");
    #endif

    // probability of recombinant haplotype
    recprob = DOrec_malX(r, s, precc_gen, precc_alpha);

    if(left == right) return log(1.0 - recprob);
    return log(recprob) - log(7.0);
}
