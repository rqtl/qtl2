// Utility functions for Diversity Outcross transition probabilities

#include <math.h>
#include <Rcpp.h>
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
 * mice are the result of random crosses among pre-CC mice with alpha_k[k]
 * denoting the proportion of those mice that are at generation k
 *
 * From eqn 6 in Broman (2011) Haplotype probabilities in advanced
 * intercross populations.  Technical report #223, Department of
 * Biostatistics and Medical Informatics, UW-Madison
 **********************************************************************/
double DOrec_auto(double r, int s)
{
, int n_k, int *k, double *alpha_k)

  double hapAA;
  int i;

  /* calculate probability of AA haplotype at generation s=1 */
  hapAA = 0.0;
  for(i=0; i<n_k; i++)
    hapAA += (alpha_k[i] * ri4way_auto_hapAA(r, k[i]+1) * (1.0-r)/2.0);

  /* later generations (using eqn (6)) */
  if(s > 1)
    hapAA = 1.0/64.0 + pow(1.0-r, (double)(s-1)) * (hapAA - 1.0/64.0);

  /* probability of recombinant, assuming random order of initial crosses */
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
double DOrec_femX_s1(double r, int n_k, int *k, double *alpha_k)
{
  double result;
  int i;

  /* calculate probability of AA haplotype at generation s=1 */
  result = 0.0;
  for(i=0; i<n_k; i++)
    result += (alpha_k[i] * (ri4way_femX_hapAA(r, k[i]+1) * (2.0-r) +
                 ri4way_femX_hapCC(r, k[i]+1) * (1.0-r)));

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
double DOrec_malX_s1(double r, int n_k, int *k, double *alpha_k)
{
  double result;
  int i;

  /* calculate probability of AA haplotype at generation s=1 */
  result = 0.0;
  for(i=0; i<n_k; i++)
    result += (alpha_k[i] * (ri4way_malX_hapAA(r, k[i]+1) * (2.0-r) +
                 ri4way_malX_hapCC(r, k[i]+1) * (1.0-r)));

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
double DOrec_femX(double r, int s, int n_k, int *k, double *alpha_k)
{
  double result, f1, m1;
  double z, ws, ys;

  if(s==1)
    result = DOrec_femX_s1(r, n_k, k, alpha_k);
  else {
    z = sqrt((1.0-r)*(9.0-r));
    ws = pow((1.0-r+z)/4.0, (double)(s-1));
    ys = pow((1.0-r-z)/4.0, (double)(s-1));

    /* calculate probability of AA haplotype at generation s=1 */
    f1 = DOrec_femX_s1(r, n_k, k, alpha_k);
    m1 = DOrec_malX_s1(r, n_k, k, alpha_k);

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
double DOrec_malX(double r, int s, int n_k, int *k, double *alpha_k)
{
  double result, f1, m1;
  double z, ws, ys;

  if(s==1)
    result = DOrec_malX_s1(r, n_k, k, alpha_k);
  else {
    z = sqrt((1.0-r)*(9.0-r));
    ws = pow((1.0-r+z)/4.0, (double)(s-1));
    ys = pow((1.0-r-z)/4.0, (double)(s-1));

    /* calculate probability of AA haplotype at generation s=1 */
    f1 = DOrec_femX_s1(r, n_k, k, alpha_k);
    m1 = DOrec_malX_s1(r, n_k, k, alpha_k);

    result = (2.0 + (64.0*m1 - 256.0*f1 + 3.0)*(1.0 - r)/z * (ys - ws) -
          (1.0 - 64*m1)*(ws+ys))/128.0;
  }

  return( 1.0 - 8.0*result );
}

/* end of DOrec.c */
/**********************************************************************
 * DOstep.c
 *
 * functions to calculate the step (aka transition) probabilities for
 * the diversity outcross population
 *
 * Karl W Broman
 * first written 23 Aug 2011
 * last modified 23 Aug 2011
 *
 **********************************************************************/

#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include "DOstep.h"

/**********************************************************************
 * transition probability for DO, autosome
 *
 * left = genotype at left locus
 * right = genotype at right locus
 * r = recombination fraction
 * s = generation of DO
 *
 * alpha_k = proportion of preCC progenitors at generation k
 * n_k = length of k and alpha_k
 *
 * This calculates Pr(right | left)
 *
 * code left, right = 0..63
 * with: first allele =  left / 8  (0..7)
 *       second allele = left % 8  (0..7)
 **********************************************************************/
double DOstep_auto(int left, int right, double r, int s,
           int n_k, int *k, double *alpha_k)
{
  int left1, left2, right1, right2;
  double recprob;

  /* pull out alleles for left and right loci */
  left1 = left / 8;
  left2 = left % 8;
  right1 = right / 8;
  right2 = right % 8;

  /* probability of recombinant haplotype */
  recprob = DOrec_auto(r, s, n_k, k, alpha_k);

  if(left1 == left2) {
    if(right1 == right2) {
      if(left1 == right1) { /* AA -> AA */
        return( (1.0 - recprob)*(1.0 - recprob) );
      }
      else { /* AA -> BB */
        return( recprob*recprob/49.0 );
      }
    }
    else {
      if(left1 == right1 || left1 == right2) { /* AA -> AB */
        return( 2.0*recprob*(1.0-recprob)/7.0 );
      }
      else { /* AA -> BC */
        return( recprob * recprob * 2.0 / 49.0 );
      }
    }
  }
  else { /* AB */
    if(right1 == right2) {
      if(left1 == right1 || left2 == right1) { /* AB -> AA */
        return( recprob * (1.0 - recprob) / 7.0 );
      }
      else { /* AB -> CC */
        return( recprob * recprob / 49.0 );
      }
    }
    else {
      if((left1==right1 && left2==right2) ||
     (left1==right2 && left2==right1)) { /* AB -> AB */
        return( recprob*recprob/49.0 + (1-recprob)*(1-recprob) );
      }
      else if(left1==right1 || left1==right2 ||
          left2==right1 || left2==right2) { /* AB -> AC */
        return( recprob*(1.0-recprob)/7.0 + recprob*recprob/49.0 );
      }
      else { /* AB -> CD */
        return( recprob*recprob*2.0/49.0 );
      }
    }
  }

} /* DOstep_auto() */

/**********************************************************************
 * transition probability for DO, female X chr
 *
 * code left, right = 0..63
 * with: first allele =  left / 8 (0..7)
 *       second allele = left % 8 (0..7)
 **********************************************************************/
double DOstep_femX(int left, int right, double r, int s,
           int n_k, int *k, double *alpha_k)
{
  int left1, left2, right1, right2;
  double recprob;

  /* pull out alleles for left and right loci */
  left1 = left / 8;
  left2 = left % 8;
  right1 = right / 8;
  right2 = right % 8;

  /* probability of recombinant haplotype */
  recprob = DOrec_femX(r, s, n_k, k, alpha_k);

  if(left1 == left2) {
    if(right1 == right2) {
      if(left1 == right1) { /* AA -> AA */
        return( (1.0 - recprob)*(1.0 - recprob) );
      }
      else { /* AA -> BB */
        return( recprob*recprob/49.0 );
      }
    }
    else {
      if(left1 == right1 || left1 == right2) { /* AA -> AB */
        return( 2.0*recprob*(1.0-recprob)/7.0 );
      }
      else { /* AA -> BC */
        return( recprob * recprob * 2.0 / 49.0 );
      }
    }
  }
  else { /* AB */
    if(right1 == right2) {
      if(left1 == right1 || left2 == right1) { /* AB -> AA */
        return( recprob * (1.0 - recprob) / 7.0 );
      }
      else { /* AB -> CC */
        return( recprob * recprob / 49.0 );
      }
    }
    else {
      if((left1==right1 && left2==right2) ||
     (left1==right2 && left2==right1)) { /* AB -> AB */
        return( recprob*recprob/49.0 + (1-recprob)*(1-recprob) );
      }
      else if(left1==right1 || left1==right2 ||
          left2==right1 || left2==right2) { /* AB -> AC */
        return( recprob*(1.0-recprob)/7.0 + recprob*recprob/49.0 );
      }
      else { /* AB -> CD */
        return( recprob*recprob*2.0/49.0 );
      }
    }
  }
} /* DOstep_femX() */

/**********************************************************************
 * transition probability for DO, male X chr
 *
 * left, right coded 1, 2, ..., 8 [just one X chr in males]
 **********************************************************************/
double DOstep_malX(int left, int right, double r, int s,
           int n_k, int *k, double *alpha_k)
{
  double recprob;

  /* probability of recombinant haplotype */
  recprob = DOrec_malX(r, s, n_k, k, alpha_k);

  if(left == right) return(1.0 - recprob);
  return(recprob/7.0);
} /* DOstep_malX() */

/* end of DOstep.c */
