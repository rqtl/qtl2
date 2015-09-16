// Utility functions for Diversity Outcross transition probabilities

#ifndef CROSS_DO_UTIL_H
#define CROSS_DO_UTIL_H

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
double ri4way_auto_hapAA(double r, int k);

/**********************************************************************
 * probability of AA haplotype on X chr in female at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_femX_hapAA(double r, int k);

/**********************************************************************
 * probability of CC haplotype on X chr in female at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_femX_hapCC(double r, int k);

/**********************************************************************
 * probability of AA haplotype on X chr in male at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_malX_hapAA(double r, int k);

/**********************************************************************
 * probability of CC haplotype on X chr in male at generation G1Fk in
 * 4-way RIL by sibling mating
 *
 * From Table 4 of Broman (2011) Genotype probabilities at intermediate
 * generations in the construction of recombinant inbred lines
 **********************************************************************/
double ri4way_malX_hapCC(double r, int k);

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
double DOrec_auto(double r, int s, IntegerVector precc_gen, NumericVector precc_alpha);

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
double DOrec_femX_s1(double r, IntegerVector precc_gen, NumericVector precc_alpha);

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
double DOrec_malX_s1(double r, IntegerVector precc_gen, NumericVector precc_alpha);

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
double DOrec_femX(double r, int s, IntegerVector precc_gen, NumericVector precc_alpha);

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
double DOrec_malX(double r, int s, IntegerVector precc_gen, NumericVector precc_alpha);

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
 * precc_alpha = proportion of preCC progenitors at generation k
 * n_precc = length of k and precc_alpha
 *
 * This calculates Pr(right | left)
 *
 * code left, right = 0..63
 * with: first allele =  left / 8  (0..7)
 *       second allele = left % 8  (0..7)
 **********************************************************************/
double DOstep_auto(int left, int right, double r, int s,
                   IntegerVector precc_gen, NumericVector precc_alpha);

/**********************************************************************
 * transition probability for DO, female X chr
 *
 * code left, right = 0..63
 * with: first allele =  left / 8 (0..7)
 *       second allele = left % 8 (0..7)
 **********************************************************************/
double DOstep_femX(int left, int right, double r, int s,
                   IntegerVector precc_gen, NumericVector precc_alpha);

/**********************************************************************
 * transition probability for DO, male X chr
 *
 * left, right coded 1, 2, ..., 8 [just one X chr in males]
 **********************************************************************/
double DOstep_malX(int left, int right, double r, int s,
                   IntegerVector precc_gen, NumericVector precc_alpha);

#endif // CROSS_DO_UTIL_H
