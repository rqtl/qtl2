#include <math.h>
#include <string>
#include <Rcpp.h>
#include "cross.h"

double RIsib::init(int true_gen, bool is_X_chr, bool ignored,
                   vector<int> cross_info)
{
    const bool forward_direction = (cross_info[0] == 0); // AA female x BB male

    if(true_gen != 1 && true_gen != 2)
        Rcpp::exception("genotype not allowed.");

    if(is_X_chr) {
        if(forward_direction) {
            if(true_gen == 1) return log(2.0)-log(3.0);
            if(true_gen == 2) return -log(3.0);
        }
        else {
            if(true_gen == 2) return log(2.0)-log(3.0);
            if(true_gen == 1) return -log(3.0);
        }
    }
    else { // autosome
        return -log(2.0);
    }

    return NA_REAL;
}

double RIsib::emit(int obs_gen, int true_gen, double error_prob, 
                   bool is_X_chr, bool ignored,
                   vector<int> cross_info)
{
    if(obs_gen==0) return 0.0; // missing

    if(obs_gen!=1 && obs_gen!=2)
        Rcpp::exception("Incorrect value for obs_gen");
    if(true_gen!=1 && true_gen!=2)
        Rcpp::exception("Incorrect value for true_gen");

    if(obs_gen == true_gen) return log(1.0 - error_prob);
    else return log(error_prob);
    
    return NA_REAL;
}

double RIsib::step(int gen_left, int gen_right, double rec_frac, 
                   bool is_X_chr, bool ignored,
                   vector<int> cross_info)
{
    if(gen_left!=1 && gen_left!=2)
        Rcpp::exception("Incorrect value for gen_left");
    if(gen_right!=1 && gen_right!=2)
        Rcpp::exception("Incorrect value for gen_right");

    if(is_X_chr)  {
        const bool forward_direction = (cross_info[0] == 0); // AA female x BB male

        if(forward_direction) {
            if(gen_left == 1) {
                if(gen_right == 1) // A -> A
                    return log(1.0 + 2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
                if(gen_right == 2) // A -> B
                    return log(2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }

            if(gen_left == 2) {
                if(gen_right == 2) // B -> B
                    return -log(1.0 + 4.0*rec_frac);
                if(gen_right == 1) // B -> A
                    return log(4.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }
        }
        else {
            if(gen_left == 1) {
                if(gen_right == 1) // A -> A
                    return -log(1.0 + 4.0*rec_frac);
                if(gen_right == 2) // A -> B
                    return log(4.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }
            
            if(gen_left == 2) {
                if(gen_right == 2) // B -> B
                    return log(1.0 + 2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
                if(gen_right == 1) // B -> A
                    return log(2.0*rec_frac) - log(1.0 + 4.0*rec_frac);
            }
        }
    }
    else { // autosome
        const double R = 4.0*rec_frac/(1+6.0*rec_frac);
        
        if(gen_left == gen_right) return log(1.0-R);
        else return log(R);
    }

    return NA_REAL; // can't get here
}

vector<int> RIsib::geno(bool is_X_chr, bool ignored,
                        vector<int> cross_info)
{
    int vals[] = {1,2};
    vector<int> result(vals, vals+2);
    return result;
}
