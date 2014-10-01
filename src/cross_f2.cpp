#include <math.h>
#include "cross.h"

double F2::init(int true_gen, 
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    if(true_gen==2) return(-log(0.5)); /* ln(0.5) */
    else return(-2.0*log(0.5)); /* ln(0.25) */
}

double F2::emit(int obs_gen, int true_gen, double error_prob, 
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    switch(obs_gen) {
    case 0: return(0.0);
    case 1: case 2: case 3:
        if(obs_gen==true_gen) return(log(1.0-error_prob));
        else return(log(error_prob)-log(0.5));
    case 4: /* AA or AB (not BB) */
        if(true_gen != 3) return(log(1.0-error_prob/2.0));
        else return(log(error_prob));
    case 5: /* AB or BB (not AA) */
        if(true_gen != 1) return(log(1.0-error_prob/2.0));
        else return(log(error_prob));
    }
    return(0.0); /* shouldn't get here */
}


double F2::step(int gen_left, int gen_right, double rf, 
                bool is_X_chr, bool is_female,
                vector<int> cross_info)
{
    switch(gen_left) {
    case 1:
        switch(gen_right) {
        case 1: return(2.0*log(1.0-rf));
        case 2: return(log(0.5)+log(1.0-rf)+log(rf));
        case 3: return(2.0*log(rf));
        }
    case 2:
        switch(gen_right) {
        case 1: case 3: return(log(rf)+log(1.0-rf));
        case 2: return(log((1.0-rf)*(1.0-rf)+rf*rf));
        }
    case 3:
        switch(gen_right) {
        case 1: return(2.0*log(rf));
        case 2: return(log(0.5)+log(1.0-rf)+log(rf));
        case 3: return(2.0*log(1.0-rf));
        }
    }
    return(log(-1.0)); /* shouldn't get here */
}
