// definition of F2 (intercross), for hmm-related stuff
#ifndef CROSS_F2_H
#define CROSS_F2_H

#include <vector>
#include <math.h>
#include "cross.h"

using namespace std;

class F2 : public Cross
{
 public:
    F2(){ type="f2"; };
    virtual ~F2(){};
    
    virtual double init(int true_gen, vector<int> cross_scheme)
    {
        if(true_gen==2) return(-log(0.5)); /* ln(0.5) */
        else return(-2.0*log(0.5)); /* ln(0.25) */
    }

    double emit(int obs_gen, int true_gen, double error_prob, vector<int> cross_scheme)
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


    double step(int gen1, int gen2, double rf, vector<int> cross_scheme)
    {
        switch(gen1) {
        case 1:
            switch(gen2) {
            case 1: return(2.0*log(1.0-rf));
            case 2: return(log(0.5)+log(1.0-rf)+log(rf));
            case 3: return(2.0*log(rf));
            }
        case 2:
            switch(gen2) {
            case 1: case 3: return(log(rf)+log(1.0-rf));
            case 2: return(log((1.0-rf)*(1.0-rf)+rf*rf));
            }
        case 3:
            switch(gen2) {
            case 1: return(2.0*log(rf));
            case 2: return(log(0.5)+log(1.0-rf)+log(rf));
            case 3: return(2.0*log(1.0-rf));
            }
        }
        return(log(-1.0)); /* shouldn't get here */
    }
    
};

#endif // CROSS_F2_H
