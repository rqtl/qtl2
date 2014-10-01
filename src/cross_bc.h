// definition of BC (backcross), for hmm-related stuff
#ifndef CROSS_BC_H
#define CROSS_BC_H

#include <vector>
#include <math.h>
#include "cross.h"

class BC : public Cross
{
 public:
    BC(){ type="bc"; };
    ~BC(){};

    double init(int true_gen, vector<int> cross_scheme)
    {
        return(-log(0.5)); /* ln(0.5) */
    }

    double emit(int obs_gen, int true_gen, double error_prob, vector<int> cross_scheme)
    {
        switch(obs_gen) {
        case 0: return(0.0);
        case 1: case 2:
            if(obs_gen==true_gen) return(log(1.0-error_prob));
            else return(log(error_prob));
        }
        return(0.0); /* shouldn't get here */
    }


    double step(int gen1, int gen2, double rf, vector<int> cross_scheme)
    {
        if(gen1==gen2) return(log(1.0-rf));
        else return(log(rf));
        return(log(-1.0)); /* shouldn't get here */
    }
    
};

#endif // CROSS_BC_H
