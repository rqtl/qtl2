// general cross class + cross factory
//
// to add a new cross type:
//     - create a file similar to cross_f2.h
//     - add include line below
//     - add if statement within Cross::Create function below
//
// to create a Cross instance using a string with cross type:
//     Cross* cross1 = Cross::Create("f2");
// or  Cross  cross1 = *(Cross::Create("f2"));

#ifndef CROSS_H
#define CROSS_H

#include <iostream>
#include <string>
#include <vector>

using namespace std;

class Cross
{
public:
    int n_gen;
    string type;

    virtual double init(int true_gen, vector<int>cross_info) {
        return 0.0;
    }

    virtual double emit(int obs_gen, int true_gen, double error_prob,
                        vector<int>cross_scheme) {
        return 0.0;
    }

    virtual double step(int gen1, int gen2, double rf,
                        vector<int>cross_scheme) {
        return 0.0;
    }

    static Cross* Create(string type);
};

#include "cross_f2.h"
#include "cross_bc.h"

Cross* Cross::Create(string type)
{
    if(type=="f2") return new F2();
    if(type=="bc") return new BC();

    return NULL;
}

#endif // CROSS.H
