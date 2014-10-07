// qtlcross "factory"
//
// to add a new cross type:
//     - create files similar to cross_f2.h and cross_f2.cpp
//     - add include line below
//     - add if statement within QTLCross::Create function below
//
// to create a QTLCross instance using a string with cross type:
//     QTLCross* cross = QTLCross::Create("f2");
// then refer to functions like cross->init()

#include "cross.h"
#include "cross_bc.h"
#include "cross_f2.h"
#include "cross_f2pk.h"
#include "cross_risib.h"
#include "cross_riself.h"
#include "cross_dh.h"
#include "cross_haploid.h"

QTLCross* QTLCross::Create(String type)
{
    if(type=="bc") return new BC();
    if(type=="f2") return new F2();
    if(type=="f2pk") return new F2PK();
    if(type=="risib") return new RISIB();
    if(type=="riself") return new RISELF();
    if(type=="dh") return new DH();
    if(type=="haploid") return new HAPLOID();

    throw std::range_error("cross type not yet supported."); 
    return NULL;
}
