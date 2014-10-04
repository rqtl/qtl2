#include "cross.h"

QTLCross* QTLCross::Create(String type)
{
    if(type=="f2") return new F2();
    if(type=="f2pk") return new F2PK();
    if(type=="bc") return new BC();
    if(type=="risib") return new RIsib();
    if(type=="riself") return new RIself();

    throw std::range_error("cross type not yet supported."); 
    return NULL;
}
