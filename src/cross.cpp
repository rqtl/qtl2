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
#include "cross_ail.h"
#include "cross_ailpk.h"
#include "cross_do.h"
#include "cross_dopk.h"
#include "cross_dof1.h"
#include "cross_hs.h"
#include "cross_hspk.h"
#include "cross_riself4.h"
#include "cross_riself8.h"
#include "cross_riself16.h"
#include "cross_risib4.h"
#include "cross_risib8.h"
#include "cross_magic19.h"
#include "cross_dh6.h"
#include "cross_ail3.h"
#include "cross_ail3pk.h"
#include "cross_genril.h"
#include "cross_genail.h"
#include <string>

QTLCross* QTLCross::Create(const String& crosstype)
{
    // first, if crosstype has length > 6 and first 6 characters are "genril" or "genail",
    // then turn remaining characters into an integer specifying the number of founders
    std::string crosstype_str = crosstype;
    if(crosstype_str.length() > 6 && crosstype_str.substr(0,6).compare("genril")==0) {
        int n_char = crosstype_str.length();
        std::string str_n_founders = crosstype_str.substr(6, n_char-6);
        int n_founders = std::atoi(str_n_founders.c_str());
        return new GENRIL(n_founders);
    }

    if(crosstype_str.length() > 6 && crosstype_str.substr(0,6).compare("genail")==0) {
        int n_char = crosstype_str.length();
        std::string str_n_founders = crosstype_str.substr(6, n_char-6);
        int n_founders = std::atoi(str_n_founders.c_str());
        return new GENAIL(n_founders);
    }

    if(crosstype=="bc")      return new BC();
    if(crosstype=="f2")      return new F2();
    if(crosstype=="f2pk")    return new F2PK();
    if(crosstype=="risib")   return new RISIB();
    if(crosstype=="riself")  return new RISELF();
    if(crosstype=="dh")      return new DH();
    if(crosstype=="haploid") return new HAPLOID();
    if(crosstype=="ail")     return new AIL();
    if(crosstype=="ailpk")   return new AILPK();
    if(crosstype=="do")      return new DO();
    if(crosstype=="dopk")    return new DOPK();
    if(crosstype=="dof1")    return new DOF1();
    if(crosstype=="hs")      return new HS();
    if(crosstype=="hspk")    return new HSPK();
    if(crosstype=="riself4") return new RISELF4();
    if(crosstype=="riself8") return new RISELF8();
    if(crosstype=="riself16") return new RISELF16();
    if(crosstype=="risib4")  return new RISIB4();
    if(crosstype=="risib8")  return new RISIB8();
    if(crosstype=="magic19") return new MAGIC19();
    if(crosstype=="dh6")     return new DH6();
    if(crosstype=="ail3")    return new AIL3();
    if(crosstype=="ail3pk")  return new AIL3PK();

    throw std::range_error("cross type not yet supported.");
    return NULL;
}
