// doubled haploid QTLCross class (for HMM)

#ifndef CROSS_DH_H
#define CROSS_DH_H

class DH : public QTLCross
{
 public:
    DH(){
        crosstype = "dh";
        phase_known_crosstype = "dh";
    };

    ~DH(){};

    // the rest of the functions match the defaults,
    // so there's no cross_dh.cpp file
};

#endif // CROSS_DH_H
