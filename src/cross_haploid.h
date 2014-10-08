#ifndef CROSS_HAPLOID_H
#define CROSS_HAPLOID_H

class HAPLOID : public QTLCross
{
 public:
    HAPLOID(){
        type = "haploid";
        phase_known_type = "haploid";
    };

    ~HAPLOID(){};

    // the rest of the functions match the defaults,
    // so there's no cross_dh.cpp file
};

#endif // CROSS_HAPLOID_H
