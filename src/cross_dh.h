#ifndef CROSS_DH_H
#define CROSS_DH_H

class DH : public QTLCross
{
 public:
    DH(){
        type = "dh";
        phase_known_type = "dh";
    };

    ~DH(){};

    // the rest of the functions match the defaults,
    // so there's no cross_dh.cpp file
};

#endif // CROSS_DH_H
