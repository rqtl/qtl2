// doubled haploid QTLCross class (for HMM)

#ifndef CROSS_DH_H
#define CROSS_DH_H

#include "cross.h"
#include "r_message.h"

class DH : public QTLCross
{
 public:
    DH(){
        crosstype = "dh";
        phase_known_crosstype = "dh";
    };

    ~DH(){};

    // check whether X chr can be handled
    const bool check_handle_x_chr(const bool any_x_chr)
    {
        if(any_x_chr) {
            r_message("X chr ignored for doubled haploids.");
            return false;
        }

        return true; // most crosses can handle the X chr
    }

    // the rest of the functions match the defaults,
    // so there's no cross_dh.cpp file
};

#endif // CROSS_DH_H
