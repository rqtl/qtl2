// haploid QTLCross class (for HMM)

#ifndef CROSS_HAPLOID_H
#define CROSS_HAPLOID_H

class HAPLOID : public QTLCross
{
 public:
    HAPLOID(){
        crosstype = "haploid";
        phase_known_crosstype = "haploid";
    };

    ~HAPLOID(){};

    // check whether X chr can be handled
    const bool check_handle_x_chr(const bool any_x_chr)
    {
        if(any_x_chr) {
            //REprintf("X chr ignored for haploids.\n");
            return false;
        }

        return true; // most crosses can handle the X chr
    }

    // the rest of the functions match the defaults,
    // so there's no cross_dh.cpp file
};

#endif // CROSS_HAPLOID_H
