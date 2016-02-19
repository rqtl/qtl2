// haploid QTLCross class (for HMM)

#ifndef CROSS_HAPLOID_H
#define CROSS_HAPLOID_H

#include "cross.h"
#include "r_message.h"

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
            r_message("X chr ignored for haploids.");
            return false;
        }

        return true; // most crosses can handle the X chr
    }

    const std::vector<std::string> geno_names(const std::vector<std::string> alleles,
                                              const bool is_x_chr);
};

#endif // CROSS_HAPLOID_H
