// print R messages or warnings from C++
#include <Rcpp.h>
#include "R_message.h"

void r_message(std::string text)
{
    Rcpp::Function msg("message");
    msg(text);
}

void r_warning(std::string text)
{
    Rcpp::Function warn("warning");
    warn(text);
}
