// print R messages or warnings from C++
#include "r_message.h"
#include <Rcpp.h>
#include <R_ext/Error.h>

void r_message(std::string text)
{
    Rcpp::Function msg("message");
    msg(text);
}

void r_warning(std::string text)
{
    const char *text_c = text.c_str();
    Rf_warning(text_c);
}
