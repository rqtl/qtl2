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
    Rf_warning("%s", text_c);
}

// Following based on code from Luke Miratrix, http://bit.ly/rcpp_assert
void __r_assert (const char *msg, const char *file, int line) {
    char buffer[100];

    snprintf(buffer, 100, "Assert failure: %s at %s line %d\n", msg, file, line);

    Rcpp::stop( buffer );
}
