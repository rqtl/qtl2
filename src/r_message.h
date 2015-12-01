// print R messages or warnings from C++
#ifndef R_MESSAGE_H
#define R_MESSAGE_H

#include <Rcpp.h>
#include <R_ext/Error.h>

void r_message(std::string text);
void r_warning(std::string text);

#undef NDEBUG // doesn't seem like it should be defined, but it is(?)

// Following based on code from Luke Miratrix, http://bit.ly/rcpp_assert
#ifdef NDEBUG
#define r_assert(EX)
#else
#define r_assert(EX) (void)((EX) || (__r_assert (#EX, __FILE__, __LINE__),0))
#endif

// r_enforce is like r_assert but always works
#define r_enforce(EX) (void)((EX) || (__r_assert (#EX, __FILE__, __LINE__),0))

#endif // R_MESSAGE_H
