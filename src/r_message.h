// print R messages or warnings from C++
#ifndef R_MESSAGE_H
#define R_MESSAGE_H

<<<<<<< HEAD
#include <string>
=======
#include <Rcpp.h>
#include <R_ext/Error.h>
>>>>>>> qtl2scan/master

#define RQTL2_NODEBUG 1 // ignore debugging code

void r_message(std::string text);
void r_warning(std::string text);

<<<<<<< HEAD
=======
// Following based on code from Luke Miratrix, http://bit.ly/rcpp_assert
#ifdef RQTL2_NODEBUG
#define r_assert(EX)
#else
#define r_assert(EX) (void)((EX) || (__r_assert (#EX, __FILE__, __LINE__),0))
#endif

// r_enforce is like r_assert but always works
#define r_enforce(EX) (void)((EX) || (__r_assert (#EX, __FILE__, __LINE__),0))
>>>>>>> qtl2scan/master

#endif // R_MESSAGE_H
