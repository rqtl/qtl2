// print R messages or warnings from C++
#ifndef R_MESSAGE_H
#define R_MESSAGE_H

#include <string>

#define RQTL2_NODEBUG 1 // ignore debugging code

void r_message(std::string text);
void r_warning(std::string text);


#endif // R_MESSAGE_H
