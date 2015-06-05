// functions to check QTL cross data/information
#ifndef CHECK_CROSS_H
#define CHECK_CROSS_H

// check if a cross type is supported
bool crosstype_supported(const String& crosstype);

// count inconsistencies in marker data
IntegerVector count_invalid_genotypes(const String& crosstype,
                                      const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                                      const bool& is_X_chr,
                                      const LogicalVector& is_female,
                                      const IntegerMatrix& cross_info); // columns are individuals

// check cross info
bool check_crossinfo(const String& crosstype,
                     const IntegerMatrix& cross_info,
                     const bool any_x_chr);

// check sex
bool check_is_female_vector(const String& crosstype,
                            const LogicalVector& is_female,
                            const bool any_x_chr);

// check if X chr can be handled
bool check_handle_x_chr(const String& crosstype,
                        const bool any_x_chr);

#endif // CHECK_CROSS_H
