// functions to check QTL cross data/information

// check if a cross type is supported
bool crosstype_supported(const String& crosstype);

// count inconsistencies in marker data
IntegerVector count_invalid_genotypes(const String& crosstype,
                                      const IntegerMatrix& genotypes, // columns are individuals, rows are markers
                                      const bool& is_X_chr,
                                      const LogicalVector& is_female,
                                      const IntegerMatrix& cross_info); // columns are individuals
