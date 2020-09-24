#ifndef _STDIO_H
    #include <stdio.h>
#endif

#ifndef _STDLIB_H
    #include <stdlib.h>
#endif

#ifndef _STRING_H
    #include <string.h>
#endif

#ifndef _STRUCT_H
    #include "struct.h"
#endif

#ifndef _HASH_H
    #include "hash.h"
#endif


#ifndef _MISC_H
    #define _MISC_H 1

    char* parse_header(char *header, long n_delimiter) ;

    long score(const char *seq, const struct matrix_rescaled *m, int cut_off_rescaled) ;

    char* get_compl_seq(char *seq, int len, int *status) ;

    long rescale_value(long value, int nfactor, const int *rescale_factors, char op) ;

#endif // _MISC_H
