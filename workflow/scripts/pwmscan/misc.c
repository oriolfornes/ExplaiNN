/*
    Contains miscellaneaous functions for matrix_scan.
*/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "errors.h"
#include "struct.h"
#include "hash.h"



char* parse_header(char *header, long n_delimiter)
{   /*  Parses a FASTA header and extract the sequence identifier. The function
        will look for token begining by '|' and ending by ' ' (space).
        The function will find a token begining after n_delimiter '|' if possible,
        otherwise it will return NULL. The end of the token will be the first
        ' ' (space) or the end of the header if no ' ' (space) is found,

        char *header     : the string containing the FASTA header.
        long n_delimiter : the number of '|' after whose the identifier is.

        Still allocated on exit :
                token

        Return :
                On success the address of token, NULL otherwise.
    */

    if(header == NULL || n_delimiter < 0)
    {   return NULL ; }

    char *token   = NULL ;
    char *token2  = NULL ;
    int n_del = 0 ;
    if(n_delimiter > 0)
    {    /* count number of '|' delimiter */
        for(size_t i=0; i<strlen(header); i++)
        {   if(header[i] == '|')
            {   n_del++ ; }
        }
        /* less '|' present than expected */
        if(n_del < n_delimiter)
        {   return NULL ; }
        /* if format is correct, n_delimiter tokens delimited by '|' */
        int i = 0 ;
        while(i<n_delimiter)
        {   if(i == 0)
            {    token = strtok(header, "|") ; }
            else
            {   token = strtok(NULL, "|") ; }
            if(token == NULL)
            {   return NULL ; }
            i++ ;
        }
        /* token points to the begining of the new token which is terminated by a '\0'. After the '\0' is
           the next token of interest with possibly the sequence identifier.  Reposition token to char just
           after '\0' to be able to search for the ' ' char. */
        token += strlen(token) + 1 ;
    }

    if(n_delimiter == 0)
    {   token = header ; }

    /* try to find ' ' to end sequence ID, otherwise go to end of header */
    token2 = strtok(token, " ") ;
    token = NULL ;
    token = (char*) calloc(strlen(token2)+1, sizeof(char)) ;
    if(token == NULL)
    {   error_calloc(stderr, __TIME__, __LINE__) ;
        return NULL ;
    }
    memcpy((void*) token, (void*) token2, strlen(token2)) ;
    return token ;
}



int score(const char *seq, const struct matrix_rescaled *m, int cut_off_rescaled)
{   /*  Computes the score of the sequence seq, given the rescaled position weight
        matrix (PWM) m.
        If the score falls under the rescaled threshold cut_off, the function
        returns 1 because the score will never be above the threshold again (max
        score per position is 0).

        char *seq                 : the sequence to score
        struct matrix_rescaled *m : the address of the rescaled PWM.
        int cut_off_rescaled      : the rescaled cut-off.

        Still allocated at exit :
                Nothing

        Returns :
                The score ( <= 0 ) if score > cut-off,
                the iteration number before returning if score < cut-off,
                1 if letter was not ATGC.

    */

    int score = 0 ;
    int ncol = m->matrix->ncol ;
    struct matrix *matrix = m->matrix ;

    int current_base ;
    for(int icol=0; icol<ncol; icol++)
    {   /* interrupt score calculation */
        if(score < cut_off_rescaled)
        {   return icol ;
        }

        current_base = seq[icol] ;

        if(current_base == HASH_N)
        {   return 1 ;
        }

        score += (matrix->index[icol][current_base]) ;
    }
    return score ;
}


char* get_compl_seq(char *seq, int len, int *status)
{   /*  Computes the reverse complement sequence of a sequence

.
        char *seq   : the sequence to reverse complement.
        int len     : the number of char from seq to revert.
        int *status : a variable to store exit informations.

        Still allocated on exit :
                seq_compl

        Status :
                 0 on success.
                -1 the reverse complement of a char was not found.
                -2 an allocation error occured.

        Returns :
                On success, the adress of the reversed sequence, NULL otherwise.
                If a base could not be converted, status is turned to -1, if
                a memory allocation error occured, status is turned to -2.

    */
    char *seq_compl = (char*) calloc(len, sizeof(char)+1) ;
    if(seq_compl == NULL)
    {   error_calloc(stderr, __TIME__, __LINE__) ;
        *status = -2 ;
        return (char*) NULL ;
    }
    int i ;
    for(i=0; i<len; i++)
    {   switch(seq[i])
        {   case 'A' :  seq_compl[len-i-1] = 'T' ;
                        break ;
            case 'C' :  seq_compl[len-i-1] = 'G' ;
                        break ;
            case 'G' :  seq_compl[len-i-1] = 'C' ;
                        break ;
            case 'T' :  seq_compl[len-i-1] = 'A' ;
                        break ;
            case 'N' :  seq_compl[len-i-1] = 'N' ;
                        break ;
            case 'W' :  seq_compl[len-i-1] = 'W' ;
                        break ;
            case 'S' :  seq_compl[len-i-1] = 'S' ;
                        break ;
            case 'M' :  seq_compl[len-i-1] = 'K' ;
                        break ;
            case 'K' :  seq_compl[len-i-1] = 'M' ;
                        break ;
            case 'R' :  seq_compl[len-i-1] = 'Y' ;
                        break ;
            case 'Y' :  seq_compl[len-i-1] = 'R' ;
                        break ;
            case 'B' :  seq_compl[len-i-1] = 'V' ;
                        break ;
            case 'V' :  seq_compl[len-i-1] = 'B' ;
                        break ;
            case 'D' :  seq_compl[len-i-1] = 'H' ;
                        break ;
            case 'H' :  seq_compl[len-i-1] = 'D' ;
                        break ;
            default :   *status = -1 ;
                        return (char*) NULL ;
        }
    }
    *status = 0 ;
    return seq_compl ;
}


long rescale_value(long value, int nfactor, const int *rescale_factors, char op)
{   /*  Rescales a value value given an array of factors rescale_factors
        containing nfactor factors.

        long value           : the value to rescale.
        int nfactor          : the number of rescaling factors.
        int *rescale_factors : the address of the rescaling factor array.
        char op              : the mathematical operator to use ('+', '-', '*', '/').
                               '+' is used by default.

        Still allocated on exit :
                    Nothing

        Returns :
                The rescaled value.
    */

    for(int i=0; i<nfactor; i++)
    {   if(op == '-')
        {   value -= rescale_factors[i] ;
            continue ;
        }
        if(op == '/')
        {   value /= rescale_factors[i] ;
            continue ;
        }
        if(op == '*')
        {   value *= rescale_factors[i] ;
            continue ;
        }
        else
        {   value += rescale_factors[i] ;
            continue ;
        }
    }
    return value ;
}


