/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      Compute the binomial coefficient N choose K = N! / (K! * (N-K)!).
 *      This represents the number of ways to choose K items out of a
 *      pool of N, such that we don't care about order.  E.g., if
 *      choosing 2 letters from the set [A B C D E], [C D] is considered
 *      the same [D C].
 *
 *      This implementation avoids overflow by alternating multiply and
 *      divide operations (rather than try to compute factorials first,
 *      which will fail for relatively small values of N or K).
 *  
 *  Arguments:
 *      n   Number of items to choose from
 *      k   Number of items chosen
 *
 *  Returns:
 *      The number of ways to choose K items from N objects.
 *
 *  Examples:
 *      #include <xtend/math.h>
 *
 *      unsigned long   n = 5, k = 2;
 *
 *      printf("Ways to choose %lu items from %lu = %lu\n",
 *              k, n, xt_n_choose_k(n, k));
 *
 *  See also:
 *      lgamma(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-12-09  Jason Bacon Begin
 ***************************************************************************/
#define _GNU_SOURCE 
#include <stdio.h>
#include <stdint.h>
#include "libxtend.h"

unsigned long   xt_n_choose_k(unsigned long n, unsigned long k)

{
    unsigned long   b, c;
    
    if ( k < 0 || k > n )
	return 0UL;
    if ( k == 0 || k == n )
	return 1UL;
    k = XT_MIN(k, n - k);
    for (b = 1, c = 0; c < k; ++c)
	b = b * (n - c) / (c + 1);
    return b;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <inttypes.h>
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      Instantaneous factorial for n in [0,20] using a lookup table.
 *
 *      Note that 21! is beyond the range of uint64_t, so programs that
 *      need it should either use a multiple precision library or rearrange
 *      the computations to avoid repeated multiplication leading to overflow.
 *      xt_n_choose_k(3) uses the latter approach to avoid computing whole
 *      factorials before finally dividing.
 *  
 *  Arguments:
 *      n   Integer in the range [0,20] (inclusive) for which n! is returned.
 *
 *  Returns:
 *      n! if n is an element of [0,20], 0 otherwise
 *
 *  Examples:
 *      #include <inttypes.h>
 *      #include <xtend/math.h>
 *
 *      printf("20! = %" PRIu64 "\n", factorial(20));
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-12-11  Jason Bacon Begin
 ***************************************************************************/

uint64_t    xt_factorial(unsigned n)

{
    static uint64_t f[] = { 1, 1, 2, 6, 24, 120, 720, 5040, 40320,
			    362880, 3628800, 39916800, 479001600,
			    6227020800, 87178291200, 1307674368000,
			    20922789888000, 355687428096000, 6402373705728000,
			    121645100408832000, 2432902008176640000 };
    
    return n <= 20 ? f[n] : 0;
}

/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      num_digits() computes the number of digits in val, assuming the
 *      given base.
 *  
 *  Arguments:
 *      val:    The number for which digits are to be counted
 *      base:   The number base, between 2 and 36
 *
 *  Returns:
 *      The number of base "base" digits in val, or -1 if base is invalid
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

int     digits(long val, unsigned base)

{
    int     d;
    
    if ( (base < 2) || (base > 36) )
	return -1;
    
    for (d=1; val != 0; ++d)
	val /= base;
    return d;
}
#ifdef __linux__
#define _GNU_SOURCE     // vasprintf()
#endif

#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <stdlib.h>
#include <stdarg.h>

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      The xt_dprintf() function, which takes a file descriptor rather
 *      than a FILE stream pointer, is provided by many systems including
 *      BSDs and Linux, but not by all.  Use of xt_xt_dprintf() from
 *      libxtend will ensure portability of code.
 *  
 *  Arguments:
 *      fd      File descriptor to which items are written
 *      format  printf-style format string
 *      ...     Additional arguments depending on format
 *
 *  Returns:
 *      The number of items written
 *
 *  Examples:
 *      int     fd;
 *
 *      if ( (fd = open(filename, O_WRONLY|O_CREAT)) != -1 )
 *      {
 *          xt_xt_dprintf(fd, "fd = %d\n", fd);
 *          ...
 *          close(fd);
 *      }
 *
 *  See also:
 *      fprintf(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-08-20  Jason Bacon Begin
 ***************************************************************************/

int     xt_dprintf(int fd, const char * restrict format, ...)

{
    char    *buff;
    int     count;
    va_list ap;
    
    va_start(ap, format);
    count = vasprintf(&buff, format, ap);
    write(fd, buff, strlen(buff));
    free(buff);
    va_end(ap);
    
    return count;
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read next delimiter-separated field from stream. The fields may be
 *      ended by any character in the string delims or by a newline ('\\\\n').
 *
 *      If the delimiter ending a field is a space, then subsequence spaces
 *      are discarded, so that multiple space characters serve as a single
 *      delimiter.
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      buff        Character buff into which field is copied
 *      buff_size   Size of the array passed to buff
 *      delims      Array of characters that may serve as delimiters
 *      len         Pointer to a variable which will receive the field length
 *
 *  Returns:
 *      Delimiter ending the field (either a member of delim or newline)
 *
 *  See also:
 *      dsv_read_field_malloc(3), dsv_skip_field(3),
 *      dsv_skip_rest_of_line(3), dsv_line_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-24  Jason Bacon Begin
 ***************************************************************************/

int     dsv_read_field(FILE *stream, char buff[], size_t buff_size,
		       const char *delims, size_t *len)

{
    size_t  c;
    char    *p;
    int     ch, ch2;
    
    for (c = 0, p = buff; (c < buff_size) && 
			  ( strchr(delims, ch = getc(stream)) == NULL) &&
			  (ch != '\n') && (ch != EOF); ++c, ++p )
	*p = ch;
    *p = '\0';
    
    if ( c == buff_size )
    {
	fprintf(stderr, "dsv_read_field(): Buffer overflow reading field.\n");
	fprintf(stderr, "Buffer size = %zu\n", buff_size);
	fputs(buff, stderr);
	// FIXME: Replace this with another sentinal value?
	// Would require all callers to handle both EOF and overflow
	exit(EX_SOFTWARE);
    }
    
    *len = c;
    
    /*
     *  Treat space specially in that multiple spaces are considered a single
     *  separator
     */
    if ( ch == ' ' )
    {
	while ( (ch2 = getc(stream)) == ch )
	    ;
	ungetc(ch2, stream);
    }
    return ch;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read next delimiter-separated field from stream, allocating a
 *      buffer to fit in the fashion of strdup(3). The fields may be
 *      ended by any character in the string delims or by a newline ('\\\\n').
 *
 *      If the delimiter ending a field is a space, then subsequence spaces
 *      are discarded, so that multiple space characters serve as a single
 *      delimiter.
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      buff        Character buffer into which field is copied
 *      buff_size   Size of the array passed to buff
 *      delims      Array of characters that may serve as delimiters
 *      len         Pointer to a variable which will receive the field length
 *
 *  Returns:
 *      Delimiter ending the field (either a member of delim or newline)
 *      or XT_MALLOC_FAILED.
 *
 *  See also:
 *      dsv_read_field(3), dsv_skip_field(3), dsv_skip_rest_of_line(3),
 *      dsv_line_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-24  Jason Bacon Begin
 ***************************************************************************/

int     dsv_read_field_malloc(FILE *stream, char **buff, size_t *buff_size,
		       const char *delims, size_t *len)

{
    size_t  c;
    int     ch, ch2;
    
    if ( *buff_size == 0 )
    {
	*buff_size = 1024;
	*buff = xt_malloc(*buff_size, sizeof(**buff));
	if ( *buff == NULL )
	    return XT_MALLOC_FAILED;
    }
    
    for (c = 0; ( ((ch = getc(stream)) != '\n') && (ch != EOF) &&
		  strchr(delims, ch) == NULL); ++c )
    {
	if ( c == *buff_size - 1 )
	{
	    *buff_size *= 2;
	    *buff = xt_realloc(*buff, *buff_size, sizeof(**buff));
	    if ( *buff == NULL )
		return XT_MALLOC_FAILED;
	}
	(*buff)[c] = ch;
    }
    (*buff)[c] = '\0';
    *len = c;

    /* Trim array */
    if ( *buff_size != c + 1 )
    {
	*buff_size = c + 1;
	*buff = xt_realloc(*buff, *buff_size, sizeof(**buff));
    }

    /*
     *  Treat space specially in that multiple spaces are considered a single
     *  separator
     */
    if ( ch == ' ' )
    {
	while ( (ch2 = getc(stream)) == ch )
	    ;
	ungetc(ch2, stream);
    }
    return ch;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read and discard next delimiter-separated field from stream. The
 *      fields may be ended by any character in the string delims or by a
 *      newline ('\\\\n').
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      delims      Array of characters that may serve as delimiters
 *      len         Length of field discarded
 *
 *  Returns:
 *      Delimiter ending the field (either a member of delim or newline)
 *
 *  See also:
 *      dsv_read_field(3), dsv_read_field_malloc(3),
 *      dsv_skip_rest_of_line(3), dsv_line_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-02-24  Jason Bacon Begin
 ***************************************************************************/

int     dsv_skip_field(FILE *stream, const char *delims, size_t *len)

{
    int     ch;
    
    for (*len = 0; (strchr(delims, ch = getc(stream)) == NULL) &&
	    (ch != '\n') && (ch != EOF); ++*len )
	;
    
    return ch;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read and discard all remaining fields in a line from stream.
 *      I.e., discard everything up to and including the next newline ('\\\\n').
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *
 *  Returns:
 *      Delimiter ending the field (should always be newline ('\\\\n'))
 *
 *  See also:
 *      dsv_read_field(3), dsv_read_field_malloc(3),
 *      dsv_skip_field(3), dsv_line_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

int     dsv_skip_rest_of_line(FILE *stream)

{
    int     ch;
    
    while ( ((ch = getc(stream)) != EOF) && (ch != '\n') )
	;
    return ch;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Read a line of an arbitrary DSV file into a dsv_line_t object.
 *      The dsv_line_t structure contains an array of strings, each
 *      holding one field from the line, and an an array of delimiters,
 *      each holding the character that ended the corresponding field.
 *      Note that each field could potentially end with a different
 *      delimiter, as multiple delimiters can be specified.
 *
 *      This function serves a purpose similar to the split() functions
 *      present in many languages.  However, it does not need to read an
 *      entire line into a character array and then split the array.
 *      Instead, it separates fields as they are read from the input stream.
 *
 *  Arguments:
 *      dsv_line    Pointer to a dsv_line_t structure to hold the fields
 *      stream      FILE stream from which the line is read
 *      delims      Array of acceptable delimiters
 *
 *  Returns:
 *      Actual delimiter of last field (should be newline)
 *
 *  Examples:
 *      dsv_line_t  line;
 *
 *      dsv_line_init(&line);
 *      while ( dsv_line_read(&line, stdin, "\\\\\t") != EOF )
 *      {
 *          dsv_line_write(line, stdout);
 *          dsv_line_free(&line);
 *      }
 *
 *  See also:
 *      dsv_line_init(3), dsv_line_free(3),
 *      dsv_line_read(3), dsv_line_write(3), dsv_line_copy(3),
 *      dsv_read_field(3), dsv_read_field_malloc(3),
 *      dsv_skip_field(3), dsv_skip_rest_of_line(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-30  Jason Bacon Begin
 ***************************************************************************/

int     dsv_line_read(dsv_line_t *dsv_line, FILE *stream, const char *delims)

{
    int     actual_delim;
    char    field[DSV_FIELD_MAX_CHARS + 1];
    size_t  actual_len;
    
    dsv_line->array_size = 32;  // Start small and double each time we run out
    dsv_line->num_fields = 0;
    
    // FIXME: Reuse previously allocated memory?
    if ( (dsv_line->fields = xt_malloc(dsv_line->array_size,
				sizeof(*dsv_line->fields))) == NULL )
    {
	fputs("dsv_line_read(): Could not allocate fields.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    
    if ( (dsv_line->delims = xt_malloc(dsv_line->array_size,
				sizeof(*dsv_line->delims))) == NULL )
    {
	fputs("dsv_line_read(): Could not allocate delims.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    
    // FIXME: Check actual_delim and/or actual_len to detect truncation
    while ( ((actual_delim = dsv_read_field(stream,
		field, DSV_FIELD_MAX_CHARS, delims, &actual_len)) != EOF) )
    {
	if ( (dsv_line->fields[dsv_line->num_fields] = strdup(field)) == NULL )
	{
	    fprintf(stderr, "dsv_line_read(): Could not strdup() field %zu.\n",
		    dsv_line->num_fields - 1);
	    exit(EX_UNAVAILABLE);
	}
	dsv_line->delims[dsv_line->num_fields++] = actual_delim;
	if ( dsv_line->num_fields == dsv_line->array_size )
	{
	    dsv_line->array_size *= 2;
	    if ( (dsv_line->fields = xt_realloc(dsv_line->fields,
		    dsv_line->array_size, sizeof(*dsv_line->fields))) == NULL )
	    {
		fputs("dsv_line_read(): Could not reallocate fields.\n", stderr);
		exit(EX_UNAVAILABLE);
	    }
	    
	    if ( (dsv_line->delims = xt_realloc(dsv_line->delims,
		    dsv_line->array_size, sizeof(*dsv_line->delims))) == NULL )
	    {
		fputs("dsv_line_read(): Could not reallocate delims.\n", stderr);
		exit(EX_UNAVAILABLE);
	    }
	}
	if ( actual_delim == '\n' )
	    break;
    }
    return actual_delim;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Write an arbitrary DSV line to the specified stream.
 *      The dsv_line_t structure contains an array of strings, each
 *      holding one field from the line, and an an array of delimiters,
 *      each holding the character that ended the corresponding field.
 *      Note that each field could potentially end with a different
 *      delimiter, as multiple delimiters can be specified.
 *
 *  Arguments:
 *      dsv_line    Pointer to dsv_line_t structure holding the fields
 *      stream      FILE stream to which fields are printed (e.g. stderr)
 *
 *  Returns:
 *      The number of fields successfully written
 *
 *  Examples:
 *      dsv_line_t  line;
 *
 *      dsv_line_init(&line);
 *      while ( dsv_line_read(&line, stdin, "\t") != EOF )
 *      {
 *          dsv_line_write(line, stdout);
 *          dsv_line_free(&line);
 *      }
 *
 *  See also:
 *      dsv_line_init(3), dsv_line_free(3),
 *      dsv_line_read(3), dsv_line_write(3), dsv_line_copy(3),
 *      dsv_read_field(3), dsv_read_field_malloc(3),
 *      dsv_skip_field(3), dsv_skip_rest_of_line(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

int     dsv_line_write(dsv_line_t *dsv_line, FILE *stream)

{
    int     c, count = 0;
    
    for (c = 0; c < dsv_line->num_fields; ++c)
    {
	if ( fprintf(stream, "%s%c", dsv_line->fields[c], dsv_line->delims[c]) == 2 )
	    ++count;
    }
    return count;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Initialize a dsv_line_t structure.
 *      The dsv_line_t structure contains an array of strings, each
 *      holding one field from the line, and an an array of delimiters,
 *      each holding the character that ended the corresponding field.
 *      Note that each field could potentially end with a different
 *      delimiter, as multiple delimiters can be specified.
 *  
 *  Arguments:
 *      dsv_line    Pointer to a dsv_lint_t object.    
 *
 *  Examples:
 *      dsv_line_t  line;
 *
 *      dsv_line_init(&line);
 *      while ( dsv_line_read(&line, stdin, "\t") != EOF )
 *      {
 *          dsv_line_write(line, stdout);
 *          dsv_line_free(&line);
 *      }
 *
 *  See also:
 *      dsv_line_init(3), dsv_line_free(3),
 *      dsv_line_read(3), dsv_line_write(3), dsv_line_copy(3),
 *      dsv_read_field(3), dsv_read_field_malloc(3),
 *      dsv_skip_field(3), dsv_skip_rest_of_line(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-11  Jason Bacon Begin
 ***************************************************************************/

void    dsv_line_init(dsv_line_t *dsv_line)

{
    dsv_line->array_size = 0;
    dsv_line->num_fields = 0;
    dsv_line->fields = NULL;
    dsv_line->delims = NULL;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Duplicate an arbitrary dsv_line_t object, allocating space for
 *      fields and delimiters as needed.
 *      The dsv_line_t structure contains an array of strings, each
 *      holding one field from the line, and an an array of delimiters,
 *      each holding the character that ended the corresponding field.
 *      Note that each field could potentially end with a different
 *      delimiter, as multiple delimiters can be specified.
 *
 *  Arguments:
 *      src     Pointer to populated dsv_line_t structure to be duplicated
 *      dest    Pointer to empty dsv_lint_t structure to receive copy
 *
 *  Returns:
 *      XT_OK or XT_MALLOC_FAILED
 *      
 *  See also:
 *      dsv_line_init(3), dsv_line_free(3),
 *      dsv_line_read(3), dsv_line_write(3), dsv_line_copy(3),
 *      dsv_read_field(3), dsv_read_field_malloc(3),
 *      dsv_skip_field(3), dsv_skip_rest_of_line(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

int     dsv_line_copy(dsv_line_t *dest, dsv_line_t *src)

{
    size_t  c;
    
    // Prune unused pointers in src
    dest->array_size = dest->num_fields = src->num_fields;
    
    // FIXME: Check malloc() success
    dest->fields = xt_malloc(dest->array_size, sizeof(*dest->fields));
    if ( dest->fields == NULL )
	return XT_MALLOC_FAILED;
    dest->delims = xt_malloc(dest->array_size, sizeof(*dest->delims));
    if ( dest->delims == NULL )
	return XT_MALLOC_FAILED;
    
    for (c = 0; c < src->num_fields; ++c)
    {
	if ( (dest->fields[c] = strdup(src->fields[c])) == NULL )
	    return XT_MALLOC_FAILED;
	dest->delims[c] = src->delims[c];
    }
    return XT_OK;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Free allocated memory for a DSV object.
 *      The dsv_line_t structure contains an array of strings, each
 *      holding one field from the line, and an an array of delimiters,
 *      each holding the character that ended the corresponding field.
 *      Note that each field could potentially end with a different
 *      delimiter, as multiple delimiters can be specified.
 *
 *  Arguments:
 *      dsv_line    Pointer to a populated dsv_line_t structure
 *
 *  Returns:
 *      The number of fields freed.  Fields set to NULL are not freed.
 *
 *  Examples:
 *      dsv_line_t  line;
 *
 *      while ( dsv_line_read(&line, stdin, "\t") != EOF )
 *      {
 *          dsv_line_write(line, stdout);
 *          dsv_line_free(&line);
 *      }
 *
 *  See also:
 *      dsv_line_init(3), dsv_line_free(3),
 *      dsv_line_read(3), dsv_line_write(3), dsv_line_copy(3),
 *      dsv_read_field(3), dsv_read_field_malloc(3),
 *      dsv_skip_field(3), dsv_skip_rest_of_line(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-01  Jason Bacon Begin
 ***************************************************************************/

int     dsv_line_free(dsv_line_t *dsv_line)

{
    int     c, count = 0;
    
    if ( dsv_line->fields != NULL )
    {
	for (c = 0; c < dsv_line->num_fields; ++c)
	    if ( dsv_line->fields[c] != NULL )
	    {
		free(dsv_line->fields[c]);
		++count;
	    }
	if ( dsv_line->fields != NULL )
	    free(dsv_line->fields);
    }
    dsv_line->num_fields = 0;
    return count;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_read_field(stream, buff, buff_size, '\\\\\t', len)
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      buff        Character buff into which field is copied
 *      buff_size   Size of the array passed to buff
 *      len         Pointer to a variable which will receive the field length
 *
 *  See also:
 *      dsv_read_field(3)
 ***************************************************************************/

int     tsv_read_field(FILE *stream, char buff[], size_t buff_size,
		       size_t *len)

{
    return dsv_read_field(stream, buff, buff_size, "\t", len);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_read_field_malloc(stream, *buff, *buff_size, '\\\\\t', len)
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      buff        Character buff into which field is copied
 *      buff_size   Size of the array passed to buff
 *      len         Pointer to a variable which will receive the field length
 *
 *  See also:
 *      dsv_read_field_malloc(3)
 ***************************************************************************/

int     tsv_read_field_malloc(FILE *stream, char **buff, size_t *buff_size,
		       size_t *len)

{
    return dsv_read_field_malloc(stream, buff, buff_size, "\t", len);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_skip_field(stream, '\\\\\t')
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      len         Length of field discarded
 *
 *  See also:
 *      dsv_skip_field(3)
 ***************************************************************************/

int     tsv_skip_field(FILE *stream, size_t *len)

{
    return dsv_skip_field(stream, "\t", len);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_skip_rest_of_line(stream)
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *
 *  See also:
 *      dsv_skip_rest_of_line(3)
 ***************************************************************************/

int     tsv_skip_rest_of_line(FILE *stream)

{
    return dsv_skip_rest_of_line(stream);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_read_field(stream, buff, buff_size, ',', len)
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      buff        Character buff into which field is copied
 *      buff_size   Size of the array passed to buff
 *      len         Pointer to a variable which will receive the field length
 *
 *  See also:
 *      dsv_read_field(3)
 ***************************************************************************/

int     csv_read_field(FILE *stream, char buff[], size_t buff_size,
		       size_t *len)

{
    return dsv_read_field(stream, buff, buff_size, ",", len);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_read_field_malloc(stream, *buff, *buff_size, ',', len)
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      buff        Character buff into which field is copied
 *      buff_size   Size of the array passed to buff
 *      len         Pointer to a variable which will receive the field length
 *
 *  See also:
 *      dsv_read_field_malloc(3)
 ***************************************************************************/

int     csv_read_field_malloc(FILE *stream, char **buff, size_t *buff_size,
		       size_t *len)

{
    return dsv_read_field_malloc(stream, buff, buff_size, ",", len);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_skip_field(stream, ',')
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *      len         Length of field discarded
 *
 *  See also:
 *      dsv_skip_field(3)
 ***************************************************************************/

int     csv_skip_field(FILE *stream, size_t *len)

{
    return dsv_skip_field(stream, ",", len);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lbiolibc
 *
 *  Description:
 *      Equivalent to dsv_skip_rest_of_line(stream)
 *
 *  Arguments:
 *      stream      FILE stream from which field is read
 *
 *  See also:
 *      dsv_skip_rest_of_line(3)
 ***************************************************************************/

int     csv_skip_rest_of_line(FILE *stream)

{
    return dsv_skip_rest_of_line(stream);
}
/***************************************************************************
 *  This file is automatically generated by gen-get-set.  Be sure to keep
 *  track of any manual changes.
 *
 *  These generated functions are not expected to be perfect.  Check and
 *  edit as needed before adding to your code.
 ***************************************************************************/

#include <string.h>
#include <ctype.h>
#include <stdbool.h>        // In case of bool
#include <stdint.h>         // In case of int64_t, etc


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for array_size member in a dsv_line_t structure.
 *      Use this function to set array_size in a dsv_line_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      new_array_size  The new value for array_size
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      size_t          new_array_size;
 *
 *      if ( dsv_line_set_array_size(&dsv_line, new_array_size)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_array_size(dsv_line_t *dsv_line_ptr, size_t new_array_size)

{
    if ( false )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	dsv_line_ptr->array_size = new_array_size;
	return BL_DSV_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for num_fields member in a dsv_line_t structure.
 *      Use this function to set num_fields in a dsv_line_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      num_fields is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      new_num_fields  The new value for num_fields
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      size_t          new_num_fields;
 *
 *      if ( dsv_line_set_num_fields(&dsv_line, new_num_fields)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_num_fields(dsv_line_t *dsv_line_ptr, size_t new_num_fields)

{
    if ( false )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	dsv_line_ptr->num_fields = new_num_fields;
	return BL_DSV_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for fields member in a dsv_line_t structure.
 *      Use this function to set fields in a dsv_line_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      fields is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      new_fields      The new value for fields
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      char **         new_fields;
 *
 *      if ( dsv_line_set_fields(&dsv_line, new_fields)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_fields(dsv_line_t *dsv_line_ptr, char ** new_fields)

{
    if ( new_fields == NULL )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	dsv_line_ptr->fields = new_fields;
	return BL_DSV_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for an array element of fields member in a dsv_line_t
 *      structure. Use this function to set dsv_line_ptr->fields[c]
 *      in a dsv_line_t object from non-member functions.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      c               Subscript to the fields array
 *      new_fields_element The new value for fields[c]
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      size_t          c;
 *      char **         new_fields_element;
 *
 *      if ( dsv_line_set_fields_ae(&dsv_line, c, new_fields_element)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      DSV_LINE_SET_FIELDS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_fields_ae(dsv_line_t *dsv_line_ptr, size_t c, char * new_fields_element)

{
    if ( new_fields_element == NULL )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	dsv_line_ptr->fields[c] = new_fields_element;
	return BL_DSV_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for fields member in a dsv_line_t structure.
 *      Use this function to set fields in a dsv_line_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_fields to dsv_line_ptr->fields.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      new_fields      The new value for fields
 *      array_size      Size of the fields array.
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      char **         new_fields;
 *      size_t          array_size;
 *
 *      if ( dsv_line_set_fields_cpy(&dsv_line, new_fields, array_size)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      DSV_LINE_SET_FIELDS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_fields_cpy(dsv_line_t *dsv_line_ptr, char ** new_fields, size_t array_size)

{
    if ( new_fields == NULL )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    dsv_line_ptr->fields[c] = new_fields[c];
	return BL_DSV_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for delims member in a dsv_line_t structure.
 *      Use this function to set delims in a dsv_line_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      delims is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      new_delims      The new value for delims
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      char *          new_delims;
 *
 *      if ( dsv_line_set_delims(&dsv_line, new_delims)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_delims(dsv_line_t *dsv_line_ptr, char * new_delims)

{
    if ( new_delims == NULL )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	dsv_line_ptr->delims = new_delims;
	return BL_DSV_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for an array element of delims member in a dsv_line_t
 *      structure. Use this function to set dsv_line_ptr->delims[c]
 *      in a dsv_line_t object from non-member functions.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      c               Subscript to the delims array
 *      new_delims_element The new value for delims[c]
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      size_t          c;
 *      char *          new_delims_element;
 *
 *      if ( dsv_line_set_delims_ae(&dsv_line, c, new_delims_element)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      DSV_LINE_SET_DELIMS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_delims_ae(dsv_line_t *dsv_line_ptr, size_t c, char  new_delims_element)

{
    if ( false )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	dsv_line_ptr->delims[c] = new_delims_element;
	return BL_DSV_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/dsv.h>
 *      -lxtend
 *
 *  Description:
 *      Mutator for delims member in a dsv_line_t structure.
 *      Use this function to set delims in a dsv_line_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_delims to dsv_line_ptr->delims.
 *
 *  Arguments:
 *      dsv_line_ptr    Pointer to the structure to set
 *      new_delims      The new value for delims
 *      array_size      Size of the delims array.
 *
 *  Returns:
 *      BL_DSV_DATA_OK if the new value is acceptable and assigned
 *      BL_DSV_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      dsv_line_t      dsv_line;
 *      char *          new_delims;
 *      size_t          array_size;
 *
 *      if ( dsv_line_set_delims_cpy(&dsv_line, new_delims, array_size)
 *              == BL_DSV_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      DSV_LINE_SET_DELIMS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  gen-get-set Auto-generated from dsv.h
 ***************************************************************************/

int     dsv_line_set_delims_cpy(dsv_line_t *dsv_line_ptr, char * new_delims, size_t array_size)

{
    if ( new_delims == NULL )
	return BL_DSV_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(dsv_line_ptr->delims, new_delims, array_size);
	return BL_DSV_DATA_OK;
    }
}
#include <stdlib.h>
#include <unistd.h>
#include <fcntl.h>
#include <sys/stat.h>
#include <sysexits.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      fast_cp() copies a file using low-level I/O with an optimal
 *      buffer size (the least common multiple of block sizes) for both
 *      source and destination filesystems.
 *  
 *  Arguments:
 *      source, dest: File names of source and destination
 *
 *  Returns:
 *      The return value of the last read(3) call.
 *
 *  See also:
 *      cp(1), read(3), write(3), fstat(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

int     xt_fast_cp(const char *source, const char *dest)

{
    int         infile,outfile;
    struct stat infile_stats,outfile_stats;
    long        x;
    long        buff_size, nbytes;
    char        *buff;
    
    /* Open source and destination files for low level transfer */
    if ( (infile = open(source,O_RDONLY)) == -1 )
	return EX_NOINPUT;
    
    if ( (outfile = open(dest,O_WRONLY|O_CREAT|O_TRUNC,0700)) == -1 )
	return EX_CANTCREAT;

    /* Create buffer of optimum size for both files, 64k max
       Optimum buffer size is any multiple of the both sector sizes
       given by struct stat member st_blksize.  COHERENT cc lacks
       st_blksize */
    fstat(infile,&infile_stats);
    fstat(outfile,&outfile_stats);
    x = lcm(infile_stats.st_blksize,outfile_stats.st_blksize);
    buff_size = XT_MIN(x,256*1024);

    if ( (buff = (char *)malloc(buff_size)) == NULL )
    {
	fputs("fast_cp(): malloc() failed.\n", stderr);
	close(infile);
	close(outfile);
	return -1;
    }
    
    /* Copy file */
    while ( (nbytes = read(infile,buff,buff_size)) > 0 )
	write(outfile,buff,nbytes);
    close(infile);
    close(outfile);
    free(buff);
    return nbytes;  /* Return 0 for success, or error code from read() */
}
#ifdef __linux__
#define _GNU_SOURCE     // vasprintf(), must come before stdio.h
#endif

#include <stdio.h>
#include <stdlib.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <sys/wait.h>
#include <string.h>
#include <stdarg.h>

/*
 *  Non-API function for completing stream initialization for ffopen()
 *  and ffdopen()
 */

ffile_t *ff_init_stream(ffile_t *stream)

{
    struct stat st;
    
    // Get optimal block size for the underlying filesystem
    if ( fstat(stream->fd, &st) != 0 )
    {
	free(stream);
	fprintf(stderr, "ffopen(): Could not stat fd %d.\n", stream->fd);
	return NULL;
    }
    stream->block_size = st.st_blksize;
    //fprintf(stderr, "Block size = %zd\n", stream->block_size);
    // Add space for a null byte
    stream->buff_size = XT_FAST_FILE_UNGETC_MAX + stream->block_size + 1;
    if ( (stream->buff = xt_malloc(1, stream->buff_size)) == NULL )
    {
	fputs("ff_init_stream(): Could not allocate buffer.\n", stderr);
	free(stream);
	return NULL;
    }
    stream->start = stream->buff + XT_FAST_FILE_UNGETC_MAX;
    stream->bytes_read = 0;
    stream->c = 0;
    return stream;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <fcntl.h>
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffopen()
 *      initializes a ffile_t stream, much as fopen() does for a FILE
 *      stream.  Unlike fopen(), ffopen() takes the same bit mask
 *      argument as open() to determine the open mode.
 *      See open(3) for details.
 *
 *      An optimally sized buffer for the underlying filesystem is allocated,
 *      along with additional space for limited ffungetc() operations.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      filename    Absolute or relative path of the file to open
 *      flags       Bit flags passed to open(3)
 *
 *  Returns:
 *      A pointer to a ffile_t object on success, NULL on failure
 *
 *  Examples:
 *      ffile_t *stream;
 *      char    *filename;
 *      
 *      // Read only
 *      stream = ffopen(filename, O_RDONLY);
 *
 *      // Overwrite
 *      stream = ffopen(filename, O_WRONLY|O_CREAT|O_TRUNC);
 *
 *      // Append
 *      stream = ffopen(filename, O_WRONLY|O_APPEND);
 *
 *  See also:
 *      open(3), ffgetc(3), ffputc(3), ffungetc(3), ffclose(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-14  Jason Bacon Begin
 ***************************************************************************/

ffile_t *ffopen(const char *filename, int flags)

{
    ffile_t     *stream;
    
    if ( (stream = xt_malloc(1, sizeof(*stream))) == NULL )
	return NULL;

    if ( flags & O_WRONLY )
	stream->fd = open(filename, flags, 0666);   // Masked by umask
    else
	stream->fd = open(filename, flags);
    if ( stream->fd == -1 )
    {
	free(stream);
	return NULL;
    }
    stream->flags = flags;

    return ff_init_stream(stream);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <fcntl.h>
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffdopen()
 *      initializes a ffile_t stream, much as fdopen() does for a FILE
 *      stream.  Unlike fdopen(), ffdopen() takes the same bit mask
 *      argument as open() to determine the open mode.
 *      See open(3) for details.
 *
 *      An optimally sized buffer for the underlying filesystem is allocated,
 *      along with additional space for limited ffungetc() operations.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      fd          Open file descriptor to which stream is attached
 *      flags       Bit flags passed to open(3)
 *
 *  Returns:
 *      A pointer to a ffile_t object on success, NULL on failure
 *
 *  Examples:
 *      ffile_t *stream;
 *      char    *filename;
 *      int     fd;
 *      
 *      fd = open(filename, O_RDONLY);
 *      stream = ffdopen(fd, O_RDONLY);
 *
 *  See also:
 *      ffopen(3), open(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-14  Jason Bacon Begin
 ***************************************************************************/

ffile_t *ffdopen(int fd, int flags)

{
    ffile_t     *stream;
    
    if ( (stream = xt_malloc(1, sizeof(*stream))) == NULL )
	return NULL;
    stream->fd = fd;
    stream->flags = flags;

    return ff_init_stream(stream);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffgetc()
 *      and the macro equivalent
 *      .B FFGETC()
 *      read a single character from a ffile_t stream opened by ffopen(3).
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      stream  Pointer to an ffile_t object
 *
 *  Returns:
 *      The character read, or EOF if no more data are available
 *
 *  Examples:
 *      ffile_t *stream;
 *      int     ch;
 *
 *      if ( (stream = ffopen(filename, O_RDONLY)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for reading.\n", filename);
 *          exit(EX_NOINPUT);
 *      }
 *      while ( (ch = FFGETC(stream)) != EOF )
 *      {
 *      }
 *      ffclose(stream);
 *
 *  See also:
 *      ffopen(3), ffputc(3), ffclose(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-14  Jason Bacon Begin
 ***************************************************************************/

int     ffgetc(ffile_t *stream)

{
    unsigned char   *start;
    
    if ( stream->c == stream->bytes_read )
    {
	/*
	 *  Move last part of buffer to ffungetc() region.  Only the last
	 *  block read should be < block_size chars, and it will never
	 *  be moved here.
	 */
	start = stream->start + stream->block_size - XT_FAST_FILE_UNGETC_MAX;
	memcpy(stream->buff, start, XT_FAST_FILE_UNGETC_MAX);
		
	if ( (stream->bytes_read =
	      read(stream->fd, stream->start, stream->block_size)) == 0 )
	    return EOF;
	stream->c = 0;
    }
    return stream->start[stream->c++];
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffputc()
 *      and the macro equivalent
 *      .B FFPUTC()
 *      write a single character to a ffile_t stream opened by ffopen(3).
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      ch      Character to write to stream
 *      stream  Pointer to an ffile_t object opened by ffopen(3)
 *
 *  Returns:
 *      The character written, or EOF if unable to write
 *
 *  Examples:
 *      char    *infilename, *outfilename;
 *      ffile_t *instream, *outstream;
 *      int     ch;
 *
 *      if ( (instream = ffopen(infilename, O_RDONLY)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for reading.\n", infilename);
 *          exit(EX_NOINPUT);
 *      }
 *      if ( (outstream = ffopen(outfilename, O_WRONLY|O_CREAT|O_TRUNC)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for writing.\n", outfilename);
 *          exit(EX_NOINPUT);
 *      }
 *      while ( (ch = FFGETC(stream)) != EOF )
 *          FFPUTC(ch, outstream);
 *      ffclose(instream);
 *      ffclose(outstream);
 *
 *  See also:
 *      ffopen(3), ffgetc(3), ffclose(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-14  Jason Bacon Begin
 ***************************************************************************/

int     ffputc(int ch, ffile_t *stream)

{
    if ( stream->c == stream->block_size )
    {
	if ( write(stream->fd, stream->start, stream->block_size) != stream->block_size )
	    return EOF;
	stream->c = 0;
    }
    //putchar(ch);
    stream->start[stream->c++] = ch;
    return ch;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffclose()
 *      closes a ffile_t stream opened by ffopen(3).  It writes out any
 *      remaining data in the output buffer, deallocates memory allocated
 *      by ffopen(3), and closes the underlying file descriptor opened by
 *      open(3).
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      stream  Pointer to an ffile_t object opened by ffopen(3)
 *
 *  Returns:
 *      The return status of the underlying close(3) call
 *
 *  Examples:
 *      char    *infilename, *outfilename;
 *      ffile_t *instream, *outstream;
 *      int     ch;
 *
 *      if ( (instream = ffopen(infilename, O_RDONLY)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for reading.\n", infilename);
 *          exit(EX_NOINPUT);
 *      }
 *      if ( (outstream = ffopen(outfilename, O_WRONLY|O_CREAT|O_TRUNC)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for writing.\n", outfilename);
 *          exit(EX_NOINPUT);
 *      }
 *      while ( (ch = FFGETC(stream)) != EOF )
 *          FFPUTC(ch, outstream);
 *      ffclose(instream);
 *      ffclose(outstream);
 *
 *  See also:
 *      ffopen(3), ffgetc(3), ffputc(3)
 *  
 *  History: 
 *  Date        Name        Modification
 *  2022-02-14  Jason Bacon Begin
 ***************************************************************************/

int     ffclose(ffile_t *stream)

{
    int     status;
    
    if ( stream->flags & O_WRONLY )
    {
	//fprintf(stderr, "ffclose() flushing output...\n");
	//stream->start[stream->c] = '\0';
	//fputs((char *)stream->start, stderr);
	write(stream->fd, stream->start, stream->c);
    }
    status = close(stream->fd);
    free(stream->buff);
    free(stream);
    return status;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffungetc()
 *      returns a single character read by ffgetc(3) to the input buffer of
 *      a stream opened by ffopen(3).  All characters from the most recently
 *      read block plus a maximum of XT_FAST_FILE_UNGETC_MAX characters
 *      from the previously read block may be returned.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      ch      Character to return to the input buffer
 *      stream  Pointer to an ffile_t object opened by ffopen(3)
 *
 *  Returns:
 *      The character written, or EOF if unable to write
 *
 *  Examples:
 *      char    *infilename;
 *      ffile_t *instream;
 *      int     ch;
 *
 *      if ( (instream = ffopen(infilename, O_RDONLY)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for reading.\n", infilename);
 *          exit(EX_NOINPUT);
 *      }
 *      if ( (ch = FFGETC(instream)) != MY_FAVORITE_CHAR )
 *          ungetc(ch, instream);
 *      ffclose(instream);
 *
 *  See also:
 *      ffopen(3), ffgetc(3), ffclose(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-18  Jason Bacon Begin
 ***************************************************************************/

int     ffungetc(int ch, ffile_t *stream)

{
    if ( stream->c > -(XT_FAST_FILE_UNGETC_MAX + 1) )
    {
	stream->start[--stream->c] = ch;
	return ch;
    }
    else
	return EOF;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffstdin()
 *      is a simple wrapper function for connecting file descriptor 0
 *      to an ffile_t object using ffdopen(3).  This is useful for
 *      high-performance filter programs, where using the traditional
 *      FILE *stdin would cause a bottleneck.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      None
 *
 *  Returns:
 *      Pointer to an ffile_t object if successful, NULL otherwise
 *
 *  Examples:
 *      ffile_t *stream;
 *
 *      // "-" as a filename argument traditionally indicates stdin
 *      if ( strcmp(argv[arg], "-") == 0 )
 *          stream = ffstdin();
 *      else
 *          stream = ffopen(argv[arg], O_RDONLY);
 *
 *  See also:
 *      ffopen(3), ffdopen(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-19  Jason Bacon Begin
 ***************************************************************************/

ffile_t *ffstdin()

{
    return ffdopen(0, O_RDONLY);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffstdout()
 *      is a simple wrapper function for connecting file descriptor 1
 *      to an ffile_t object using ffdopen(3).  This is useful for
 *      high-performance filter programs, where using the traditional
 *      FILE *stdout would cause a bottleneck.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      None
 *
 *  Returns:
 *      Pointer to an ffile_t object if successful, NULL otherwise
 *
 *  Examples:
 *      ffile_t *stream;
 *
 *      // "-" as a filename argument traditionally indicates stdout
 *      if ( strcmp(argv[arg], "-") == 0 )
 *          stream = ffstdout();
 *      else
 *          stream = ffopen(argv[arg], O_WRONLY|O_CREAT|O_TRUNC);
 *
 *  See also:
 *      ffopen(3), ffdopen(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-19  Jason Bacon Begin
 ***************************************************************************/

ffile_t *ffstdout()

{
    return ffdopen(1, O_WRONLY|O_APPEND);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffpopen(3)
 *      creates a pipe for interprocess communication, runs the specified
 *      command, connecting the command's standard input or standard
 *      output to the pipe, and returning a pointer to a ffile_t object
 *      connected to the other end.
 *
 *      It behaves much like popen(3), except that it returns a fast-file
 *      fffile_t pointer rather than a standard I/O FILE pointer, and
 *      accepts a full set of open(3) flags rather than the fopen(3)
 *      type strings "r", "w", etc.
 *
 *      This allows the calling program to spawn a child process
 *      and read its standard output or write to its standard input as
 *      easily as reading or writing a file.
 *
 *      The stream should be closed with ffpclose(3) rather than ffclose(3)
 *      in order to wait for the child process to complete and return its
 *      exit status.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      cmd     Full command to execute as the child, passed to sh(1)
 *      flags   Open mode flags passed to open(3)
 *
 *  Returns:
 *      Pointer to a ffile_t object on success, NULL otherwise
 *
 *  Examples:
 *      ffile_t *instream;
 *
 *      if ( (instream = ffpopen("xzcat file.xz", O_RDONLY)) == NULL )
 *      {
 *          fprintf(stderr, "Failed to read xzcat file.xz.\n");
 *          exit(EX_NOINPUT);
 *      }
 *
 *      ffpclose(instream);
 *
 *  See also:
 *      ffopen(3), ffpclose(3), popen(3), open(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-19  Jason Bacon Begin
 ***************************************************************************/

ffile_t *ffpopen(const char *cmd, int flags)

{
    pid_t   pid;
    int     fd[2];
    ffile_t *stream = NULL;
    char    *argv[XT_FAST_FILE_MAX_ARGS];
    
    if ( pipe(fd) == 0 )
    {
	if ( (pid = fork()) == 0 )  // Child process
	{
	    // Use shell to process redirection, etc.
	    argv[0] = "sh";
	    argv[1] = "-c";
	    argv[2] = (char *)cmd;
	    argv[3] = NULL;

	    if ( flags == O_RDONLY )    // O_RDONLY = 0x0, not bits
	    {
		// Child runs command and writes standard output to pipe
		// Readers won't get EOF until last descriptor is closed
		// so don't leave this lying around
		close(fd[0]);   // Not used by child
		close(1);
		if ( dup(fd[1]) != 1 )
		{
		    fprintf(stderr, "%s: dup() failed to return 1.\n",
			    __FUNCTION__);
		    return NULL;
		}
		execvp("/bin/sh", argv);
		return NULL;    // Should not be reached
	    }
	    else
	    {
		// Child runs command and reads standard input from pipe
		// Readers won't get EOF until last descriptor is closed
		// so don't leave this lying around
		close(fd[1]);   // Not used by child
		close(0);
		if ( dup(fd[0]) != 0 )
		{
		    fprintf(stderr, "%s: dup() failed to return 0.\n",
			    __FUNCTION__);
		    return NULL;
		}
		execvp("/bin/sh", argv);
		return NULL;    // Should not be reached
	    }
	}
	else
	{
	    if ( flags == O_RDONLY )    // O_RDONLY = 0x0, no bits
	    {
		// Parent reads from child via pipe
		// Readers won't get EOF until last descriptor is closed
		// so don't leave this lying around
		close(fd[1]);   // Not used by parent
		if ( (stream = ffdopen(fd[0], O_RDONLY)) == NULL )
		    return NULL;
	    }
	    else
	    {
		// Parent writes to child via pipe
		// Readers won't get EOF until last descriptor is closed
		// so don't leave this lying around
		close(fd[0]);   // Not used by parent
		if ( (stream = ffdopen(fd[1], O_WRONLY)) == NULL )
		    return NULL;
	    }
    
	    // Set pid in ffile_t stream for waitpid() in ffpclose()
	    stream->child_pid = pid;
	    return stream;
	}
    }
    return stream;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffpclose(3)
 *      closes a stream opened by ffpopen(3), and
 *      waits for the child process to complete and returns its
 *      exit status.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      stream  ffile_t stream opened by ffpopen(3)
 *
 *  Returns:
 *      Exit status of the child process spawned by ffpopen(3), or -1 on error
 *
 *  Examples:
 *      ffile_t *instream;
 *
 *      if ( (instream = ffpopen("xzcat file.xz", O_RDONLY)) == NULL )
 *      {
 *          fprintf(stderr, "Failed to read xzcat file.xz.\n");
 *          exit(EX_NOINPUT);
 *      }
 *
 *      ffpclose(instream);
 *
 *  See also:
 *      ffopen(3), ffpclose(3), popen(3), open(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-19  Jason Bacon Begin
 ***************************************************************************/

int     ffpclose(ffile_t *stream)

{
    int     status = 0;
    pid_t   pid = stream->child_pid;
    
    if ( pid == 0 )
    {
	fprintf(stderr, "%s(): No child PID available.  Was the stream opened with ffpopen()?\n",
		__FUNCTION__);
	return -1;
    }
    
    ffclose(stream);
    
    // Compatibility with pclose()
    waitpid(pid, &status, 0);
    //fprintf(stderr, "Back from waitpid().\n");
    
    return status;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      .B xt_ffopen(3)
 *      opens a raw data file using ffopen() or a gzipped, bzipped, or
 *      xzipped file using ffpopen(), returning a pointer to a ffile_t
 *      stream.  Must be used in conjunction with
 *      xt_ffclose() to ensure that ffclose() or ffpclose() is called where
 *      appropriate.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *
 *  Arguments:
 *      filename:   Name of the file to be opened
 *      mode:       Bit mask as used by open()
 *
 *  Returns:
 *      A pointer to the FILE structure or NULL if open failed
 *
 *  See also:
 *      fopen(3), popen(3), gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

ffile_t *xt_ffopen(const char *filename, int flags)

{
    char    *ext = strrchr(filename, '.'),
	    cmd[XT_CMD_MAX_CHARS + 1];
    
    if ( ext == NULL )
    {
	// FIXME: Use __FUNCTION__ in all such messages
	fprintf(stderr, "%s(): No filename extension on %s.\n",
		__FUNCTION__, filename);
	return NULL;
    }

    //fprintf(stderr, "flags = %x\n", flags);
    if ( flags == O_RDONLY )    // O_RDONLY = 0x0, no bits set
    {
	//fprintf(stderr, "Reading from %s...\n", filename);
	if ( strcmp(ext, ".gz") == 0 )
	{
// Big Sur zcat requires a .Z extension and CentOS 7 lacks gzcat
#ifdef __APPLE__
	    snprintf(cmd, XT_CMD_MAX_CHARS, "gzcat %s", filename);
#else
	    snprintf(cmd, XT_CMD_MAX_CHARS, "zcat %s", filename);
#endif
	    return ffpopen(cmd, flags);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "bzcat %s", filename);
	    return ffpopen(cmd, flags);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "xzcat %s", filename);
	    return ffpopen(cmd, flags);
	}
	else
	    return ffopen(filename, flags);
    }
    else    // O_WRONLY
    {
	//fprintf(stderr, "Writing to %s...\n", filename);
	if ( strcmp(ext, ".gz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "gzip -c > %s", filename);
	    return ffpopen(cmd, flags);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "bzip2 -c > %s", filename);
	    return ffpopen(cmd, flags);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "xz -c > %s", filename);
	    return ffpopen(cmd, flags);
	}
	else
	    return ffopen(filename, flags);
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      .B xt_ffclose(3)
 *      closes a ffile_t stream with ffclose() or ffpclose() as appropriate.
 *      Automatically determines the proper close function to call using
 *      S_ISFIFO on the stream stat structure.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *
 *  Arguments:
 *      stream: Pointer to the ffile_t structure to be closed
 *
 *  Returns:
 *      The value returned by ffclose() or ffpclose()
 *
 *  See also:
 *      ffpopen(3), xt_ffopen(3), gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-10  Jason Bacon Begin
 ***************************************************************************/

int     xt_ffclose(ffile_t *stream)

{
    struct stat stat;
    
    fstat(stream->fd, &stat);
    if ( S_ISFIFO(stat.st_mode) )
	return ffpclose(stream);
    else
	return ffclose(stream);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffprintf(3)
 *      writes formatted data to a ffile_t stream the same was as
 *      fprintf(3) writes to a FILE stream.
 *
 *      The ffile_t system is simpler than and several times as
 *      fast as FILE on typical systems.  It is intended for processing
 *      large files character-by-character, where low-level block I/O
 *      is not convenient, but FILE I/O causes a bottleneck.
 *  
 *  Arguments:
 *      stream  Pointer to an ffile_t object opened by ffopen(3)
 *      format  Format string indicating how remaining arguments are printed
 *
 *  Returns:
 *      The number of characters written
 *
 *  Examples:
 *      ffile_t *stream;
 *      int     count = 1;
 *
 *      if ( (stream = ffopen(filename, O_WRONLY|O_CREAT|O_TRUNC)) == NULL )
 *      {
 *          fprintf(stderr, "Could not open %s.\n", filename);
 *          exit(EX_CANTCREAT);
 *      }
 *      ffprintf(stream, "%d\n", count);
 *      ffclose(stream);
 *
 *  See also:
 *      fprintf(3), ffopen(3), ffclose(3), ffputc(3), ffputs(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-19  Jason Bacon Begin
 ***************************************************************************/

int     ffprintf(ffile_t *stream, const char *format, ...)

{
    va_list ap;
    char    *buff;
    int     chars_printed, c;
    
    va_start(ap, format);
    chars_printed = vasprintf(&buff, format, ap);
    for (c = 0; buff[c] != '\0'; ++c)
	FFPUTC(buff[c], stream);
    free(buff);
    return chars_printed;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      ffputs() writes a null-terminated string to the given ffile_t
 *      stream.  It is fnuctionally equivalent to fputs() with FILE.
 *  
 *  Arguments:
 *      string      A null-terminated string
 *      stream      Pointer to an ffile_t structure opened with ffopen()
 *
 *  Returns:
 *      A non-negative integer on success, EOF on failure
 *
 *  Examples:
 *      ffile_t *outstream;
 *      char    *buff;
 *
 *      if ( (outstream = ffopen(outfilename, O_WRONLY|O_CREAT|O_TRUNC)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for writing.\n", outfilename);
 *          exit(EX_NOINPUT);
 *      }
 *      ffputs("Hello, world!\n", outstream);
 *      ffclose(outstream);
 *
 *  See also:
 *      fputs(3), ffgets(3), ffopen(3), ffclose(3), ffputc(3), ffprintf(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-07-29  Jason Bacon Begin
 ***************************************************************************/

int     ffputs(const char *string, ffile_t *stream)

{
    size_t  c;
    int     status = 0;
    
    for (c = 0; (status >= 0) && (string[c] != '\0'); ++c)
	status = FFPUTC(string[c], stream);
    return status;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      ffgets() writes a line of text from the given ffile_t
 *      stream.  It is fnuctionally equivalent to fgets() with FILE.
 *      The maximum number of characters read is size - 1, to allow
 *      for a null-terminator byte.
 *  
 *  Arguments:
 *      string      A character array into which the line is read
 *      size        Size of the character array
 *      stream      Pointer to an ffile_t structure opened with ffopen()
 *
 *  Returns:
 *      A non-negative integer on success, EOF on failure
 *
 *  Examples:
 *      ffile_t *instream;
 *      char    buff[BUFF_SIZE];
 *
 *      if ( (instream = ffopen(outfilename, O_RDONLY)) == NULL )
 *      {
 *          fprintf(stderr, "Cannot open %s for writing.\n", outfilename);
 *          exit(EX_NOINPUT);
 *      }
 *      ffgets(buff, BUFF_SIZE, instream);
 *      ffclose(instream);
 *
 *  See also:
 *      fgets(3), ffputs(3), ffopen(3), ffclose(3), ffputc(3), ffprintf(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-07-29  Jason Bacon Begin
 ***************************************************************************/

char    *ffgets(char *string, size_t size, ffile_t *stream)

{
    size_t  c;
    int     ch;
    
    c = 0;
    while ( (c < size - 1) && ((ch = FFGETC(stream)) != '\n') && (ch != EOF) )
	string[c++] = ch;
    if ( (c == 0) && (ch == EOF) )
	return NULL;
    else
	return string;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/fast-file.h>
 *      -lxtend
 *
 *  Description:
 *      .B ffread_line_malloc()
 *      reads a single line of text (up to the next newline or EOF)
 *      from stream, allocating and/or extending the provided buffer if
 *      needed.
 *  
 *  Arguments:
 *      stream:     ffile_t stream from which field is read
 *      buff:       Character buffer into which field is copied
 *      buff_size:  Size of the array passed to buff
 *      len:        Pointer to a variable which will receive the field length
 *
 *  Returns:
 *      Delimiter ending the read: either newline or EOF
 *
 *  Examples:
 *      ffile_t *stream;
 *      char    *buff;
 *      size_t  buff_len, len;
 *
 *      while ( ffile_read_line_malloc(stream, buff, &buff_len, &len) != EOF )
 *      {
 *      }
 *
 *  See also:
 *      dsv_read_field_malloc(3), ffgetc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-20  Jason Bacon Begin
 ***************************************************************************/

int     ffread_line_malloc(ffile_t *stream, char **buff, size_t *buff_size,
			   size_t *len)

{
    size_t  c;
    int     ch;
    
    if ( *buff_size == 0 )
    {
	*buff_size = 1024;
	*buff = xt_malloc(*buff_size, sizeof(**buff));
	if ( *buff == NULL )
	    return XT_MALLOC_FAILED;
    }
    
    for (c = 0; ( ((ch = FFGETC(stream)) != '\n') && (ch != EOF) ); ++c)
    {
	if ( c == *buff_size - 1 )
	{
	    *buff_size *= 2;
	    *buff = xt_realloc(*buff, *buff_size, sizeof(**buff));
	    if ( *buff == NULL )
		return XT_MALLOC_FAILED;
	}
	(*buff)[c] = ch;
    }
    (*buff)[c] = '\0';
    *len = c;

    /* Trim array */
    if ( *buff_size != c + 1 )
    {
	*buff_size = c + 1;
	*buff = xt_realloc(*buff, *buff_size, sizeof(**buff));
    }
    return ch;
}
#include <unistd.h>
#include <fcntl.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      fd_purge() reads and discards unwanted input data (such as leftover
 *      input from a keyboard or mouse) from the input file descriptor fd.
 *  
 *  Arguments:
 *      fd: File descriptor to purge
 *
 *  See also:
 *      read(3), fcntl(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

void    xt_fd_purge(int fd)

{
    char    buff[128];
    int     old_flags;
    
    old_flags = fcntl(fd,F_GETFL,0);
    
    /* Prevent processing from waiting for new input */
    fcntl(fd,F_SETFL,O_NONBLOCK);
    
    /* Remove pending characters */
    while ( read(fd,buff,128) != -1 )
	;
    
    fcntl(fd,F_SETFL,old_flags);
}
#include <stdio.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      fgetline() reads a line of text from a FILE stream.  Input is
 *      terminated when a newline or end of file is encountered,
 *      or when maxlen characters have been read.  Note that up to maxlen
 *      characters may be stored, NOT INCLUDING THE NULL TERMINATOR BYTE,
 *      hence the buffer should be at least maxlen+1 bytes long. Unlike
 *      fgets(3), fgetline() does not store the trailing newline character
 *      in the string.
 *  
 *  Arguments:
 *      fp:     Input stream from which to read
 *      buff:   Character array into which line is read
 *      maxlen: Size of array buff, not counting null byte
 *
 *  Returns:
 *      The number of bytes read, or EOF if EOF is encountered before a newline
 *
 *  See also:
 *      fgets(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/
 
size_t  xt_fgetline(FILE *fp, char *buff, size_t maxlen)

{
    char    *p = buff, 
	    *end = buff+maxlen;
    int     ch;

    /* Read to end of line, end of file, or maxlen characters */
    while ( ((ch = getc(fp)) != EOF) && (ch != '\n') && (p < end) )
	*p++ = ch;
    *p = '\0';      /* Replace \n with \0, or add \0 after last char */
    
    if (ch == EOF)
	return EOF;
    else
	return p - buff;  /* Return string length */
}
#include <errno.h>
#include <sys/stat.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      file_mod_cmp() compares the modification times of file1 and file2.
 *      It compares modification times on two files using the same rules
 *      as "make", and returns an strcmp(3) compatible status value indicating
 *      which is older.  A file that doesn't exist is considered
 *      older than the big bang.
 *  
 *  Arguments:
 *      file1, file2: Names of two filesystem objects whose time stamps
 *      are to be compared
 *
 *  Returns:
 *      A value < 0 if file1 is older or does not exist.
 *      A value > 0 if file2 is older or does not exist.
 *      0 if the files have identical modification times, or neither exists.
 *
 *  See also:
 *      make(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

int     xt_file_mod_cmp(const char *file1, const char *file2)

{
    struct stat stats1, stats2;
    int     rc1,rc2;
    extern int errno;

    rc1 = stat(file1, &stats1);
    rc2 = stat(file2, &stats2);
    
    /* Both or neither exist */
    if ( rc1 == rc2 )
    {
	if (rc1 == 0)   /* Both files exist */
	    return(stats1.st_mtime - stats2.st_mtime);
	else
	    return(0);  /* Neither file exists */
    }
    else    /* One of the two files exists */
    {
	if ((rc1 == -1) && (errno == ENOENT))   /* file1 doesn't exist */
	    return (-1);
	else if ((rc2 == -1) && (errno == ENOENT))   /* file2 dne */
	    return (1);
    }
    return(0);
}

/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      Computes the greatest common divisor of two natural
 *      numbers a and b.
 *  
 *  Arguments:
 *      a, b: Numbers for which to find GCD
 *
 *  Returns:
 *      The greatest common divisor of a and b.
 *
 *  See also:
 *      lcm(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

unsigned long   gcd(unsigned long a, unsigned long b)

{
    if ( a == 0 ) return b;
    if ( b == 0 ) return a;
    if ( a < b )
	return gcd(b,a);
    else
	return gcd(b, a % b);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      Computes the least common multiple of two natural
 *      numbers a and b.  Note that this function may fail for relatively
 *      small values, as their LCM may be beyond the range of a 32-bit
 *      integer.
 *  
 *  Arguments:
 *      a, b: Numbers for which to find LCM
 *
 *  Returns:
 *      The least common multiple of a and b.
 *
 *  See also:
 *      gcd(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

unsigned long   lcm(unsigned long a,unsigned long b)

{
    return a * b / gcd(a,b);
}
#include <string.h>
#include <unistd.h>
#include <sys/types.h>
#include <pwd.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      get_home_dir() determines the full pathname of the process owner's
 *      home directory.  The information is retrieved using a call to
 *      getpwuid(3), and copied to the argument "dir".
 *   
 *      The name is stored in dir up to maxlen characters.
 *      Note that up to maxlen characters are stored, not including the 
 *      null terminator, hence the buffer should be at least maxlen+1
 *      bytes long.
 *  
 *  Arguments:
 *      dir:    Character buffer to receive home directory path
 *      maxlen: Max characters to copy to dir, not including null byte
 *
 *  Returns:
 *      A pointer to dir, or NULL upon failure.
 *
 *  See also:
 *      getuid(3), getpwuid(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

char   *xt_get_home_dir(char *dir, size_t maxlen)

{
    int     user;
    struct passwd *pwentry;

    /* Determine who the user is */
    user = getuid();

    /* Get password file entry */
    if ((pwentry = getpwuid(user)) == NULL)
	return (NULL);
 
    strlcpy(dir, pwentry->pw_dir,maxlen);
    return dir;
}

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      Move file from pathname src to pathname dest. First attempt to
 *      rename using rename(3).  This will fail if src and dest are
 *      in different filesystems.  Then attempt to copy the file using
 *      fast_cp(3), an optimized cross-filesystem file copy routine.
 *  
 *  Arguments:
 *      src     Original filename
 *      dest    New filename
 *
 *  Returns:
 *      0 on success, otherwise error code from fastcp(3)
 *
 *  Examples:
 *      char    *old_name, char *new_name;
 *
 *      if ( mv(old_name, new_name) != 0 )
 *      {
 *          fprintf(stderr, "Failed to move %s to %s.\n", old_name, new_name);
 *          ...
 *      }
 *
 *  See also:
 *      fast_cp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-10  Jason Bacon Begin
 ***************************************************************************/

#include <unistd.h>

int     mv(const char *src, const char *dest)

{
    int     status = 0;
    
    if ( rename(src,dest) != 0 )
    {
	if ( (status = xt_fast_cp(src,dest)) == 0 )
	    unlink(src);
	else
	    unlink(dest);
    }
    return status;
}

/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      This is a function that compares two doubles as a service to
 *      polymorphic functions such as qsort(3), bsearch(3), etc.  The
 *      address of double_cmp() is passed as an argument to perform the
 *      data type specific comparison on behalf of the sort of search function.
 *  
 *  Arguments:
 *      n1, n2  Pointers to two double values
 *
 *  Returns:
 *      A value > 0 if *n1 is greater than *n2
 *      A value < 0 if *n1 is less than *n2
 *      0 if the values are equal
 *
 *  Examples:
 *      double  list[LIST_SIZE];
 *
 *      // sizeof(*list) will continue to work if we change the data type
 *      // We'll still need to change the cmp function, though
 *      qsort(list, LIST_SIZE, sizeof(*list),
 *            (int (*)(const void *, const void *))double_cmp);
 *      
 *  See also:
 *      qsort(3), heapsort(3), mergesort(3), bsearch(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-07  Jason Bacon Begin
 ***************************************************************************/

int     double_cmp(const double *n1, const double *n2)

{
    /*
     *  Don't just return *n1 - *n2, since it will truncate to 0 in some
     *  cases where one is actually greater
     */
    
    if ( *n1 > *n2 )
	return 1;
    else if ( *n1 < *n2 )
	return -1;
    else
	return 0;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      This is a function that compares two floats as a service to
 *      polymorphic functions such as qsort(3), bsearch(3), etc.  The
 *      address of float_cmp() is passed as an argument to perform the
 *      data type specific comparison on behalf of the sort of search function.
 *  
 *  Arguments:
 *      n1, n2  Pointers to two float values
 *
 *  Returns:
 *      A value > 0 if *n1 is greater than *n2
 *      A value < 0 if *n1 is less than *n2
 *      0 if the values are equal
 *
 *  Examples:
 *      float  list[LIST_SIZE];
 *
 *      // sizeof(*list) will continue to work if we change the data type
 *      // We'll still need to change the cmp function, though
 *      qsort(list, LIST_SIZE, sizeof(*list),
 *            (int (*)(const void *, const void *))float_cmp);
 *      
 *  See also:
 *      qsort(3), heapsort(3), mergesort(3), bsearch(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-07  Jason Bacon Begin
 ***************************************************************************/

int     float_cmp(const float *n1, const float *n2)

{
    /*
     *  Don't just return *n1 - *n2, since it will truncate to 0 in some
     *  cases where one is actually greater
     */
    
    if ( *n1 > *n2 )
	return 1;
    else if ( *n1 < *n2 )
	return -1;
    else
	return 0;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      This is a function that compares two long longs as a service to
 *      polymorphic functions such as qsort(3), bsearch(3), etc.  The
 *      address of long_long_cmp() is passed as an argument to perform the
 *      data type specific comparison on behalf of the sort of search function.
 *  
 *  Arguments:
 *      n1, n2  Pointers to two long long values
 *
 *  Returns:
 *      A value > 0 if *n1 is greater than *n2
 *      A value < 0 if *n1 is less than *n2
 *      0 if the values are equal
 *
 *  Examples:
 *      long long  list[LIST_SIZE];
 *
 *      // sizeof(*list) will continue to work if we change the data type
 *      // We'll still need to change the cmp function, though
 *      qsort(list, LIST_SIZE, sizeof(*list),
 *            (int (*)(const void *, const void *))long_long_cmp);
 *      
 *  See also:
 *      qsort(3), heapsort(3), mergesort(3), bsearch(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-07  Jason Bacon Begin
 ***************************************************************************/

int     long_long_cmp(const long long *n1, const long long *n2)

{
    /*
     *  Don't just return *n1 - *n2, since it might exceed the range of
     *  an int
     */
    
    if ( *n1 > *n2 )
	return 1;
    else if ( *n1 < *n2 )
	return -1;
    else
	return 0;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      This is a function that compares two longs as a service to
 *      polymorphic functions such as qsort(3), bsearch(3), etc.  The
 *      address of long_cmp() is passed as an argument to perform the
 *      data type specific comparison on behalf of the sort of search function.
 *  
 *  Arguments:
 *      n1, n2  Pointers to two long values
 *
 *  Returns:
 *      A value > 0 if *n1 is greater than *n2
 *      A value < 0 if *n1 is less than *n2
 *      0 if the values are equal
 *
 *  Examples:
 *      long  list[LIST_SIZE];
 *
 *      // sizeof(*list) will continue to work if we change the data type
 *      // We'll still need to change the cmp function, though
 *      qsort(list, LIST_SIZE, sizeof(*list),
 *            (int (*)(const void *, const void *))long_cmp);
 *      
 *  See also:
 *      qsort(3), heapsort(3), mergesort(3), bsearch(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-07  Jason Bacon Begin
 ***************************************************************************/

int     long_cmp(const long *n1, const long *n2)

{
    /*
     *  Don't just return *n1 - *n2, since it might exceed the range of
     *  an int
     */
    
    if ( *n1 > *n2 )
	return 1;
    else if ( *n1 < *n2 )
	return -1;
    else
	return 0;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      This is a function that compares two ints as a service to
 *      polymorphic functions such as qsort(3), bsearch(3), etc.  The
 *      address of int_cmp() is passed as an argument to perform the
 *      data type specific comparison on behalf of the sort of search function.
 *  
 *  Arguments:
 *      n1, n2  Pointers to two int values
 *
 *  Returns:
 *      A value > 0 if *n1 is greater than *n2
 *      A value < 0 if *n1 is less than *n2
 *      0 if the values are equal
 *
 *  Examples:
 *      int  list[LIST_SIZE];
 *
 *      // sizeof(*list) will continue to work if we change the data type
 *      // We'll still need to change the cmp function, though
 *      qsort(list, LIST_SIZE, sizeof(*list),
 *            (int (*)(const void *, const void *))int_cmp);
 *      
 *  See also:
 *      qsort(3), heapsort(3), mergesort(3), bsearch(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-07  Jason Bacon Begin
 ***************************************************************************/

int     int_cmp(const int *n1, const int *n2)

{
    return *n1 - *n2;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/math.h>
 *      -lxtend
 *
 *  Description:
 *      This is a function that compares two shorts as a service to
 *      polymorphic functions such as qsort(3), bsearch(3), etc.  The
 *      address of short_cmp() is passed as an argument to perform the
 *      data type specific comparison on behalf of the sort of search function.
 *  
 *  Arguments:
 *      n1, n2  Pointers to two short values
 *
 *  Returns:
 *      A value > 0 if *n1 is greater than *n2
 *      A value < 0 if *n1 is less than *n2
 *      0 if the values are equal
 *
 *  Examples:
 *      short  list[LIST_SIZE];
 *
 *      // sizeof(*list) will continue to work if we change the data type
 *      // We'll still need to change the cmp function, though
 *      qsort(list, LIST_SIZE, sizeof(*list),
 *            (int (*)(const void *, const void *))short_cmp);
 *      
 *  See also:
 *      qsort(3), heapsort(3), mergesort(3), bsearch(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-07  Jason Bacon Begin
 ***************************************************************************/

int     short_cmp(const short *n1, const short *n2)

{
    return *n1 - *n2;
}
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <sysexits.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/proc.h>
 *      -lxtend
 *
 *  Description:
 *      Breaks a shell command into an argv[] style array suitable
 *      for spawnvp() or execv*().  A copy of cmd is created using
 *      strshellcpy(), which expands certain shell features such as
 *      variables and paths starting with '~'.  The copy is then
 *      modified by replacing separators with '\0' and the argv[] array
 *      is populated with pointers to each token in the copy.
 *  
 *  Arguments:
 *      argv:   Pointer array to be filled with command tokens
 *      cmd:    Raw command string with limited meta-character support
 *              from strshellcpy(3)
 *
 *  Returns:
 *      Pointer to strdup() copy of cmd, which should be freed as soon
 *      as possible when argv[] is no longer needed.
 *
 *
 *  Examples:
 *      char *cmd, *argv[], *expanded_cmd;
 *
 *      expanded_cmd = parse_cmd(argv, cmd);
 *      spawnvp(P_WAIT, P_NOECHO, argv, NULL, NULL, NULL);
 *      free(expanded_cmd);
 *
 *  See also:
 *      spawnvp(3), spawnlp(3), exec(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

char    *parse_cmd(char *argv[], int max_args, const char *cmd)

{
    char    *cmd_copy;
    int     c;

    if ( (cmd_copy = malloc(XT_CMD_MAX_CHARS)) == NULL )
    {
	fprintf(stderr, "parse_cmd(): malloc failed.\n");
	exit(EX_UNAVAILABLE);
    }
    
    /* Expand shell meta-characters */
    // FIXME: Make sure strshellcpy() stops 1 short of max and remove -1
    strshellcpy(cmd_copy, cmd, XT_CMD_MAX_CHARS - 1);
    
    /* Break command into tokens for argv[] */
    // FIXME: Replace deprecated strtok() with strsep()
    argv[0] = strtok(cmd_copy, " \t");
    for (c = 1; (c < max_args) && (argv[c] = strtok(NULL, " \t")) != NULL; ++c)
	;
    return cmd_copy;
}
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <string.h>
#include <netdb.h>
#include <arpa/inet.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Resolve a host name to an IP address.
 *  
 *  Arguments:
 *      hostname    Name of the host to be resolved
 *      ip          Character array to receive IP address
 *      ip_buff_len Size of ip array including null byte
 *
 *  Returns:
 *      XT_OK on success, XT_FAIL otherwise
 *
 *  Examples:
 *      #define IP_MAX_CHARS    64
 *
 *      char    *hostname = "my.site.edu",
 *              ip[IP_MAX_CHARS + 1];
 *
 *      if ( resolve_hostname(hostname, ip, IP_MAX_CHARS + 1) == XT_OK )
 *      {
 *      }
 *
 *  See also:
 *      gethostbyname(3), getaddrinfo(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-09-28  Jason Bacon Begin
 ***************************************************************************/

int     resolve_hostname(const char *hostname, char *ip, size_t ip_buff_len)

{
    struct hostent  *ent;
    struct in_addr  **address_list;

    /*
     *  FIXME: Reimplement with getaddrinfo() to better support IPv6
     *  gethostbyname() is simpler and will suffice for now
     */
    
    if ( (ent = gethostbyname(hostname)) == NULL )
    {
	herror("resolve_hostname(): gethostbyname() failed");
	fprintf(stderr, "hostname = %s\n", hostname);
	fputs("Check /etc/hosts and /etc/resolv.conf.\n", stderr);
	return XT_FAIL;
    }

    // Just take first address
    address_list = (struct in_addr **)ent->h_addr_list;
    strlcpy(ip, inet_ntoa(*address_list[0]), ip_buff_len);
    
    return XT_OK;
}
#include <string.h>
#include <unistd.h>
#include <sys/stat.h>
#include <sys/types.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      rmkdir() recursively creates a directory from within a compiled
 *      program in the same way as "mkdir -r" from shell.
 *  
 *  Arguments:
 *      path:   Absolute or relative pathname of the directory to create
 *      mode:   Permissions and other bits passed to mkdir(2)
 *
 *  Returns:
 *      0 on success, -1 on failure
 *
 *  See also:
 *      mkdir(1), mkdir(2)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-10  Jason Bacon Begin
 ***************************************************************************/

int     xt_rmkdir(const char *path, mode_t mode)

{
    char    *parent_end;
    
    /* First try given path */
    if ( mkdir(path,mode) == 0 )
	return 0;
    else
    {
	/* Recursively attempt to make parent directories */
	parent_end = strrchr(path,'/');
	if ( parent_end == NULL )
	{
	    /* Ran out of ancestors - give it up */
	    return -1;
	}
	else
	{
	    /* Try to make parent with recursive call */
	    *parent_end = '\0';
	    if ( xt_rmkdir(path,mode) == 0 )
	    {
		/* If parent successfully made, try again to make current */
		*parent_end = '/';
		return mkdir(path,mode);
	    }
	    else
	    {
		*parent_end = '/';
		return -1;
	    }
	}
    }
}
#include <stdio.h>
#include <sysexits.h>
#include <stdlib.h>
#include <ctype.h>


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/stdlib.h>
 *      -lxtend
 *
 *  Description:
 *      The 
 *      .B romantoi() function converts a string containing a valid
 *      Roman numeral to an integer, much like strtol().  It rejects
 *      non-normalized values, such as IIIII, XXXXX, or CCCCC, which
 *      should be written as V, L, and D, respectively.  IIII, XXXX, and
 *      CCCC are accepted in place of IV, XL, and CD.  Any number of
 *      consecutive Ms (1000s) are accepted, since there is no larger digit.
 *      Like strtol(), it returns the
 *      address of the first character not converted as part of the
 *      number.  This can be used to verify that the number ended as
 *      it should have, perhaps with a '\0' byte.
 *  
 *  Arguments:
 *      nptr:   Pointer to the first character of the Roman numeral
 *      endptr: Address of a pointer variable to receive the end of the string
 *
 *  Returns:
 *      The integer value of the Roman numeral if converted, or 0 if
 *      an invalid numeral is detected.
 *
 *  Examples:
 *      char    string[] = "XIV", *end;
 *      int     n;
 *
 *      n = romantoi(string, &end);
 *      if ( *end == '\0' )
 *          printf("%d\n", n);
 *      else
 *          fprintf(stderr, "Error converting %s.\n", string);
 *
 *  See also:
 *      strtol(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-12-14  Jason Bacon Begin
 ***************************************************************************/

int     romantoi(const char *nptr, char **endptr)

{
    int     digit, next_digit, previous_digit, val, consecutive;
    const char    *p;
    
    // Array of values using subscripts from 'C' to 'X'
    const static int  digits[] = 
    {
	// I = 1, V = 5, X = 10, L = 50, C = 100, D = 500, M = 1000
	100, 500, 0, 0, 0, 0, 1, 0, 0, 50, 1000, 0, 0, 0, 0, 0,
	0, 0, 0, 5, 0, 10
    };

    // FIXME: Check for more than 3 consecutive identical digits
    val = 0;
    previous_digit = 0;
    p = nptr;
    while ( isalpha(*p) )
    {
	digit = digits[toupper(*p) - 'C'];
	// fprintf(stderr, "digit = %d\n", digit);
	
	// Can't have more than 3 I's in a row
	if ( digit == previous_digit )
	{
	    ++consecutive;
	    
	    // IIIII should be V, XXXXX should be L, etc.
	    if ( ((consecutive > 4) && (digit != 1000)) ||
		 ((consecutive > 1) &&
		  ((digit == 5) || (digit == 50) || (digit == 500))) )
	    {
		fprintf(stderr, "romantoi(): Invalid Roman numeral: %s.\n",
			nptr);
		return 0;
	    }
	}
	else
	    consecutive = 1;
	
	if ( digit != 0 )
	{
	    if ( ! isalpha(p[1]) )
		next_digit = 0;
	    else
		next_digit = digits[toupper(*(p + 1)) - 'C'];
	    if ( next_digit > digit )
	    {
		// Only 1 lesser digit allowed before a greater one.
		// E.g. IV is valid, IIV is not.
		if ( consecutive > 1 )
		{
		    fprintf(stderr, "romantoi(): Invalid Roman numeral: %s.\n",
			    nptr);
		    return 0;
		}
		val += next_digit - digit;  // IV, IX, XL, XC, DC, CM
		++p;
	    }
	    else
		val += digit;
	}
	previous_digit = digit;
	++p;
    }
    
    *endptr = (char *)p;
    return val;
}

#include <stdio.h>
#include <stdarg.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/proc.h>
 *      -lxtend
 *
 *  Description:
 *      spawnvp() and spawnlp() are wrappers around fork(2) and exec(3)
 *      which make it easy to run a child process without an intermediate
 *      shell process as is used by system(3).  The spawnlp() function
 *      spawns a child process using a variable argument list.  The 6th
 *      argument is passed to argv[0] of the child, the 7th to argv[1], etc.
 *
 *      The spawnvp() function spawns a process using the command contained
 *      in an argv[] array constructed by the caller.  spawnlp() automatically
 *      constructs such an argv[] array and calls spawnvp().
 *
 *      The calling process waits for the child to complete if P_WAIT is
 *      passed to parent_action, or continues immediately if P_NOWAIT
 *      is passed.  If P_ECHO is passed as the echo argument, the command
 *      is echoed, the command is echoed to the parent's stdout.
 *
 *      If infile, outfile, or errfile are not NULL, then the corresponding
 *      file streams stdin, stdout, or stderr are redirected to the filename
 *      provided.
 *  
 *  Arguments:
 *      parent_action:  P_WAIT or P_NOWAIT
 *      echo:           P_ECHO or P_NOECHO
 *      infile:         File to which stdin of child is redirected or NULL
 *      outfile:        File to which stdout of child is redirected or NULL
 *      errfile:        File to which stderr of child is redirected or NULL
 *
 *  Returns:
 *      The exit status of the child process if P_WAIT is passed
 *      The PID of the child process if P_NOWAIT is passed
 *
 *  See also:
 *      spawnvp(3), fork(2), exec(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

int     spawnlp(int parent_action, int echo,
		char *infile, char *outfile, char *errfile,
		char *arg0, ...)

{
    va_list list;
    char    *argv[100];
    int     c;
    
    va_start(list,arg0);
    argv[0] = arg0;
    for (c=1; (argv[c] = (char *)va_arg(list,char *)) != NULL; ++c)
	;
    return(spawnvp(parent_action,echo,argv,infile,outfile,errfile));
}
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <fcntl.h>
#include <unistd.h>
#include <errno.h>
#include <signal.h>
#include <sys/types.h>
#include <sys/wait.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/proc.h>
 *      -lxtend
 *
 *  Description:
 *      spawnvp() and spawnlp() are wrappers around fork(2) and exec(3)
 *      which make it easy to run a child process without an intermediate
 *      shell process as is used by system(3).  The spawnlp() function
 *      spawns a child process using a variable argument list.  The 6th
 *      argument is passed to argv[0] of the child, the 7th to argv[1], etc.
 *
 *      The spawnvp() function spawns a process using the command contained
 *      in an argv[] array constructed by the caller.  spawnlp() automatically
 *      constructs such an argv[] array and calls spawnvp().
 *
 *      The calling process waits for the child to complete if P_WAIT is
 *      passed to parent_action, or continues immediately if P_NOWAIT
 *      is passed.  If P_ECHO is passed as the echo argument, the command
 *      is echoed, the command is echoed to the parent's stdout.
 *
 *      If infile, outfile, or errfile are not NULL, then the corresponding
 *      file streams stdin, stdout, or stderr are redirected to the filename
 *      provided.
 *  
 *  Arguments:
 *      parent_action:  P_WAIT or P_NOWAIT
 *      echo:           P_ECHO or P_NOECHO
 *      infile:         File to which stdin of child is redirected or NULL
 *      outfile:        File to which stdout of child is redirected or NULL
 *      errfile:        File to which stderr of child is redirected or NULL
 *
 *  Returns:
 *      The exit status of the child process if P_WAIT is passed
 *      The PID of the child process if P_NOWAIT is passed
 *
 *  See also:
 *      spawnlp(3), fork(2), exec(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

// FIXME: Ugly stop-gap for now.  sig_t is not portable, so don't use it.
#if defined(__sun__)
typedef void (*sig_t)(int);
#endif

int     spawnvp(int parent_action, int echo, char *argv[],
		char *infile, char *outfile, char *errfile)

{
    int     stat = 0;
    pid_t   pid;
    char    **p;
    extern int  errno;
    sig_t   oldsig;
    
    switch(echo)
    {
	case    P_ECHO:     /* Echo command */
	    for (p = argv; *p != NULL; ++p)
		printf("%s ",*p);
	    putchar('\n');
	    fflush(stdout);
	case    P_NOECHO:
	    break;
	default:
	    fprintf(stderr,
		"spawnvp(): Invalid echo flag: must be ECHO or NO_ECHO.\n");
	    exit(1);
    }
    
    /* If in child process, exec the new program */
    if ((pid = fork()) == 0)
    {
	redirect(infile,outfile,errfile);
	signal(SIGINT,SIG_DFL); /* Allow child process to be interrupted */
	execvp(argv[0], argv);
	exit(errno|0x80);   /* Return errno - all I could think of */
    }
    else    /* If parent, wait for child to croak */
    {
	switch ( parent_action )
	{   
	    case    P_WAIT:
		/* wait() may fail is SIGCHLD isn't SIG_DFL */
		oldsig = signal(SIGCHLD,SIG_DFL);
		waitpid(pid,&stat,0);
		signal(SIGCHLD,oldsig);
		return stat;
	    case    P_NOWAIT:
		return pid;
	    default:
		fprintf(stderr,"spawnvp(): Invalid parent action.\n");
		exit(1);
	}
    }
    /* Dummy return for some compilers (IRIX) that think the code
       can actually get here */
    return 0;
}


/*************************************************************************
 * Name:
 *  Redirect stdin, stdout and stderr if corresponding argument isn't NULL
 *
 * Description: 
 *  This function redirects the stdin, stdout, and stderr of the current
 *  process to the files named by the corresponding arguments.  The original
 *  file streams are not preserved.  If you need to restore any of these
 *  streams to their original state, they must be saved (e.g. using dup(),
 *  dup2(), or ttyname()) prior to calling redirect().
 * 
 * Author: 
 *  Jason W. Bacon
 ****************************************************************************/
 
void    redirect(
    char    *infile,    /* If not NULL, stdin is redirected from this file */
    char    *outfile,   /* If not NULL, stdout is redirected to this file */
    char    *errfile    /* If not NULL, stderr is redirected to this file */
    )

{
    if (infile != NULL)
    {
	close(0);
	if ( open(infile, O_RDONLY) == -1 )
	    fprintf(stderr,"redirect(): Cannot open infile: %s.\n",infile);
    }
    if (outfile != NULL)
    {
	close(1);
	if ( open(outfile, O_WRONLY | O_CREAT | O_TRUNC, 0600) == -1 )
	    fprintf(stderr,"redirect(): Cannot open outfile: %s.\n",outfile);
    }
    if (errfile != NULL)
    {
	close(2);
	if ( strcmp(errfile,outfile) == 0 )
	{
	    if ( dup(1) == -1 )
		fprintf(stderr,"redirect(): Cannot open errfile: %s.\n",errfile);
	}
	else
	{
	    if ( open(errfile, O_WRONLY | O_CREAT | O_TRUNC, 0600) == -1 )
		fprintf(stderr,"redirect(): Cannot open errfile: %s.\n",errfile);
	}
    }
}

#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strlupper(3) copies a string from src to dest, up to a maximum of
 *      dest_size - 1 characters.
 *      It behaves exactly like strlcpy(3), except that any lower
 *      case characters in the string are converted to upper case.
 *  
 *  Arguments:
 *      src         Pointer to null-terminated string to be copied
 *      dest        Pointer to a character array to receive the copy
 *      dest_size   Size of the destination array
 *
 *  Returns:
 *      Size of the src string.  If this differs from dest_size, then
 *      we knoiw the copy is truncated.
 *
 *  Examples:
 *      char    src[] = "Some text",
 *      dest    [DEST_SIZE + 1];
 *
 *      if ( strlupper(dest, src, DEST_SIZE + 1) != DEST_SIZE + 1 )
 *          fputs("Warning: String truncated.\n", stderr);
 *
 *  See also:
 *      strllower(3), strlcpy(3), strlcat(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-04  Jason Bacon Begin
 ***************************************************************************/

size_t  strlupper(char *dest, const char *src, size_t dest_size)

{
    size_t  c;
    
    for (c = 0; (src[c] != '\0') && (c < dest_size - 1); ++c)
	dest[c] = toupper(src[c]);
    dest[c] = '\0';
    while ( src[c] != '\0' )
	++c;
    return c;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strupper(3) converts lower case characters in a string to upper
 *      case, overwriting the original.  It is functionally equivalent to
 *      strlupper(str, str, strlen(str)).  It is implemented separately for
 *      efficiency, since using strlupper(3) for this purpose requires
 *      knowing or computing the length of the string and passing three
 *      arguments instead of one.
 *  
 *  Arguments:
 *      str         Pointer to null-terminated string to be copied
 *
 *  Returns:
 *      Size of the src string.
 *
 *  Examples:
 *      char    src[] = "Some text",
 *      dest    [DEST_SIZE + 1];
 *
 *      if ( strlupper(dest, src, DEST_SIZE + 1) != DEST_SIZE + 1 )
 *          fputs("Warning: String truncated.\n", stderr);
 *
 *  See also:
 *      strllower(3), strlcpy(3), strlcat(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-04  Jason Bacon Begin
 ***************************************************************************/

size_t  strupper(char *str)

{
    size_t  c;
    
    for (c = 0; str[c] != '\0'; ++c)
	str[c] = toupper(str[c]);
    return c;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strllower(3) copies a string from src to dest, up to a maximum of
 *      dest_size - 1 characters.
 *      It behaves exactly like strlcpy(3), except that any upper
 *      case characters in the string are converted to lower case.
 *  
 *  Arguments:
 *      src         Pointer to null-terminated string to be copied
 *      dest        Pointer to a character array to receive the copy
 *      dest_size   Size of the destination array
 *
 *  Returns:
 *      Size of the src string.  If this differs from dest_size, then
 *      we knoiw the copy is truncated.
 *
 *  Examples:
 *      char    src[] = "Some text",
 *      dest    [DEST_SIZE + 1];
 *
 *      if ( strllower(dest, src, DEST_SIZE + 1) != DEST_SIZE + 1 )
 *          fputs("Warning: String truncated.\n", stderr);
 *
 *  See also:
 *      strllower(3), strlcpy(3), strlcat(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-04  Jason Bacon Begin
 ***************************************************************************/

size_t  strllower(char *dest, const char *src, size_t dest_size)

{
    size_t  c;
    
    for (c = 0; (src[c] != '\0') && (c < dest_size - 1); ++c)
	dest[c] = tolower(src[c]);
    dest[c] = '\0';
    while ( src[c] != '\0' )
	++c;
    return c;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strlower(3) converts upper case characters in a string to lower
 *      case, overwriting the original.  It is functionally equivalent to
 *      strllower(str, str, strlen(str)).  It is implemented separately for
 *      efficiency, since using strllower(3) for this purpose requires
 *      knowing or computing the length of the string and passing three
 *      arguments instead of one.
 *  
 *  Arguments:
 *      str         Pointer to null-terminated string to be copied
 *
 *  Returns:
 *      Size of the src string.
 *
 *  Examples:
 *      char    src[] = "Some text",
 *      dest    [DEST_SIZE + 1];
 *
 *      if ( strllower(dest, src, DEST_SIZE + 1) != DEST_SIZE + 1 )
 *          fputs("Warning: String truncated.\n", stderr);
 *
 *  See also:
 *      strllower(3), strlcpy(3), strlcat(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-04  Jason Bacon Begin
 ***************************************************************************/

size_t  strlower(char *str)

{
    size_t  c;
    
    for (c = 0; str[c] != '\0'; ++c)
	str[c] = tolower(str[c]);
    return c;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Append an argv style list of arguments to a string.  This is
 *      useful for constructing a command to be passed to a shell via
 *      system() or similar methods.
 *
 *  Arguments:
 *      string              String to which argv elements are appended
 *      argv                Character pointer array to a list of elements
 *      first_arg           Index of first argument to append
 *      string_buff_size    Size of string array including null byte
 *
 *  Returns:
 *      Length of string + all argv elements.  If this is greater than
 *      string_buff_size, then the string has been truncated.
 *
 *  Examples:
 *      char    cmd[CMD_MAX + 1] = "ls",
 *              *argv[] = { "-l", NULL };
 *
 *      if ( str_argv_cat(cmd, argv, 0, CMD_MAX + 1) > CMD_MAX + 1 )
 *          fputs("string is truncated.\n", stderr);
 *      else
 *          system(cmd);
 *
 *  See also:
 *      strlcpy(3), strlcat(3), snprintf(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-09-30  Jason Bacon Begin
 ***************************************************************************/

size_t  str_argv_cat(char *string, char *argv[], size_t first_arg,
		     size_t string_buff_size)

{
    size_t  c,
	    len;
    
    len = strlen(string);
    for (c = first_arg; argv[c] != NULL; ++c)
    {
	len += strlen(argv[c]);
	strlcat(string, argv[c], string_buff_size);
	strlcat(string, " ", string_buff_size);
    }
    return len;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strblank() returns true if the null-terminated string contains only
 *      whitespace or nothing at all.  It is the null-terminated string
 *      equivalent of isblank(3), which tests a single character.
 *  
 *  Arguments:
 *      string: A null-terminated string
 *
 *  Returns:
 *      true is string contains only whitespace, or nothing
 *      false if any non-whitespace characters are present
 *
 *  See also:
 *      isblank(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

int     strblank(const char *string)

{
    while ( *string != '\0' )
    {
	if ( !isspace((size_t)(*string)) )   /* Not blank? */
	    return 0;
	++string;
    }
    return 1;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Determine whether a string is a valid integer by attempting to
 *      convert it using strtoll().
 *
 *  Arguments:
 *      string: The string to be tested
 *      base:   The expected base of the integer (usually 8, 10, or 16)
 *
 *  Returns:
 *      Non-zero value if the string represents an integer, zero otherwise
 *
 *  See also:
 *      strtoll(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-24  Jason Bacon Begin
 ***************************************************************************/

int     strisint(const char *string, int base)

{
    char    *end;
    
    strtoll(string, &end, base);
    return *end == '\0';
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Determine whether a string is a valid real number by attempting to
 *      convert it using strtod().
 *
 *  Arguments:
 *      string: The string to be tested
 *
 *  Returns:
 *      Non-zero value if the string represents a real number, zero otherwise
 *
 *  See also:
 *      strtod(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-24  Jason Bacon Begin
 ***************************************************************************/

int     strisreal(const char *string)

{
    char    *end;
    
    strtod(string, &end);
    return *end == '\0';
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strlbasecpy() is a convenience for copying a string to a non-zero
 *      starting position in another string.  The caller provides the address
 *      of the destination to which string should be copied, the base address
 *      of the array containing the destination, and the TOTAL length of the
 *      destination array.  strlbasecpy() then computes the distance from
 *      the destination address to the end of the array and like strlcpy(),
 *      prevents overrun from occurring.
 *  
 *  Arguments:
 *      dest:       Address to which src is copied
 *      dest_base:  Base address of the array containing dest address
 *      src:        Address of null-terminated string to be copied
 *      dstsize:    Size of the dest array
 *
 *  Returns:
 *      The original value of dest
 *
 *  See also:
 *      strlcpy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

char   *strlbasecpy(char *dest, const char *dest_base, const char *src,
		    size_t dstsize)

{
    char        *save_dest;
    const char  *end;

    save_dest = dest;
    dstsize -= dest-dest_base;
    end = src + dstsize;
    while ((*src != '\0') && (src < end - 1))
	*dest++ = *src++;
    *dest = '\0';
    return (save_dest);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Compare two strings via indirect pointers.  This can be used by
 *      qsort(), heapsort(), etc. to sort an argv-style pointer array,
 *      swapping only the pointers rather than the string contents.
 *
 *  Arguments:
 *      p1, p2: Pointers to pointers to the strings to compare
 *
 *  Returns:
 *      0 if the strings are the same, a value < 0 if the string at p1
 *      is lexically < that at p2, a value > 0 otherwise.
 *
 *  See also:
 *      strcmp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-04  Jason Bacon Begin
 ***************************************************************************/

int     strptrcmp(const char **p1, const char **p2)

{
    return strcmp(*p1, *p2);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Compare two strings without regard for case
 *      via indirect pointers.  This can be used by
 *      qsort(), heapsort(), etc. to sort an argv-style pointer array,
 *      swapping only the pointers rather than the string contents.
 *
 *  Arguments:
 *      p1, p2: Pointers to pointers to the strings to compare
 *
 *  Returns:
 *      0 if the strings are the same, a value < 0 if the string at p1
 *      is lexically < that at p2, a value > 0 otherwise.
 *
 *  See also:
 *      strcmp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-04  Jason Bacon Begin
 ***************************************************************************/

int     strptrcasecmp(const char **p1, const char **p2)

{
    return strcmp(*p1, *p2);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strshellcpy() expands a string containing shell meta-characters,
 *      usually a shell command, as a Unix shell would do before execution.
 *      This is useful if you want to avoid spawning an unnecessary shell
 *      process, such as when using fork(2) and exec(3) directly instead
 *      of using system(3).
 *
 *      Currently supports:
 *          ~/: process owner's home directory
 *          $: environment variable
 *  
 *  Arguments:
 *      src:        String containing meta-characters
 *      dest:       Expanded string
 *      dest_len:   Size of destination array
 *
 *  Returns:
 *      0 on success, -1 if expansion did not fit in dest_len
 *
 *  See also:
 *      sh(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

int     strshellcpy(char *dest, const char *src, size_t dest_len)

{
    char    home[PATH_MAX + 1],*p,*val,var[PATH_MAX + 1];
    int     c;
    
    while ( dest_len && (*src != '\0') )
    {
	switch(*src)
	{
	    case    '~':
		++src;
		if ( *src == '/' )  /* Process owner's home dir */
		{
		    xt_get_home_dir(home,PATH_MAX);
		    for (p=home; dest_len-- && (*p != '\0'); )
			*dest++ = *p++;
		}
		else
		{
		    /* Get somebody else's home dir */
		}
		break;
	    case    '$':
		++src;
		/* Get ENV variable name */
		for (p=var, c=0; (c<PATH_MAX) && ISIDENT((size_t)(*src)); ++c)
		    *p++ = *src++;
		*p = '\0';
		
		/* Get value and copy to dest command */
		if ( (val = getenv(var)) != NULL )
		    while ( dest_len-- && (*val != '\0') )
			*dest++ = *val++;
		break;
	    default:
		*dest++ = *src++;
		--dest_len;
	}
    }
    *dest = '\0';
    if ( (dest_len == 0) && (*src != '\0') )
	return -1;
    else
	return 0;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Copy string src to dest, reducing to a maximum length of dstsize if
 *      necessary by replacing the center portion of the string with "...".
 *  
 *  Arguments:
 *      dest:   Pointer to a character array of at least len dstsize, where
 *              the possibly modified string is stored.
 *      src:    Point to a character array containing the string to compact
 *      dstsize: Length of the dest array
 *
 *  Returns:
 *      The length of the original string (like strlcpy(3))
 *
 *  Examples:
 *      #define MAXLEN  127
 *      char    limited_str[MAXLEN + 1],
 *              original_str[SOME_OTHER_LEN + 1];
 *      
 *      strsqueeze(limited_str, original_str, MAXLEN + 1);
 *
 *  See also:
 *      strlcpy(3), strlcat(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-29  Jason Bacon Begin
 ***************************************************************************/

size_t  strsqueeze(char *dest, const char *src, size_t dstsize)

{
    size_t  len = strlen(src),
	    left_len,
	    right_len;
    
    if ( len <= dstsize )
	strlcpy(dest, src, dstsize);
    else
    {
	left_len = (dstsize - 3) / 2;
	right_len = dstsize - left_len - 3;
	memcpy(dest, src, left_len);
	strlcat(dest, "...", dstsize);
	strlcat(dest, src + len - right_len + 1, dstsize);
    }
    return len;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Translate characters in a string similar to tr(1), where each
 *      occurrence of from[i] is replaced with to[i] in string.
 *
 *      Currently only replacement is supported and flags is ignored.
 *      In the future, additional features such as deletion and
 *      compaction may be supported.
 *  
 *  Arguments:
 *      string  The string to be transformed
 *      from    Characters to be replaced with those in to
 *      to      Characters to replace those in from
 *      flags   Bit mask to enable optional features simlar to tr(1)
 *
 *  Examples:
 *      char    string[] = "Hello";
 *
 *      // Convert string to "HELLO"
 *      strtr(string, "elo", "ELO", 0);
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-10  Jason Bacon Begin
 ***************************************************************************/

void    strtr(char *string, const char *from, const char *to, int flags)

{
    char    *p,
	    *i;
    
    for (p = string; *p != '\0'; ++p)
    {
	i = strchr(from, *p);
	if ( i != NULL )
	{
	    //fprintf(stderr, "Replacing %c with %c\n", *p, to[i - from]);
	    //fflush(stderr);
	    *p = to[i - from];
	}
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      Trim unwanted characters off both ends of a string.  Typically
 *      this is whitespace added to a comma-separated file or similar.
 *
 *  Arguments:
 *      string  The string to be trimmed
 *      fat     A string containing a list of characters to be removed
 *
 *  See also:
 *      strsep(3)
 *
 *  Examples:
 *      char    string[] = "  Alfred E. Neumann."
 *
 *      strtrim(string, " .");
 *      puts(string);
 *
 *      Output is "Aldred E. Neumann"
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-09-24  Jason Bacon Begin
 ***************************************************************************/

void    strtrim(char *string, const char *fat)

{
    char    *start, *end;
    
    for (start = string; (*start != '\0') && strchr(fat, *start); ++start)
	;
    for (end = start; *end != '\0'; ++end)
	;
    while ( (end >= string) && strchr(fat, *end) )
	--end;
    end[1] = '\0';
    if ( (start > string) && (end > start) )
	memmove(string, start, end - start + 2);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      strviscpy() copies a string from src to dest, converting invisible
 *      characters to a visible format much like the vis command or cat -v.
 *  
 *  Arguments:
 *      src:    Source string containing invisible characters
 *      dest:   Destination array to receive modified string
 *      maxlen: Maximum number of characters in dest, not including null byte
 *
 *  Returns:
 *      A pointer to dest
 *
 *  See also:
 *      vis(1), cat(1), strlcpy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

char    *strviscpy(unsigned char *dest, const unsigned char *src,
		size_t maxlen)

{
    char    *d = (char *)dest;

    if ( (src == NULL) || (dest == NULL) )
	return NULL;
    
    while ( (*src != '\0') && (maxlen > 0) )
    {
	if ( ((unsigned char)*src < 128) && isgraph(*src) )
	{
	    *d++ = *src++;
	    --maxlen;
	}
	else
	{
	    if ( maxlen > 4 )
	    {
		snprintf(d,maxlen,"\\%03o",(unsigned char)*src);
		++src;
		d+=4;
		maxlen-=4;
	    }
	}
    }
    *d = '\0';
    return (char *)dest;
}


/***************************************************************************
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      ltostrn() is a small, fast integer to string converter that can
 *      convert using any base from 2 to 36.  It is the converse of strtol(3).
 *      The size of the char buffer passed should be 1 more than maxlen to
 *      allow for a null byte.
 *  
 *  Arguments:
 *      string: char array to receive ascii text
 *      val:    value to convert to ascii
 *      base:   number base for conversion, must be between 2 and 36
 *      maxlen: size of 'string' array - 1 (to account for null byte)
 *
 *  Returns:
 *      A pointer to the converted string, or NULL if the string buffer was
 *      not big enough for all the digits
 *
 *  See also:
 *      strtol(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

char    *ltostrn(char string[], long val, unsigned base, size_t maxlen)

{
    char    temp[maxlen+1], *p, *s = string;
    int     digit;
    
    if ( base > 36 )
	return NULL;
    
    /* Tack on - sign if negative */
    if ( val < 0 )
    {
	*s++ = '-';
	val = -val;
    }
    
    /* Convert val to ascii digits (in reverse order) */
    for (p=temp; (val > 0) && (maxlen > 0); val /= base, --maxlen)
    {
	digit = val % base;
	if ( digit < 10 )
	    *p++ = digit+'0';
	else
	    *p++ = digit-10+'a';
    }
    
    /* Reverse digits */
    while ( p > temp )
	*s++ = *(--p);
    *s = '\0';
    if ( val > 0 )
	return NULL;
    else
	return string;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      .B str2u64()
 *      is a super-fast hash function that converts a string of 8 or fewer
 *      characters to a 64-bit integer.  Strings of more than 8 characters may
 *      also be hashed, though collisions will occur (same hash value for more
 *      than one string) if the first 8 characters are the same.
 *
 *      This sort of hashing is useful for storing lists
 *      of very short strings, as it eliminates the need to use strdup(),
 *      strlcpy(), and strcmp() for processing.  Strings can be compared
 *      for equality using a straight integer comparison.  Strings of 7
 *      or fewer characters can still be accessed as a string by simply
 *      casting to char * for output, then using lexical comparison with strcmp(),
 *      etc.  A string of 8 characters will not have a null-terminator.
 *
 *      The value returned varies depending on endianness.  Hence, hash
 *      values generated on one architecture will need to be byte swapped
 *      before comparison to values generated under a different endianness.
 *  
 *  Arguments:
 *      str     String to convert
 *
 *  Returns:
 *      uint64_t integer containing the characters in str
 *
 *  Examples:
 *      char        *s1 = "hello!", s2 = "Hello!";
 *      uint64_t    v1, v2;
 *      
 *      v1 = str2u64(s1);
 *      v2 = str2u64(s2);
 *      if ( v1 != v2 )
 *          printf("%s and %s are different.\n", (char *)&v1, (char *)&v2);
 *
 *  See also:
 *      strdup(3), strcmp(3), strlcpy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-02  Jason Bacon Begin
 ***************************************************************************/

uint64_t    str2u64(const char *str)

{
    size_t      c;
    char        *p;
    uint64_t    v = 0;
    
    for (c = 0, p = (char *)&v; (c < sizeof(v)) && (str[c] != '\0'); ++c, ++p)
	*p = str[c];
    return v;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/string.h>
 *      -lxtend
 *
 *  Description:
 *      .B strsplit()
 *      splits a string into tokens separated by any character
 *      in the string argument sep.
 *      The function interface is similar to split() in awk, except that
 *      sep is a simple list of characters rather than a regular expression.
 *
 *      The array argument should be the address of a char ** variable.
 *      strsplit() allocates memory for the pointers as needed and
 *      assigns one token to each pointer.
 *
 *      strsplit() should only be used when an array of strings
 *      representing the tokens is actually needed.  In cases where each
 *      token can be immediately processed and forgotten, use a loop with
 *      strsep().  Introducing arrays into a program unnecessarily should
 *      avoided as a habit to maximize speed and minimize memory use.
 *
 *      Caution: strsplit() is destructive: It replaces the separators
 *      in string with null bytes.  To preserve the original string,
 *      duplicate it with strdup() first and pass the copy to strsplit().
 *  
 *  Arguments:
 *      string  String to be parsed for tokens
 *      array   Pointer array to be filled with tokens
 *      sep     Character string listing all recognized separators
 *
 *  Returns:
 *      The number of tokens into which string is separated, or 0 if
 *      a memory allocation or other failure occurred.
 *
 *  Examples:
 *      char    *string = "1,2,3,4,5", *copy, **array;
 *      size_t  c, tokens;
 *
 *      copy = strdup(string);
 *      tokens = strsplit(copy, &array, ",");
 *      for (int c = 0; c < tokens; ++c)
 *          puts(array[c]);
 *
 *  See also:
 *      strsep(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-12  Jason Bacon Begin
 ***************************************************************************/

int     strsplit(char *string, char ***array, const char *sep)

{
    size_t  c,
	    array_size = 64;

    if ((*array = xt_malloc(array_size, sizeof(*array))) == NULL )
    {
	fprintf(stderr, "strsplit(): malloc() failed.\n");
	return 0;
    }
    
    for (c = 0; ((*array)[c] = strsep(&string, sep)) != NULL; )
    {
	if ( ++c == array_size )
	{
	    *array = xt_realloc(*array, array_size *= 2, sizeof(*array));
	    if ( *array == NULL )
	    {
		fprintf(stderr, "strsplit(): malloc() failed.\n");
		return 0;
	    }
	}
    }
    
    // Trim
    *array = xt_realloc(*array, c, sizeof(*array));
    return c;
}

/***************************************************************************
 *  Library:
 *      #include <xtend/time.h>
 *      -lextend
 *
 *  Description:
 *      difftimeofday() returns the difference, in microseconds, between two 
 *      time values returned by gettimeofday(3).  This function can be used
 *      to get a good estimate of the real time elapsed in a process between
 *      any two points (where calls to gettimeofday(3) are strategically
 *      placed.)
 *
 *      Use of these functions should have minimal impact on run time,
 *      unless called many times to measure time of a function with a very
 *      short run time.
 *  
 *  Arguments:
 *      later, earlier: timeval structures populated by gettimeofday(3)
 *
 *  Returns:
 *      The difference between the two times in microseconds
 *
 *  See also:
 *      gettimeofday(2), xt_tic(3), xt_toc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

time_t  xt_difftimeofday(struct timeval *later, struct timeval *earlier)

{
    return 1000000 * (later->tv_sec - earlier->tv_sec) +
	    (later->tv_usec - earlier->tv_usec);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/time.h>
 *      -lxtend
 *
 *  Description:
 *      xt_tic() records the current time in a struct timeval structure.
 *      It is a simple wrapper around gettimeofday(2) meant for use with
 *      xt_toc(3), which reports elapsed time since the xt_tic() call.
 *
 *      The xt_tic() and xt_toc() functions are used to accurately determine
 *      the elapsed time of a segment of code, such as a loop that is
 *      suspected to be costly.  xt_tic() is inserted into the program just
 *      before the code and xt_toc() immediately after.
 *  
 *  Arguments:
 *      start_time  A struct timeval structure populated by xt_tic()
 *
 *  Returns:
 *      The exit status of gettimeofday(2)
 *
 *  Examples:
 *      struct timeval  start_time;
 *      struct rusage   start_usage;
 *
 *      xt_tic(&start_time, &start_usage);
 *      // Code for which elapsed time is to be measured
 *      for (c = 0; c < bignum; ++c)
 *      {
 *          ...
 *      }
 *      xt_toc(stderr, "Elapsed time for loop:\n", &start_time, &start_usage);
 *
 *  See also:
 *      xt_toc(3), difftimeofday(3), gettimeofday(2)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-20  Jason Bacon Begin
 ***************************************************************************/

int     xt_tic(struct timeval *start_time, struct rusage *start_usage)

{
    getrusage(RUSAGE_SELF, start_usage);
    return gettimeofday(start_time, NULL);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/time.h>
 *      -lxtend
 *
 *  Description:
 *      xt_toc() reports the elapsed time, user time, and system time
 *      since start_time and start_usage, which should have been populated
 *      by xt_tic(3) at the beginning of the interval being measured.
 *      Time is reported in microseconds, and if elapsed time is greater
 *      than one second, days, hours, and seconds are also reported.
 *
 *      The xt_tic() and xt_toc() functions are used to accurately determine
 *      the elapsed time of a segment of code, such as a loop that is
 *      suspected to be costly.  xt_tic() is inserted into the program just
 *      before the code and xt_toc() immediately after.
 *  
 *  Arguments:
 *      stream      FILE stream to which output is sent
 *      message     Optional message to print before time stats, or NULL
 *      start_time  A struct timeval structure populated by xt_tic()
 *
 *  Returns:
 *      The time difference in microseconds
 *
 *  Examples:
 *      struct timeval  start_time;
 *      struct rusage   start_usage;
 *
 *      xt_tic(&start_time, &start_usage);
 *      // Code for which elapsed time is to be measured
 *      for (c = 0; c < bignum; ++c)
 *      {
 *          ...
 *      }
 *      xt_toc(stderr, "Elapsed time for loop:\n", &start_time, &start_usage);
 *
 *  See also:
 *      xt_tic(3), difftimeofday(3), gettimeofday(2)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-20  Jason Bacon Begin
 ***************************************************************************/

unsigned long xt_toc(FILE *stream, const char *message,
		     struct timeval *start_time, struct rusage *start_usage)

{
    struct timeval  end_time;
    struct rusage   end_usage;
    unsigned long   diff, hours, minutes, seconds;
    
    if ( message != NULL )
	fputs(message, stream);
    gettimeofday(&end_time, NULL);
    diff = xt_difftimeofday(&end_time, start_time);
    fprintf(stream, "Elapsed time     = %10lu microseconds", diff);
    if ( diff >= 1000000 )
    {
	seconds = diff / 1000000;
	minutes = seconds / 60;
	hours = minutes / 24;
	fprintf(stream, " (%lu hours, %lu minutes, %lu seconds)",
		hours, minutes, seconds);
    }
    putc('\n', stream);
    getrusage(RUSAGE_SELF, &end_usage);
    fprintf(stream, "User time        = %10lu microseconds\n",
	    end_usage.ru_utime.tv_sec * 1000000 + end_usage.ru_utime.tv_usec -
	    (start_usage->ru_utime.tv_sec * 1000000 + start_usage->ru_utime.tv_usec));
    fprintf(stream, "Sys time         = %10lu microseconds\n",
	    end_usage.ru_stime.tv_sec * 1000000 + end_usage.ru_stime.tv_usec -
	    (start_usage->ru_stime.tv_sec * 1000000 + start_usage->ru_stime.tv_usec));
    return diff;
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      Verify that the filename extension on filename is either the
 *      valid extension provided or that extension followed by a
 *      compression extension, e.g. .gz, .bz2, .xz.
 *
 *  Arguments:
 *      filename:   Name of the file to be checked
 *      valid_ext:  Valid extension, not include .gz, .bz2, or .xz
 *
 *  Returns:
 *      true (defined in stdbool.h) if filename has a valid extension
 *      false if not
 *
 *  See also:
 *      gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-05-04  Jason Bacon Begin
 ***************************************************************************/

bool    xt_valid_extension(const char *filename, const char *valid_ext)

{
    char    *zip_exts[] = { ".gz", ".bz2", ".xz" },
	    *ext,
	    *compressed;
    size_t  c;

    if ( (ext = strrchr(filename, '.')) != NULL )
    {
	if ( strcmp(ext, valid_ext) == 0 )
	    return true;
	for (c = 0; c < sizeof(zip_exts) / sizeof(*zip_exts); ++c)
	{
	    if ( strcmp(ext, zip_exts[c]) == 0 )
	    {
		compressed = strdup(filename);
		// Already confirmed there's an extension, clip it
		*strrchr(compressed, '.') = '\0';
		if ( ((ext = strrchr(compressed, '.')) != NULL) &&
		     (strcmp(ext, valid_ext) == 0) )
		{
		    free(compressed);
		    return true;
		}
		free(compressed);
		break;
	    }
	}
    }
    fprintf(stderr, "Error: %s should have a %s[.%s] extension\n",
	    filename, valid_ext, "gz|bz2|xz");
    return false;
}
#include <stdio.h>
#include <stdarg.h>
#include <stdlib.h>
#include <sysexits.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/proc.h>
 *      -lxtend
 *
 *  Description:
 *      va_usage() is a simple convenience function that takes a
 *      printf-style variable argument list, prints a message to stderr,
 *      and terminates the proces with an exit status of EX_USAGE.
 *      The message should indicate correct command-line usage of the
 *      calling program as conventional for Unix commands.
 *  
 *  Arguments:
 *      format_string:  printf-style format string
 *      Additional arguments to match placeholders in format_string
 *
 *  Returns:
 *      Does not return, terminates the calling process
 *
 *  See also:
 *      printf(3), vfprintf(3), exit(3), sysexits(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

void    va_usage(const char *format_string, ...)

{
    va_list list;
    char    new_format[XT_FORMAT_MAX_CHARS + 1];
    
    va_start(list,format_string);
    snprintf(new_format, XT_FORMAT_MAX_CHARS, "Usage: %s", format_string);
    vfprintf(stderr,new_format,list);
    exit(EX_USAGE);
}
#include <string.h>
#include <sysexits.h>
#include <sys/stat.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      Open a raw data file using fopen() or a gzipped, bzipped, or
 *      xzipped file using popen().  Must be used in conjunction with
 *      xt_fclose() to ensure that fclose() or pclose() is called where
 *      appropriate.
 *
 *  Arguments:
 *      filename:   Name of the file to be opened
 *      mode:       "r" or "w", passed to fopen() or popen()
 *
 *  Returns:
 *      A pointer to the FILE structure or NULL if open failed
 *
 *  See also:
 *      fopen(3), popen(3), gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

FILE    *xt_fopen(const char *filename, const char *mode)

{
    char    *ext = strrchr(filename, '.'),
	    cmd[XT_CMD_MAX_CHARS + 1];
    
    if ( (strcmp(mode, "r") != 0 ) && (strcmp(mode, "w") != 0) )
    {
	fprintf(stderr, "xt_fopen(): Only \"r\" and \"w\" modes supported.\n");
	return NULL;
    }
    
    if ( ext == NULL )
    {
	fprintf(stderr, "xt_fopen(): No filename extension on %s.\n", filename);
	return NULL;
    }

    if ( *mode == 'r' )
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
// Big Sur zcat requires a .Z extension and CentOS 7 lacks gzcat
#ifdef __APPLE__
	    snprintf(cmd, XT_CMD_MAX_CHARS, "gzcat %s", filename);
#else
	    snprintf(cmd, XT_CMD_MAX_CHARS, "zcat %s", filename);
#endif
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "bzcat %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "xzcat %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
    else    // "w"
    {
	if ( strcmp(ext, ".gz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "gzip -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".bz2") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "bzip2 -c > %s", filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".xz") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "xz -c > %s", filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
}


/***************************************************************************
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      Close a FILE stream with fclose() or pclose() as appropriate.
 *      Automatically determines the proper close function to call using
 *      S_ISFIFO on the stream stat structure.
 *
 *  Arguments:
 *      stream: The FILE structure to be closed
 *
 *  Returns:
 *      The value returned by fclose() or pclose()
 *
 *  See also:
 *      fopen(3), popen(3), gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-10  Jason Bacon Begin
 ***************************************************************************/

int     xt_fclose(FILE *stream)

{
    struct stat stat;
    
    fstat(fileno(stream), &stat);
    if ( S_ISFIFO(stat.st_mode) )
	return pclose(stream);
    else
	return fclose(stream);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      .B xt_inhale_strings()
 *      reads a list of strings from a file, one per line, into a pointer
 *      array.  Memory is allocated for the pointer array and for each
 *      string.
 *
 *      Memory should be freed using xt_free_strings(3) as soon as the
 *      strings are no longer needed.
 *
 *      Inhaling large amounts of data into arrays should generally be
 *      avoided in favor of more memory-efficient use-once-and-discard
 *      strategies, but may be advantageous for small lists of strings
 *      accessed repeatedly, or necessary for a few tasks such as sorting.
 *  
 *  Arguments:
 *      stream  FILE * from which strings are read, one per line
 *      list    Pointer to a char ** (poiner array), populated with strings
 *
 *  Returns:
 *      The number of strings read, XT_READ_IO_ERR on read error
 *
 *  Examples:
 *      FILE    *instream;
 *      char    **strings;
 *      ssize_t string_count;
 *
 *      string_count = xt_inhale_strings(instream, &strings);
 *      ...
 *      xt_free_strings(strings);
 *
 *  See also:
 *      xt_free_strings(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-21  Jason Bacon Begin
 ***************************************************************************/

ssize_t xt_inhale_strings(FILE *stream, char ***list)

{
    size_t  list_size = 1024,
	    c,
	    buff_size,
	    len;
    char    *temp;
    
    if ( (*list = (char **)xt_malloc(list_size, sizeof(*list))) == NULL )
    {
	fprintf(stderr, "load_strings(): Unable to allocate list.\n");
	return EX_UNAVAILABLE;
    }
    
    buff_size = 0;  // Make xt_read_line_malloc() allocate a new string
    for (c = 0; xt_read_line_malloc(stream, &temp, &buff_size, &len) != EOF; ++c)
    {
	if ( c == list_size - 1 )
	{
	    list_size *= 2;
	    if ( (*list = (char **)xt_realloc(*list, list_size, sizeof(*list))) == NULL )
	    {
		fprintf(stderr, "load_strings(): Unable to reallocate list.\n");
		return EX_UNAVAILABLE;
	    }
	}
	(*list)[c] = temp;
	buff_size = 0;  // Make xt_read_line_malloc() allocate a new string
    }
    (*list)[c] = NULL;  // So xt_free_strings() doesn't need count
    return c;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/file.h>
 *      -lxtend
 *
 *  Description:
 *      .B xt_read_line_malloc()
 *      reads a single line of text (up to the next newline or EOF)
 *      from stream, allocating and/or extending the provided buffer if
 *      needed.
 *
 *      The buff_size argument must be initilized to 0 if buff has
 *      not been previously allocated.  This will cause an initial
 *      allocation to occur.  If buff has been previously allocated,
 *      the buff_size must accurately reflect the allocated memory size.
 *      This will happen naturally when reusing buff in a loop, as shown
 *      in the example below.
 *  
 *  Arguments:
 *      stream:     FILE stream from which field is read
 *      buff:       Character buffer into which field is copied
 *      buff_size:  Size of the array passed to buff
 *      len:        Pointer to a variable which will receive the field length
 *
 *  Returns:
 *      Delimiter ending the read: either newline or EOF
 *
 *  Examples:
 *      FILE    *stream;
 *      char    *buff;
 *      size_t  buff_len, len;
 *
 *      // Reuse buff to minimize malloc() calls.  buff will be extended
 *      // as needed when longer strings are read.  Initialize buff here
 *      // rather than above for the most cohesive code.
 *      buff_len = 0;
 *      while ( ffile_read_line_malloc(stream, buff, &buff_len, &len) != EOF )
 *      {
 *      }
 *
 *  See also:
 *      dsv_read_field_malloc(3), ffgetc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-20  Jason Bacon Begin
 ***************************************************************************/

int     xt_read_line_malloc(FILE *stream, char **buff, size_t *buff_size,
			    size_t *len)

{
    size_t  c;
    int     ch;
    
    if ( *buff_size == 0 )
    {
	*buff_size = 1024;
	*buff = xt_malloc(*buff_size, sizeof(**buff));
	if ( *buff == NULL )
	    return XT_MALLOC_FAILED;
    }
    
    for (c = 0; ( ((ch = getc(stream)) != '\n') && (ch != EOF) ); ++c)
    {
	if ( c == *buff_size - 1 )
	{
	    *buff_size *= 2;
	    *buff = xt_realloc(*buff, *buff_size, sizeof(**buff));
	    if ( *buff == NULL )
		return XT_MALLOC_FAILED;
	}
	(*buff)[c] = ch;
    }
    (*buff)[c] = '\0';
    //fprintf(stderr, "buff = %s\n", *buff);
    *len = c;

    /* Trim array */
    if ( *buff_size != c + 1 )
    {
	*buff_size = c + 1;
	*buff = xt_realloc(*buff, *buff_size, sizeof(**buff));
    }
    //fprintf(stderr, "Returning %d\n", ch);
    return ch;
}
#include <stdlib.h>

/***************************************************************************
 *  Library:
 *      #include <xtend/mem.h>
 *      -lxtend
 *
 *  Description:
 *      xt_malloc() is a simple wrapper around malloc(3) that requires two
 *      arguments representing the number of objects to allocate and the
 *      size of an element.  This prevents the very common mistake with
 *      malloc(3) of forgetting to multiply by the size of an element.
 *
 *      Specifying the size using sizeof(*variable) has the advantage of
 *      being type-independent.  I.e. if you change the type of the variable,
 *      this code need not be updated.  Simply add one * to whatever
 *      the return value is assigned to.
 *  
 *  Arguments:
 *      nelem:  Number of objects to allocate
 *      size:   Size of a single object
 *
 *  Examples:
 *      size_t      widget_list_size = 1024;
 *      widget_t    *widgets;
 *
 *      widgets = xt_malloc(widget_list_size, sizeof(*widgets));
 *
 *  Returns:
 *      Address of the newly allocated array, or NULL if allocation failed
 *
 *  See also:
 *      malloc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

void    *xt_malloc(size_t nelem, size_t size)

{
    return  malloc(nelem * size);
}


/***************************************************************************
 *  Library:
 *      #include <xtend/mem.h>
 *      -lxtend
 *
 *  Description:
 *      xt_realloc() is a simple wrapper around realloc(3) that requires three
 *      arguments representing the original array, the new number of objects
 *      to allocate and the size of an element.  This prevents the very
 *      common mistake with realloc(3) of forgetting to multiply by the size
 *      of an element.
 *
 *      Specifying the size using sizeof(*variable) has the advantage of
 *      being type-independent.  I.e. if you change the type of the variable,
 *      this code need not be updated.  Simply add one * to whatever
 *      the return value is assigned to.
 *  
 *  Arguments:
 *      array:  Address of the previously allocated array
 *      nelem:  Number of objects to allocate
 *      size:   Size of a single object
 *
 *  Examples:
 *      size_t      widget_list_size = 1024;
 *      widget_t    *widgets;
 *
 *      widgets = xt_malloc(widget_list_size, sizeof(*widgets));
 *      ...
 *      widget_list_size *= 2;
 *      widgets = xt_realloc(widgets, widget_list_size, sizeof(*widgets));
 *
 *  Returns:
 *      Address of the newly allocated array, or NULL if allocation failed
 *
 *  See also:
 *      realloc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  Circa 1990  Jason Bacon Begin
 ***************************************************************************/

void    *xt_realloc(void *array, size_t nelem, size_t size)

{
    return  realloc(array, nelem * size);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/mem.h>
 *      -lxtend
 *
 *  Description:
 *      .B xt_free_strings(3)
 *      frees all memory for an argv-style pointer array of strings,
 *      such as those allocated by xt_inhale_strings(3).  The pointer
 *      array itself and each string it points to must have been
 *      dynamically allocated with malloc(3), strdup(3), or similar.
 *  
 *  Arguments:
 *      list    A dynamically allocated char ** array
 *
 *  Examples:
 *      FILE    *instream;
 *      char    **strings;
 *      ssize_t string_count;
 *
 *      string_count = xt_inhale_strings(instream, &strings);
 *      ...
 *      xt_free_strings(strings);
 *
 *  See also:
 *      xt_inhale_strings(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-21  Jason Bacon Begin
 ***************************************************************************/

void    xt_free_strings(char **list)

{
    size_t  c;
    
    for (c = 0; list[c] != NULL; ++c)
	free(list[c]);
    free(list);
}
#include <string.h>
#include <stdlib.h>

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <xtend/stdlib.h>
 *      -lxtend
 *
 *  Description:
 *      Shuffle an array of objects using the Fisher-Yates method, which
 *      ensures equal probability of all arrangements.
 *  
 *  Arguments:
 *      base    Base address of the array (address of the first element)
 *      nelem   Number of elements in the array
 *      size    Size of one element
 *
 *  Examples:
 *      type_t  *list;
 *      size_t  list_size = 100;
 *
 *      if ( (list = xt_malloc(list_size, sizeof(*list))) == NULL)
 *      {
 *          fprintf(stderr, "xt_malloc() failed.\n");
 *          exit(EX_UNAVAILABLE);
 *      }
 *      ...
 *      xt_shuffle(list, list_size, sizeof(*list));
 *
 *  See also:
 *      qsort(3), heapsort(3), mergesort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-10-23  Jason Bacon Begin
 ***************************************************************************/

void    xt_shuffle(void *base, size_t nelem, size_t size)

{
    size_t  c, c1;
    char    temp[size];
    
    // Fisher-Yates shuffle
    for (c = 0; c < nelem - 1; ++c)
    {
	// c <= c1 < pair_nelem
	// Yes, we want the possibility of swapping an element with itself
	// Leaving it in place should have the same probability as
	// every other possibility
	c1 = c + random() % (nelem - c);
	memcpy((void *)temp, base + c * size, size);
	memcpy(base + c * size, base + c1 * size, size);
	memcpy(base+ c1 * size, (void *)temp, size);
    }
}
