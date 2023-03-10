#ifndef _XTEND_COMMON_H_
#define _XTEND_COMMON_H_

#define XT_FORMAT_MAX_CHARS 4096
#define XT_CMD_MAX_CHARS    4096

#define XT_OK                   0
// FIXME: Return this instead of EOF in dsv_read*()
// Don't trust that EOF is -1 on all platforms
#define XT_READ_EOF             -1
#define XT_READ_BUFF_OVERFLOW   -2
#define XT_READ_IO_ERR          -3

#define XT_FAIL                 -4
#define XT_MALLOC_FAILED        -5

#endif // _XTEND_COMMON_H
#ifndef _XTEND_CTYPE_H_
#define _XTEND_CTYPE_H_

#define ISIDENT(c)  ( isalnum(c) | ((c)=='_') )

#endif  // _XTEND_CTYPE_H_
    
/*
 *  Generated by /usr/local/bin/auto-gen-get-set
 *
 *  Accessor macros.  Use these to access structure members from functions
 *  outside the dsv_line_t class.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
 */

#define DSV_LINE_ARRAY_SIZE(ptr)        ((ptr)->array_size)
#define DSV_LINE_NUM_FIELDS(ptr)        ((ptr)->num_fields)
#define DSV_LINE_FIELDS(ptr)            ((ptr)->fields)
#define DSV_LINE_FIELDS_AE(ptr,c)       ((ptr)->fields[c])
#define DSV_LINE_DELIMS(ptr)            ((ptr)->delims)
#define DSV_LINE_DELIMS_AE(ptr,c)       ((ptr)->delims[c])
#ifndef _XTEND_DSV_H_
#define _XTEND_DSV_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#define DSV_DATA_OK             0
#define DSV_DATA_INVALID        -1      // Catch-all for non-specific error
#define DSV_DATA_OUT_OF_RANGE   -2

#define DSV_INIT                { 0, 0, NULL, NULL }
#define DSV_FIELD_MAX_CHARS     32767

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/dsv.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      .B dsv_line_t
 *      is a generic structure for holding a line of data from a delimiter
 *      separated file such as CSV or TSV.
 *  
 *  Examples:
 *
 *  See also:
 *      dsv_line_read(3), dsv_line_write(3), dsv_line_copy(3), dsv_line_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-08  Jason Bacon Begin
 ***************************************************************************/

typedef struct
{
    size_t      array_size,
		num_fields;
    char        **fields,
		*delims;
}   dsv_line_t;

#define DSV_LINE_ARRAY_SIZE(ptr)        ((ptr)->array_size)
#define DSV_LINE_NUM_FIELDS(ptr)        ((ptr)->num_fields)
#define DSV_LINE_FIELDS(ptr)            ((ptr)->fields)
#define DSV_LINE_FIELDS_AE(ptr,c)       ((ptr)->fields[c])
#define DSV_LINE_DELIMS(ptr)            ((ptr)->delims)
#define DSV_LINE_DELIMS_AE(ptr,c)       ((ptr)->delims[c])

/* dsv.c */
int dsv_read_field(FILE *stream, char buff[], size_t buff_size, const char *delims, size_t *len);
int dsv_read_field_malloc(FILE *stream, char **buff, size_t *buff_size, const char *delims, size_t *len);
int dsv_skip_field(FILE *stream, const char *delims, size_t *len);
int dsv_skip_rest_of_line(FILE *stream);
void dsv_line_init(dsv_line_t *dsv_line);
int dsv_line_read(dsv_line_t *dsv_line, FILE *stream, const char *delims);
int dsv_line_write(dsv_line_t *dsv_line, FILE *stream);
int dsv_line_copy(dsv_line_t *dest, dsv_line_t *src);
int dsv_line_free(dsv_line_t *dsv_line);
int tsv_read_field(FILE *stream, char buff[], size_t buff_size, size_t *len);
int tsv_read_field_malloc(FILE *stream, char **buff, size_t *buff_size, size_t *len);
int tsv_skip_field(FILE *stream, size_t *len);
int tsv_skip_rest_of_line(FILE *stream);
int csv_read_field(FILE *stream, char buff[], size_t buff_size, size_t *len);
int csv_read_field_malloc(FILE *stream, char **buff, size_t *buff_size, size_t *len);
int csv_skip_field(FILE *stream, size_t *len);
int csv_skip_rest_of_line(FILE *stream);

#endif  // _XTEND_DSV_H_

/*
 *  Generated by /usr/local/bin/auto-gen-get-set
 *
 *  Mutator functions for setting with no sanity checking.  Use these to
 *  set structure members from functions outside the dsv_line_t
 *  class.  These macros perform no data validation.  Hence, they achieve
 *  maximum performance where data are guaranteed correct by other means.
 *  Use the mutator functions (same name as the macro, but lower case)
 *  for more robust code with a small performance penalty.
 *
 *  These generated macros are not expected to be perfect.  Check and edit
 *  as needed before adding to your code.
 */

/* temp-dsv-mutators.c */
int dsv_line_set_array_size(dsv_line_t *dsv_line_ptr, size_t new_array_size);
int dsv_line_set_num_fields(dsv_line_t *dsv_line_ptr, size_t new_num_fields);
int dsv_line_set_fields(dsv_line_t *dsv_line_ptr, char **new_fields);
int dsv_line_set_fields_ae(dsv_line_t *dsv_line_ptr, size_t c, char *new_fields_element);
int dsv_line_set_fields_cpy(dsv_line_t *dsv_line_ptr, char **new_fields, size_t array_size);
int dsv_line_set_delims(dsv_line_t *dsv_line_ptr, char *new_delims);
int dsv_line_set_delims_ae(dsv_line_t *dsv_line_ptr, size_t c, char new_delims_element);
int dsv_line_set_delims_cpy(dsv_line_t *dsv_line_ptr, char *new_delims, size_t array_size);

/* Return values for mutator functions */
#define BL_DSV_DATA_OK              0
#define BL_DSV_DATA_INVALID         -1      // Catch-all for non-specific error
#define BL_DSV_DATA_OUT_OF_RANGE    -2

#ifndef _XTEND_FAST_FILE_H_
#define _XTEND_FAST_FILE_H_

#ifndef _FCNTL_H_
#include <fcntl.h>
#endif

#ifndef _UNISTD_H_
#include <unistd.h>
#endif

#ifndef _XT_COMMON_H_
#endif

/*
 *  These macro implementations of ffgetc() and ffputc() show significantly
 *  lower CPU usage.
 */

#define FFGETC(st) \
    ((st)->c == (st)->bytes_read ? \
	((st)->bytes_read = read((st)->fd, (st)->start, (st)->block_size)) == 0 ? \
	    EOF \
	: ((st)->c = 0, (st)->start[(st)->c++]) \
    : (st)->start[(st)->c++])

#define FFPUTC(ch, st) \
    ((st)->c == (st)->block_size ? \
	write((st)->fd, (st)->start, (st)->block_size) != (st)->block_size ? \
	    EOF \
	: ((st)->c = 0, (st)->start[(st)->c++] = ch) \
    : ((st)->start[(st)->c++] = ch))

#define XT_FAST_FILE_UNGETC_MAX 64L
#define XT_FAST_FILE_MAX_ARGS   128

#define FFILE_INIT  { NULL, NULL, 0, 0, 0, 0, 0, 0, 0 }

typedef struct
{
    unsigned char   *buff;
    unsigned char   *start;
    ssize_t         bytes_read;
    ssize_t         c;
    ssize_t         block_size;
    ssize_t         buff_size;
    int             fd;
    int             flags;
    pid_t           child_pid;
}   ffile_t;

/* fast-file.c */
ffile_t *ff_init_stream(ffile_t *stream);
ffile_t *ffopen(const char *filename, int flags);
ffile_t *ffdopen(int fd, int flags);
int ffgetc(ffile_t *stream);
int ffputc(int ch, ffile_t *stream);
int ffclose(ffile_t *stream);
int ffungetc(int ch, ffile_t *stream);
ffile_t *ffstdin(void);
ffile_t *ffstdout(void);
ffile_t *ffpopen(const char *cmd, int flags);
int ffpclose(ffile_t *stream);
ffile_t *xt_ffopen(const char *filename, int flags);
int xt_ffclose(ffile_t *stream);
int ffprintf(ffile_t *stream, const char *format, ...);
int ffread_line_malloc(ffile_t *stream, char **buff, size_t *buff_size, size_t *len);
int ffputs(const char *string, ffile_t *stream);
char *ffgets(char *string, size_t size, ffile_t *stream);

#endif  // _XTEND_FAST_FILE_H_
#ifndef _XTEND_FILE_H_
#define _XTEND_FILE_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _SYS_STAT_H_
#include <sys/stat.h>   // mode_t on Darwin
#endif

#ifndef __bool_true_false_are_defined
#include <stdbool.h>
#endif

#ifndef _XT_COMMON_H_
#endif

// Added 2022-02-03
#define fgetline(fp, buff, maxlen) \
	    _Pragma("message(\"fgetline() is deprecated.  Use xt_fgetline().\")") \
	    xt_fgetline(fp, buff, maxlen)

#define valid_extension(filename, valid_ext) \
	    _Pragma("message(\"valid_extension() is deprecated.  Use xt_valid_extension().\")") \
	    xt_valid_extension(filename, valid_ext)

#define fast_cp(source, dest) \
	    _Pragma("message(\"fast_cp() is deprecated.  Use xt_fast_cp().\")") \
	    xt_fast_cp(source, dest)

#define file_mod_cmp(file1, file2) \
	    _Pragma("message(\"file_mod_cmp() is deprecated.  Use xt_file_mod_cmp().\")") \
	    xt_file_mod_cmp(file1, file2)

#define fd_purge(fd) \
	    _Pragma("message(\"fd_purge() is deprecated.  Use xt_fd_purge().\")") \
	    xt_fd_purge(fd)

#define get_home_dir(dir, maxlen) \
	    _Pragma("message(\"get_home_dir() is deprecated.  Use xt_get_home_dir().\")") \
	    xt_get_home_dir(dir, maxlen)

#define rmkdir(path, mode) \
	    _Pragma("message(\"rmkdir() is deprecated.  Use xt_rmkdir().\")") \
	    xt_rmkdir(path, mode)

/* valid-extension.c */
bool xt_valid_extension(const char *filename, const char *valid_ext);

/* fast-cp.c */
int xt_fast_cp(const char *source, const char *dest); //__attribute__((deprecated("Use xt_fast_cp()")));

/* file-mod-cmp.c */
int xt_file_mod_cmp(const char *file1, const char *file2);

/* fd-purge.c */
void xt_fd_purge(int fd);

/* fgetline.c */
size_t xt_fgetline(FILE *fp, char *buff, size_t maxlen);

/* get-home-dir.c */
char *xt_get_home_dir(char *dir, size_t maxlen);

/* rmkdir.c */
int xt_rmkdir(const char *path, mode_t mode);

/* xt-file.c */
FILE *xt_fopen(const char *filename, const char *mode);
int xt_fclose(FILE *stream);
ssize_t xt_inhale_strings(FILE *stream, char ***list);
int xt_read_line_malloc(FILE *stream, char **buff, size_t *buff_size, size_t *len);

/* dprintf.c */
int xt_dprintf(int fd, const char * restrict format, ...);

#endif // _XTEND_FILE_H_
#ifndef _XTEND_MATH_H_
#define _XTEND_MATH_H_

#define XT_MIN(a,b) ((a) < (b) ? (a) : (b))
#define XT_MAX(a,b) ((a) > (b) ? (a) : (b))

#ifndef _INTTYPES_H_
#include <inttypes.h>
#endif

/* gcd.c */
unsigned long gcd(unsigned long a, unsigned long b);
unsigned long lcm(unsigned long a, unsigned long b);

/* digits.c */
int digits(long val, unsigned base);

/* numeric_cmp.c */
int double_cmp(const double *d1, const double *d2);
int float_cmp(const float *d1, const float *d2);
int long_long_cmp(const long long *d1, const long long *d2);
int long_cmp(const long *d1, const long *d2);
int int_cmp(const int *d1, const int *d2);
int short_cmp(const short *d1, const short *d2);

/* combinatorics.c */
unsigned long xt_n_choose_k(unsigned long n, unsigned long k);
uint64_t xt_factorial(unsigned n);

#endif  // _XTEND_MATH_H_
#ifndef _XTEND_MEM_H_
#define _XTEND_MEM_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

/* xt-malloc.c */
void *xt_malloc(size_t nelem, size_t size);
void *xt_realloc(void *array, size_t nelem, size_t size);
void xt_free_strings(char **list);

#endif // _XTEND_MEM_H_
#ifndef _XTEND_NET_H_
#define _XTEND_NET_H_

int resolve_hostname(const char *hostname, char *ip, size_t ip_buff_len);

#ifndef _XTEND_COMMON_H_
#endif

#endif  // _XTEND_NET_H_
#ifndef _XTEND_PROC_H_
#define _XTEND_PROC_H_

/*
 *  Process control
 */

/* spawn*() parent_action */
#define P_NOWAIT  0
#define P_WAIT    1

/* spawn*() echo */
#define P_NOECHO  0
#define P_ECHO    1

#define P_TERM_STATUS(s)    ((s) & 0xff)
#define P_EXIT_CODE(s)      (((s) & 0x0000ff00) >> 8)
#define P_EXEC_FAILED(s)    ((s) & 0x8000)

#ifndef _XTEND_COMMON_H_
#endif

/* parse-cmd.c */
char *parse_cmd(char *argv[], int max_args, const char *cmd);

/* spawnlp.c */
int spawnlp(int parent_action, int echo, char *infile, char *outfile, char *errfile, char *arg0, ...);

/* spawnvp.c */
int spawnvp(int parent_action, int echo, char *argv[], char *infile, char *outfile, char *errfile);
void redirect(char *infile, char *outfile, char *errfile);

/* va-usage.c */
void va_usage(const char *format_string, ...);

#endif  // _XTEND_PROC_H_
#ifndef _XTEND_ARRAY_H_
#define _XTEND_ARRAY_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

/* xt-shuffle.c */
void xt_shuffle(void *base, size_t nelem, size_t size);

/* roman.c */
int romantoi(const char *nptr, char **endptr);

#endif // _XTEND_ARRAY_H_
#ifndef _XTEND_STRING_H_
#define _XTEND_STRING_H_

#ifdef __linux__
#define strlcpy(dest,src,len)   strcpy(dest,src)
#define strlcat(dest,src,len)   strcat(dest,src)
#endif

#ifndef _STDIO_H_
#include <stdio.h>  // size_t
#endif

#ifndef _INTTYPES_H_
#include <inttypes.h>
#endif

/* string.c */
size_t strlupper(char *dest, const char *src, size_t dest_size);
size_t strupper(char *str);
size_t strllower(char *dest, const char *src, size_t dest_size);
size_t strlower(char *str);
size_t str_argv_cat(char *string, char *argv[], size_t first_arg, size_t string_buff_size);
int strblank(const char *string);
int strisint(const char *string, int base);
int strisreal(const char *string);
char *strlbasecpy(char *dest, const char *dest_base, const char *src, size_t dstsize);
int strptrcmp(const char **p1, const char **p2);
int strptrcasecmp(const char **p1, const char **p2);
int strshellcpy(char *dest, const char *src, size_t dest_len);
size_t strsqueeze(char *dest, const char *src, size_t dstsize);
void strtr(char *string, const char *from, const char *to, int flags);
void strtrim(char *string, const char *fat);
char *strviscpy(unsigned char *dest, const unsigned char *src, size_t maxlen);
char *ltostrn(char string[], long val, unsigned base, size_t maxlen);
uint64_t str2u64(const char *str);
int strsplit(char *string, char ***array, const char *sep);

#endif  // _XTEND_STRING_H_
#ifndef _XTEND_TIME_H_
#define _XTEND_TIME_H_

#ifndef _STDIO_H_
#include <stdio.h>
#endif

#ifndef _SYS_TIME_H_
#include <sys/time.h>
#endif

#ifndef _SYS_RESOURCE_H
#include <sys/resource.h>
#endif

#define difftimeofday(later, earlier) \
	    _Pragma("message(\"difftimeofday() is deprecated.  Use xt_difftimeofday().\")") \
	    xt_difftimeofday(later, earlier)

/* difftimeofday.c */
time_t xt_difftimeofday(struct timeval *later, struct timeval *earlier);
int xt_tic(struct timeval *start_time, struct rusage *start_usage);
unsigned long xt_toc(FILE *stream, const char *message, struct timeval *start_time, struct rusage *start_usage);

#endif  // _XTEND_TIME_H_
