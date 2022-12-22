#include <stdio.h>
#include <ctype.h>


#include "libxtend.h"
#include "biolibc.h"

/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/align.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate the leftmost (farthest 5') match for sequence little within
 *      sequence big, tolerating the given percentage of mismatched bases.
 *
 *      The content of little is assumed to be all upper case.  This
 *      improves speed by avoiding numerous redundant toupper()
 *      conversions on the same string, assuming multiple big strings will
 *      be searched for little, as in adapter removal and read mapping.
 *      Use strlupper(3) or strupper(3) before calling this function if
 *      necessary.
 *
 *      A minimum of min_match bases must match between little and
 *      big.  This mainly matters near the end of big, where
 *      remaining bases are fewer than the length of little.
 *
 *      A maximum of max_mismatch_percent mismatched bases are tolerated
 *      to allow for read errors. This is taken as a percent of little, or
 *      the same percent of remaining bases in big, whichever is smaller.
 *      Note that the NUMBER of allowed mismatched bases tolerated is
 *      truncated from the percent calculation.  E.g. using 10% tolerance,
 *      0 mismatched bases are tolerated among 9 total bases, or 1 mismatch
 *      among 10 total.
 *
 *      Higher values of max_mismatch_percent will results in slightly
 *      longer run times, more alignments detected, and a higher risk of
 *      false-positives (falsely identifying other big sequences as matching
 *      little.
 *
 *      Indels (insertions and deletions) are not currently handled.
 *
 *      Note that alignment is not an exact science.  We cannot detect every
 *      true little sequence without falsely detecting other sequences, since
 *      it is impossible to know whether any given sequence is really from
 *      the source of interest (e.g. an adapter) or naturally
 *      occurring from another source.  The best we can do is guestimate
 *      what will provide the most true positives (best statistical power)
 *      and fewest false positives.
 *
 *      In the case of adapter removal,
 *      it is also not usually important to remove every adapter, but only to
 *      minimize adapter contamination.  Failing to align a small percentage
 *      of sequences due to adapter contamination will not change the story
 *      told by the downstream analysis.  Nor will erroneously trimming off
 *      the 3' end of a small percentage of reads containing natural
 *      sequences resembling adapters.  Just trimming exact matches of
 *      the adapter sequence will generally remove 99% or more of the
 *      adapter contamination and minimize false-positives.  Tolerating
 *      1 or 2 differences has been shown to do slightly better overall.
 *      Modern read mapping software is also tolerant of adapter
 *      contamination and can clip adapters as needed.
 *
 *  Arguments:
 *      params      bl_align_t parameters.  Only min_match and
 *                  max_mismatch_percent are used.
 *      big         Sequence to be searched for matches to little
 *      little      Sequence to be located within big
 *
 *  Returns:
 *      Index of little sequence within big if found, index of null
 *      terminator of big otherwise
 *
 *  Examples:
 *      bl_param_t  params;
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      bl_align_set_min_match(&params, 3);
 *      bl_align_set_max_mismatch_percent(&params, 10);
 *      index = bl_align_map_seq_sub(&params,
 *          BL_FASTQ_SEQ(&read), BL_FASTQ_SEQ_LEN(&read),
 *          little, strlen(adapter)3, 10);
 *      if ( index != BL_FASTQ_SEQ_LEN(&read) )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_align_map_seq_exact(3), bl_align_set_min_match(3),
 *      bl_align_set_max_mismatch_percent(3), bl_fastq_3p_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_align_map_seq_sub(const bl_align_t *params,
	    const char *big, size_t big_len,
	    const char *little, size_t little_len)

{
    // The strlen() looks expensive, but tests show that eliminating it
    // doesn't reduce run time measurably
    size_t      mismatch, max_mismatch,
		start, bc, lc,
		md, little_mm,
		min_match = params->min_match;
    
    // Start at 5' end assuming 5' littles already removed
    // Cutadapt uses a semiglobal alignment algorithm to find littles.
    // Not sure what the benefit of this is over exact matching. I would
    // assume that errors in little sequences are extremely rare.
    // https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm

    // Convert max mismatch percentage to a divisor for the string len
    md = 100 / params->max_mismatch_percent;
    little_mm = little_len / md;  // Max mismatch based on little len
    // Could stop at big_len - min_match, but the extra math
    // outweights the few iterations saved
    for (start = 0; start < big_len; ++start)
    {
	// Terminate loop as soon as max_mismatch is reached, before
	// checking other conditions
	max_mismatch = XT_MIN((big_len - start) / md, little_mm);
	for (bc = start, lc = 0, mismatch = 0;
	     (mismatch <= max_mismatch) &&
	     (lc < little_len) && (bc < big_len); ++bc, ++lc)
	{
	    if ( toupper(big[bc]) != little[lc] )
		++mismatch;
	}
	if ( mismatch <= max_mismatch )
	{
	    if ( lc - mismatch >= min_match )
		return start;
	}
    }
    return big_len;   // Location of '\0' terminator
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/align.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate the leftmost (farthest 5') match for sequence little within
 *      sequence big, using exact matching only.
 *
 *      The content of little is assumed to be all upper case.  This
 *      improves speed by avoiding numerous redundant toupper()
 *      conversions on the same string, assuming multiple big strings will
 *      be searched for little, as in adapter removal and read mapping.
 *      Use strlupper(3) or strupper(3) before calling this function if
 *      necessary.
 *
 *      A minimum of min_match bases must match between little and
 *      big.  This mainly matters near the end of big, where
 *      remaining bases are fewer than the length of little.
 *
 *      Note that alignment is not an exact science.  We cannot detect every
 *      true little sequence without falsely detecting other sequences, since
 *      it is impossible to know whether any given sequence is really from
 *      the source of interest (e.g. an adapter) or naturally
 *      occurring from another source.  The best we can do is guestimate
 *      what will provide the most true positives (best statistical power)
 *      and fewest false positives.
 *
 *      In the case of adapter removal,
 *      it is also not usually important to remove every adapter, but only to
 *      minimize adapter contamination.  Failing to align a small percentage
 *      of sequences due to adapter contamination will not change the story
 *      told by the downstream analysis.  Nor will erroneously trimming off
 *      the 3' end of a small percentage of reads containing natural
 *      sequences resembling adapters.  Just trimming exact matches of
 *      the adapter sequence will generally remove 99% or more of the
 *      adapter contamination and minimize false-positives.  Tolerating
 *      1 or 2 differences has been shown to do slightly better overall.
 *      Modern read mapping software is also tolerant of adapter
 *      contamination and can clip adapters as needed.
 *
 *  Arguments:
 *      params      bl_align_t parameters.  Only min_match is used.
 *      big         Sequence to be searched for matches to little
 *      little      Sequence to be located within big
 *
 *  Returns:
 *      Index of little sequence within big if found, index of null
 *      terminator of big otherwise
 *
 *  Examples:
 *      bl_param_t  params;
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      bl_align_set_min_match(&params, 3);
 *      index = bl_align_map_seq_exact(&params,
 *          BL_FASTQ_SEQ(&read), BL_FASTQ_SEQ_LEN(&read),
 *          little, strlen(adapter)3, 10);
 *      if ( index != BL_FASTQ_SEQ_LEN(&read) )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_align_map_seq_sub(3), bl_align_set_min_match(3),
 *      bl_fastq_3p_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_align_map_seq_exact(const bl_align_t *params,
	    const char *big, size_t big_len,
	    const char *little, size_t little_len)

{
    size_t  start, bc, lc;
    
    // Start at 5' end assuming 5' adapters already removed
    // Cutadapt uses a semiglobal alignment algorithm to find adapters.
    // Not sure what the benefit of this is over exact matching. I would
    // assume that errors in adapter sequences are extremely rare.
    // https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm

    // Could stop at big_len - min_match, but the extra math
    // outweights the few iterations saved
    for (start = 0; start < big_len; ++start)
    {
	for (bc = start, lc = 0; (toupper(big[bc]) == little[lc]) &&
	     (lc < little_len); ++bc, ++lc)
	    ;
	if ( (lc == little_len) || ((bc == big_len) &&
	     (lc >= params->min_match)) )
	    return start;
    }
    return big_len;   // Location of '\0' terminator
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
 *      #include <biolibc/align.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for min_match member in a bl_align_t structure.
 *      Use this function to set min_match in a bl_align_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      min_match is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_align_ptr    Pointer to the structure to set
 *      new_min_match   The new value for min_match
 *
 *  Returns:
 *      BL_ALIGN_DATA_OK if the new value is acceptable and assigned
 *      BL_ALIGN_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_align_t      bl_align;
 *      size_t          new_min_match;
 *
 *      if ( bl_align_set_min_match(&bl_align, new_min_match)
 *              == BL_ALIGN_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from align.h
 ***************************************************************************/

int     bl_align_set_min_match(
	    bl_align_t *bl_align_ptr,
	    size_t new_min_match
	)

{
    if ( false )
	return BL_ALIGN_DATA_OUT_OF_RANGE;
    else
    {
	bl_align_ptr->min_match = new_min_match;
	return BL_ALIGN_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/align.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for max_mismatch_percent member in a bl_align_t structure.
 *      Use this function to set max_mismatch_percent in a bl_align_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      max_mismatch_percent is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_align_ptr    Pointer to the structure to set
 *      new_max_mismatch_percent The new value for max_mismatch_percent
 *
 *  Returns:
 *      BL_ALIGN_DATA_OK if the new value is acceptable and assigned
 *      BL_ALIGN_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_align_t      bl_align;
 *      unsigned        new_max_mismatch_percent;
 *
 *      if ( bl_align_set_max_mismatch_percent(&bl_align, new_max_mismatch_percent)
 *              == BL_ALIGN_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from align.h
 ***************************************************************************/

int     bl_align_set_max_mismatch_percent(
	    bl_align_t *bl_align_ptr,
	    unsigned new_max_mismatch_percent
	)

{
    if ( false )
	return BL_ALIGN_DATA_OUT_OF_RANGE;
    else
    {
	bl_align_ptr->max_mismatch_percent = new_max_mismatch_percent;
	return BL_ALIGN_DATA_OK;
    }
}
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>   // PRId64




/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over header lines in bed input stream, leaving the FILE
 *      structure pointing to the first character in the first line of data.
 *      The header is copied to a temporary file whose FILE pointer 
 *      is returned.
 *
 *  Arguments:
 *      stream  Pointer to the FILE structure for reading the BED stream
 *
 *  Returns:
 *      Pointer to the FILE structure of the temporary file.
 *
 *  Examples:
 *      FILE    *header, *bed_stream;
 *      ...
 *      header = bl_bed_skip_header(bed_stream);
 *
 *  See also:
 *      bl_bed_read(3), xt_fopen(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_bed_skip_header(FILE *bed_stream)

{
    char    start[7] = "xxxxxx";
    size_t  count;
    int     ch, c;
    FILE    *header_stream = tmpfile();

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like peak-classifier to replicate the
     *  header in output files.
     */
    
    while ( ((count=fread(start, 1, 7, bed_stream)) == 7) && 
	    ((memcmp(start, "browser", 7) == 0) ||
	    (memcmp(start, "track", 5) == 0) ||
	    (*start == '#')) )
    {
	fwrite(start, count, 1, header_stream);
	do
	{
	    ch = getc(bed_stream);
	    putc(ch, header_stream);
	}   while ( (ch != '\n') && (ch != EOF) );
    }
    
    // Rewind to start of first non-header line
    for (c = count - 1; c >= 0; --c)
	ungetc(start[c], bed_stream);
    rewind(header_stream);
    return header_stream;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read next entry (line) from a BED file.  The line must have at
 *      least the first 3 fields (chrom, start, and end).  It may
 *      have up to 12 fields, all of which must be in the correct order
 *      according to the BED specification.
 *
 *      If field_mask is not BL_BED_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in bed_feature.
 *      Possible mask values are:
 *
 *      BL_BED_FIELD_ALL
 *      BL_BED_FIELD_NAME
 *      BL_BED_FIELD_SCORE
 *      BL_BED_FIELD_STRAND
 *      BL_BED_FIELD_THICK
 *      BL_BED_FIELD_RGB
 *      BL_BED_FIELD_BLOCK
 *
 *      The chrom, start, and end fields are required and therefore have
 *      no corresponding mask bits. The thickStart and thickEnd fields must
 *      occur together or not at all, so only a single bit BL_BED_FIELD_THICK
 *      selects both of them.  Likewise, blockCount, blockSizes and
 *      blockStarts must all be present or omitted, so BL_BED_FIELD_BLOCK
 *      masks all three.
 *
 *  Arguments:
 *      bed_feature     Pointer to a bl_bed_t structure
 *      bed_stream      A FILE stream from which to read the line
 *      field_mask      Bit mask indicating which fields to store in bed_feature
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered at the start of a line
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_bed_read(stdin, &bed_feature, BL_BED_FIELD_ALL);
 *      bl_bed_read(bed_stream, &bed_feature,
 *                       BL_BED_FIELD_NAME|BL_BED_FIELD_SCORE);
 *
 *  See also:
 *      bl_bed_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_bed_read(bl_bed_t *bed_feature, FILE *bed_stream,
	    bed_field_mask_t field_mask)

{
    char    *end,
	    strand[BL_BED_STRAND_MAX_CHARS + 1],
	    block_count_str[BL_BED_BLOCK_COUNT_MAX_DIGITS + 1],
	    block_size_str[BL_BED_BLOCK_SIZE_MAX_DIGITS + 1],
	    block_start_str[BL_BED_BLOCK_START_MAX_DIGITS + 1],
	    chrom_start_str[BL_POSITION_MAX_DIGITS + 1],
	    chrom_end_str[BL_POSITION_MAX_DIGITS + 1],
	    score_str[BL_BED_SCORE_MAX_DIGITS + 1],
	    thick_start_str[BL_POSITION_MAX_DIGITS + 1],
	    thick_end_str[BL_POSITION_MAX_DIGITS + 1];
    size_t  len;
    int     delim;
    unsigned long   block_count;
    unsigned    c;
    
    // FIXME: Respect field_mask
    
    // Chromosome
    if ( tsv_read_field(bed_stream, bed_feature->chrom,
			BL_CHROM_MAX_CHARS, &len) == EOF )
    {
	// fputs("bl_bed_read(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BL_READ_EOF;
    }
    
    // Feature start position
    if ( tsv_read_field(bed_stream, chrom_start_str,
			BL_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_bed_read(): Got EOF reading start position: %s.\n",
		chrom_start_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	bed_feature->chrom_start = strtoul(chrom_start_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_bed_read(): Invalid start position: %s\n",
		    chrom_start_str);
	    return BL_READ_TRUNCATED;
	}
    }
    
    // Feature end position
    // FIXME: Check for > or < start if strand + or -
    if ( (delim = tsv_read_field(bed_stream, chrom_end_str,
			BL_POSITION_MAX_DIGITS, &len)) == EOF )
    {
	fprintf(stderr, "bl_bed_read(): Got EOF reading end position: %s.\n",
		chrom_end_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	bed_feature->chrom_end = strtoul(chrom_end_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_bed_read(): Invalid end position: %s\n",
		    chrom_end_str);
	    return BL_READ_TRUNCATED;
	}
    }

    bed_feature->fields = 3;
    
    // Read NAME field if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, bed_feature->name,
			    BL_BED_NAME_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading name: %s.\n",
		    bed_feature->name);
	    return BL_READ_TRUNCATED;
	}
	++bed_feature->fields;
    }
    
    // Read SCORE if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, score_str,
			    BL_POSITION_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading score: %s.\n",
		    score_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->score = strtoul(score_str, &end, 10);
	    if ( (*end != '\0') || (bed_feature->score > 1000) )
	    {
		fprintf(stderr,
			"bl_bed_read(): Invalid feature score: %s\n",
			score_str);
		return BL_READ_TRUNCATED;
	    }
	}
	++bed_feature->fields;
    }
    
    // Read strand if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, strand,
			    BL_BED_STRAND_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading strand: %s.\n",
		    bed_feature->name);
	    return BL_READ_TRUNCATED;
	}
	if ( (len != 1) || ((*strand != '+') && (*strand != '-') && (*strand != '.')) )
	{
	    fprintf(stderr, "bl_bed_read(): Strand must be + or - or .: %s\n",
		    strand);
	    return BL_READ_TRUNCATED;
	}
	bed_feature->strand = *strand;
	++bed_feature->fields;
    }
    
    // Read thick start position if present
    // Must be followed by thick end position, > or < for + or - strand
    // Feature start position
    if ( delim != '\n' )
    {
	if ( tsv_read_field(bed_stream, thick_start_str,
			    BL_POSITION_MAX_DIGITS, &len) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading thick start "
		    "POS: %s.\n", thick_start_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->thick_start =
		strtoul(thick_start_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid thick start "
				"position: %s\n",
				thick_start_str);
		return BL_READ_TRUNCATED;
	    }
	}
	
	if ( delim == '\n' )
	{
	    fprintf(stderr, "bl_bed_read(): Found thick start, but no thick end.\n");
	    return BL_READ_TRUNCATED;
	}
    
	if ( tsv_read_field(bed_stream, thick_end_str,
			    BL_POSITION_MAX_DIGITS, &len) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading thick end "
		    "POS: %s.\n", thick_end_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    bed_feature->thick_end =
		strtoul(thick_end_str, &end, 10);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid thick end "
				"position: %s\n",
				thick_end_str);
		return BL_READ_TRUNCATED;
	    }
	}
	bed_feature->fields += 2;
    }

    // Read RGB string field if present
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, bed_feature->item_rgb,
			    BL_BED_ITEM_RGB_MAX_CHARS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading RGB: %s.\n",
		    bed_feature->name);
	    return BL_READ_TRUNCATED;
	}
	++bed_feature->fields;
    }

    /*
     *  Read block count if present
     *  Must be followed by comma-separated list of sizes
     *  and comma-separated list of start positions
     */
    if ( delim != '\n' )
    {
	if ( (delim = tsv_read_field(bed_stream, block_count_str,
			    BL_BED_BLOCK_COUNT_MAX_DIGITS, &len)) == EOF )
	{
	    fprintf(stderr, "bl_bed_read(): Got EOF reading block count: %s.\n",
		    score_str);
	    return BL_READ_TRUNCATED;
	}
	else
	{
	    block_count = strtoul(block_count_str, &end, 10);
	    if ( (*end != '\0') || (block_count > 65535) )
	    {
		fprintf(stderr,
			"bl_bed_read(): Invalid block count: %s\n",
			score_str);
		return BL_READ_TRUNCATED;
	    }
	    bed_feature->block_count = block_count;
	}
	bed_feature->block_sizes = xt_malloc(bed_feature->block_count,
					sizeof(*bed_feature->block_sizes));
	if ( bed_feature->block_sizes == NULL )
	{
	    fputs("bl_bed_read(): Cannot allocate block_sizes.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
	bed_feature->block_starts = xt_malloc(bed_feature->block_count,
					sizeof(*bed_feature->block_starts));
	if ( bed_feature->block_starts == NULL )
	{
	    fputs("bl_bed_read(): Cannot allocate block_starts.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
	if ( delim == '\n' )
	{
	    fputs("bl_bed_read(): Found block count, but no sizes.\n", stderr);
	    return BL_READ_TRUNCATED;
	}
	
	// Read comma-separated sizes
	c = 0;
	do
	{
	    delim = dsv_read_field(bed_stream, block_size_str,
			    BL_BED_BLOCK_SIZE_MAX_DIGITS, ",\t", &len);
	    bed_feature->block_sizes[c++] = strtoul(block_size_str, &end, 10);
	    //fprintf(stderr, "Block size[%u] = %s\n", c-1, block_size_str);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid block size: %s\n",
			block_size_str);
		return BL_READ_TRUNCATED;
	    }
	}   while ( delim == ',' );
	if ( c != bed_feature->block_count )
	{
	    fprintf(stderr, "bl_bed_read(): Block count = %u  Sizes = %u\n",
		    bed_feature->block_count, c);
	    return BL_READ_MISMATCH;
	}
	if ( delim == '\n' )
	{
	    fputs("bl_bed_read(): Found block sizes, but no starts.\n", stderr);
	    return BL_READ_TRUNCATED;
	}
	
	// Read comma-separated starts
	c = 0;
	do
	{
	    delim = dsv_read_field(bed_stream, block_start_str,
			    BL_BED_BLOCK_START_MAX_DIGITS, ",\t", &len);
	    bed_feature->block_starts[c++] = strtoul(block_start_str, &end, 10);
	    //fprintf(stderr, "Block start[%u] = %s\n", c-1, block_start_str);
	    if ( *end != '\0' )
	    {
		fprintf(stderr, "bl_bed_read(): Invalid block start: %s\n",
			block_start_str);
		return BL_READ_TRUNCATED;
	    }
	}   while ( delim == ',' );
	if ( c != bed_feature->block_count )
	{
	    fprintf(stderr, "bl_bed_read(): Block count = %u  Sizes = %u\n",
		    bed_feature->block_count, c);
	    return BL_READ_MISMATCH;
	}
	bed_feature->fields += 3;
    }

    //fprintf(stderr, "Bed fields = %u\n", bed_feature->fields);
    /*
     *  There shouldn't be anything left at this point.  Once block reads
     *  are implemented, we should error out of delim != '\n'
     */
    
    if ( delim != '\n' )
    {
	fputs("bl_bed_read(): Extra columns found.\n", stderr);
	return BL_READ_EXTRA_COLS;
    }
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write fields from one line of a bed file to the specified FILE
 *      stream.  If field_mask is not BL_BED_FIELD_ALL, only selected fields
 *      are written.
 *
 *      If field_mask is not BL_BED_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate marker for that field,
 *      such as a '.', rather than writing the real data.
 *      Possible mask values are:
 *
 *      BL_BED_FIELD_NAME
 *      BL_BED_FIELD_SCORE
 *      BL_BED_FIELD_STRAND
 *      BL_BED_FIELD_THICK
 *      BL_BED_FIELD_RGB
 *      BL_BED_FIELD_BLOCK
 *
 *      The chrom, start, and end fields are required and therefore have
 *      no corresponding mask bits. The thickStart and thickEnd fields must
 *      occur together or not at all, so only a single bit BL_BED_FIELD_THICK
 *      selects both of them.  Likewise, blockCount, blockSizes and
 *      blockStarts must all be present or omitted, so BL_BED_FIELD_BLOCK
 *      masks all three.
 *
 *  Arguments:
 *      bed_feature     Pointer to the bl_bed_t structure to output
 *      bed_stream      FILE stream to which TSV bed line is written
 *      field_mask      Bit mask indicating which fields to output
 *
 *  Returns:
 *      BL_WRITE_OK on success
 *      BL_WRITE_ERROR on failure (errno may provide more information)
 *
 *  Examples:
 *      bl_bed_write(stdout, &bed_feature, BL_BED_FIELD_ALL);
 *      bl_bed_write(bed_stream, &bed_feature,
 *                        BL_BED_FIELD_NAME|BL_BED_FIELD_SCORE);
 *
 *  See also:
 *      bl_bed_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_bed_write(bl_bed_t *bed_feature, FILE *bed_stream,
	    bed_field_mask_t field_mask)

{
    unsigned    c;
    
    // FIXME: Respect field_mask
    // FIXME: Check fprintf() return codes
    fprintf(bed_stream, "%s\t%" PRId64 "\t%" PRId64,
	    bed_feature->chrom,
	    bed_feature->chrom_start, bed_feature->chrom_end);
    if ( bed_feature->fields > 3 )
	fprintf(bed_stream, "\t%s", bed_feature->name);
    if ( bed_feature->fields > 4 )
	fprintf(bed_stream, "\t%u", bed_feature->score);
    if ( bed_feature->fields > 5 )
	fprintf(bed_stream, "\t%c", bed_feature->strand);
    if ( bed_feature->fields > 6 )
	fprintf(bed_stream, "\t%" PRId64 "\t%" PRId64,
		bed_feature->thick_start, bed_feature->thick_end);
    if ( bed_feature->fields > 8 )
	fprintf(bed_stream, "\t%s", bed_feature->item_rgb);
    if ( bed_feature->fields > 9 )
    {
	fprintf(bed_stream, "\t%u\t", bed_feature->block_count);
	for (c = 0; c < bed_feature->block_count - 1; ++c)
	    fprintf(bed_stream, "%" PRId64 ",", bed_feature->block_sizes[c]);
	fprintf(bed_stream, "%" PRId64 "\t", bed_feature->block_sizes[c]);
	for (c = 0; c < bed_feature->block_count - 1; ++c)
	    fprintf(bed_stream, "%" PRId64 ",", bed_feature->block_starts[c]);
	fprintf(bed_stream, "%" PRId64, bed_feature->block_starts[c]);
    }
    putc('\n', bed_stream);
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Make sure the BED input is sorted by chrom and start position.
 *
 *  Arguments:
 *      bed_feature     Pointer to BED structure containing current entry
 *      last_chrom      Chromosome of the previous BED entry
 *      last_start      Start position of the previous BED entry
 *
 *  Returns:
 *      Nothing: Terminates process if input is out of order
 *
 *  See also:
 *      bl_bed_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-08  Jason Bacon Begin
 ***************************************************************************/

void    bl_bed_check_order(bl_bed_t *bed_feature, char last_chrom[],
			   int64_t last_start)

{
    if ( bl_chrom_name_cmp(bed_feature->chrom, last_chrom) == 0 )
    {
	if ( bed_feature->chrom_start < last_start )
	{
	    fprintf(stderr, "peak-classifier: BED file not sorted by start position.\n");
	    exit(EX_DATAERR);
	}
    }
    else if ( bl_chrom_name_cmp(bed_feature->chrom, last_chrom) < 0 )
    {
	fprintf(stderr, "peak-classifier: BED file not sorted by chrom.\n");
	fprintf(stderr, "%s, %s\n", bed_feature->chrom, last_chrom);
	exit(EX_DATAERR);
    }
}

/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Compare the position of a BED feature to that of a GFF feature.
 *      Return 0 if the features overlap, < 0 if the BED feature is upstream
 *      of the GFF feature, > 0 if the BED feature is downstream of the GFF
 *      feature.
 *
 *      If the features overlap, populate the bl_overlap_t structure
 *      pointed to by overlap.  The structure contains the lengths of the
 *      two features, the start and end positions of the overlapping region,
 *      and the length of the overlap.  Positions in overlap are 1-based and
 *      inclusive at both ends (like most bioinformatics formats and unlike
 *      BED).
 *
 *  Arguments:
 *      bed_feature     Pointer to the bl_bed_t structure to compare
 *      gff_feature     Pointer to the bl_gff_t structure to compare
 *      overlap         Pointer to the bl_overlap_t structure to receive
 *                      comparison results
 *
 *  Returns:
 *      A value < 0 if the BED feature is upstream of the GFF feature
 *      A value > 0 if the BED feature is downstream of the GFF feature
 *      0 if the BED feature overlaps the GFF feature
 *
 *  Examples:
 *      if ( bl_bed_gff_cmp(&bed_feature, &gff_feature, *overlap) == 0 )
 *
 *  See also:
 *      bl_bed_read(3), bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_bed_gff_cmp(bl_bed_t *bed_feature, bl_gff_t *gff_feature,
		       bl_overlap_t *overlap)

{
    int         chrom_cmp;
    int64_t    bed_start, bed_end, bed_len,
		gff_start, gff_end, gff_len;
    
    chrom_cmp = bl_chrom_name_cmp(BL_BED_CHROM(bed_feature),
					 BL_GFF_SEQID(gff_feature));
    if ( chrom_cmp == 0 )
    {
	/*
	 *  BED positions are 0-based, with end non-inclusive, which can
	 *  also be viewed as an inclusive 1-based coordinate
	 *  GFF is 1-based, both ends inclusive
	 */
	
	if ( BL_BED_CHROM_END(bed_feature) < BL_GFF_START(gff_feature) )
	{
	    bl_overlap_set_all(overlap, 0, 0, 0, 0);
	    return -1;
	}
	else if ( BL_BED_CHROM_START(bed_feature) + 1 > BL_GFF_END(gff_feature) )
	{
	    bl_overlap_set_all(overlap, 0, 0, 0, 0);
	    return 1;
	}
	else
	{
	    bed_start = BL_BED_CHROM_START(bed_feature);
	    bed_end = BL_BED_CHROM_END(bed_feature);
	    gff_start = BL_GFF_START(gff_feature);
	    gff_end = BL_GFF_END(gff_feature);
	    bed_len = bed_end - bed_start;
	    gff_len = gff_end - gff_start + 1;
	    bl_overlap_set_all(overlap, bed_len, gff_len,
			    XT_MAX(bed_start+1, gff_start),
			    XT_MIN(bed_end, gff_end));
	    return 0;
	}
    }
    return chrom_cmp;
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
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of chrom member in a bl_bed_t
 *      structure. Use this function to set bl_bed_ptr->chrom[c]
 *      in a bl_bed_t object from non-member functions.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      c               Subscript to the chrom array
 *      new_chrom_element The new value for chrom[c]
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      size_t          c;
 *      char            new_chrom_element;
 *
 *      if ( bl_bed_set_chrom_ae(&bl_bed, c, new_chrom_element)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_CHROM_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_chrom_ae(
	    bl_bed_t *bl_bed_ptr,
	    size_t c,
	    char new_chrom_element
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->chrom[c] = new_chrom_element;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for chrom member in a bl_bed_t structure.
 *      Use this function to set chrom in a bl_bed_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_chrom to bl_bed_ptr->chrom.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_chrom       The new value for chrom
 *      array_size      Size of the chrom array.
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      char            new_chrom;
 *      size_t          array_size;
 *
 *      if ( bl_bed_set_chrom_cpy(&bl_bed, new_chrom, array_size)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_CHROM(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_chrom_cpy(
	    bl_bed_t *bl_bed_ptr,
	    char new_chrom[],
	    size_t array_size
	)

{
    if ( new_chrom == NULL )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_bed_ptr->chrom, new_chrom, array_size);
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for chrom_start member in a bl_bed_t structure.
 *      Use this function to set chrom_start in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      chrom_start is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_chrom_start The new value for chrom_start
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t        new_chrom_start;
 *
 *      if ( bl_bed_set_chrom_start(&bl_bed, new_chrom_start)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_chrom_start(
	    bl_bed_t *bl_bed_ptr,
	    int64_t new_chrom_start
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->chrom_start = new_chrom_start;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for chrom_end member in a bl_bed_t structure.
 *      Use this function to set chrom_end in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      chrom_end is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_chrom_end   The new value for chrom_end
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t        new_chrom_end;
 *
 *      if ( bl_bed_set_chrom_end(&bl_bed, new_chrom_end)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_chrom_end(
	    bl_bed_t *bl_bed_ptr,
	    int64_t new_chrom_end
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->chrom_end = new_chrom_end;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of name member in a bl_bed_t
 *      structure. Use this function to set bl_bed_ptr->name[c]
 *      in a bl_bed_t object from non-member functions.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      c               Subscript to the name array
 *      new_name_element The new value for name[c]
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      size_t          c;
 *      char            new_name_element;
 *
 *      if ( bl_bed_set_name_ae(&bl_bed, c, new_name_element)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_NAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_name_ae(
	    bl_bed_t *bl_bed_ptr,
	    size_t c,
	    char new_name_element
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->name[c] = new_name_element;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for name member in a bl_bed_t structure.
 *      Use this function to set name in a bl_bed_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_name to bl_bed_ptr->name.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_name        The new value for name
 *      array_size      Size of the name array.
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      char            new_name;
 *      size_t          array_size;
 *
 *      if ( bl_bed_set_name_cpy(&bl_bed, new_name, array_size)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_NAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_name_cpy(
	    bl_bed_t *bl_bed_ptr,
	    char new_name[],
	    size_t array_size
	)

{
    if ( new_name == NULL )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_bed_ptr->name, new_name, array_size);
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for score member in a bl_bed_t structure.
 *      Use this function to set score in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      score is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_score       The new value for score
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      unsigned short  new_score;
 *
 *      if ( bl_bed_set_score(&bl_bed, new_score)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_score(
	    bl_bed_t *bl_bed_ptr,
	    unsigned short new_score
	)

{
    if ( new_score > 1000 )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->score = new_score;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for strand member in a bl_bed_t structure.
 *      Use this function to set strand in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      strand is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_strand      The new value for strand
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      char            new_strand;
 *
 *      if ( bl_bed_set_strand(&bl_bed, new_strand)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_strand(
	    bl_bed_t *bl_bed_ptr,
	    char new_strand
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->strand = new_strand;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for thick_start member in a bl_bed_t structure.
 *      Use this function to set thick_start in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      thick_start is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_thick_start The new value for thick_start
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t        new_thick_start;
 *
 *      if ( bl_bed_set_thick_start(&bl_bed, new_thick_start)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_thick_start(
	    bl_bed_t *bl_bed_ptr,
	    int64_t new_thick_start
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->thick_start = new_thick_start;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for thick_end member in a bl_bed_t structure.
 *      Use this function to set thick_end in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      thick_end is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_thick_end   The new value for thick_end
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t        new_thick_end;
 *
 *      if ( bl_bed_set_thick_end(&bl_bed, new_thick_end)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_thick_end(
	    bl_bed_t *bl_bed_ptr,
	    int64_t new_thick_end
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->thick_end = new_thick_end;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of item_rgb member in a bl_bed_t
 *      structure. Use this function to set bl_bed_ptr->item_rgb[c]
 *      in a bl_bed_t object from non-member functions.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      c               Subscript to the item_rgb array
 *      new_item_rgb_element The new value for item_rgb[c]
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      size_t          c;
 *      char            new_item_rgb_element;
 *
 *      if ( bl_bed_set_item_rgb_ae(&bl_bed, c, new_item_rgb_element)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_ITEM_RGB_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_item_rgb_ae(
	    bl_bed_t *bl_bed_ptr,
	    size_t c,
	    char new_item_rgb_element
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->item_rgb[c] = new_item_rgb_element;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for item_rgb member in a bl_bed_t structure.
 *      Use this function to set item_rgb in a bl_bed_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_item_rgb to bl_bed_ptr->item_rgb.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_item_rgb    The new value for item_rgb
 *      array_size      Size of the item_rgb array.
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      char            new_item_rgb;
 *      size_t          array_size;
 *
 *      if ( bl_bed_set_item_rgb_cpy(&bl_bed, new_item_rgb, array_size)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_ITEM_RGB(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_item_rgb_cpy(
	    bl_bed_t *bl_bed_ptr,
	    char new_item_rgb[],
	    size_t array_size
	)

{
    if ( new_item_rgb == NULL )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_bed_ptr->item_rgb, new_item_rgb, array_size);
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for block_count member in a bl_bed_t structure.
 *      Use this function to set block_count in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      block_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_block_count The new value for block_count
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      unsigned short  new_block_count;
 *
 *      if ( bl_bed_set_block_count(&bl_bed, new_block_count)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_block_count(
	    bl_bed_t *bl_bed_ptr,
	    unsigned short new_block_count
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->block_count = new_block_count;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for block_sizes member in a bl_bed_t structure.
 *      Use this function to set block_sizes in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      block_sizes is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_block_sizes The new value for block_sizes
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t *      new_block_sizes;
 *
 *      if ( bl_bed_set_block_sizes(&bl_bed, new_block_sizes)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_block_sizes(
	    bl_bed_t *bl_bed_ptr,
	    int64_t * new_block_sizes
	)

{
    if ( new_block_sizes == NULL )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->block_sizes = new_block_sizes;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of block_sizes member in a bl_bed_t
 *      structure. Use this function to set bl_bed_ptr->block_sizes[c]
 *      in a bl_bed_t object from non-member functions.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      c               Subscript to the block_sizes array
 *      new_block_sizes_element The new value for block_sizes[c]
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      size_t          c;
 *      int64_t *      new_block_sizes_element;
 *
 *      if ( bl_bed_set_block_sizes_ae(&bl_bed, c, new_block_sizes_element)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_BLOCK_SIZES_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_block_sizes_ae(
	    bl_bed_t *bl_bed_ptr,
	    size_t c,
	    int64_t  new_block_sizes_element
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->block_sizes[c] = new_block_sizes_element;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for block_sizes member in a bl_bed_t structure.
 *      Use this function to set block_sizes in a bl_bed_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_block_sizes to bl_bed_ptr->block_sizes.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_block_sizes The new value for block_sizes
 *      array_size      Size of the block_sizes array.
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t *      new_block_sizes;
 *      size_t          array_size;
 *
 *      if ( bl_bed_set_block_sizes_cpy(&bl_bed, new_block_sizes, array_size)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_BLOCK_SIZES(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_block_sizes_cpy(
	    bl_bed_t *bl_bed_ptr,
	    int64_t * new_block_sizes,
	    size_t array_size
	)

{
    if ( new_block_sizes == NULL )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_bed_ptr->block_sizes[c] = new_block_sizes[c];
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for block_starts member in a bl_bed_t structure.
 *      Use this function to set block_starts in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      block_starts is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_block_starts The new value for block_starts
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t *      new_block_starts;
 *
 *      if ( bl_bed_set_block_starts(&bl_bed, new_block_starts)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_block_starts(
	    bl_bed_t *bl_bed_ptr,
	    int64_t * new_block_starts
	)

{
    if ( new_block_starts == NULL )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->block_starts = new_block_starts;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of block_starts member in a bl_bed_t
 *      structure. Use this function to set bl_bed_ptr->block_starts[c]
 *      in a bl_bed_t object from non-member functions.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      c               Subscript to the block_starts array
 *      new_block_starts_element The new value for block_starts[c]
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      size_t          c;
 *      int64_t *      new_block_starts_element;
 *
 *      if ( bl_bed_set_block_starts_ae(&bl_bed, c, new_block_starts_element)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_BLOCK_STARTS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_block_starts_ae(
	    bl_bed_t *bl_bed_ptr,
	    size_t c,
	    int64_t  new_block_starts_element
	)

{
    if ( false )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->block_starts[c] = new_block_starts_element;
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for block_starts member in a bl_bed_t structure.
 *      Use this function to set block_starts in a bl_bed_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_block_starts to bl_bed_ptr->block_starts.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_block_starts The new value for block_starts
 *      array_size      Size of the block_starts array.
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      int64_t *      new_block_starts;
 *      size_t          array_size;
 *
 *      if ( bl_bed_set_block_starts_cpy(&bl_bed, new_block_starts, array_size)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_BED_SET_BLOCK_STARTS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_block_starts_cpy(
	    bl_bed_t *bl_bed_ptr,
	    int64_t * new_block_starts,
	    size_t array_size
	)

{
    if ( new_block_starts == NULL )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_bed_ptr->block_starts[c] = new_block_starts[c];
	return BL_BED_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/bed.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for fields member in a bl_bed_t structure.
 *      Use this function to set fields in a bl_bed_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      fields is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_bed_ptr      Pointer to the structure to set
 *      new_fields      The new value for fields
 *
 *  Returns:
 *      BL_BED_DATA_OK if the new value is acceptable and assigned
 *      BL_BED_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_bed_t        bl_bed;
 *      unsigned short  new_fields;
 *
 *      if ( bl_bed_set_fields(&bl_bed, new_fields)
 *              == BL_BED_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from bed.h
 ***************************************************************************/

int     bl_bed_set_fields(
	    bl_bed_t *bl_bed_ptr,
	    unsigned short new_fields
	)

{
    if ( (new_fields < 3) || (new_fields > 9) )
	return BL_BED_DATA_OUT_OF_RANGE;
    else
    {
	bl_bed_ptr->fields = new_fields;
	return BL_BED_DATA_OK;
    }
}
#include <stdio.h>
#include <sysexits.h>
#include <ctype.h>
#include <stdlib.h>
#include <limits.h>

/***************************************************************************
 *  Library:
 *      #include <biolibc/biostring.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Perform a numeric comparison of two chrom names.
 *
 *      The names may contain a prefix of non-digits, such as "chr".
 *      Characters that follow must be a chrom number or letter.
 *      Numbers are considered less than letters (e.g. 22 < X).  As such,
 *      if either is a letter, they are compared lexically.  If both are
 *      numbers, they are converted to integers and compared numerically.
 *
 *      Use bl_chrom_name_cmp() only if you need to know which string is
 *      < or >.  If only checking for equality/inequality, strcmp() will be
 *      faster.
 *
 *  Arguments:
 *      name1, name2    Names of two chroms
 *
 *  Returns:
 *      A value < 1 if name1 is numerically < name2
 *      A value > 1 if name1 is numerically > name2
 *      0 if name1 == name2
 *
 *  See also:
 *      strcmp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-07  Jason Bacon Begin
 ***************************************************************************/

int     bl_chrom_name_cmp(const char *name1, const char *name2)

{
    const char      *p1 = name1, *p2 = name2;
    char            *end;
    unsigned long   c1, c2;

    /* Skip identical portions of strings, e.g. "chr" prefix */
    while ( (*p1 == *p2) && (*p1 != '\0') )
	++p1, ++p2;
    
    /*
     *  Next should be a number or letter ID such as X or Y.  If either
     *  ID is not a number, simply compare lexically.
     *  ISO character order will take care of it since letters come after
     *  digits (chrX > chr22) and everything comes after null
     *  (chr22 > chr2).  This also handles the case where the names are
     *  the same (both *p1 and *p2 are '\0') or we reached the end of one
     *  of them (chr2, chr22).
     */
    if ( !isdigit(*p1) || !isdigit(*p2) )
	return *p1 - *p2;

    /* Both IDs are numeric, so perform an integer compare */
    c1 = strtoul(p1, &end, 10);
    if ( *end != '\0' )
    {
	fprintf(stderr,
		"bl_chrom_name_cmp(): Invalid chrom ID: %s\n", name1);
	exit(EX_DATAERR);
    }
    c2 = strtoul(p2, &end, 10);
    if ( *end != '\0' )
    {
	fprintf(stderr,
		"bl_chrom_name_cmp(): Invalid chrom ID: %s\n", name2);
	exit(EX_DATAERR);
    }
    return c1 - c2;
}
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>



/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a FASTA record from a FILE stream.  Each record must begin
 *      with a description line (beginning with '>'), which is then
 *      followed by one or more lines of sequence data.  The end of the
 *      sequence is marked either by the next description line or EOF.
 *      If desc_len and seq_len are 0 (e.g. the structure is initialized
 *      with BL_FASTA_INIT or bl_fasta_init(3), or has been freed with
 *      bl_fasta_free(3), then
 *      memory is allocated for the description and sequence.
 *
 *      Otherwise, the existing allocated buffers are reused.  Hence, when
 *      reading many FASTA records of the same length, only one allocation
 *      is needed.  In any case, the buffers are automatically enlarged if
 *      they become full and automatically trimmed to the actual data size
 *      after reading is complete.
 *
 *      Buffer memory should be freed as soon as possible by calling
 *      bl_fasta_free(3).
 *  
 *  Arguments:
 *      fasta_stream    FILE stream from which FASTA data are read
 *      record          Pointer to a bl_fasta_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fasta_t  rec = BL_FASTA_INIT;
 *
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *      bl_fasta_free(&rec);
 *
 *  See also:
 *      bl_fasta_write(3), bl_fastq_read(3), bl_fastq_write(3),
 *      bl_fasta_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fasta_read(bl_fasta_t *record, FILE *fasta_stream)

{
    int     ch,
	    last_ch;
    size_t  len;
    
    /* Skip comment lines */
    while ( ((ch = getc(fasta_stream)) == ';') && (ch != EOF) )
	while ( ((ch = getc(fasta_stream)) != '\n') && (ch != EOF) )
	    ;
    
    if ( ch == EOF )
	return BL_READ_EOF;
    
    /* Every record should begin with a '>' */
    if ( ch == '>' )    // Desc
    {
	ungetc(ch, fasta_stream);
	ch = dsv_read_field_malloc(fasta_stream, &record->desc,
			    &record->desc_array_size, "", &record->desc_len);
	if ( record->desc == NULL )
	{
	    fprintf(stderr, "bl_fasta_read(): Could not allocate desc.\n");
	    exit(EX_UNAVAILABLE);
	}
	
	/* Should not encounter EOF while reading description line */
	/* Every description should be followed by at least one seq line */
	if ( ch == EOF )
	{
	    fprintf(stderr, "bl_fasta_read(): Record truncated in desc %s.\n",
		    record->desc);
	    return BL_READ_TRUNCATED;
	}
	
	/*
	 *  Read sequence lines
	 */
	
	if ( record->seq_array_size == 0 )
	{
	    // 128 MiB will hold many chromosomes and minimize reallocs
	    record->seq_array_size = 128 * 1024 * 1024;
	    //fprintf(stderr, "Allocating initial array of %zu\n", record->seq_array_size);
	    record->seq = xt_malloc(record->seq_array_size, sizeof(*record->seq));
	    if ( record->seq == NULL )
	    {
		fprintf(stderr, "bl_fasta_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	
	len = 0;
	do
	{
	    if ( ch != '\n' )
		record->seq[len++] = ch;
	    if ( len == record->seq_array_size - 1 )
	    {
		record->seq_array_size *= 2;
		record->seq = xt_realloc(record->seq, record->seq_array_size,
		    sizeof(*record->seq));
		if ( record->seq == NULL )
		{
		    fprintf(stderr, "bl_fasta_read(): Could not reallocate seq.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	    last_ch = ch;
	}   while ( ((ch = getc(fasta_stream)) != '>') && (ch != EOF) );
	record->seq[len] = '\0';
	record->seq_len = len;
	
	if ( last_ch != '\n' )
	    fprintf(stderr, "bl_fasta_read(): Missing newline at end of seq %s.\n",
		    record->seq);

	/* Trim array */
	if ( record->seq_array_size != record->seq_len + 1 )
	{
	    record->seq_array_size = record->seq_len + 1;
	    record->seq = xt_realloc(record->seq, record->seq_array_size,
		sizeof(*record->desc));
	}
	if ( ch == '>' )
	    ungetc(ch, fasta_stream);
	return BL_READ_OK;
    }
    else
	return BL_READ_BAD_DATA;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a FASTA record to the specified FILE stream, writing at most
 *      max_line_len sequence characters per line.  The special value
 *      BL_FASTA_LINE_UNLIMITED indicates no line length limit.
 *  
 *  Arguments:
 *      fasta_stream    FILE stream to which data are written
 *      record          Pointer to a bl_fasta_t structure to be written
 *      max_line_len    Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      bl_fasta_t  rec = BL_FASTA_INIT;
 *
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *      bl_fasta_free(&rec);
 *
 *  See also:
 *      bl_fasta_read(3), bl_fastq_read(3), bl_fastq_write(3),
 *
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fasta_write(bl_fasta_t *record, FILE *fasta_stream,
	    size_t max_line_len)

{
    size_t  c;
    int     save_ch = ' ';  // Silence false uninit warning on CentOS
    
    if ( fprintf(fasta_stream, "%s\n", record->desc) < 0 )
	return BL_WRITE_FAILURE;
    
    for (c = 0; c < record->seq_len; c += max_line_len)
    {
	// Temporarily null-terminate segment of string to be printed
	if ( record->seq_len - c > max_line_len )
	{
	    save_ch = record->seq[c + max_line_len];
	    record->seq[c + max_line_len] = '\0';
	}
	
	// Print segment
	if ( fprintf(fasta_stream, "%s\n", record->seq + c) < 0 )
	    return BL_WRITE_FAILURE;
	
	// Remove temporary null-termination
	if ( record->seq_len - c > max_line_len )
	    record->seq[c + max_line_len] = save_ch;
    }
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fast.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_fasta_read()
 *  
 *  Arguments:
 *      record  Pointer to a previously populated bl_fasta_t structure
 *
 *  Examples:
 *      bl_fasta_t  rec = BL_FASTA_INIT;
 *
 *      while ( bl_fasta_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fasta_write(stdout, &rec, BL_FASTA_LINE_UNLIMITED);
 *      bl_fasta_free(&rec);
 *
 *  See also:
 *      bl_fasta_read(3), bl_fasta_write(3)
 *      bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_fasta_free(bl_fasta_t *record)

{
    free(record->seq);
    free(record->desc);
    record->desc = record->seq = NULL;
    record->desc_array_size = record->seq_array_size = 0;
    record->desc_len = record->seq_len = 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_fasta_t structure.  This must be done before
 *      passing it to bl_fasta_read() for the first time, so that
 *      bl_fasta_read() will know to allocate memory for the fields.
 *  
 *  Arguments:
 *      record  Pointer to the bl_fasta_t structure to initialize.
 *
 *  Examples:
 *      bl_fasta_t  rec;
 *
 *      bl_fasta_init(&rec);
 *      bl_fasta_read(stdin, &rec);
 *
 *  See also:
 *      bl_fasta_read(3), bl_fasta_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_fasta_init(bl_fasta_t *record)

{
    record->desc = record->seq = NULL;
    record->desc_array_size = record->seq_array_size = 0;
    record->desc_len = record->seq_len = 0;
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
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc member in a bl_fasta_t structure.
 *      Use this function to set desc in a bl_fasta_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      desc is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_desc        The new value for desc
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      char *          new_desc;
 *
 *      if ( bl_fasta_set_desc(&bl_fasta, new_desc)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_desc(
	    bl_fasta_t *bl_fasta_ptr,
	    char * new_desc
	)

{
    if ( new_desc == NULL )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->desc = new_desc;
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of desc member in a bl_fasta_t
 *      structure. Use this function to set bl_fasta_ptr->desc[c]
 *      in a bl_fasta_t object from non-member functions.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      c               Subscript to the desc array
 *      new_desc_element The new value for desc[c]
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      size_t          c;
 *      char *          new_desc_element;
 *
 *      if ( bl_fasta_set_desc_ae(&bl_fasta, c, new_desc_element)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTA_SET_DESC_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_desc_ae(
	    bl_fasta_t *bl_fasta_ptr,
	    size_t c,
	    char  new_desc_element
	)

{
    if ( false )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->desc[c] = new_desc_element;
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc member in a bl_fasta_t structure.
 *      Use this function to set desc in a bl_fasta_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_desc to bl_fasta_ptr->desc.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_desc        The new value for desc
 *      array_size      Size of the desc array.
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      char *          new_desc;
 *      size_t          array_size;
 *
 *      if ( bl_fasta_set_desc_cpy(&bl_fasta, new_desc, array_size)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTA_SET_DESC(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_desc_cpy(
	    bl_fasta_t *bl_fasta_ptr,
	    char * new_desc,
	    size_t array_size
	)

{
    if ( new_desc == NULL )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_fasta_ptr->desc, new_desc, array_size);
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_fasta_t structure.
 *      Use this function to set seq in a bl_fasta_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_seq         The new value for seq
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      char *          new_seq;
 *
 *      if ( bl_fasta_set_seq(&bl_fasta, new_seq)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_seq(
	    bl_fasta_t *bl_fasta_ptr,
	    char * new_seq
	)

{
    if ( new_seq == NULL )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->seq = new_seq;
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seq member in a bl_fasta_t
 *      structure. Use this function to set bl_fasta_ptr->seq[c]
 *      in a bl_fasta_t object from non-member functions.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      c               Subscript to the seq array
 *      new_seq_element The new value for seq[c]
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      size_t          c;
 *      char *          new_seq_element;
 *
 *      if ( bl_fasta_set_seq_ae(&bl_fasta, c, new_seq_element)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTA_SET_SEQ_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_seq_ae(
	    bl_fasta_t *bl_fasta_ptr,
	    size_t c,
	    char  new_seq_element
	)

{
    if ( false )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->seq[c] = new_seq_element;
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_fasta_t structure.
 *      Use this function to set seq in a bl_fasta_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seq to bl_fasta_ptr->seq.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_seq         The new value for seq
 *      array_size      Size of the seq array.
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      char *          new_seq;
 *      size_t          array_size;
 *
 *      if ( bl_fasta_set_seq_cpy(&bl_fasta, new_seq, array_size)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTA_SET_SEQ(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_seq_cpy(
	    bl_fasta_t *bl_fasta_ptr,
	    char * new_seq,
	    size_t array_size
	)

{
    if ( new_seq == NULL )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_fasta_ptr->seq, new_seq, array_size);
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc_array_size member in a bl_fasta_t structure.
 *      Use this function to set desc_array_size in a bl_fasta_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      desc_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_desc_array_size The new value for desc_array_size
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      size_t          new_desc_array_size;
 *
 *      if ( bl_fasta_set_desc_array_size(&bl_fasta, new_desc_array_size)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_desc_array_size(
	    bl_fasta_t *bl_fasta_ptr,
	    size_t new_desc_array_size
	)

{
    if ( false )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->desc_array_size = new_desc_array_size;
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_array_size member in a bl_fasta_t structure.
 *      Use this function to set seq_array_size in a bl_fasta_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_seq_array_size The new value for seq_array_size
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      size_t          new_seq_array_size;
 *
 *      if ( bl_fasta_set_seq_array_size(&bl_fasta, new_seq_array_size)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_seq_array_size(
	    bl_fasta_t *bl_fasta_ptr,
	    size_t new_seq_array_size
	)

{
    if ( false )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->seq_array_size = new_seq_array_size;
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc_len member in a bl_fasta_t structure.
 *      Use this function to set desc_len in a bl_fasta_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      desc_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_desc_len    The new value for desc_len
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      size_t          new_desc_len;
 *
 *      if ( bl_fasta_set_desc_len(&bl_fasta, new_desc_len)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_desc_len(
	    bl_fasta_t *bl_fasta_ptr,
	    size_t new_desc_len
	)

{
    if ( false )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->desc_len = new_desc_len;
	return BL_FASTA_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_len member in a bl_fasta_t structure.
 *      Use this function to set seq_len in a bl_fasta_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fasta_ptr    Pointer to the structure to set
 *      new_seq_len     The new value for seq_len
 *
 *  Returns:
 *      BL_FASTA_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTA_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fasta_t      bl_fasta;
 *      size_t          new_seq_len;
 *
 *      if ( bl_fasta_set_seq_len(&bl_fasta, new_seq_len)
 *              == BL_FASTA_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fasta.h
 ***************************************************************************/

int     bl_fasta_set_seq_len(
	    bl_fasta_t *bl_fasta_ptr,
	    size_t new_seq_len
	)

{
    if ( false )
	return BL_FASTA_DATA_OUT_OF_RANGE;
    else
    {
	bl_fasta_ptr->seq_len = new_seq_len;
	return BL_FASTA_DATA_OK;
    }
}
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>



/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a FASTQ record from a FILE stream.  Each record must begin
 *      with a description line (beginning with '@'), which is then
 *      followed by one or more lines of sequence data, a separator line
 *      beginning with '+', and a line of quality scores.  The end of the
 *      sequence is marked either by the next description line or EOF.
 *      If desc_len and seq_len are 0 (e.g. the structure is initialized
 *      with BL_FASTQ_INIT or bl_fastq_init(3), has been freed with
 *      bl_fastq_free(3), then memory is allocated for each line.
 *
 *      Otherwise, the existing allocated buffers are reused.  Hence, when
 *      reading many FASTQ records of the same length, only one allocation
 *      is needed.  In any case, the buffers are automatically enlarged if
 *      they become full and automatically trimmed to the actual data size
 *      after reading is complete.
 *
 *      Buffer memory should be freed as soon as possible by calling
 *      bl_fastq_free(3).
 *  
 *  Arguments:
 *      fastq_stream    FILE stream from which FASTQ data are read
 *      record          Pointer to a bl_fastq_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_write(3), bl_fastq_read(3), bl_fastq_write(3),
 *      bl_fastq_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastq_read(bl_fastq_t *record, FILE *fastq_stream)

{
    int     ch,
	    last_ch;
    size_t  len;
    
    /* Skip comment lines */
    while ( ((ch = getc(fastq_stream)) == ';') && (ch != EOF) )
	while ( ((ch = getc(fastq_stream)) != '\n') && (ch != EOF) )
	    ;
    
    if ( ch == EOF )
	return BL_READ_EOF;
    
    /* Every record should begin with a '@' */
    if ( ch == '@' )    // Desc
    {
	/*
	 *  Read description
	 */

	ungetc(ch, fastq_stream);
	ch = dsv_read_field_malloc(fastq_stream, &record->desc,
			    &record->desc_array_size, "", &record->desc_len);
	if ( record->desc == NULL )
	{
	    fprintf(stderr, "bl_fastq_read(): Could not allocate desc.\n");
	    exit(EX_UNAVAILABLE);
	}
	
	/* Should not encounter EOF while reading description line */
	/* Every description should be followed by at least one seq line */
	if ( ch == EOF )
	{
	    fprintf(stderr, "bl_fastq_read(): Record truncated in desc %s.\n",
		    record->desc);
	    return BL_READ_TRUNCATED;
	}
	else if ( ch != '\n' )
	{
	    fprintf(stderr, "bl_fastq_read(): Bad data after desc %s\n", record->desc);
	    return BL_READ_BAD_DATA;
	}

	/*
	 *  Read sequence lines.  May span multiple lines so can't use
	 *  dsv_read_field_malloc().
	 */
	
	if ( record->seq_array_size == 0 )
	{
	    // Easily hold a short read sequence and minimize reallocs for long
	    record->seq_array_size = 1024;
	    record->seq = xt_malloc(record->seq_array_size,
				    sizeof(*record->seq));
	    if ( record->seq == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	
	len = 0;
	do
	{
	    if ( ch != '\n' )
		record->seq[len++] = ch;
	    if ( len == record->seq_array_size - 1 )
	    {
		record->seq_array_size *= 2;
		record->seq = xt_realloc(record->seq, record->seq_array_size,
		    sizeof(*record->seq));
		if ( record->seq == NULL )
		{
		    fprintf(stderr, "bl_fastq_read(): Could not reallocate seq.\n");
		    exit(EX_UNAVAILABLE);
		}
	    }
	    last_ch = ch;
	}   while ( ((ch = getc(fastq_stream)) != '+') && (ch != EOF) );
	record->seq[len] = '\0';
	record->seq_len = len;

	if ( last_ch != '\n' )
	    fprintf(stderr, "bl_fasta_read(): Missing newline at end of seq %s.\n",
		    record->qual);

	/* 
	 * Trim array.  realloc() can carry a significant cost, but it does
	 * not affect overall performance here, probably because I/O is
	 * the major bottleneck.
	 */
	if ( record->seq_array_size != record->seq_len + 1 )
	{
	    record->seq_array_size = record->seq_len + 1;
	    record->seq = xt_realloc(record->seq, record->seq_array_size,
		sizeof(*record->seq));
	}

	/* Should not encounter EOF while reading sequence lines */
	/* Every sequence should be followed by a + separator line */
	if ( ch == EOF )
	{
	    fprintf(stderr, "bl_fastq_read(): Record truncated in seq %s.\n",
		    record->seq);
	    return BL_READ_TRUNCATED;
	}
	else if (ch != '+')
	{
	    fprintf(stderr, "bl_fasq_read(): Bad data after seq %s\n", record->seq);
	    return BL_READ_BAD_DATA;
	}
	// Put '+' back so it's read into plus field
	ungetc(ch, fastq_stream);
	    
	/*
	 *  Read + separator
	 */
	
	ch = dsv_read_field_malloc(fastq_stream, &record->plus,
			    &record->plus_array_size, "", &record->plus_len);
	if ( record->plus == NULL )
	{
	    fprintf(stderr, "bl_fastq_read(): Could not allocate plus.\n");
	    exit(EX_UNAVAILABLE);
	}
	
	/* Should not encounter EOF while reading plus line */
	/* Every plus should be followed by at least one qual line */
	if ( ch == EOF )
	{
	    fprintf(stderr, "bl_fastq_read(): Record truncated in plus %s.\n",
		    record->plus);
	    return BL_READ_TRUNCATED;
	}
	else if ( ch != '\n' )
	{
	    fprintf(stderr, "bl_fasq_read(): Bad data after plus %s\n",
		    record->plus);
	    return BL_READ_BAD_DATA;
	}

	/*
	 *  Read quality string.  May span multiple lines so can't use
	 *  dsv_read_field_malloc().
	 */

	// FIXME: This could be problematic with bad data where qual len
	// doesn't match seq len
	if ( record->qual_array_size == 0 )
	{
	    // Must match sequence len and no need to trim
	    record->qual_array_size = record->seq_array_size;
	    record->qual = xt_malloc(record->qual_array_size, sizeof(*record->qual));
	    if ( record->qual == NULL )
	    {
		fprintf(stderr, "bl_fastq_read(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
	
	len = 0;
	do
	{
	    /* Read at least one full line, since '@' can be a quality score */
	    while ( ((ch = getc(fastq_stream)) != '\n') && (ch != EOF) )
	    {
		record->qual[len++] = ch;
		if ( len == record->qual_array_size - 1 )
		{
		    record->qual_array_size *= 2;
		    record->qual = xt_realloc(record->qual, record->qual_array_size,
			sizeof(*record->qual));
		    if ( record->qual == NULL )
		    {
			fprintf(stderr, "bl_fastq_read(): Could not reallocate qual.\n");
			exit(EX_UNAVAILABLE);
		    }
		}
	    }
	    last_ch = ch;
	}   while ( ((ch = getc(fastq_stream)) != '@') && (ch != EOF) );
	record->qual[len] = '\0';
	record->qual_len = len;
	
	if ( last_ch != '\n' )
	    fprintf(stderr, "bl_fasta_read(): Missing newline at end of qual %s.\n",
		    record->qual);

	/*
	 *  This is where EOF should occur since we read past newlines
	 *  No need to trim since qual must be the same size as seq
	 */

	// Put '@' back so it's read into next desc
	if ( ch == '@' )
	    ungetc(ch, fastq_stream);

	return BL_READ_OK;
    }
    else
	return BL_READ_BAD_DATA;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a FASTQ record to the specified FILE stream, writing at most
 *      max_line_len sequence characters per line.  The special value
 *      BL_FASTQ_LINE_UNLIMITED indicates no line length limit.
 *  
 *  Arguments:
 *      fastq_stream    FILE stream to which data are written
 *      record          Pointer to a bl_fastq_t structure to be written
 *      max_line_len    Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_read(3), bl_fastq_write(3),
 *
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastq_write(bl_fastq_t *record, FILE *fastq_stream,
		       size_t max_line_len)

{
    size_t  c;
    int     save_ch = ' ';  // Silence false uninit warning on CentOS
    
    if ( fprintf(fastq_stream, "%s\n", record->desc) < 0 )
	return BL_WRITE_FAILURE;
    
    if ( max_line_len == BL_FASTQ_LINE_UNLIMITED )
    {
	if ( fprintf(fastq_stream, "%s\n", record->seq) < 0 )
	    return BL_WRITE_FAILURE;
    }
    else
    {
	for (c = 0; c < record->seq_len; c += max_line_len)
	{
	    // Temporarily null-terminate segment of string to be printed
	    if ( record->seq_len - c > max_line_len )
	    {
		save_ch = record->seq[c + max_line_len];
		record->seq[c + max_line_len] = '\0';
	    }
	    
	    // Print segment
	    if ( fprintf(fastq_stream, "%s\n", record->seq + c) < 0 )
		return BL_WRITE_FAILURE;
    
	    // Remove temporary null-termination
	    if ( record->seq_len - c > max_line_len )
		record->seq[c + max_line_len] = save_ch;
	}
    }
    
    if ( fprintf(fastq_stream, "%s\n", record->plus) < 0 )
	return BL_WRITE_FAILURE;
    
    if ( max_line_len == BL_FASTQ_LINE_UNLIMITED )
    {
	if ( fprintf(fastq_stream, "%s\n", record->qual) < 0 )
	    return BL_WRITE_FAILURE;
    }
    else
    {
	for (c = 0; c < record->qual_len; c += max_line_len)
	{
	    // Temporarily null-terminate segment of string to be printed
	    if ( record->qual_len - c > max_line_len )
	    {
		save_ch = record->qual[c + max_line_len];
		record->qual[c + max_line_len] = '\0';
	    }
    
	    // Print segment
	    if ( fprintf(fastq_stream, "%s\n", record->qual + c) < 0 )
		return BL_WRITE_FAILURE;
    
	    // Remove temporary null-termination
	    if ( record->qual_len - c > max_line_len )
		record->qual[c + max_line_len] = save_ch;
	}
    }
    
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_fastq_read()
 *  
 *  Arguments:
 *      record  Pointer to a previously populated bl_fastq_t structure
 *
 *  Examples:
 *      bl_fastq_t  rec = BL_FASTQ_INIT;
 *
 *      while ( bl_fastq_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastq_write(stdout, &rec, BL_FASTQ_LINE_UNLIMITED);
 *      bl_fastq_free(&rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_write(3)
 *      bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-28  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastq_free(bl_fastq_t *record)

{
    free(record->seq);
    free(record->desc);
    free(record->plus);
    free(record->qual);
    record->desc = record->seq = record->plus = record->qual = NULL;
    record->desc_array_size = record->seq_array_size = 
	record->plus_array_size = record->qual_array_size = 0;
    record->desc_len = record->seq_len = 
	record->plus_len = record->qual_len = 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_fastq_t structure.  This must be done before
 *      passing it to bl_fastq_read() for the first time, so that
 *      bl_fastq_read() will know to allocate memory for the fields.
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastq_t structure to initialize.
 *
 *  Examples:
 *      bl_fastq_t  rec;
 *
 *      bl_fastq_init(&rec);
 *      bl_fastq_read(stdin, &rec);
 *
 *  See also:
 *      bl_fastq_read(3), bl_fastq_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastq_init(bl_fastq_t *record)

{
    record->desc = record->seq = record->plus = record->qual = NULL;
    record->desc_array_size = record->seq_array_size = 
	record->plus_array_size = record->qual_array_size = 0;
    record->desc_len = record->seq_len = 
	record->plus_len = record->qual_len = 0;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Trim the 3' end of a FASTQ sequence and qualit string at location
 *      new_len.
 *  
 *  Arguments:
 *      read        FASTQ read to be trimmed
 *      new_len     New length and location of the null terminators
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if new_len is between 0 and original length,
 *      BL_FASTQ_DATA_INVALID otherwise.
 *
 *  Examples:
 *      bl_fastq_t  read;
 *      char        *adapter;
 *      size_t      index;
 *
 *      index = bl_fastq_find_adapter_smart(&read, adapter, 3, 10);
 *      if ( BL_FASTQ_SEQ_AE(&read, index) != '\0' )
 *          bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_smart(3), bl_fastq_find_adapter_exact(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_3p_trim(bl_fastq_t *read, size_t new_len)

{
    if ( (new_len >= 0) && (new_len <= read->seq_len) )
    {
	read->seq_len = read->qual_len = new_len;
	read->seq[new_len] = read->qual[new_len] = '\0';
	// FIXME: realloc?
	return BL_FASTQ_DATA_OK;
    }
    else
	return BL_FASTQ_DATA_INVALID;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Locate start of a low-quality 3' end in a FASTQ read.  This
 *      function uses the same algorithm as fastq and cutadapt as of the
 *      time of writing.  Namely, it starts at the 3' end of the quality
 *      string and sums (base quality - minimum quality) while moving in
 *      the 5' direction.  This sum will be < 0 as long as the average
 *      base quality is < minimum quality.  It also keeps track of where
 *      the minimum of this sum occurs.  When the sum become > 0, we have
 *      reached a point where the average quality of the 3' end is
 *      satisfactory, and it is assumed it will remain that way if we
 *      continue in the 5' direction.  ( Illumina reads tend to drop in
 *      quality near the 3' end. )  The location of the minimum sum is
 *      then returned, since the average quality of everything in the 5'
 *      direction must be satisfactory.
 *  
 *  Arguments:
 *      read        FASTQ read to be searched
 *      min_qual    Minimum quality of bases to keep
 *      phred_base  Offset into the ISO character set used by PHRED scores
 *                  (usually 33 for modern data)
 *
 *  Returns:
 *      Index of first low-quality base at the 3' end if found,
 *      index of NULL terminator otherwise
 *
 *  Examples:
 *      bl_fastq_t  read;
 *      
 *      ...
 *      index = bl_fastq_find_3p_low_qual(&read, 20, 33);
 *      bl_fastq_3p_trim(&read, index);
 *
 *  See also:
 *      bl_fastq_find_adapter_smart(3), bl_fastq_find_adapter_exact(3),
 *      bl_fastq_trim(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_find_3p_low_qual(const bl_fastq_t *read, unsigned min_qual,
			unsigned phred_base)

{
    ssize_t     c,
		cut_pos;
    long        sum,
		min_sum;
    
    /*
     *  Use same algorithm as BWA/cutadapt
     *  https://cutadapt.readthedocs.io/en/stable/algorithms.html#quality-trimming-algorithm
     *  score-minimum will be negative for bases we want to remove
     *  Sum score-minimum for each base starting at end until sum > 0
     *  Use the position of the minimum sum as the trim point
     *  Verified using 42, 40, 26, 27, 8, 7, 11, 4, 2, 3 example from link
     */

    if ( read->seq_len != read->qual_len )
    {
	fprintf(stderr, "bl_fastq_find_3p_low_qual(): qual_len != seq_len.\n");
	exit(EX_DATAERR);
    }
    
    sum = min_sum = 0;
    c = read->qual_len - 1;
    cut_pos = read->seq_len;
    while ( (c >= 0) && (sum <= 0) )
    {
	// Revert promotions to unsigned
	sum = (long)sum + read->qual[c] - phred_base - min_qual;
	if ( sum < min_sum )
	{
	    // fprintf(stderr, "%zu %c %c %ld\n", c, read->seq[c], read->qual[c], sum);
	    min_sum = sum;
	    cut_pos = c;
	}
	--c;
    }
    // fprintf(stderr, "Returning %zd\n", cut_pos);
    return cut_pos;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Compare the read names of two FASTQ reads.  This is useful when
 *      processing paired-end data, which must be kept in-sync.  I.e.
 *      if a sequence if removed from a 5' file, the same sequence should
 *      be removed from the 3' file whether or not it meets quality
 *      minimums.
 *  
 *  Arguments:
 *      read1, read2    FASTQ reads to compare   
 *
 *  Returns:
 *      0 if read1 and read2 have the same name
 *      < 0 if read1 name is lexically less than read2 name
 *      > 0 if read1 name is lexically greater than read2 name
 *
 *  Examples:
 *      s1 = bl_fastq_read(&fastq_rec[0], tp->instream1);
 *      s2 = bl_fastq_read(&fastq_rec[1], tp->instream2);
 *      if ( (s1 == BL_READ_OK) && (s2 == BL_READ_OK) )
 *      {
 *          if ( bl_fastq_name_cmp(&fastq_rec[0], &fastq_rec[1]) != 0 )
 *          {
 *              fprintf(stderr, "fastq-trim: Paired files out of sync.\n");
 *              trim_close_files(tp);
 *              exit(EX_DATAERR);
 *          }
 *          ...
 *      }
 *
 *  See also:
 *      bl_fastq_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-01-02  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastq_name_cmp(bl_fastq_t *read1, bl_fastq_t *read2)

{
    // FIXME: This is a hack based on test data.  Find out how to 
    // properly compare names in arbitrary FASTQ files
    // Description up to first space char is the same for R1 and R2 files
    char    *p1 = strchr(BL_FASTQ_DESC(read1), ' ');
    char    *p2 = strchr(BL_FASTQ_DESC(read2), ' ');
    int     save_p1, save_p2, status;
    
    // Temporarily null-terminate descriptions at first space char
    // Not thread-safe
    save_p1 = *p1;
    save_p2 = *p2;
    *p1 = *p2 = '\0';
    status = strcmp(BL_FASTQ_DESC(read1), BL_FASTQ_DESC(read2));
    *p1 = save_p1;
    *p2 = save_p2;
    return status;
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
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc member in a bl_fastq_t structure.
 *      Use this function to set desc in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      desc is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_desc        The new value for desc
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_desc;
 *
 *      if ( bl_fastq_set_desc(&bl_fastq, new_desc)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_desc(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_desc
	)

{
    if ( new_desc == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->desc = new_desc;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of desc member in a bl_fastq_t
 *      structure. Use this function to set bl_fastq_ptr->desc[c]
 *      in a bl_fastq_t object from non-member functions.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      c               Subscript to the desc array
 *      new_desc_element The new value for desc[c]
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          c;
 *      char *          new_desc_element;
 *
 *      if ( bl_fastq_set_desc_ae(&bl_fastq, c, new_desc_element)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_DESC_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_desc_ae(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t c,
	    char  new_desc_element
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->desc[c] = new_desc_element;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc member in a bl_fastq_t structure.
 *      Use this function to set desc in a bl_fastq_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_desc to bl_fastq_ptr->desc.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_desc        The new value for desc
 *      array_size      Size of the desc array.
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_desc;
 *      size_t          array_size;
 *
 *      if ( bl_fastq_set_desc_cpy(&bl_fastq, new_desc, array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_DESC(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_desc_cpy(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_desc,
	    size_t array_size
	)

{
    if ( new_desc == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_fastq_ptr->desc, new_desc, array_size);
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_fastq_t structure.
 *      Use this function to set seq in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_seq         The new value for seq
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_seq;
 *
 *      if ( bl_fastq_set_seq(&bl_fastq, new_seq)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_seq(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_seq
	)

{
    if ( new_seq == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->seq = new_seq;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seq member in a bl_fastq_t
 *      structure. Use this function to set bl_fastq_ptr->seq[c]
 *      in a bl_fastq_t object from non-member functions.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      c               Subscript to the seq array
 *      new_seq_element The new value for seq[c]
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          c;
 *      char *          new_seq_element;
 *
 *      if ( bl_fastq_set_seq_ae(&bl_fastq, c, new_seq_element)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_SEQ_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_seq_ae(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t c,
	    char  new_seq_element
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->seq[c] = new_seq_element;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_fastq_t structure.
 *      Use this function to set seq in a bl_fastq_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seq to bl_fastq_ptr->seq.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_seq         The new value for seq
 *      array_size      Size of the seq array.
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_seq;
 *      size_t          array_size;
 *
 *      if ( bl_fastq_set_seq_cpy(&bl_fastq, new_seq, array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_SEQ(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_seq_cpy(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_seq,
	    size_t array_size
	)

{
    if ( new_seq == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_fastq_ptr->seq, new_seq, array_size);
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for plus member in a bl_fastq_t structure.
 *      Use this function to set plus in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      plus is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_plus        The new value for plus
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_plus;
 *
 *      if ( bl_fastq_set_plus(&bl_fastq, new_plus)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_plus(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_plus
	)

{
    if ( new_plus == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->plus = new_plus;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of plus member in a bl_fastq_t
 *      structure. Use this function to set bl_fastq_ptr->plus[c]
 *      in a bl_fastq_t object from non-member functions.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      c               Subscript to the plus array
 *      new_plus_element The new value for plus[c]
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          c;
 *      char *          new_plus_element;
 *
 *      if ( bl_fastq_set_plus_ae(&bl_fastq, c, new_plus_element)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_PLUS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_plus_ae(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t c,
	    char  new_plus_element
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->plus[c] = new_plus_element;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for plus member in a bl_fastq_t structure.
 *      Use this function to set plus in a bl_fastq_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_plus to bl_fastq_ptr->plus.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_plus        The new value for plus
 *      array_size      Size of the plus array.
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_plus;
 *      size_t          array_size;
 *
 *      if ( bl_fastq_set_plus_cpy(&bl_fastq, new_plus, array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_PLUS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_plus_cpy(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_plus,
	    size_t array_size
	)

{
    if ( new_plus == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_fastq_ptr->plus, new_plus, array_size);
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual member in a bl_fastq_t structure.
 *      Use this function to set qual in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_qual        The new value for qual
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_qual;
 *
 *      if ( bl_fastq_set_qual(&bl_fastq, new_qual)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_qual(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_qual
	)

{
    if ( new_qual == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->qual = new_qual;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of qual member in a bl_fastq_t
 *      structure. Use this function to set bl_fastq_ptr->qual[c]
 *      in a bl_fastq_t object from non-member functions.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      c               Subscript to the qual array
 *      new_qual_element The new value for qual[c]
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          c;
 *      char *          new_qual_element;
 *
 *      if ( bl_fastq_set_qual_ae(&bl_fastq, c, new_qual_element)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_QUAL_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_qual_ae(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t c,
	    char  new_qual_element
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->qual[c] = new_qual_element;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual member in a bl_fastq_t structure.
 *      Use this function to set qual in a bl_fastq_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_qual to bl_fastq_ptr->qual.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_qual        The new value for qual
 *      array_size      Size of the qual array.
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      char *          new_qual;
 *      size_t          array_size;
 *
 *      if ( bl_fastq_set_qual_cpy(&bl_fastq, new_qual, array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_FASTQ_SET_QUAL(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_qual_cpy(
	    bl_fastq_t *bl_fastq_ptr,
	    char * new_qual,
	    size_t array_size
	)

{
    if ( new_qual == NULL )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_fastq_ptr->qual, new_qual, array_size);
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc_array_size member in a bl_fastq_t structure.
 *      Use this function to set desc_array_size in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      desc_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_desc_array_size The new value for desc_array_size
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_desc_array_size;
 *
 *      if ( bl_fastq_set_desc_array_size(&bl_fastq, new_desc_array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_desc_array_size(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_desc_array_size
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->desc_array_size = new_desc_array_size;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_array_size member in a bl_fastq_t structure.
 *      Use this function to set seq_array_size in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_seq_array_size The new value for seq_array_size
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_seq_array_size;
 *
 *      if ( bl_fastq_set_seq_array_size(&bl_fastq, new_seq_array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_seq_array_size(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_seq_array_size
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->seq_array_size = new_seq_array_size;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for plus_array_size member in a bl_fastq_t structure.
 *      Use this function to set plus_array_size in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      plus_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_plus_array_size The new value for plus_array_size
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_plus_array_size;
 *
 *      if ( bl_fastq_set_plus_array_size(&bl_fastq, new_plus_array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_plus_array_size(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_plus_array_size
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->plus_array_size = new_plus_array_size;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual_array_size member in a bl_fastq_t structure.
 *      Use this function to set qual_array_size in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_qual_array_size The new value for qual_array_size
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_qual_array_size;
 *
 *      if ( bl_fastq_set_qual_array_size(&bl_fastq, new_qual_array_size)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_qual_array_size(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_qual_array_size
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->qual_array_size = new_qual_array_size;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for desc_len member in a bl_fastq_t structure.
 *      Use this function to set desc_len in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      desc_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_desc_len    The new value for desc_len
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_desc_len;
 *
 *      if ( bl_fastq_set_desc_len(&bl_fastq, new_desc_len)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_desc_len(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_desc_len
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->desc_len = new_desc_len;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_len member in a bl_fastq_t structure.
 *      Use this function to set seq_len in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_seq_len     The new value for seq_len
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_seq_len;
 *
 *      if ( bl_fastq_set_seq_len(&bl_fastq, new_seq_len)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_seq_len(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_seq_len
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->seq_len = new_seq_len;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for plus_len member in a bl_fastq_t structure.
 *      Use this function to set plus_len in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      plus_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_plus_len    The new value for plus_len
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_plus_len;
 *
 *      if ( bl_fastq_set_plus_len(&bl_fastq, new_plus_len)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_plus_len(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_plus_len
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->plus_len = new_plus_len;
	return BL_FASTQ_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastq.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual_len member in a bl_fastq_t structure.
 *      Use this function to set qual_len in a bl_fastq_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastq_ptr    Pointer to the structure to set
 *      new_qual_len    The new value for qual_len
 *
 *  Returns:
 *      BL_FASTQ_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTQ_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastq_t      bl_fastq;
 *      size_t          new_qual_len;
 *
 *      if ( bl_fastq_set_qual_len(&bl_fastq, new_qual_len)
 *              == BL_FASTQ_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastq.h
 ***************************************************************************/

int     bl_fastq_set_qual_len(
	    bl_fastq_t *bl_fastq_ptr,
	    size_t new_qual_len
	)

{
    if ( false )
	return BL_FASTQ_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastq_ptr->qual_len = new_qual_len;
	return BL_FASTQ_DATA_OK;
    }
}
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>

/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a FASTA or FASTQ record from a FILE stream by calling
 *      bl_read_fasta(3) or bl_read_fastq(3).  The bl_fastx_t structure
 *      must first be initialized by assigning it BL_FASTX_INIT and
 *      calling bl_fastx_init(3).
 *      See bl_fasta_read(3) and bl_fastq_read(3) for further details.
 *
 *  Arguments:
 *      fastx_stream    FILE stream from which FASTA data are read
 *      record          Pointer to a bl_fastx_t structure to receive data
 *
 *  Returns:
 *      BL_READ_OK upon successful read of description and sequence
 *      BL_READ_BAD_DATA if something is amiss with input format
 *      BL_READ_EOF if no more data are available
 *
 *  Examples:
 *      bl_fastx_t  rec = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &rec);
 *      while ( bl_fastx_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastx_write(stdout, &rec, BL_FASTX_LINE_UNLIMITED);
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_write(3), bl_fastq_read(3), bl_fastq_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastx_read(bl_fastx_t *record, FILE *fastx_stream)

{
    switch(record->format)
    {
	case BL_FASTX_FORMAT_FASTA:
	    return bl_fasta_read(&record->fasta, fastx_stream);
	case BL_FASTX_FORMAT_FASTQ:
	    return bl_fastq_read(&record->fastq, fastx_stream);
    }
    fprintf(stderr, "bl_fastx_read(): Input format is unknown.  Call bl_fastx_init() first.\n");
    return BL_READ_UNKNOWN_FORMAT;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a FASTA or FASTQ record from a FILE stream by calling
 *      bl_fasta_write(3) or bl_fastq_write(3).  The bl_fastx_t structure
 *      must first be initialized by assigning it BL_FASTX_INIT and
 *      calling bl_fastx_init(3), and then populated by bl_fastx_read(3)
 *      or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *      See bl_fasta_write(3) and bl_fastq_write(3) for further details.
 *  
 *  Arguments:
 *      fastx_stream    FILE stream to which data are written
 *      record          Pointer to a bl_fastx_t structure to be written
 *      max_line_len    Maximum length of a sequence line in output
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE if a write error occurs.
 *
 *  Examples:
 *      bl_fastx_t  rec = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &rec);
 *      while ( bl_fastx_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastx_write(stdout, &rec, BL_FASTX_LINE_UNLIMITED);
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_read(3), bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_fastx_write(bl_fastx_t *record, FILE *fastx_stream,
	    size_t max_line_len)

{
    switch(BL_FASTX_FORMAT(record))
    {
	case BL_FASTX_FORMAT_FASTA:
	    return bl_fasta_write(&record->fasta, fastx_stream, max_line_len);
	    break;
	case BL_FASTX_FORMAT_FASTQ:
	    return bl_fastq_write(&record->fastq, fastx_stream, max_line_len);
	    break;
    }
    fprintf(stderr, "bl_fasta_write(): File format is unknown.\n");
    return BL_WRITE_FAILURE;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fast.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_fastx_read()
 *  
 *  Arguments:
 *      record  Pointer to a previously populated bl_fastx_t structure
 *
 *  Examples:
 *      bl_fastx_t  rec = BL_FASTX_INIT;
 *
 *      while ( bl_fastx_read(stdin, &rec) != BL_READ_EOF )
 *          bl_fastx_write(stdout, &rec, BL_FASTX_LINE_UNLIMITED);
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_read(3), bl_fastx_write(3)
 *      bl_fastq_read(3), bl_fastq_write(3),
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-07-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastx_free(bl_fastx_t *record)

{
    switch(BL_FASTX_FORMAT(record))
    {
	case BL_FASTX_FORMAT_FASTA:
	    bl_fasta_free(&record->fasta);
	    break;
	case BL_FASTX_FORMAT_FASTQ:
	    bl_fastq_free(&record->fastq);
	    break;
    }
    record->format = BL_FASTX_FORMAT_UNKNOWN;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fasta.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_fastx_t structure by peaking at the first character
 *      of the description string to determine whether the stream is FASTA
 *      or FASTQ, and then initializing the appropriate structure within
 *      the bl_fastx_t structure.  This must be done before
 *      passing it to bl_fastx_read() for the first time, so that
 *      bl_fastx_read() will know to allocate memory for the fields.
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure to initialize.
 *
 *  Examples:
 *      bl_fastx_t  rec = BL_FASTX_INIT;
 *
 *      bl_fastx_init(&rec);
 *      bl_fastx_read(stdin, &rec);
 *      ...
 *      bl_fastx_free(&rec);
 *
 *  See also:
 *      bl_fastx_read(3), bl_fastx_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_fastx_init(bl_fastx_t *record, FILE *fastx_stream)

{
    int     ch;
    
    if ( record->format != BL_FASTX_FORMAT_UNKNOWN )
    {
	fputs("bl_fastx_init(): Warning: format should be unknown.\n"
	      "bl_fastx_t variables should be initialized with BL_FASTX_INIT.\n"
	      "This could also indicate a previously used structure that has\n"
	      "not been freed, which is a memory leak.\n", stderr);
    }
    
    /* Skip comment lines */
    while ( ((ch = getc(fastx_stream)) == ';') && (ch != EOF) )
	while ( ((ch = getc(fastx_stream)) != '\n') && (ch != EOF) )
	    ;
    
    if ( ch == EOF )
    {
	fputs("bl_fastq_init(): EOF encountered.\n", stderr);
	exit(EX_DATAERR);
    }
    ungetc(ch, fastx_stream);
    switch(ch)
    {
	case '>':
	    record->format = BL_FASTX_FORMAT_FASTA;
	    bl_fasta_init(&record->fasta);
	    break;
	case '@':
	    record->format = BL_FASTX_FORMAT_FASTQ;
	    bl_fastq_init(&record->fastq);
	    break;
	default:
	    fprintf(stderr, "bl_fastx_init(): Unexpected first char: %c\n", ch);
	    fputs("Should be '>' or '@'.\n", stderr);
	    exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the description string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the description string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Desc string = %s\n", bl_fastx_desc(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_desc(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_DESC(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_DESC(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_desc(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the description string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the description string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Desc string length = %zu\n", bl_fastx_desc_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_desc_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_DESC_LEN(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_DESC_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_desc_len(): File format is unknown.\n");
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the seq string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the seq string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Seq string = %s\n", bl_fastx_seq(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_seq(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_SEQ(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_SEQ(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_seq(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the seq string of a FASTA or FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the seq string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Seq string length = %zu\n", bl_fastx_seq_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_seq_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    return BL_FASTA_SEQ_LEN(&record->fasta);
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_SEQ_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_seq_len(): File format is unknown.\n");
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the + string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and NULL
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the + string, or NULL if FASTA
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("+ string = %s\n", bl_fastx_plus(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_plus(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_plus(): Warning: Attempt to access + field in a FASTA stream.\n", stderr);
	    return NULL;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_PLUS(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_plus(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the + string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and 0
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the + string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("+ string length = %zu\n", bl_fastx_plus_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_plus_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_plus_len(): Warning: Attempt to access + length field in a FASTA stream.\n", stderr);
	    return 0;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_PLUS_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_plus_len(): File format is unknown.\n");
    return 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return a pointer to the qual string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and NULL
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Pointer to the qual string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Qual string = %s\n", bl_fastx_qual(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

char    *bl_fastx_qual(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_qual(): Warning: Attempt to access + field in a FASTA stream.\n", stderr);
	    return NULL;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_QUAL(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_qual(): File format is unknown.\n");
    return NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the length of the qual string of a FASTQ
 *      record.  The record must be initialized with bl_fastx_init(3)
 *      and populated by bl_fastx_read(3) or other means.  If the format
 *      if the fastx stream is FASTA, a warning is generated and 0
 *      is returned.  Previously used
 *      variables may be reused to process another record in the same
 *      format (FASTA or FASTQ) or reinitialized by bl_fastx_free(3);
 *  
 *  Arguments:
 *      record  Pointer to the bl_fastx_t structure
 *
 *  Returns:
 *      Length of the qual string
 *
 *  Examples:
 *      bl_fastx_t  record = BL_FASTX_INIT;
 *
 *      bl_fastx_init(stdin, &record);
 *      bl_fastx_read(stdin, &record);
 *      printf("Qual string length = %zu\n", bl_fastx_qual_len(&record));
 *      bl_fastx_free(&record);
 *
 *  See also:
 *      bl_fastx_init(3), bl_fastx_read(3), bl_fastx_write(3),
 *      bl_fastx_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-08-25  Jason Bacon Begin
 ***************************************************************************/

size_t  bl_fastx_qual_len(bl_fastx_t *record)

{
    switch(record->format)
    {
	case    BL_FASTX_FORMAT_FASTA:
	    fputs("bl_fastx_qual_len(): Warning: Attempt to access + length field in a FASTA stream.\n", stderr);
	    return 0;
	case    BL_FASTX_FORMAT_FASTQ:
	    return BL_FASTQ_QUAL_LEN(&record->fastq);
    }
    fprintf(stderr, "bl_fasta_qual_len(): File format is unknown.\n");
    return 0;
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
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for format member in a bl_fastx_t structure.
 *      Use this function to set format in a bl_fastx_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      format is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastx_ptr    Pointer to the structure to set
 *      new_format      The new value for format
 *
 *  Returns:
 *      BL_FASTX_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastx_t      bl_fastx;
 *      int             new_format;
 *
 *      if ( bl_fastx_set_format(&bl_fastx, new_format)
 *              == BL_FASTX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastx.h
 ***************************************************************************/

int     bl_fastx_set_format(
	    bl_fastx_t *bl_fastx_ptr,
	    int new_format
	)

{
    if ( false )
	return BL_FASTX_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastx_ptr->format = new_format;
	return BL_FASTX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for fasta member in a bl_fastx_t structure.
 *      Use this function to set fasta in a bl_fastx_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      fasta is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastx_ptr    Pointer to the structure to set
 *      new_fasta       The new value for fasta
 *
 *  Returns:
 *      BL_FASTX_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastx_t      bl_fastx;
 *      bl_fasta_t      new_fasta;
 *
 *      if ( bl_fastx_set_fasta(&bl_fastx, new_fasta)
 *              == BL_FASTX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastx.h
 ***************************************************************************/

int     bl_fastx_set_fasta(
	    bl_fastx_t *bl_fastx_ptr,
	    bl_fasta_t new_fasta
	)

{
    if ( false )
	return BL_FASTX_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastx_ptr->fasta = new_fasta;
	return BL_FASTX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/fastx.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for fastq member in a bl_fastx_t structure.
 *      Use this function to set fastq in a bl_fastx_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      fastq is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_fastx_ptr    Pointer to the structure to set
 *      new_fastq       The new value for fastq
 *
 *  Returns:
 *      BL_FASTX_DATA_OK if the new value is acceptable and assigned
 *      BL_FASTX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_fastx_t      bl_fastx;
 *      bl_fastq_t      new_fastq;
 *
 *      if ( bl_fastx_set_fastq(&bl_fastx, new_fastq)
 *              == BL_FASTX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from fastx.h
 ***************************************************************************/

int     bl_fastx_set_fastq(
	    bl_fastx_t *bl_fastx_ptr,
	    bl_fastq_t new_fastq
	)

{
    if ( false )
	return BL_FASTX_DATA_OUT_OF_RANGE;
    else
    {
	bl_fastx_ptr->fastq = new_fastq;
	return BL_FASTX_DATA_OK;
    }
}
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>
#include <inttypes.h>       // PRId64





/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over header lines in GFF input stream.  The FILE pointer
 *      gff_stream is advanced to the first character of the first line
 *      after the header.  The header is copied to a temporary file and and
 *      the function returns a FILE pointer to the header stream.
 *
 *  Arguments:
 *      gff_stream  FILE pointer to the open GFF file
 *
 *  Returns:
 *      A FILE pointer to a temporary file containing a copy of the header
 *
 *  See also:
 *      bl_gff_read(3), bl_gff_copy_header(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-08  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_gff_skip_header(FILE *gff_stream)

{
    int     ch;
    FILE    *header_stream = tmpfile();

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like peak-classifier to replicate the
     *  header in output files.
     */
    
    while ( (ch = getc(gff_stream)) == '#' )
    {
	putc(ch, header_stream);
	do
	{
	    ch = getc(gff_stream);
	    putc(ch, header_stream);
	}   while ( (ch != '\n') && (ch != EOF) );
    }
    
    // Rewind to start of first non-header line
    if ( ch != EOF )
	ungetc(ch, gff_stream);
    rewind(header_stream);
    return header_stream;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy GFF header from one FILE stream to another.  This is meant to
 *      be used in conjunction with bl_gff_skip_header(), which stores the
 *      header in a temporary file.
 *
 *  Arguments:
 *      header_stream   Open FILE stream of GFF header
 *      gff_stream      FILE stream to which header is copied
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE or BL_READ_* on failure
 *
 *  See also:
 *      bl_gff_skip_header(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-03  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_copy_header(FILE *header_stream, FILE *gff_stream)

{
    int     ch;
    
    rewind(header_stream);
    while ( (ch = getc(header_stream)) != EOF )
	if ( putc(ch, gff_stream) == EOF )
	    return BL_WRITE_FAILURE;
    rewind(header_stream);
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read next feature (line) from a GFF file.
 *
 *      feature must be initialized using BL_GFF_INIT or bl_gff_init()
 *      before being passed to this function.
 *
 *      bl_gff_read() will allocate memory for string fields as needed.
 *      The object should be passed to bl_gff_free() as soon as possible
 *      after the data are no longer needed.
 *
 *      If passed an object that is not in an initialized state,
 *      bl_gff_read() will free and initialize it before repopulating it
 *      with a new feature.
 *
 *      If field_mask is not BL_GFF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in feature.
 *      That field in the structure is then populated with an appropriate
 *      marker, such as '.'.  Possible mask values are:
 *
 *      BL_GFF_FIELD_ALL
 *      BL_GFF_FIELD_SEQID
 *      BL_GFF_FIELD_SOURCE
 *      BL_GFF_FIELD_TYPE
 *      BL_GFF_FIELD_START
 *      BL_GFF_FIELD_END
 *      BL_GFF_FIELD_SCORE
 *      BL_GFF_FIELD_STRAND
 *      BL_GFF_FIELD_PHASE
 *      BL_GFF_FIELD_ATTRIBUTES
 *
 *  Arguments:
 *      feature         Pointer to a bl_gff_t structure
 *      gff_stream      A FILE stream from which to read the line
 *      field_mask      Bit mask indicating which fields to store in feature
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered after a complete feature
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_gff_skip_header(stdin);
 *      bl_gff_init(&feature);
 *
 *      bl_gff_read(&feature, stdin, BL_GFF_FIELD_ALL);
 *      bl_gff_read(&feature, gff_stream,
 *          BL_GFF_FIELD_SEQID|BL_GFF_FIELD_START|BL_GFF_FIELD_END);
 *
 *  See also:
 *      bl_gff_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_read(bl_gff_t *feature, FILE *gff_stream,
	    gff_field_mask_t field_mask)

{
    char    *end,
	    line[BL_GFF_LINE_MAX_CHARS + 1],
	    strand_str[BL_GFF_STRAND_MAX_CHARS + 1],
	    phase_str[BL_GFF_PHASE_MAX_DIGITS + 1],
	    start_str[BL_POSITION_MAX_DIGITS + 1],
	    end_str[BL_POSITION_MAX_DIGITS + 1],
	    score_str[BL_GFF_SCORE_MAX_DIGITS + 1];
    size_t  len;
    int     delim,
	    ch;
    
    // Use this as a model for other _read() functions?
    // Makes reusing a structure easy without risk of memory leaks
    if ( feature->attributes != NULL )
	bl_gff_free(feature);
    
    // Check for group terminators (Line with just ###)
    // FIXME: Rely on parent ID instead of ###?
    if ( (ch = getc(gff_stream)) == '#' )
    {
	fgets(line, BL_GFF_LINE_MAX_CHARS, gff_stream);
	if ( strcmp(line, "##\n") == 0 )
	{
	    strlcpy(feature->type, "###", BL_GFF_TYPE_MAX_CHARS);
	    return BL_READ_OK;
	}
    }
    else if ( ch != EOF )
	ungetc(ch, gff_stream);

    feature->file_pos = ftell(gff_stream);
    
    // FIXME: Respect field_mask
    
    // 1 Chromosome
    if ( tsv_read_field(gff_stream, feature->seqid,
			BL_CHROM_MAX_CHARS, &len) == EOF )
    {
	return BL_READ_EOF;
    }
    
    // 2 Source
    if ( tsv_read_field(gff_stream, feature->source,
			BL_GFF_SOURCE_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading SOURCE: %s.\n",
		feature->source);
	return BL_READ_TRUNCATED;
    }

    // 3 Feature
    if ( tsv_read_field(gff_stream, feature->type,
			BL_GFF_TYPE_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading feature: %s.\n",
		feature->type);
	return BL_READ_TRUNCATED;
    }
    
    // 4 Feature start position
    if ( tsv_read_field(gff_stream, start_str,
			BL_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading start POS: %s.\n",
		start_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	feature->start = strtoul(start_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_gff_read(): Invalid feature position: %s\n",
		    start_str);
	    return BL_READ_TRUNCATED;
	}
    }
    
    // 5 Feature end position
    if ( tsv_read_field(gff_stream, end_str,
			BL_POSITION_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading end POS: %s.\n",
		end_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	feature->end = strtoul(end_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_gff_read(): Invalid feature position: %s\n",
		    end_str);
	    return BL_READ_TRUNCATED;
	}
    }

    // 6 Score
    if ( tsv_read_field(gff_stream, score_str,
			BL_GFF_SCORE_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading SCORE: %s.\n",
		score_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	feature->score = strtod(score_str, &end);
	if ( *end != '\0' )
	    feature->score = BL_GFF_SCORE_UNAVAILABLE;
    }
    
    
    // 7 Strand
    if ( tsv_read_field(gff_stream, strand_str,
			BL_GFF_STRAND_MAX_CHARS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading STRAND: %s.\n",
		strand_str);
	return BL_READ_TRUNCATED;
    }
    else
	feature->strand = *strand_str;
    
    // 8 Phase (bases to start of next codon: 0, 1, or 2. "." if unavailable)
    if ( tsv_read_field(gff_stream, phase_str,
			BL_GFF_PHASE_MAX_DIGITS, &len) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading PHASE: %s.\n",
		phase_str);
	return BL_READ_TRUNCATED;
    }
    else
	feature->phase = *phase_str;

    // 9 Attributes
    if ( (delim = tsv_read_field_malloc(gff_stream, &feature->attributes,
			&feature->attributes_array_size,
			&feature->attributes_len)) == EOF )
    {
	fprintf(stderr, "bl_gff_read(): Got EOF reading ATTRIBUTES: %s.\n",
		feature->attributes);
	return BL_READ_TRUNCATED;
    }
    //fprintf(stderr, "%s %zu\n", feature->attributes,
    //        strlen(feature->attributes));
    
    // printf("delim = %u\n", delim);
    if ( delim != '\n' )
	dsv_skip_rest_of_line(gff_stream);

    // Extract feature ID from attributes
    feature->feature_id = bl_gff_extract_attribute(feature, "ID");

    // Extract feature name from attributes
    feature->feature_name = bl_gff_extract_attribute(feature, "Name");
    if ( feature->feature_name == NULL )
    {
	if ( (feature->feature_name = strdup("unnamed")) == NULL )
	    fprintf(stderr, "bl_gff_read(): Could not strdup() feature_name.\n");
    }

    // Extract feature parent from attributes
    feature->feature_parent = bl_gff_extract_attribute(feature, "Parent");
    if ( feature->feature_parent == NULL )
    {
	if ( (feature->feature_parent = strdup("noparent")) == NULL )
	    fprintf(stderr, "bl_gff_read(): Could not strdup() feature_parent.\n");
    }
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write fields from a GFF feature to the specified FILE
 *      stream.
 *
 *      If field_mask is not BL_GFF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate marker for that field,
 *      such as a '.', rather than writing the real data.
 *      Possible mask values are:
 *
 *      BL_GFF_FIELD_ALL
 *      BL_GFF_FIELD_SEQID
 *      BL_GFF_FIELD_SOURCE
 *      BL_GFF_FIELD_TYPE
 *      BL_GFF_FIELD_START
 *      BL_GFF_FIELD_END
 *      BL_GFF_FIELD_SCORE
 *      BL_GFF_FIELD_STRAND
 *      BL_GFF_FIELD_PHASE
 *      BL_GFF_FIELD_ATTRIBUTES
 *
 *  Arguments:
 *      feature     Pointer to the bl_gff_t structure to output
 *      gff_stream  FILE stream to which TSV gff line is written
 *      field_mask  Bit mask indicating which fields to output
 *
 *  Returns:
 *      BL_WRITE_OK on success
 *      BL_WRITE_ERROR on failure (errno may provide more information)
 *
 *  Examples:
 *      bl_gff_write(&feature, stdout, BL_GFF_FIELD_ALL);
 *      bl_gff_write(&feature, gff_stream,
 *          BL_GFF_FIELD_SEQID|BL_GFF_FIELD_START|BL_GFF_FIELD_END);
 *
 *  See also:
 *      bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-05  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_write(bl_gff_t *feature, FILE *gff_stream,
	    gff_field_mask_t field_mask)

{
    int     printed = 0;
    
    /* FIXME: Fully test and enable this
    if ( field_mask & BL_GFF_FIELD_SEQID )
	printed += fprintf(gff_stream, "%s", feature->seqid);
    else
	printed += putc('.', gff_stream);
	
    if ( field_mask & BL_GFF_FIELD_SOURCE )
	printed += fprintf(gff_stream, "\t%s", feature->source);
    else
	printed += fprintf(gff_stream, "\t.");

    if ( field_mask & BL_GFF_FIELD_TYPE )
	printed += fprintf(gff_stream, "\t%s", feature->type);
    else
	printed += fprintf(gff_stream, "\t.");

    if ( field_mask & BL_GFF_FIELD_START )
	printed += fprintf(gff_stream, "\t%" PRId64, feature->start);
    else
	printed += fprintf(gff_stream, "\t-1");
    
    if ( field_mask & BL_GFF_FIELD_END )
	printed += fprintf(gff_stream, "\t%" PRId64, feature->end);
    else
	printed += fprintf(gff_stream, "\t-1");
    
    if ( field_mask & BL_GFF_FIELD_SCORE )
	printed += fprintf(gff_stream, "\t%f", feature->score);
    else
	printed += fprintf(gff_stream, "\t.");
    
    if ( field_mask & BL_GFF_FIELD_STRAND )
	printed += fprintf(gff_stream, "\t%c", feature->strand);
    else
	printed += fprintf(gff_stream, "\t.");
    
    if ( field_mask & BL_GFF_FIELD_PHASE )
	printed += fprintf(gff_stream, "\t%c", feature->phase);
    else
	printed += fprintf(gff_stream, "\t.");
    
    if ( field_mask & BL_GFF_FIELD_ATTRIBUTES )
	printed += fprintf(gff_stream, "\t%s", feature->attributes);
    else
	printed += fprintf(gff_stream, "\t.");
    putc('\n', gff_stream);
    */
    fprintf(gff_stream,
	"%s\t%s\t%s\t%" PRId64 "\t%" PRId64 "\t%f\t%c\t%c\t%s\n",
	feature->seqid, feature->source, feature->type,
	feature->start, feature->end, feature->score,
	feature->strand, feature->phase, feature->attributes);
    return printed;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy GFF fields to a BED structure to the extent possible.  Since
 *      GFF and BED files do not necessarily contain the same information,
 *      some information may be lost or filled in with appropriate markers.
 *
 *  Arguments:
 *      gff_feature  Pointer to the bl_gff_t structure to copy
 *      bed_feature  Pointer to the bl_bed_t structure to receive data
 *
 *  See also:
 *      bl_bed_read(3), bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-19  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_to_bed(bl_gff_t *gff_feature, bl_bed_t *bed_feature)

{
    char    name[BL_BED_NAME_MAX_CHARS + 1],
	    strand = BL_GFF_STRAND(gff_feature);
    
    // Update this if/when more fields are converted
    bl_bed_set_fields(bed_feature, 6);
    bl_bed_set_score(bed_feature, 0);
    
    bl_bed_set_chrom_cpy(bed_feature, BL_GFF_SEQID(gff_feature), BL_CHROM_MAX_CHARS + 1);
    /*
     *  BED start is 0-based and inclusive
     *  GFF is 1-based and inclusive
     */
    bl_bed_set_chrom_start(bed_feature, BL_GFF_START(gff_feature) - 1);
    /*
     *  BED end is 0-base and inclusive (or 1-based and non-inclusive)
     *  GFF is the same
     */
    bl_bed_set_chrom_end(bed_feature, BL_GFF_END(gff_feature));
    snprintf(name, BL_BED_NAME_MAX_CHARS + 1, "%s", BL_GFF_TYPE(gff_feature));
    bl_bed_set_name_cpy(bed_feature, name, BL_BED_NAME_MAX_CHARS + 1);
    bl_bed_set_score(bed_feature, 0);  // FIXME: Take as arg?
    if ( bl_bed_set_strand(bed_feature, strand) != BL_BED_DATA_OK )
    {
	fputs("bl_gff_to_bed(): bl_bed_set_strand() failed.\n", stderr);
	exit(EX_DATAERR);
    }
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated for a bl_gff_t object
 *  
 *  Arguments:
 *      feature     Pointer to the bl_gff_t object
 *
 *  Examples:
 *      bl_gff_t    feature;
 *
 *      bl_gff_free(&feature);
 *
 *  See also:
 *      bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_free(bl_gff_t *feature)

{
    if ( feature->attributes != NULL )
    {
	/*fprintf(stderr, "Freeing %s %p %zu %s\n",
		feature->type, feature->attributes,
		strlen(feature->attributes), feature->attributes);
	fflush(stderr);
	*/
	free(feature->attributes);
    }
    if ( feature->feature_id != NULL )
	free(feature->feature_id);
    if ( feature->feature_name != NULL )
	free(feature->feature_name);
    bl_gff_init(feature);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Extract an attribute value of a feature given the attribute name.
 *      Common attribute names include "ID" and "Name".  Attributes are
 *      embedded in the GFF attributes field in the form name=value;, e.g.
 *      ID=gene:ENSDARG00000029944;Name=parpbp.
 *  
 *  Arguments:
 *      Attribute name, such as "ID" or "Name"
 *
 *  Returns:
 *      Attribute value (text after '='), or NULL if name is not found
 *
 *  Examples:
 *      bl_gff_t    feature;
 *
 *      if ( bl_gff_extract_attribute(&feature, "Name") != NULL )
 *      {
 *      }
 *
 *  See also:
 *      bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-05  Jason Bacon Begin
 ***************************************************************************/

char    *bl_gff_extract_attribute(bl_gff_t *feature, const char *attr_name)

{
    char    *attribute = NULL,
	    *start,
	    *val_start,
	    *end;
    size_t  len = strlen(attr_name);

    // Find attribute beginning with "attr_name="
    for (start = feature->attributes; (*start != '\0'); ++start)
    {
	if ( (memcmp(start, attr_name, len) == 0) && (start[len] == '=') )
	{
	    val_start = start + len + 1;
	    // ; separates attributes, last one terminated by null byte
	    if ( (end = strchr(val_start, ';')) != NULL )
		*end = '\0';    // Not thread safe
	    if ( (attribute = strdup(val_start)) == NULL )
		fprintf(stderr, "%s: strdup() failed.\n", __FUNCTION__);
	    if ( end != NULL )
		*end = ';';
	}
    }
    return attribute;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_gff_t object, setting all fields to sentinel
 *      values such as 0, NULL, or '.' as appropriate.  Note that bl_gff_t
 *      objects defined as structures, not pointers to structures, can
 *      also be initialized with the BL_GFF_INIT macro.
 *  
 *  Arguments:
 *      feature     Address of a bl_gff_t structure
 *
 *  Examples:
 *      bl_gff_t    feature1 = BL_GFF_INIT,
 *                  *feature2;
 *
 *      if ( (feature2 = xt_malloc(1, sizeof(*feature2))) != NULL )
 *          bl_gff_init(feature2);
 *
 *  See also:
 *      bl_gff_read(3), bl_gff_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-16  Jason Bacon Begin
 ***************************************************************************/

void    bl_gff_init(bl_gff_t *feature)

{
    feature->seqid[0] = feature->source[0] = feature->type[0] = '.';
    feature->seqid[1] = feature->source[1] = feature->type[1] = '\0';
    feature->start = feature->end = 0;
    feature->score = 0.0;
    feature->strand = feature->phase = '.';
    feature->attributes = feature->feature_id = feature->feature_name = NULL;
    feature->attributes_array_size = feature->attributes_len = 0;
    feature->file_pos = 0;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-23  Jason Bacon Begin
 ***************************************************************************/

bl_gff_t    *bl_gff_dup(bl_gff_t *feature)

{
    bl_gff_t    *copy;
    
    if ( (copy = xt_malloc(1, sizeof(bl_gff_t))) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate new bl_gff_t object.\n",
		__FUNCTION__);
	return NULL;
    }
    bl_gff_init(copy);
    return bl_gff_copy(copy, feature);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <>
 *      -l
 *
 *  Description:
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-23  Jason Bacon Begin
 ***************************************************************************/

bl_gff_t    *bl_gff_copy(bl_gff_t *copy, bl_gff_t *feature)

{
    strlcpy(copy->seqid, feature->seqid, BL_CHROM_MAX_CHARS + 1);
    strlcpy(copy->source, feature->source, BL_GFF_SOURCE_MAX_CHARS + 1);
    strlcpy(copy->type, feature->type, BL_GFF_TYPE_MAX_CHARS + 1);
    copy->start = feature->start;
    copy->end = feature->end;
    copy->score = feature->score;
    copy->strand = feature->strand;
    copy->phase = feature->phase = '.';
    
    if ( (copy->attributes = strdup(feature->attributes)) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate attributes.\n", __FUNCTION__);
	free(copy);
	return NULL;
    }
    
    if ( feature->feature_id == NULL )
	copy->feature_id = NULL;
    else if ( (copy->feature_id = strdup(feature->feature_id)) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate attributes.\n", __FUNCTION__);
	free(copy->attributes);
	free(copy);
	return NULL;
    }

    if ( feature->feature_name == NULL )
	copy->feature_name = NULL;
    else if ( (copy->feature_name = strdup(feature->feature_name)) == NULL )
    {
	fprintf(stderr, "%s: Could not allocate attributes.\n", __FUNCTION__);
	free(copy->attributes);
	free(copy->feature_id);
	free(copy);
	return NULL;
    }
    
    copy->file_pos = feature->file_pos;
    
    return copy;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Compare the positions of a GFF feature and a SAM alignment and
 *      return a status value much like strcmp().  0 is returned if the
 *      feature and alignment overlap.  A value < 0 is returned if the
 *      feature is entirely "before" the alignment, i.e. it is on an
 *      earlier chromosome according to bl_chrom_name_cmp(3), or on the
 *      same chromosome at a lower position.  A value > 0 is returned
 *      if the feature is entirely "after" the alignment, i.e. on a later
 *      chromosome or same chromosome and higher position.
 *
 *      This function is mainly intended for programs that sweep properly
 *      sorted GFF and SAM files locating overlaps in a single pass.
 *
 *      A converse function, bl_sam_gff_cmp(3) is also provided so that
 *      the programmer can choose the more intuitive interface.
 *  
 *  Arguments:
 *      feature     Pointer to a bl_gff_t object
 *      alignment   Pointer to a bl_sam_t object
 *
 *  Returns:
 *      A value < 0 if the the feature is entirely before the alignment
 *      A value > 0 if the the feature is entirely after the alignment
 *      0 if the feature and the alignment overlap
 *
 *  See also:
 *      bl_gff_sam_cmp(3), bl_chrom_name_cmp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-06  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_sam_cmp(bl_gff_t *feature, bl_sam_t *alignment)

{
    return -bl_sam_gff_cmp(alignment, feature);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the amount of overlap between a GFF feature and a SAM
 *      alignment.
 *  
 *  Arguments:
 *      feature     Pointer to a bl_gff_t object
 *      alignment   Pointer to a bl_sam_t object
 *
 *  Returns:
 *      The number of bases of overlap between the feature and alignment.
 *      A zero or negative return value indicates no overlap.
 *
 *  See also:
 *      bl_sam_gff_overlap(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-07  Jason Bacon Begin
 ***************************************************************************/

int64_t bl_gff_sam_overlap(bl_gff_t *feature, bl_sam_t *alignment)

{
    int64_t alignment_end = BL_SAM_POS(alignment) + BL_SAM_SEQ_LEN(alignment),
	    overlap_start = XT_MAX(BL_GFF_START(feature), BL_SAM_POS(alignment)),
	    overlap_end = XT_MIN(BL_GFF_END(feature), alignment_end);
    
    //fprintf(stderr, "%" PRId64 " %" PRId64 "\n", overlap_start, overlap_end);
    //fprintf(stderr, "Coverage = %" PRId64 "\n", overlap_end - overlap_start + 1);
    return overlap_end - overlap_start + 1;
}

#include <stdio.h>
#include <stdlib.h>
#include <string.h>



/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      The gff_index_t class maintains an in-memory index of GFF
 *      features, containing the GFF fields SEQ_ID, START, and END,
 *      and the offset into the file as reported by ftell(3), or by
 *      bl_gff_read(3), which records the file position of each GFF
 *      feature it reads.
 *
 *      .B bl_gff_index_add_pos(3)
 *      adds a GFF feature with file position file_pos to the index.
 *      Features of interest, perhaps only genes or only exons, can
 *      be added to the index on-the fly while reading through a GFF
 *      file with bl_gff_read(3).
 *
 *      The index can later be searched or traversed forward or backward
 *      to quickly find
 *      the location of a feature and reposition a FILE pointer to it
 *      using fseek(3).  This system eliminates the need to inhale
 *      large numbers of GFF features into memory.
 *  
 *  Arguments:
 *      gi      Pointer to gff_index_t object to which a record will be added
 *      feature Pointer to GFF feature to be indexed
 *
 *  Returns:
 *      BL_GFF_INDEX_OK on success, BL_GFF_MALLOC_FAILED if memory could
 *      not be allocated
 *
 *  Examples:
 *      bl_gff_index_t  gi;
 *      bl_gff_t        feature;
 *
 *      if ( bl_gff_read(&feature, gff_stream, BL_GFF_FIELD_ALL) == BL_READ_OK )
 *      {
 *          if ( bl_gff_index_add(&gi, &feature) != BL_GFF_INDEX_OK )
 *              fprintf(stderr, "Error addind to GFF index.\n");
 *      }
 *
 *  See also:
 *      bl_gff_index_seek_reverse(3), bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_index_add(bl_gff_index_t *gi, bl_gff_t *feature)

{
    if ( gi->count == gi->array_size )
    {
	gi->array_size += 65536;
	gi->file_pos = xt_realloc(gi->file_pos, gi->array_size, sizeof(*gi->file_pos));
	if ( gi->file_pos == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->start = xt_realloc(gi->start, gi->array_size, sizeof(*gi->start));
	if ( gi->start == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->end = xt_realloc(gi->end, gi->array_size, sizeof(*gi->end));
	if ( gi->end == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
	gi->seqid = xt_realloc(gi->seqid, gi->array_size, sizeof(*gi->seqid));
	if ( gi->seqid == NULL )
	    return BL_GFF_INDEX_MALLOC_FAILED;
    }
    gi->file_pos[gi->count] = BL_GFF_FILE_POS(feature);
    gi->start[gi->count] = BL_GFF_START(feature);
    gi->end[gi->count] = BL_GFF_END(feature);
    
    // Hash chr to a 64-bit integer, padding with 0s beyond the end
    if ( (gi->seqid[gi->count] = strdup(BL_GFF_SEQID(feature))) == NULL )
	return BL_GFF_INDEX_MALLOC_FAILED;
    ++gi->count;
    return BL_GFF_INDEX_OK;
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      The gff_index_t class maintains an in-memory index of GFF
 *      features, containing the GFF fields SEQ_ID, START, and END,
 *      and the offset into the file as reported by ftell(3), or by
 *      bl_gff_read(3), which records the file position of each GFF
 *      feature it reads.
 *
 *      .B bl_gff_index_seek_reverse(3)
 *      moved the FILE pointer stream to feature_count indexed features
 *      upstream of feature, or to the most upstream feature within
 *      max_nt of feature.  A max_nt of 0 indicates no maximum distance,
 *      i.e. the search will proceed to feature_count features behind
 *      feature or to the beginning of the file, whichever is encontered
 *      first.
 *
 *      The max_nt parameter refers to the END of a feature, i.e.
 *      .B bl_gff_index_seek_reverse()
 *      will back up to a feature that overlaps the position of feature
 *      minus max_nt.  The START position of the feature moved to could
 *      be more than max_nt nucleotides behind the START of feature.
 *
 *      Note that this function counts only *indexed* features, i.e. those
 *      added to gi by bl_gff_index_add(3), not all features in the GFF
 *      file.  An application may only add genes to the index, for example,
 *      ignoring exons, etc.
 *  
 *  Arguments:
 *      gi              Pointer to the gff_index_t object used to search
 *      feature         Feature from which search starts
 *      feature_count   Number of indexed features to back up from feature
 *      max_nt          Maximum number of nucleotides to back up
 *
 *  Returns:
 *      The return value of fseek(), i.e. 0 upon success, -1 on error
 *
 *  Examples:
 *      bl_gff_index_t  gi;
 *      bl_gff_t        feature;
 *      long            new_pos;
 *
 *      new_pos = bl_gff_index_seek_reverse(&gi, &feature, 4, 200000);
 *
 *  See also:
 *      bl_gff_index_add(3), fseek(3), bl_gff_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-01  Jason Bacon Begin
 ***************************************************************************/

int     bl_gff_index_seek_reverse(bl_gff_index_t *gi, FILE *stream,
	    bl_gff_t *feature, int64_t feature_count, int64_t max_nt)

{
    ssize_t     c;
    char        *ref_seqid = BL_GFF_SEQID(feature);
    int64_t     ref_start = BL_GFF_START(feature),
		end = XT_MAX(ref_start - max_nt, 0),
		f;

    // First find the reference feature where the search begins
    for (c = gi->count - 1; (c >= 0) &&
			    (gi->start[c] != ref_start) &&
			    (strcmp(gi->seqid[c], ref_seqid) != 0); --c)
	;
    
    // Now back up feature_count features or to the leftmost feature
    // overlapping with the ref feature start - max_nt
    for (f = feature_count; (f > 0) &&
			    (c > 0) &&
			    (strcmp(gi->seqid[c], ref_seqid) == 0) &&
			    (gi->end[c] > end); --f, --c)
	;
    if ( gi->end[c] < end )
	++c;

    return fseek(stream, gi->file_pos[c], SEEK_SET);
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
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for array_size member in a bl_gff_index_t structure.
 *      Use this function to set array_size in a bl_gff_index_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_array_size  The new value for array_size
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      size_t          new_array_size;
 *
 *      if ( bl_gff_index_set_array_size(&bl_gff_index, new_array_size)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_array_size(
	    bl_gff_index_t *bl_gff_index_ptr,
	    size_t new_array_size
	)

{
    if ( false )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->array_size = new_array_size;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for count member in a bl_gff_index_t structure.
 *      Use this function to set count in a bl_gff_index_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_count       The new value for count
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      size_t          new_count;
 *
 *      if ( bl_gff_index_set_count(&bl_gff_index, new_count)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_count(
	    bl_gff_index_t *bl_gff_index_ptr,
	    size_t new_count
	)

{
    if ( false )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->count = new_count;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for file_pos member in a bl_gff_index_t structure.
 *      Use this function to set file_pos in a bl_gff_index_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      file_pos is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_file_pos    The new value for file_pos
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      long *          new_file_pos;
 *
 *      if ( bl_gff_index_set_file_pos(&bl_gff_index, new_file_pos)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_file_pos(
	    bl_gff_index_t *bl_gff_index_ptr,
	    long * new_file_pos
	)

{
    if ( new_file_pos == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->file_pos = new_file_pos;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of file_pos member in a bl_gff_index_t
 *      structure. Use this function to set bl_gff_index_ptr->file_pos[c]
 *      in a bl_gff_index_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      c               Subscript to the file_pos array
 *      new_file_pos_element The new value for file_pos[c]
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      size_t          c;
 *      long *          new_file_pos_element;
 *
 *      if ( bl_gff_index_set_file_pos_ae(&bl_gff_index, c, new_file_pos_element)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_FILE_POS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_file_pos_ae(
	    bl_gff_index_t *bl_gff_index_ptr,
	    size_t c,
	    long  new_file_pos_element
	)

{
    if ( false )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->file_pos[c] = new_file_pos_element;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for file_pos member in a bl_gff_index_t structure.
 *      Use this function to set file_pos in a bl_gff_index_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_file_pos to bl_gff_index_ptr->file_pos.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_file_pos    The new value for file_pos
 *      array_size      Size of the file_pos array.
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      long *          new_file_pos;
 *      size_t          array_size;
 *
 *      if ( bl_gff_index_set_file_pos_cpy(&bl_gff_index, new_file_pos, array_size)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_FILE_POS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_file_pos_cpy(
	    bl_gff_index_t *bl_gff_index_ptr,
	    long * new_file_pos,
	    size_t array_size
	)

{
    if ( new_file_pos == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_gff_index_ptr->file_pos[c] = new_file_pos[c];
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seqid member in a bl_gff_index_t structure.
 *      Use this function to set seqid in a bl_gff_index_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seqid is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_seqid       The new value for seqid
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      char **         new_seqid;
 *
 *      if ( bl_gff_index_set_seqid(&bl_gff_index, new_seqid)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_seqid(
	    bl_gff_index_t *bl_gff_index_ptr,
	    char ** new_seqid
	)

{
    if ( new_seqid == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->seqid = new_seqid;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seqid member in a bl_gff_index_t
 *      structure. Use this function to set bl_gff_index_ptr->seqid[c]
 *      in a bl_gff_index_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      c               Subscript to the seqid array
 *      new_seqid_element The new value for seqid[c]
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      size_t          c;
 *      char **         new_seqid_element;
 *
 *      if ( bl_gff_index_set_seqid_ae(&bl_gff_index, c, new_seqid_element)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_SEQID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_seqid_ae(
	    bl_gff_index_t *bl_gff_index_ptr,
	    size_t c,
	    char * new_seqid_element
	)

{
    if ( new_seqid_element == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->seqid[c] = new_seqid_element;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seqid member in a bl_gff_index_t structure.
 *      Use this function to set seqid in a bl_gff_index_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seqid to bl_gff_index_ptr->seqid.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_seqid       The new value for seqid
 *      array_size      Size of the seqid array.
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      char **         new_seqid;
 *      size_t          array_size;
 *
 *      if ( bl_gff_index_set_seqid_cpy(&bl_gff_index, new_seqid, array_size)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_SEQID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_seqid_cpy(
	    bl_gff_index_t *bl_gff_index_ptr,
	    char ** new_seqid,
	    size_t array_size
	)

{
    if ( new_seqid == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_gff_index_ptr->seqid[c] = new_seqid[c];
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for start member in a bl_gff_index_t structure.
 *      Use this function to set start in a bl_gff_index_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      start is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_start       The new value for start
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      int64_t *      new_start;
 *
 *      if ( bl_gff_index_set_start(&bl_gff_index, new_start)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_start(
	    bl_gff_index_t *bl_gff_index_ptr,
	    int64_t * new_start
	)

{
    if ( new_start == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->start = new_start;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of start member in a bl_gff_index_t
 *      structure. Use this function to set bl_gff_index_ptr->start[c]
 *      in a bl_gff_index_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      c               Subscript to the start array
 *      new_start_element The new value for start[c]
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      size_t          c;
 *      int64_t *      new_start_element;
 *
 *      if ( bl_gff_index_set_start_ae(&bl_gff_index, c, new_start_element)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_START_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_start_ae(
	    bl_gff_index_t *bl_gff_index_ptr,
	    size_t c,
	    int64_t  new_start_element
	)

{
    if ( false )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->start[c] = new_start_element;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for start member in a bl_gff_index_t structure.
 *      Use this function to set start in a bl_gff_index_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_start to bl_gff_index_ptr->start.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_start       The new value for start
 *      array_size      Size of the start array.
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      int64_t *      new_start;
 *      size_t          array_size;
 *
 *      if ( bl_gff_index_set_start_cpy(&bl_gff_index, new_start, array_size)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_START(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_start_cpy(
	    bl_gff_index_t *bl_gff_index_ptr,
	    int64_t * new_start,
	    size_t array_size
	)

{
    if ( new_start == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_gff_index_ptr->start[c] = new_start[c];
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for end member in a bl_gff_index_t structure.
 *      Use this function to set end in a bl_gff_index_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      end is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_end         The new value for end
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      int64_t *      new_end;
 *
 *      if ( bl_gff_index_set_end(&bl_gff_index, new_end)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_end(
	    bl_gff_index_t *bl_gff_index_ptr,
	    int64_t * new_end
	)

{
    if ( new_end == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->end = new_end;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of end member in a bl_gff_index_t
 *      structure. Use this function to set bl_gff_index_ptr->end[c]
 *      in a bl_gff_index_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      c               Subscript to the end array
 *      new_end_element The new value for end[c]
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      size_t          c;
 *      int64_t *      new_end_element;
 *
 *      if ( bl_gff_index_set_end_ae(&bl_gff_index, c, new_end_element)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_END_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_end_ae(
	    bl_gff_index_t *bl_gff_index_ptr,
	    size_t c,
	    int64_t  new_end_element
	)

{
    if ( false )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_index_ptr->end[c] = new_end_element;
	return BL_GFF_INDEX_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff-index.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for end member in a bl_gff_index_t structure.
 *      Use this function to set end in a bl_gff_index_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_end to bl_gff_index_ptr->end.
 *
 *  Arguments:
 *      bl_gff_index_ptr Pointer to the structure to set
 *      new_end         The new value for end
 *      array_size      Size of the end array.
 *
 *  Returns:
 *      BL_GFF_INDEX_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_INDEX_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_index_t  bl_gff_index;
 *      int64_t *      new_end;
 *      size_t          array_size;
 *
 *      if ( bl_gff_index_set_end_cpy(&bl_gff_index, new_end, array_size)
 *              == BL_GFF_INDEX_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_INDEX_SET_END(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from gff-index.h
 ***************************************************************************/

int     bl_gff_index_set_end_cpy(
	    bl_gff_index_t *bl_gff_index_ptr,
	    int64_t * new_end,
	    size_t array_size
	)

{
    if ( new_end == NULL )
	return BL_GFF_INDEX_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_gff_index_ptr->end[c] = new_end[c];
	return BL_GFF_INDEX_DATA_OK;
    }
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
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seqid member in a bl_gff_t
 *      structure. Use this function to set bl_gff_ptr->seqid[c]
 *      in a bl_gff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the seqid array
 *      new_seqid_element The new value for seqid[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char            new_seqid_element;
 *
 *      if ( bl_gff_set_seqid_ae(&bl_gff, c, new_seqid_element)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SEQID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_seqid_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_seqid_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->seqid[c] = new_seqid_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seqid member in a bl_gff_t structure.
 *      Use this function to set seqid in a bl_gff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seqid to bl_gff_ptr->seqid.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_seqid       The new value for seqid
 *      array_size      Size of the seqid array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_seqid;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_seqid_cpy(&bl_gff, new_seqid, array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SEQID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_seqid_cpy(bl_gff_t *bl_gff_ptr, char new_seqid[], size_t array_size)

{
    if ( new_seqid == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->seqid, new_seqid, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of source member in a bl_gff_t
 *      structure. Use this function to set bl_gff_ptr->source[c]
 *      in a bl_gff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the source array
 *      new_source_element The new value for source[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char            new_source_element;
 *
 *      if ( bl_gff_set_source_ae(&bl_gff, c, new_source_element)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SOURCE_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_source_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_source_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->source[c] = new_source_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for source member in a bl_gff_t structure.
 *      Use this function to set source in a bl_gff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_source to bl_gff_ptr->source.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_source      The new value for source
 *      array_size      Size of the source array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_source;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_source_cpy(&bl_gff, new_source, array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_SOURCE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_source_cpy(bl_gff_t *bl_gff_ptr, char new_source[], size_t array_size)

{
    if ( new_source == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->source, new_source, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of type member in a bl_gff_t
 *      structure. Use this function to set bl_gff_ptr->type[c]
 *      in a bl_gff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the type array
 *      new_type_element The new value for type[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char            new_type_element;
 *
 *      if ( bl_gff_set_type_ae(&bl_gff, c, new_type_element)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_TYPE_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_type_ae(bl_gff_t *bl_gff_ptr, size_t c, char new_type_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->type[c] = new_type_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for type member in a bl_gff_t structure.
 *      Use this function to set type in a bl_gff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_type to bl_gff_ptr->type.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_type        The new value for type
 *      array_size      Size of the type array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_type;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_type_cpy(&bl_gff, new_type, array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_TYPE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_type_cpy(bl_gff_t *bl_gff_ptr, char new_type[], size_t array_size)

{
    if ( new_type == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->type, new_type, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for start member in a bl_gff_t structure.
 *      Use this function to set start in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      start is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_start       The new value for start
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      int64_t         new_start;
 *
 *      if ( bl_gff_set_start(&bl_gff, new_start)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_start(bl_gff_t *bl_gff_ptr, int64_t new_start)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->start = new_start;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for end member in a bl_gff_t structure.
 *      Use this function to set end in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      end is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_end         The new value for end
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      int64_t         new_end;
 *
 *      if ( bl_gff_set_end(&bl_gff, new_end)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_end(bl_gff_t *bl_gff_ptr, int64_t new_end)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->end = new_end;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for score member in a bl_gff_t structure.
 *      Use this function to set score in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      score is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_score       The new value for score
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      double          new_score;
 *
 *      if ( bl_gff_set_score(&bl_gff, new_score)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_score(bl_gff_t *bl_gff_ptr, double new_score)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->score = new_score;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for strand member in a bl_gff_t structure.
 *      Use this function to set strand in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      strand is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_strand      The new value for strand
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_strand;
 *
 *      if ( bl_gff_set_strand(&bl_gff, new_strand)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_strand(bl_gff_t *bl_gff_ptr, char new_strand)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->strand = new_strand;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for phase member in a bl_gff_t structure.
 *      Use this function to set phase in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      phase is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_phase       The new value for phase
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char            new_phase;
 *
 *      if ( bl_gff_set_phase(&bl_gff, new_phase)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_phase(bl_gff_t *bl_gff_ptr, char new_phase)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->phase = new_phase;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes member in a bl_gff_t structure.
 *      Use this function to set attributes in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      attributes is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_attributes  The new value for attributes
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_attributes;
 *
 *      if ( bl_gff_set_attributes(&bl_gff, new_attributes)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes(bl_gff_t *bl_gff_ptr, char * new_attributes)

{
    if ( new_attributes == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->attributes = new_attributes;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of attributes member in a bl_gff_t
 *      structure. Use this function to set bl_gff_ptr->attributes[c]
 *      in a bl_gff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the attributes array
 *      new_attributes_element The new value for attributes[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char *          new_attributes_element;
 *
 *      if ( bl_gff_set_attributes_ae(&bl_gff, c, new_attributes_element)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_ATTRIBUTES_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes_ae(bl_gff_t *bl_gff_ptr, size_t c, char  new_attributes_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->attributes[c] = new_attributes_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes member in a bl_gff_t structure.
 *      Use this function to set attributes in a bl_gff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_attributes to bl_gff_ptr->attributes.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_attributes  The new value for attributes
 *      array_size      Size of the attributes array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_attributes;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_attributes_cpy(&bl_gff, new_attributes, array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_ATTRIBUTES(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes_cpy(bl_gff_t *bl_gff_ptr, char * new_attributes, size_t array_size)

{
    if ( new_attributes == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->attributes, new_attributes, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes_array_size member in a bl_gff_t structure.
 *      Use this function to set attributes_array_size in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      attributes_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_attributes_array_size The new value for attributes_array_size
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          new_attributes_array_size;
 *
 *      if ( bl_gff_set_attributes_array_size(&bl_gff, new_attributes_array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes_array_size(bl_gff_t *bl_gff_ptr, size_t new_attributes_array_size)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->attributes_array_size = new_attributes_array_size;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for attributes_len member in a bl_gff_t structure.
 *      Use this function to set attributes_len in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      attributes_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_attributes_len The new value for attributes_len
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          new_attributes_len;
 *
 *      if ( bl_gff_set_attributes_len(&bl_gff, new_attributes_len)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_attributes_len(bl_gff_t *bl_gff_ptr, size_t new_attributes_len)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->attributes_len = new_attributes_len;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_id member in a bl_gff_t structure.
 *      Use this function to set feature_id in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature_id is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_id  The new value for feature_id
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_id;
 *
 *      if ( bl_gff_set_feature_id(&bl_gff, new_feature_id)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_id(bl_gff_t *bl_gff_ptr, char * new_feature_id)

{
    if ( new_feature_id == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_id = new_feature_id;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of feature_id member in a bl_gff_t
 *      structure. Use this function to set bl_gff_ptr->feature_id[c]
 *      in a bl_gff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the feature_id array
 *      new_feature_id_element The new value for feature_id[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char *          new_feature_id_element;
 *
 *      if ( bl_gff_set_feature_id_ae(&bl_gff, c, new_feature_id_element)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_ID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_id_ae(bl_gff_t *bl_gff_ptr, size_t c, char  new_feature_id_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_id[c] = new_feature_id_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_id member in a bl_gff_t structure.
 *      Use this function to set feature_id in a bl_gff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_feature_id to bl_gff_ptr->feature_id.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_id  The new value for feature_id
 *      array_size      Size of the feature_id array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_id;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_feature_id_cpy(&bl_gff, new_feature_id, array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_ID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_id_cpy(bl_gff_t *bl_gff_ptr, char * new_feature_id, size_t array_size)

{
    if ( new_feature_id == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->feature_id, new_feature_id, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_name member in a bl_gff_t structure.
 *      Use this function to set feature_name in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature_name is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_name The new value for feature_name
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_name;
 *
 *      if ( bl_gff_set_feature_name(&bl_gff, new_feature_name)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_name(bl_gff_t *bl_gff_ptr, char * new_feature_name)

{
    if ( new_feature_name == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_name = new_feature_name;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of feature_name member in a bl_gff_t
 *      structure. Use this function to set bl_gff_ptr->feature_name[c]
 *      in a bl_gff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the feature_name array
 *      new_feature_name_element The new value for feature_name[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char *          new_feature_name_element;
 *
 *      if ( bl_gff_set_feature_name_ae(&bl_gff, c, new_feature_name_element)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_NAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_name_ae(bl_gff_t *bl_gff_ptr, size_t c, char  new_feature_name_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_name[c] = new_feature_name_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_name member in a bl_gff_t structure.
 *      Use this function to set feature_name in a bl_gff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_feature_name to bl_gff_ptr->feature_name.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_name The new value for feature_name
 *      array_size      Size of the feature_name array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_name;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_feature_name_cpy(&bl_gff, new_feature_name, array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_NAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_name_cpy(bl_gff_t *bl_gff_ptr, char * new_feature_name, size_t array_size)

{
    if ( new_feature_name == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->feature_name, new_feature_name, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_parent member in a bl_gff_t structure.
 *      Use this function to set feature_parent in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature_parent is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_parent The new value for feature_parent
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_parent;
 *
 *      if ( bl_gff_set_feature_parent(&bl_gff, new_feature_parent)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_parent(bl_gff_t *bl_gff_ptr, char * new_feature_parent)

{
    if ( new_feature_parent == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_parent = new_feature_parent;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of feature_parent member in a bl_gff_t
 *      structure. Use this function to set bl_gff_ptr->feature_parent[c]
 *      in a bl_gff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      c               Subscript to the feature_parent array
 *      new_feature_parent_element The new value for feature_parent[c]
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      size_t          c;
 *      char *          new_feature_parent_element;
 *
 *      if ( bl_gff_set_feature_parent_ae(&bl_gff, c, new_feature_parent_element)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_PARENT_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_parent_ae(bl_gff_t *bl_gff_ptr, size_t c, char  new_feature_parent_element)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->feature_parent[c] = new_feature_parent_element;
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature_parent member in a bl_gff_t structure.
 *      Use this function to set feature_parent in a bl_gff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_feature_parent to bl_gff_ptr->feature_parent.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_feature_parent The new value for feature_parent
 *      array_size      Size of the feature_parent array.
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      char *          new_feature_parent;
 *      size_t          array_size;
 *
 *      if ( bl_gff_set_feature_parent_cpy(&bl_gff, new_feature_parent, array_size)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_GFF_SET_FEATURE_PARENT(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_feature_parent_cpy(bl_gff_t *bl_gff_ptr, char * new_feature_parent, size_t array_size)

{
    if ( new_feature_parent == NULL )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_gff_ptr->feature_parent, new_feature_parent, array_size);
	return BL_GFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for file_pos member in a bl_gff_t structure.
 *      Use this function to set file_pos in a bl_gff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      file_pos is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_gff_ptr      Pointer to the structure to set
 *      new_file_pos    The new value for file_pos
 *
 *  Returns:
 *      BL_GFF_DATA_OK if the new value is acceptable and assigned
 *      BL_GFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_gff_t        bl_gff;
 *      long            new_file_pos;
 *
 *      if ( bl_gff_set_file_pos(&bl_gff, new_file_pos)
 *              == BL_GFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-29  gen-get-set Auto-generated from gff.h
 ***************************************************************************/

int     bl_gff_set_file_pos(bl_gff_t *bl_gff_ptr, long new_file_pos)

{
    if ( false )
	return BL_GFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_gff_ptr->file_pos = new_file_pos;
	return BL_GFF_DATA_OK;
    }
}
#include <stdio.h>
#include <ctype.h>

/***************************************************************************
 *  Library:
 *      #include <biolibc/translate.h>
 *      -lbiolibc
 *
 *  Description:
 *      Locate the next start codon in stream and report its position in
 *      the input.  Reported positions are 0-based offsets from the file
 *      position at the time of the call.  Hence, to find the absolute
 *      positions of codons within a file stream across multiple calls to
 *      next_start_codon() or next_stop_codon(), their return values should
 *      added to the previous return value + 3 (the size of the previous
 *      codon).  For example, givent the following input:
 *
 *      acaucauguucguggugacc
 *
 *      A call to next_start_codon() will return 5 (off set of AUG from the
 *      beginning of the sequence and a subsequent call to next_stop_codon()
 *      will return 7 (offset of UGA from the first base after the AUG).
 *  
 *  Arguments:
 *      rna_stream  FILE stream containing RNA sequence data
 *      codon       4-character buffer to receive codon sequence
 *                  Set to "" if no codon found
 *
 *  Returns:
 *      1-based position of the next start codon within the stream
 *      or EOF if no codon is found.
 *
 *  Examples:
 *      unsigned long   pos;
 *      char            codon[4];
 *      
 *      if ( (pos = next_start_codon(stdin)) != EOF )
 *      {
 *          printf("Start codon at %lu.\n", pos);
 *          if ( (pos = next_stop_codon(stdin, codon)) != EOF )
 *              printf("Stop codon %s at %lu.\n", codon, pos);
 *          else
 *              puts("No stop codon found.");
 *      }
 *      else
 *          puts("No start codon found.");
 *
 *  See also:
 *      next_stop_codon(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-17  Jason Bacon Begin
 ***************************************************************************/

long    bl_next_start_codon(FILE *rna_stream, char codon[4])

{
    int     ch1,
	    ch2,
	    ch3;
    long    pos = 0;
    
    // Blank and null-terminate
    codon[0] = codon[3] = '\0';
    
    /*
     *  Not actually scanning to EOF.  Returns from inside loop if
     *  start codon is found.
     *  Start codon is AUG.
     */
    
    while ( !feof(rna_stream) )
    {
	while ( ((ch1 = toupper(getc(rna_stream))) != EOF) && (ch1 != 'A') )
	    ++pos;
	if ( ch1 != EOF )
	{
	    ++pos;  // Count the 'A'
	    if ( ((ch2 = toupper(getc(rna_stream))) == 'U') || (ch2 == 'T') )
	    {
		if ( (ch3 = toupper(getc(rna_stream))) == 'G' )
		{
		    codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		    return pos - 1;
		}
		else if ( ch3 != EOF )
		{
		    ungetc(ch3, rna_stream);
		    ungetc(ch2, rna_stream);
		}
	    }
	    else if ( ch2 != EOF )
		ungetc(ch2, rna_stream);
	}
    }
    return EOF;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/translate.h>
 *      -lbiolibc
 *
 *  Description:
 *      Locate the next stop codon in stream and report its position in
 *      the input.  Reported positions are 0-based offsets from the file
 *      position at the time of the call.  Hence, to find the absolute
 *      positions of codons within a file stream across multiple calls to
 *      next_start_codon() or next_stop_codon(), their return values should
 *      added to the previous return value + 3 (the size of the previous
 *      codon).  For example, givent the following input:
 *
 *      acaucauguucguggugacc
 *
 *      A call to next_start_codon() will return 5 (off set of AUG from the
 *      beginning of the sequence and a subsequent call to next_stop_codon()
 *      will return 7 (offset of UGA from the first base after the AUG).
 *  
 *  Arguments:
 *      rna_stream  FILE stream containing RNA sequence data
 *      codon       4-character buffer to receive codon sequence
 *                  Set to "" if no codon found
 *
 *  Returns:
 *      1-based position of the next stop codon within the stream
 *      or EOF if no codon is found.
 *
 *  Examples:
 *      unsigned long   pos;
 *      char            codon[4];
 *      
 *      if ( (pos = next_start_codon(stdin)) != EOF )
 *      {
 *          printf("Start codon at %lu.\n", pos);
 *          if ( (pos = next_stop_codon(stdin, codon)) != EOF )
 *              printf("Stop codon %s at %lu.\n", codon, pos);
 *          else
 *              puts("No stop codon found.");
 *      }
 *      else
 *          puts("No start codon found.");
 *
 *  See also:
 *      next_stop_codon(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-10-17  Jason Bacon Begin
 ***************************************************************************/

long    bl_next_stop_codon(FILE *rna_stream, char codon[4])

{
    int     ch1,
	    ch2,
	    ch3;
    long    pos = 0;
    
    // Blank and null-terminate
    codon[0] = codon[3] = '\0';
    
    /*
     *  Not actually scanning to EOF.  Returns from inside loop if
     *  stop codon is found.
     *  Valid codons are UAG, UAA, UGA.
     */
    
    while ( !feof(rna_stream) )
    {
	while ( ((ch1 = toupper(getc(rna_stream))) != EOF) &&
		(ch1 != 'U') && (ch1 != 'T') )
	    ++pos;
	if ( ch1 != EOF )
	{
	    ++pos;  // Count the 'U' or 'T'
	    if ( (ch2 = toupper(getc(rna_stream))) == 'A' )
	    {
		if ( (ch3 = toupper(getc(rna_stream))) == 'G' || (ch3 == 'A') )
		{
		    codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		    return pos - 1;
		}
		else if ( ch3 != EOF )
		{
		    ungetc(ch3, rna_stream);
		    ungetc(ch2, rna_stream);
		}
	    }
	    else if ( ch2 == 'G' )
	    {
		if ( (ch3 = toupper(getc(rna_stream))) == 'A' )
		{
		    codon[0] = ch1; codon[1] = ch2; codon[2] = ch3;
		    return pos - 1;
		}
		else if ( ch3 != EOF )
		{
		    ungetc(ch3, rna_stream);
		    ungetc(ch2, rna_stream);
		}
	    }
	    else if ( ch2 != EOF )
		ungetc(ch2, rna_stream);
	}
    }
    return EOF;
}
#include <string.h>
#include <sys/stat.h>


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Set all fields in a bl_overlap_t structure.  Start and end
 *      positions are 1-based regardless of the feature type.  (BED file
 *      positions are 0-based and must be adjusted before passed to this
 *      function.)
 *
 *  Arguments:
 *      feature1_len      Length of feature 1
 *      feature2_len      Length of feature 2
 *      overlap_start    Start position of overlap relative to start of feature 1
 *      overlap_end      End position of overlap relative to start of feature 1
 *
 *  Returns:
 *      BL_DATA_OK upon success.
 *      BL_DATA_INVALID is arguments don't make sense.
 *
 *  Examples:
 *          bed_start = BL_BED_CHROM_START(bed_feature);
 *          bed_end = BL_BED_CHROM_END(bed_feature);
 *          gff_start = BL_GFF_CHROM_START(gff_feature);
 *          gff_end = BL_GFF_CHROM_END(gff_feature);
 *          bed_len = bed_end - bed_start;
 *          gff_len = gff_end - gff_start + 1;
 *          bl_overlap_set_all(overlap, bed_len, gff_len,
 *                          XT_MAX(bed_start+1, gff_start),
 *                          XT_MIN(bed_end, gff_end));
 *
 *  See also:
 *      bl_overlap_print(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_overlap_set_all(bl_overlap_t *overlap,
			int64_t feature1_len, int64_t feature2_len,
			int64_t overlap_start, int64_t overlap_end)

{
    overlap->feature1_len = feature1_len;
    overlap->feature2_len = feature2_len;
    overlap->overlap_start = overlap_start;
    overlap->overlap_end = overlap_end;
    overlap->overlap_len = overlap_end - overlap_start + 1;
    
    // FIXME: Return BL_DATA_INVALID if sanity checks fail
    return BL_OVERLAP_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Print all fields in a bl_overlap_t structure for debugging.
 *
 *  Arguments:
 *      stream      FILE stream to which data are printed (e.g. stderr)
 *      overlap     Address of a bl_overlap_t structure
 *      feature1_name     Name of field 1 to print with data
 *      feature2_name     Name of field 2 to print with data
 *
 *  Returns:
 *      Return status from fprintf(3)
 *
 *  See also:
 *      bl_overlap_set_all(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_overlap_print(bl_overlap_t *overlap, FILE *stream,
			  char *feature1_name, char *feature2_name)

{
    char    feature1_len[16], feature2_len[16];
    
    strlcpy(feature1_len, feature1_name, 12);
    strlcat(feature1_len, " len", 16);
    strlcpy(feature2_len, feature2_name, 12);
    strlcat(feature2_len, " len", 16);
    return fprintf(stream, "%-16s: %" PRId64 "\n"
	   "%-16s: %" PRId64 "\n"
	   "Overlap start   : %" PRId64 "\n"
	   "Overlap end     : %" PRId64 "\n"
	   "Overlap length  : %" PRId64 "\n",
	   feature1_len, overlap->feature1_len,
	   feature2_len, overlap->feature2_len,
	   overlap->overlap_start, overlap->overlap_end, overlap->overlap_len);
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
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature1_len member in a bl_overlap_t structure.
 *      Use this function to set feature1_len in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature1_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_feature1_len The new value for feature1_len
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_feature1_len;
 *
 *      if ( bl_overlap_set_feature1_len(&bl_overlap, new_feature1_len)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_feature1_len(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_feature1_len
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->feature1_len = new_feature1_len;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for feature2_len member in a bl_overlap_t structure.
 *      Use this function to set feature2_len in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      feature2_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_feature2_len The new value for feature2_len
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_feature2_len;
 *
 *      if ( bl_overlap_set_feature2_len(&bl_overlap, new_feature2_len)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_feature2_len(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_feature2_len
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->feature2_len = new_feature2_len;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for overlap_start member in a bl_overlap_t structure.
 *      Use this function to set overlap_start in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      overlap_start is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_overlap_start The new value for overlap_start
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_overlap_start;
 *
 *      if ( bl_overlap_set_overlap_start(&bl_overlap, new_overlap_start)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_overlap_start(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_overlap_start
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->overlap_start = new_overlap_start;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for overlap_end member in a bl_overlap_t structure.
 *      Use this function to set overlap_end in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      overlap_end is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_overlap_end The new value for overlap_end
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_overlap_end;
 *
 *      if ( bl_overlap_set_overlap_end(&bl_overlap, new_overlap_end)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_overlap_end(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_overlap_end
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->overlap_end = new_overlap_end;
	return BL_OVERLAP_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/overlap.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for overlap_len member in a bl_overlap_t structure.
 *      Use this function to set overlap_len in a bl_overlap_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      overlap_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_overlap_ptr  Pointer to the structure to set
 *      new_overlap_len The new value for overlap_len
 *
 *  Returns:
 *      BL_OVERLAP_DATA_OK if the new value is acceptable and assigned
 *      BL_OVERLAP_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_overlap_t    bl_overlap;
 *      int64_t        new_overlap_len;
 *
 *      if ( bl_overlap_set_overlap_len(&bl_overlap, new_overlap_len)
 *              == BL_OVERLAP_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from overlap.h
 ***************************************************************************/

int     bl_overlap_set_overlap_len(
	    bl_overlap_t *bl_overlap_ptr,
	    int64_t new_overlap_len
	)

{
    if ( false )
	return BL_OVERLAP_DATA_OUT_OF_RANGE;
    else
    {
	bl_overlap_ptr->overlap_len = new_overlap_len;
	return BL_OVERLAP_DATA_OK;
    }
}
#include <stdio.h>
#include <stdlib.h>
#include <sysexits.h>
#include <string.h>


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a position list with an initial array size of
 *      array_size.  The size will be increased by bl_pos_list_add_position()
 *      if necessary.
 *
 *  Arguments:
 *      pos_list    Pointer to the bl_pos_list_t structure to initialize
 *      array_size  Initial size of the array of positions
 *
 *  See also:
 *      bl_pos_list_add_position(3), bl_pos_list_free(3), bl_pos_list_from_csv(3),
 *      bl_pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_pos_list_allocate(bl_pos_list_t *pos_list, size_t array_size)

{
    if ( (pos_list->count != 0) || (pos_list->array_size != 0) ||
	 (pos_list->positions != NULL) )
    {
	fputs("bl_pos_list_allocate(): List is not blank.\n", stderr);
	fputs("Was it previously allocated?\n", stderr);
	fputs("Did you forget to initialize it with POS_LIST_INIT?\n", stderr);
	exit(EX_SOFTWARE);
    }
    pos_list->positions = xt_malloc(array_size, sizeof(*pos_list->positions));
    if ( pos_list->positions == NULL )
    {
	fputs("bl_pos_list_allocate(): Could not allocate positions.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    pos_list->array_size = array_size;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free the array of positions in position list and reset array_size
 *      and position count to 0.
 *
 *  Arguments:
 *      pos_list    Pointer to the position_list_t structure to reset
 *
 *  See also:
 *      bl_pos_list_allocate(3), bl_pos_list_add_position(3), bl_pos_list_from_csv(3),
 *      bl_pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_pos_list_free(bl_pos_list_t *pos_list)

{
    if ( pos_list->positions == NULL )
    {
	fputs("bl_pos_list_free(): List pointer is NULL.\n", stderr);
	fputs("Was it previously allocated?\n", stderr);
	exit(EX_SOFTWARE);
    }
    pos_list->count = 0;
    pos_list->array_size = 0;
    free(pos_list->positions);
    pos_list->positions = NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Add another position to a pos_list, expanding the array if needed.
 *
 *  Arguments:
 *      pos_list    Pointer to the bl_pos_list_t structure to add to
 *      position    New position to be added to the list
 *
 *  See also:
 *      bl_pos_list_allocate(3), bl_pos_list_free(3), bl_pos_list_from_csv(3),
 *      bl_pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     bl_pos_list_add_position(bl_pos_list_t *pos_list, int64_t position)

{
    if ( pos_list->count == pos_list->array_size )
    {
	pos_list->array_size *= 2;
	pos_list->positions = xt_realloc(pos_list->positions,
	    pos_list->array_size, sizeof(*pos_list->positions));
	if ( pos_list == NULL )
	{
	    fputs("bl_pos_list_add_position(): Could not reallocate positions.\n", stderr);
	    exit(EX_UNAVAILABLE);
	}
    }
    pos_list->positions[pos_list->count++] = position;
    return BL_POS_LIST_DATA_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Convert a comma-separated list of positions to a bl_pos_list_t list.
 *      The array_size argument should be your best guess at the final size
 *      of the list.  Choosing a large enough value is not critical since
 *      it will be extended by bl_pos_list_add_position() if needed.
 *
 *  Arguments:
 *      pos_list    Pointer to the bl_pos_list_t to receive the list
 *      list_str    Character string containing comma-separated list of positions
 *      array_size  Initial array size
 *
 *  Returns:
 *      The number of positions added to the list
 *
 *  See also:
 *      bl_pos_list_allocate(3), bl_pos_list_add_position(3), bl_pos_list_free(3),
 *      bl_pos_list_sort(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     bl_pos_list_from_csv(bl_pos_list_t *pos_list, const char *bounds_str,
		       size_t array_size)

{
    char        *copy, *p, *token, *end;
    size_t      c;
    int64_t    position;
    
    if ( (copy = strdup(bounds_str)) == NULL )
    {
	fputs("peak-classifier: Cannot allocate temporary bounds string.\n", stderr);
	exit(EX_UNAVAILABLE);
    }
    bl_pos_list_allocate(pos_list, array_size);
    for (p = copy, c = 0; (c < BL_POS_LIST_ARRAY_SIZE(pos_list)) &&
			  ((token = strsep(&p, ",")) != NULL); ++c)
    {
	position = strtoull(token, &end, 10);
	if ( *end != '\0' )
	    return BL_POS_LIST_DATA_INVALID;
	else
	    bl_pos_list_add_position(pos_list, position);
    }
    return c;
}


/***************************************************************************
 *  Compare two int64_t values for sort functions.  Difference may
 *  exceed the range of an int, so don't just subtract.
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

int     position_cmp_ascending(const int64_t *pos1, const int64_t *pos2)

{
    if ( *pos1 == *pos2 )
	return 0;
    else if ( *pos1 > *pos2 )
	return 1;
    else
	return -1;
}


int     position_cmp_descending(const int64_t *pos1, const int64_t *pos2)

{
    return -position_cmp_ascending(pos1, pos2);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Sort a position list in either ascending or descending order.
 *
 *  Arguments:
 *      pos_list    Pointer to the position_list_t structure to sort
 *      order       POS_LIST_ASCENDING or POS_LIST_DESCENDING
 *
 *  See also:
 *      bl_pos_list_allocate(3), bl_pos_list_add_position(3), bl_pos_list_free(3),
 *      bl_pos_list_from_csv(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2021-04-17  Jason Bacon Begin
 ***************************************************************************/

void    bl_pos_list_sort(bl_pos_list_t *pos_list, bl_pos_list_sort_order_t order)

{
    if ( order == BL_POS_LIST_ASCENDING )
	qsort(pos_list->positions, pos_list->count, sizeof(pos_list->positions[0]),
	     (int (*)(const void *,const void *))position_cmp_ascending);
    else
	qsort(pos_list->positions, pos_list->count, sizeof(pos_list->positions[0]),
	     (int (*)(const void *,const void *))position_cmp_descending);
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
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for array_size member in a bl_pos_list_t structure.
 *      Use this function to set array_size in a bl_pos_list_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the structure to set
 *      new_array_size  The new value for array_size
 *
 *  Returns:
 *      BL_POS_LIST_DATA_OK if the new value is acceptable and assigned
 *      BL_POS_LIST_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      size_t          new_array_size;
 *
 *      if ( bl_pos_list_set_array_size(&bl_pos_list, new_array_size)
 *              == BL_POS_LIST_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_array_size(
	    bl_pos_list_t *bl_pos_list_ptr,
	    size_t new_array_size
	)

{
    if ( false )
	return BL_POS_LIST_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->array_size = new_array_size;
	return BL_POS_LIST_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for count member in a bl_pos_list_t structure.
 *      Use this function to set count in a bl_pos_list_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the structure to set
 *      new_count       The new value for count
 *
 *  Returns:
 *      BL_POS_LIST_DATA_OK if the new value is acceptable and assigned
 *      BL_POS_LIST_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      size_t          new_count;
 *
 *      if ( bl_pos_list_set_count(&bl_pos_list, new_count)
 *              == BL_POS_LIST_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_count(
	    bl_pos_list_t *bl_pos_list_ptr,
	    size_t new_count
	)

{
    if ( false )
	return BL_POS_LIST_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->count = new_count;
	return BL_POS_LIST_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for positions member in a bl_pos_list_t structure.
 *      Use this function to set positions in a bl_pos_list_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      positions is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the structure to set
 *      new_positions   The new value for positions
 *
 *  Returns:
 *      BL_POS_LIST_DATA_OK if the new value is acceptable and assigned
 *      BL_POS_LIST_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      int64_t *      new_positions;
 *
 *      if ( bl_pos_list_set_positions(&bl_pos_list, new_positions)
 *              == BL_POS_LIST_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_positions(
	    bl_pos_list_t *bl_pos_list_ptr,
	    int64_t * new_positions
	)

{
    if ( new_positions == NULL )
	return BL_POS_LIST_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->positions = new_positions;
	return BL_POS_LIST_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of positions member in a bl_pos_list_t
 *      structure. Use this function to set bl_pos_list_ptr->positions[c]
 *      in a bl_pos_list_t object from non-member functions.
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the structure to set
 *      c               Subscript to the positions array
 *      new_positions_element The new value for positions[c]
 *
 *  Returns:
 *      BL_POS_LIST_DATA_OK if the new value is acceptable and assigned
 *      BL_POS_LIST_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      size_t          c;
 *      int64_t *      new_positions_element;
 *
 *      if ( bl_pos_list_set_positions_ae(&bl_pos_list, c, new_positions_element)
 *              == BL_POS_LIST_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_POS_LIST_SET_POSITIONS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_positions_ae(
	    bl_pos_list_t *bl_pos_list_ptr,
	    size_t c,
	    int64_t  new_positions_element
	)

{
    if ( false )
	return BL_POS_LIST_DATA_OUT_OF_RANGE;
    else
    {
	bl_pos_list_ptr->positions[c] = new_positions_element;
	return BL_POS_LIST_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/pos-list.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for positions member in a bl_pos_list_t structure.
 *      Use this function to set positions in a bl_pos_list_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_positions to bl_pos_list_ptr->positions.
 *
 *  Arguments:
 *      bl_pos_list_ptr Pointer to the structure to set
 *      new_positions   The new value for positions
 *      array_size      Size of the positions array.
 *
 *  Returns:
 *      BL_POS_LIST_DATA_OK if the new value is acceptable and assigned
 *      BL_POS_LIST_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_pos_list_t   bl_pos_list;
 *      int64_t *      new_positions;
 *      size_t          array_size;
 *
 *      if ( bl_pos_list_set_positions_cpy(&bl_pos_list, new_positions, array_size)
 *              == BL_POS_LIST_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_POS_LIST_SET_POSITIONS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from pos-list.h
 ***************************************************************************/

int     bl_pos_list_set_positions_cpy(
	    bl_pos_list_t *bl_pos_list_ptr,
	    int64_t * new_positions,
	    size_t array_size
	)

{
    if ( new_positions == NULL )
	return BL_POS_LIST_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_pos_list_ptr->positions[c] = new_positions[c];
	return BL_POS_LIST_DATA_OK;
    }
}
#include <sysexits.h>
#include <stdlib.h>
#include <string.h>
#include <inttypes.h>




/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Check that the newly read SAM alignment comes after the previous
 *      one, assuming the input is sorted first by chrom and then
 *      position.  The previous chrom and position are stored in
 *      sam_buff (and initialized so that the first SAM alignment read is
 *      always OK).
 *  
 *  Arguments:
 *      sam_buff        Pointer to a SAM buffer with recent alignments
 *      sam_alignment   Pointer to the most recently read alignment
 *
 *  See also:
 *      bl_sam_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_check_order(bl_sam_buff_t *sam_buff,
			     bl_sam_t *sam_alignment)

{
    /*fprintf(stderr, "Previous SAM: %s %zu, Current SAM: %s %zu\n",
	    sam_buff->previous_rname, sam_buff->previous_pos,
	    sam_alignment->rname, sam_alignment->pos);*/
    if ( strcmp(sam_alignment->rname, sam_buff->previous_rname) == 0 )
    {
	// Silly to assign when already ==, but sillier to add another check
	if (sam_alignment->pos < sam_buff->previous_pos )
	    bl_sam_buff_out_of_order(sam_buff, sam_alignment);
	else
	    sam_buff->previous_pos = sam_alignment->pos;
    }
    else if ( bl_chrom_name_cmp(sam_alignment->rname,
				  sam_buff->previous_rname) < 0 )
	bl_sam_buff_out_of_order(sam_buff, sam_alignment);
    else
    {
	strlcpy(sam_buff->previous_rname, sam_alignment->rname, BL_SAM_RNAME_MAX_CHARS);
	sam_buff->previous_pos = sam_alignment->pos;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a SAM alignment buffer for holding recently read SAM
 *      alignments.  This is useful, for example, when scanning a SAM
 *      stream for alignments overlapping a certain region or position.
 *      The buffer array is set to a
 *      reasonable initial size and extended as far as BL_SAM_BUFF_MAX_SIZE
 *      by bl_sam_buff_add_alignment(3) if needed.  A minimum MAPQ value
 *      is stored in the bl_sam_buff_t structure for filtering with
 *      bl_sam_buff_alignment_ok(3).
 *  
 *  Arguments:
 *      sam_buff    Pointer to a the bl_sam_buff_t structure to initialize
 *      mapq_min    User-selected minimum MAPQ value
 *
 *  See also:
 *      bl_sam_buff_check_order(3), bl_sam_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_init(bl_sam_buff_t *sam_buff, unsigned int mapq_min,
			 size_t max_alignments)

{
    size_t  c;
    
    sam_buff->buff_size = BL_SAM_BUFF_START_SIZE;
    sam_buff->max_alignments = max_alignments;
    sam_buff->buffered_count = 0;
    sam_buff->max_count = 0;
    sam_buff->previous_pos = 0;
    *sam_buff->previous_rname = '\0';
    
    sam_buff->mapq_min = mapq_min;
    sam_buff->mapq_low = UINT64_MAX;
    sam_buff->mapq_high = 0;
    sam_buff->mapq_sum = 0;
    sam_buff->reads_used = 0;

    sam_buff->total_alignments = 0;
    sam_buff->trailing_alignments = 0;
    sam_buff->discarded_alignments = 0;
    sam_buff->discarded_score_sum = 0;
    sam_buff->min_discarded_score = SIZE_MAX;
    sam_buff->max_discarded_score = 0;
    sam_buff->discarded_trailing = 0;
    sam_buff->unmapped_alignments = 0;
    
    /*
     *  Dynamically allocating the pointers is probably senseless since they
     *  take very little space compared to the alignment data.  By the time
     *  the pointer array takes a significant amount of RAM, you're probably
     *  already thrashing to accomodate the sequence data.  The pointer array
     *  size is capped by BL_SAM_BUFF_MAX_SIZE to prevent memory exhaustion.
     *  We may save a few megabytes with this, though.
     */
    sam_buff->alignments =
	(bl_sam_t **)xt_malloc(sam_buff->buff_size,
				   sizeof(bl_sam_t **));
    for (c = 0; c < sam_buff->buff_size; ++c)
	sam_buff->alignments[c] = NULL;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Add a new alignment to the buffer, expanding the array as needed
 *      up to BL_SAM_BUFF_MAX_SIZE.
 *  
 *  Arguments:
 *      sam_buff    Pointer to bl_sam_buff_t structure where alignments are buffered
 *      sam_alignment   New SAM alignment to add to buffer
 *
 *  See also:
 *      bl_sam_buff_init(3), bl_sam_buff_check_order(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_buff_add_alignment(bl_sam_buff_t *sam_buff,
			       bl_sam_t *sam_alignment)

{
    size_t  old_buff_size,
	    c;

    bl_sam_buff_check_order(sam_buff, sam_alignment);
    
    sam_buff->mapq_low = XT_MIN(sam_buff->mapq_low, BL_SAM_MAPQ(sam_alignment));
    sam_buff->mapq_high = XT_MAX(sam_buff->mapq_high, BL_SAM_MAPQ(sam_alignment));
    sam_buff->mapq_sum += BL_SAM_MAPQ(sam_alignment);
    ++sam_buff->reads_used;

    // Just allocate the static fields, bl_sam_copy() does the rest
    if ( sam_buff->alignments[sam_buff->buffered_count] == NULL )
    {
	//fprintf(stderr, "Allocating alignment #%zu\n", sam_buff->buffered_count);
	sam_buff->alignments[sam_buff->buffered_count] = 
	    xt_malloc(1, sizeof(bl_sam_t));
	if ( sam_buff->alignments[sam_buff->buffered_count] == NULL )
	{
	    fprintf(stderr, "bl_sam_buff_add_alignment(): Could not allocate alignments.\n");
	    exit(EX_UNAVAILABLE);
	}
	// Redundant to bl_sam_copy()
	// bl_sam_init(sam_buff->alignments[sam_buff->buffered_count], 0);
    }
    else
	bl_sam_free(sam_buff->alignments[sam_buff->buffered_count]);

    //fprintf(stderr, "Adding alignment #%zu...\n", sam_buff->buffered_count);
    //fprintf(stderr, "buff_size = %zu\n", sam_buff->buff_size);
    bl_sam_copy(sam_buff->alignments[sam_buff->buffered_count], sam_alignment);
    
    ++sam_buff->buffered_count;

    if ( sam_buff->buffered_count > sam_buff->max_count )
    {
	sam_buff->max_count = sam_buff->buffered_count;
	// fprintf(stderr, "sam_buff->max_count = %zu\n", sam_buff->max_count);
    }
    
    if ( sam_buff->buffered_count == sam_buff->max_alignments )
    {
	fprintf(stderr,
		"bl_sam_buff_add_alignment(): Hit maximum alignments=%zu.\n",
		sam_buff->max_alignments);
	fprintf(stderr, "Aborting add to prevent runaway memory use.\n");
	fprintf(stderr, "Check your SAM input.\n");
	return BL_SAM_BUFF_ADD_FAILED;
    }
    
    if ( sam_buff->buffered_count == sam_buff->buff_size )
    {
	fprintf(stderr,
		"bl_sam_buff_add_alignment(): Hit buff_size=%zu, doubling buffer size.\n",
		sam_buff->buff_size);
	fprintf(stderr, "RNAME: %s  POS: %" PRId64 " LEN: %zu\n",
		BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
		BL_SAM_SEQ_LEN(sam_alignment));
	old_buff_size = sam_buff->buff_size;
	sam_buff->buff_size *= 2;
	sam_buff->alignments =
	    (bl_sam_t **)xt_realloc(sam_buff->alignments,
					sam_buff->buff_size,
					sizeof(bl_sam_t **));
	for (c = old_buff_size; c < sam_buff->buff_size; ++c)
	    sam_buff->alignments[c] = NULL;
    }
    return BL_SAM_BUFF_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Report SAM input sort error and terminate process.
 *  
 *  Arguments:
 *      sam_buff        Pointer to bl_sam_buff_t structure
 *      sam_alignment   Offending alignment out of order with previous
 *
 *  Returns:
 *      Does not return.
 *
 *  See also:
 *      bl_sam_buff_alignment_ok(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_out_of_order(bl_sam_buff_t *sam_buff, bl_sam_t *sam_alignment)

{
    fprintf(stderr, "Error: SAM input must be sorted by chrom and then position.\n");
    fprintf(stderr, "Found %s,%" PRId64 " after %s,%" PRId64 ".\n",
	    BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
	    sam_buff->previous_rname, sam_buff->previous_pos);
    exit(EX_DATAERR);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free an element of the SAM alignment array by first freeing all
 *      memory allocated by the bl_sam_t structure and then freeing
 *      memory allocated for the structure itself.
 *  
 *  Arguments:
 *      sam_buff    Pointer to the bl_sam_buff_t structure holding alignments
 *      c           Index of the alignment to be freed (0-based)
 *
 *  See also:
 *      bl_sam_buff_init(3), bl_sam_buff_add_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_free_alignment(bl_sam_buff_t *sam_buff, size_t c)

{
    bl_sam_free(sam_buff->alignments[c]);
    bl_sam_init(sam_buff->alignments[c]);
    if ( sam_buff->alignments[c] != NULL )
    {
	free(sam_buff->alignments[c]);
	sam_buff->alignments[c] = NULL;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free nelem SAM alignments at the head of the queue and shift
 *      remaining elements forward nelem positions.
 *  
 *  Arguments:
 *      sam_buff    Pointer to bl_sam_buff_t structure holding alignments
 *      nelem       Number of alignments to free
 *
 *  See also:
 *      bl_sam_buff_free_alignment(3)
 *
 *  FIXME: Use circular queuing for better efficiency
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_buff_shift(bl_sam_buff_t *sam_buff, size_t nelem)

{
    size_t  c;

    /*
     *  FIXME: A circular queue would be more efficient, but won't matter
     *  much to the bottom line to avoid shifting a few pointers on top
     *  of everything else going on here
     */
    
    /* Make susre elements to be removed are freed */
    for (c = 0; c < nelem; ++c)
	bl_sam_buff_free_alignment(sam_buff, c);

    /* Shift elements */
    for (c = 0; c < sam_buff->buffered_count - nelem; ++c)
	sam_buff->alignments[c] = sam_buff->alignments[c + nelem];
    
    /* Clear vacated elements */
    while ( c < sam_buff->buffered_count )
	sam_buff->alignments[c++] = NULL;
    
    sam_buff->buffered_count -= nelem;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Verify that an alignment meets the quality requirements for the
 *      given SAM buffer.  Each bl_sam_buff_t structure contains
 *      specifications for minimum quality scores as well as statistics
 *      on alignments checked, including the number of discarded alignments
 *      and a sum of the alignment scores (for calculating the average
 *      score of discarded alignments).
 *  
 *  Arguments:
 *      sam_buff    Pointer to bl_sam_buff_t structure with quality specs
 *      sam_alignment   Pointer to new alignment
 *
 *  Returns:
 *      true if the alignment meets the quality requirements
 *      false otherwise
 *
 *  See also:
 *      bl_sam_buff_check_order(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/
 
bool    bl_sam_buff_alignment_ok(bl_sam_buff_t *sam_buff,
			      bl_sam_t *sam_alignment)

{
    if ( sam_alignment->flag & BAM_FUNMAP )
    {
	++sam_buff->unmapped_alignments;
#ifdef DEBUG
	fprintf(stderr, "Discarding unmapped read: %s,%zu,0x%04x\n",
		BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
		BL_SAM_FLAG(sam_alignment));
#endif
	return false;
    }
    else if ( BL_SAM_MAPQ(sam_alignment) < BL_SAM_BUFF_MAPQ_MIN(sam_buff) )
    {
	++sam_buff->discarded_alignments;
	sam_buff->discarded_score_sum += BL_SAM_MAPQ(sam_alignment);
	if ( BL_SAM_MAPQ(sam_alignment) < sam_buff->min_discarded_score )
	    sam_buff->min_discarded_score = BL_SAM_MAPQ(sam_alignment);
	if ( BL_SAM_MAPQ(sam_alignment) > sam_buff->max_discarded_score )
	    sam_buff->max_discarded_score = BL_SAM_MAPQ(sam_alignment);

#ifdef DEBUG
	fprintf(stderr, "bl_sam_buff_alignment_ok(): Discarding low quality alignment: %s,%zu MAPQ=%u\n",
		BL_SAM_RNAME(sam_alignment), BL_SAM_POS(sam_alignment),
		BL_SAM_MAPQ(sam_alignment));
#endif
	return false;
    }
    else
	return true;
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
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for buff_size member in a bl_sam_buff_t structure.
 *      Use this function to set buff_size in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      buff_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_buff_size   The new value for buff_size
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      size_t          new_buff_size;
 *
 *      if ( bl_sam_buff_set_buff_size(&bl_sam_buff, new_buff_size)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_buff_size(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    size_t new_buff_size
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->buff_size = new_buff_size;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for max_alignments member in a bl_sam_buff_t structure.
 *      Use this function to set max_alignments in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      max_alignments is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_max_alignments The new value for max_alignments
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      size_t          new_max_alignments;
 *
 *      if ( bl_sam_buff_set_max_alignments(&bl_sam_buff, new_max_alignments)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_max_alignments(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    size_t new_max_alignments
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->max_alignments = new_max_alignments;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for alignments member in a bl_sam_buff_t structure.
 *      Use this function to set alignments in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      alignments is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_alignments  The new value for alignments
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      bl_sam_t **     new_alignments;
 *
 *      if ( bl_sam_buff_set_alignments(&bl_sam_buff, new_alignments)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_alignments(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    bl_sam_t ** new_alignments
	)

{
    if ( new_alignments == NULL )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->alignments = new_alignments;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of alignments member in a bl_sam_buff_t
 *      structure. Use this function to set bl_sam_buff_ptr->alignments[c]
 *      in a bl_sam_buff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      c               Subscript to the alignments array
 *      new_alignments_element The new value for alignments[c]
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      size_t          c;
 *      bl_sam_t **     new_alignments_element;
 *
 *      if ( bl_sam_buff_set_alignments_ae(&bl_sam_buff, c, new_alignments_element)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_BUFF_SET_ALIGNMENTS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_alignments_ae(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    size_t c,
	    bl_sam_t * new_alignments_element
	)

{
    if ( new_alignments_element == NULL )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->alignments[c] = new_alignments_element;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for alignments member in a bl_sam_buff_t structure.
 *      Use this function to set alignments in a bl_sam_buff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_alignments to bl_sam_buff_ptr->alignments.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_alignments  The new value for alignments
 *      array_size      Size of the alignments array.
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      bl_sam_t **     new_alignments;
 *      size_t          array_size;
 *
 *      if ( bl_sam_buff_set_alignments_cpy(&bl_sam_buff, new_alignments, array_size)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_BUFF_SET_ALIGNMENTS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_alignments_cpy(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    bl_sam_t ** new_alignments,
	    size_t array_size
	)

{
    if ( new_alignments == NULL )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_sam_buff_ptr->alignments[c] = new_alignments[c];
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for buffered_count member in a bl_sam_buff_t structure.
 *      Use this function to set buffered_count in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      buffered_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_buffered_count The new value for buffered_count
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      size_t          new_buffered_count;
 *
 *      if ( bl_sam_buff_set_buffered_count(&bl_sam_buff, new_buffered_count)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_buffered_count(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    size_t new_buffered_count
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->buffered_count = new_buffered_count;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for max_count member in a bl_sam_buff_t structure.
 *      Use this function to set max_count in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      max_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_max_count   The new value for max_count
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      size_t          new_max_count;
 *
 *      if ( bl_sam_buff_set_max_count(&bl_sam_buff, new_max_count)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_max_count(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    size_t new_max_count
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->max_count = new_max_count;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for previous_pos member in a bl_sam_buff_t structure.
 *      Use this function to set previous_pos in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      previous_pos is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_previous_pos The new value for previous_pos
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_previous_pos;
 *
 *      if ( bl_sam_buff_set_previous_pos(&bl_sam_buff, new_previous_pos)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_previous_pos(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_previous_pos
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->previous_pos = new_previous_pos;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of previous_rname member in a bl_sam_buff_t
 *      structure. Use this function to set bl_sam_buff_ptr->previous_rname[c]
 *      in a bl_sam_buff_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      c               Subscript to the previous_rname array
 *      new_previous_rname_element The new value for previous_rname[c]
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      size_t          c;
 *      char            new_previous_rname_element;
 *
 *      if ( bl_sam_buff_set_previous_rname_ae(&bl_sam_buff, c, new_previous_rname_element)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_BUFF_SET_PREVIOUS_RNAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_previous_rname_ae(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    size_t c,
	    char new_previous_rname_element
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->previous_rname[c] = new_previous_rname_element;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for previous_rname member in a bl_sam_buff_t structure.
 *      Use this function to set previous_rname in a bl_sam_buff_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_previous_rname to bl_sam_buff_ptr->previous_rname.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_previous_rname The new value for previous_rname
 *      array_size      Size of the previous_rname array.
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      char            new_previous_rname;
 *      size_t          array_size;
 *
 *      if ( bl_sam_buff_set_previous_rname_cpy(&bl_sam_buff, new_previous_rname, array_size)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_BUFF_SET_PREVIOUS_RNAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_previous_rname_cpy(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    char new_previous_rname[],
	    size_t array_size
	)

{
    if ( new_previous_rname == NULL )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_buff_ptr->previous_rname, new_previous_rname, array_size);
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for mapq_min member in a bl_sam_buff_t structure.
 *      Use this function to set mapq_min in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      mapq_min is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_mapq_min    The new value for mapq_min
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_mapq_min;
 *
 *      if ( bl_sam_buff_set_mapq_min(&bl_sam_buff, new_mapq_min)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_mapq_min(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_mapq_min
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->mapq_min = new_mapq_min;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for mapq_low member in a bl_sam_buff_t structure.
 *      Use this function to set mapq_low in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      mapq_low is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_mapq_low    The new value for mapq_low
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_mapq_low;
 *
 *      if ( bl_sam_buff_set_mapq_low(&bl_sam_buff, new_mapq_low)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_mapq_low(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_mapq_low
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->mapq_low = new_mapq_low;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for mapq_high member in a bl_sam_buff_t structure.
 *      Use this function to set mapq_high in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      mapq_high is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_mapq_high   The new value for mapq_high
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_mapq_high;
 *
 *      if ( bl_sam_buff_set_mapq_high(&bl_sam_buff, new_mapq_high)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_mapq_high(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_mapq_high
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->mapq_high = new_mapq_high;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for mapq_sum member in a bl_sam_buff_t structure.
 *      Use this function to set mapq_sum in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      mapq_sum is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_mapq_sum    The new value for mapq_sum
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_mapq_sum;
 *
 *      if ( bl_sam_buff_set_mapq_sum(&bl_sam_buff, new_mapq_sum)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_mapq_sum(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_mapq_sum
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->mapq_sum = new_mapq_sum;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for reads_used member in a bl_sam_buff_t structure.
 *      Use this function to set reads_used in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      reads_used is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_reads_used  The new value for reads_used
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_reads_used;
 *
 *      if ( bl_sam_buff_set_reads_used(&bl_sam_buff, new_reads_used)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_reads_used(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_reads_used
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->reads_used = new_reads_used;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for total_alignments member in a bl_sam_buff_t structure.
 *      Use this function to set total_alignments in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      total_alignments is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_total_alignments The new value for total_alignments
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_total_alignments;
 *
 *      if ( bl_sam_buff_set_total_alignments(&bl_sam_buff, new_total_alignments)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_total_alignments(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_total_alignments
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->total_alignments = new_total_alignments;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for trailing_alignments member in a bl_sam_buff_t structure.
 *      Use this function to set trailing_alignments in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      trailing_alignments is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_trailing_alignments The new value for trailing_alignments
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_trailing_alignments;
 *
 *      if ( bl_sam_buff_set_trailing_alignments(&bl_sam_buff, new_trailing_alignments)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_trailing_alignments(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_trailing_alignments
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->trailing_alignments = new_trailing_alignments;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for discarded_alignments member in a bl_sam_buff_t structure.
 *      Use this function to set discarded_alignments in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      discarded_alignments is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_discarded_alignments The new value for discarded_alignments
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_discarded_alignments;
 *
 *      if ( bl_sam_buff_set_discarded_alignments(&bl_sam_buff, new_discarded_alignments)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_discarded_alignments(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_discarded_alignments
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->discarded_alignments = new_discarded_alignments;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for discarded_score_sum member in a bl_sam_buff_t structure.
 *      Use this function to set discarded_score_sum in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      discarded_score_sum is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_discarded_score_sum The new value for discarded_score_sum
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_discarded_score_sum;
 *
 *      if ( bl_sam_buff_set_discarded_score_sum(&bl_sam_buff, new_discarded_score_sum)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_discarded_score_sum(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_discarded_score_sum
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->discarded_score_sum = new_discarded_score_sum;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for discarded_trailing member in a bl_sam_buff_t structure.
 *      Use this function to set discarded_trailing in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      discarded_trailing is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_discarded_trailing The new value for discarded_trailing
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_discarded_trailing;
 *
 *      if ( bl_sam_buff_set_discarded_trailing(&bl_sam_buff, new_discarded_trailing)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_discarded_trailing(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_discarded_trailing
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->discarded_trailing = new_discarded_trailing;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for min_discarded_score member in a bl_sam_buff_t structure.
 *      Use this function to set min_discarded_score in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      min_discarded_score is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_min_discarded_score The new value for min_discarded_score
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_min_discarded_score;
 *
 *      if ( bl_sam_buff_set_min_discarded_score(&bl_sam_buff, new_min_discarded_score)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_min_discarded_score(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_min_discarded_score
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->min_discarded_score = new_min_discarded_score;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for max_discarded_score member in a bl_sam_buff_t structure.
 *      Use this function to set max_discarded_score in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      max_discarded_score is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_max_discarded_score The new value for max_discarded_score
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_max_discarded_score;
 *
 *      if ( bl_sam_buff_set_max_discarded_score(&bl_sam_buff, new_max_discarded_score)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_max_discarded_score(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_max_discarded_score
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->max_discarded_score = new_max_discarded_score;
	return BL_SAM_BUFF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam-buff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for unmapped_alignments member in a bl_sam_buff_t structure.
 *      Use this function to set unmapped_alignments in a bl_sam_buff_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      unmapped_alignments is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_buff_ptr Pointer to the structure to set
 *      new_unmapped_alignments The new value for unmapped_alignments
 *
 *  Returns:
 *      BL_SAM_BUFF_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_BUFF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_buff_t   bl_sam_buff;
 *      int64_t        new_unmapped_alignments;
 *
 *      if ( bl_sam_buff_set_unmapped_alignments(&bl_sam_buff, new_unmapped_alignments)
 *              == BL_SAM_BUFF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-07  gen-get-set Auto-generated from sam-buff.h
 ***************************************************************************/

int     bl_sam_buff_set_unmapped_alignments(
	    bl_sam_buff_t *bl_sam_buff_ptr,
	    int64_t new_unmapped_alignments
	)

{
    if ( false )
	return BL_SAM_BUFF_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_buff_ptr->unmapped_alignments = new_unmapped_alignments;
	return BL_SAM_BUFF_DATA_OK;
    }
}
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sysexits.h>





/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over header lines in SAM input stream.  The FILE pointer
 *      sam_stream is advanced to the first character of the first line
 *      after the header.  The header is copied to a temporary file and and
 *      the function returns a FILE pointer to the header stream.
 *
 *  Arguments:
 *      sam_stream  FILE pointer to the open sam file
 *
 *  Returns:
 *      A FILE pointer to a temporary file containing a copy of the header
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_copy_header(3)
 *  
 *  Arguments:
 *
 *  Returns:
 *
 *  Examples:
 *
 *  Files:
 *
 *  Environment
 *
 *  See also:
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-08  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_sam_skip_header(FILE *sam_stream)

{
    int     ch;
    FILE    *header_stream = tmpfile();

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like peak-classifier to replicate the
     *  header in output files.
     */
    
    while ( (ch = getc(sam_stream)) == '@' )
    {
	putc(ch, header_stream);
	do
	{
	    ch = getc(sam_stream);
	    putc(ch, header_stream);
	}   while ( (ch != '\n') && (ch != EOF) );
    }
    
    // Rewind to start of first non-header line
    if ( ch != EOF )
	ungetc(ch, sam_stream);
    rewind(header_stream);
    return header_stream;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy SAM header from one FILE stream to another.  This is meant to
 *      be used in conjunction with bl_sam_skip_header(), which stores the
 *      header in a temporary file.
 *
 *  Arguments:
 *      header_stream   Open FILE stream of SAM header
 *      sam_stream      FILE stream to which header is copied
 *
 *  Returns:
 *      BL_WRITE_OK upon success, BL_WRITE_FAILURE or BL_READ_* on failure
 *
 *  See also:
 *      bl_sam_skip_header(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-02-03  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_copy_header(FILE *header_stream, FILE *sam_stream)

{
    int     ch;
    
    rewind(header_stream);
    while ( (ch = getc(header_stream)) != EOF )
	if ( putc(ch, sam_stream) == EOF )
	    return BL_WRITE_FAILURE;
    rewind(header_stream);
    return BL_WRITE_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read next alignment (line) from a SAM stream.
 *
 *      If field_mask is not BL_SAM_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in alignment.
 *      That field in the structure is then populated with an appropriate
 *      marker, such as '.'.  Possible mask values are:
 *
 *      BL_SAM_FIELD_ALL
 *      BL_SAM_FIELD_QNAME
 *      BL_SAM_FIELD_FLAG
 *      BL_SAM_FIELD_RNAME
 *      BL_SAM_FIELD_POS
 *      BL_SAM_FIELD_MAPQ
 *      BL_SAM_FIELD_CIGAR
 *      BL_SAM_FIELD_RNEXT
 *      BL_SAM_FIELD_PNEXT
 *      BL_SAM_FIELD_TLEN
 *      BL_SAM_FIELD_SEQ
 *      BL_SAM_FIELD_QUAL
 *
 *  Arguments:
 *      sam_stream  A FILE stream from which to read the line
 *      alignment   Pointer to a bl_sam_t structure
 *      field_mask  Bit mask indicating which fields to store in alignment
 *
 *  Returns:
 *      BL_READ_OK on successful read
 *      BL_READ_EOF if EOF is encountered after a complete feature
 *      BL_READ_TRUNCATED if EOF or bad data is encountered elsewhere
 *
 *  Examples:
 *      bl_sam_read(stdin, &alignment, BL_SAM_FIELD_ALL);
 *      bl_sam_read(sam_stream, &alignment,
 *                         BL_SAM_FIELD_QNAME|BL_SAM_FIELD_POS|BL_SAM_FIELD_TLEN);
 *
 *  See also:
 *      bl_sam_write(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_read(bl_sam_t *alignment, FILE *sam_stream,
			   sam_field_mask_t field_mask)

{
    char    mapq_str[BL_SAM_MAPQ_MAX_CHARS + 1],
	    pos_str[BL_POSITION_MAX_DIGITS + 1],
	    flag_str[BL_SAM_FLAG_MAX_DIGITS + 1],
	    *end;
    size_t  len;
    static int64_t   previous_pos = 0;
    int     delim;
    
    if ( field_mask & BL_SAM_FIELD_QNAME )
	delim = tsv_read_field(sam_stream, alignment->qname,
			BL_SAM_QNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	*alignment->qname = '\0';
    }
    if ( delim == EOF )
	return BL_READ_EOF;

    // 2 Flag
    if ( field_mask & BL_SAM_FIELD_FLAG )
	delim = tsv_read_field(sam_stream, flag_str, BL_SAM_FLAG_MAX_DIGITS, &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading flag: %s.\n",
		flag_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_FLAG )
    {
	alignment->flag = strtoul(flag_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid position: %s\n",
		    flag_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->flag = 0;    // FIXME: Is there a better choice?
    
    // 3 RNAME
    if ( field_mask & BL_SAM_FIELD_RNAME )
	delim = tsv_read_field(sam_stream, alignment->rname,
			       BL_SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	*alignment->rname = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading rname: %s.\n",
		alignment->rname);
	return BL_READ_TRUNCATED;
    }
    
    // 4 POS
    if ( field_mask & BL_SAM_FIELD_POS )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading pos: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_POS )
    {
	alignment->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid position: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
	previous_pos = alignment->pos;
    }
    else
	alignment->pos = 0;
    
    // 5 MAPQ
    if ( field_mask & BL_SAM_FIELD_MAPQ )
	delim = tsv_read_field(sam_stream, mapq_str, BL_SAM_MAPQ_MAX_CHARS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading mapq: %s.\n",
		mapq_str);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_MAPQ )
    {
	alignment->mapq = strtoul(mapq_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid mapq: %s\n",
		    mapq_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->mapq = 0;
    
    // 6 CIGAR
    if ( field_mask & BL_SAM_FIELD_CIGAR )
	delim = tsv_read_field_malloc(sam_stream, &alignment->cigar,
			       &alignment->cigar_array_size,
			       &alignment->cigar_len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	alignment->cigar_len = 0;
	// Do not set to NULL or set array_size to 0.  Leave buffer
	// allocated for reuse.
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading cigar: %s.\n",
		alignment->cigar);
	return BL_READ_TRUNCATED;
    }
    
    // 7 RNEXT
    if ( field_mask & BL_SAM_FIELD_RNEXT )
	delim = tsv_read_field(sam_stream, alignment->rnext,
			       BL_SAM_RNAME_MAX_CHARS, &len);
    else
    {
	delim = tsv_skip_field(sam_stream, &len);
	*alignment->rnext = '\0';
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading rnext: %s.\n",
		alignment->rnext);
	return BL_READ_TRUNCATED;
    }
    
    // 8 PNEXT
    if ( field_mask & BL_SAM_FIELD_PNEXT )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading pnext: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_PNEXT )
    {
	alignment->pnext = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid pnext: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->pnext = 0;
    
    // 9 TLEN
    if ( field_mask & BL_SAM_FIELD_TLEN )
	delim = tsv_read_field(sam_stream, pos_str, BL_POSITION_MAX_DIGITS,
			       &len);
    else
	delim = tsv_skip_field(sam_stream, &len);
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading tlen: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    if ( field_mask & BL_SAM_FIELD_TLEN )
    {
	alignment->tlen = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr, "bl_sam_read(): Invalid tlen: %s\n",
		    pos_str);
	    fprintf(stderr, "qname = %s rname = %s\n",
		    alignment->qname, alignment->rname);
	    fprintf(stderr, "previous_pos = %" PRId64 "\n", previous_pos);
	    exit(EX_DATAERR);
	}
    }
    else
	alignment->tlen = 0;
    
    // 10 SEQ
    if ( field_mask & BL_SAM_FIELD_SEQ )
	delim = tsv_read_field_malloc(sam_stream, &alignment->seq,
		    &alignment->seq_array_size, &alignment->seq_len);
    else
    {
	delim = tsv_skip_field(sam_stream, &alignment->seq_len);
	alignment->seq_len = 0;
	// Do not set to NULL or set array_size to 0.  Leave buffer
	// allocated for reuse.
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading seq: %s.\n",
		alignment->seq);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_SEQ )
    {
	// May be allocated by bl_sam_init() or bl_sam_copy()
	if ( alignment->seq == NULL )
	{
	    if ( (alignment->seq = xt_malloc(alignment->seq_len + 1,
		    sizeof(*alignment->seq))) == NULL )
	    {
		fprintf(stderr, "bl_sam_read(): Could not allocate seq.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
    }
    
    // 11 QUAL, should be last field
    if ( field_mask & BL_SAM_FIELD_QUAL )
    {
	delim = tsv_read_field_malloc(sam_stream, &alignment->qual,
		    &alignment->qual_array_size,
		    &alignment->qual_len);
    }
    else
    {
	delim = tsv_skip_field(sam_stream, &alignment->qual_len);
	alignment->qual_len = 0;
	// Do not set to NULL or set array_size to 0.  Leave buffer
	// allocated for reuse.
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_sam_read(): Got EOF reading qual: %s.\n",
		alignment->qual);
	return BL_READ_TRUNCATED;
    }

    if ( field_mask & BL_SAM_FIELD_QUAL )
    {
	// May be allocated by bl_sam_init() or bl_sam_copy()
	if ( alignment->qual == NULL )
	{
	    if ( (alignment->qual = xt_malloc(alignment->qual_len + 1,
		    sizeof(*alignment->qual))) == NULL )
	    {
		fprintf(stderr, "bl_sam_read(): Could not allocate qual.\n");
		exit(EX_UNAVAILABLE);
	    }
	}
    
	if ( (alignment->qual_len != 1) &&
	     (alignment->seq_len != alignment->qual_len) )
	    fprintf(stderr, "bl_sam_read(): Warning: qual_len != seq_len for %s,%" PRId64 "\n",
		    alignment->rname, alignment->pos);
    }
    
    // Some SRA CRAMs have 11 fields, most have 12
    // Discard everything after the 11th
    if ( delim == '\t' )
	while ( getc(sam_stream) != '\n' )
	    ;

    /*fprintf(stderr,"bl_sam_read(): %s,%" PRId64 ",%zu\n",
	    BL_SAM_RNAME(alignment), BL_SAM_POS(alignment),
	    BL_SAM_SEQ_LEN(alignment));*/
    
    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Copy a SAM alignment as efficiently as possible, allocating memory
 *      as needed.
 *
 *  Arguments:
 *      dest    Pointer to bl_sam_t structure to receive copy
 *      src     Pointer to bl_sam_t structure to be copied
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_init(3), bl_sam_free(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_copy(bl_sam_t *dest, bl_sam_t *src)

{
    strlcpy(dest->qname, src->qname, BL_SAM_QNAME_MAX_CHARS + 1);
    dest->flag = src->flag;
    strlcpy(dest->rname, src->rname, BL_SAM_RNAME_MAX_CHARS + 1);
    dest->pos = src->pos;
    dest->mapq = src->mapq;

    if ( src->cigar != NULL )
    {
	dest->cigar = strdup(src->cigar);
	if ( dest->cigar == NULL )
	{
	    fprintf(stderr, "bl_sam_copy(): Could not allocate cigar.\n");
	    exit(EX_UNAVAILABLE);
	}
	dest->cigar_array_size = src->cigar_len + 1;
	dest->cigar_len = src->cigar_len;
    }
    dest->cigar_array_size = src->cigar_array_size;
    dest->cigar_len = src->cigar_len;
    
    strlcpy(dest->rnext, src->rnext, BL_SAM_RNAME_MAX_CHARS + 1);
    dest->pnext = src->pnext;
    dest->tlen = src->tlen;

    if ( src->seq != NULL )
    {
	if ( (dest->seq = strdup(src->seq)) == NULL )
	{
	    fprintf(stderr, "bl_sam_copy(): Could not allocate seq.\n");
	    exit(EX_UNAVAILABLE);
	}
	dest->seq_array_size = src->seq_len + 1;
	dest->seq_len = src->seq_len;
    }
    dest->seq_array_size = src->seq_array_size;
    dest->seq_len = src->seq_len;
    
    //fprintf(stderr, "src->seq = %s %zu\n", src->seq, src->seq_len);
    //fprintf(stderr, "src->qual = %s %zu\n", src->qual, src->qual_len);
    
    /*
     *  qual is an optional field.  If skipped using tsv_skip_field()
     *  qual_len will be non-zero, but qual will be NULL.
     */
    if ( src->qual != NULL )
    {
	if ( (dest->qual = strdup(src->qual)) == NULL )
	{
	    fprintf(stderr, "bl_sam_copy(): Could not allocate qual.\n");
	    exit(EX_UNAVAILABLE);
	}
	dest->qual_array_size = src->qual_len + 1;
	dest->qual_len = src->qual_len;
    }
    dest->qual_array_size = src->qual_array_size;
    dest->qual_len = src->qual_len;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free memory allocated by bl_sam_read() or
 *      bl_sam_init().
 *
 *  Arguments:
 *      alignment   Pointer to bl_sam_t structure to be freed.
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_init(3), bl_sam_copy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_free(bl_sam_t *alignment)

{
    if ( alignment->cigar != NULL )
	free(alignment->cigar);
    if ( alignment->seq != NULL )
	free(alignment->seq);
    if ( alignment->qual != NULL )
	free(alignment->qual);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_sam_t structure, allocating memory for
 *      sequence and quality strings according to seq_len.  Passing a
 *      seq_len of 0 prevents memory allocation from occurring.
 *
 *      Only BL_SAM_FIELD_SEQ and BL_SAM_FIELD_QUAL are meaningful bits in
 *      field_mask, as they determine whether memory is allocated.  All
 *      other fields are unconditionally initialized to 0, NULL, or blank.
 *
 *  Arguments:
 *      alignment   Pointer to bl_sam_t structure to initialize
 *      seq_len     Length of sequence and quality strings
 *      field_mask  Bit mask indicating which fields will be used
 *
 *  See also:
 *      bl_sam_read(3), bl_sam_free(3), bl_sam_copy(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-29  Jason Bacon Begin
 ***************************************************************************/

void    bl_sam_init(bl_sam_t *alignment)

{
    *alignment->qname = '\0';
    alignment->flag = 0;
    *alignment->rname = '\0';
    alignment->pos = 0;
    alignment->mapq = 0;
    alignment->cigar = NULL;
    *alignment->rnext = '\0';
    alignment->pnext = 0;
    alignment->tlen = 0;
    alignment->seq = NULL;
    alignment->qual = NULL;
    alignment->seq_array_size = 0;
    alignment->seq_len = 0;
    alignment->qual_array_size = 0;
    alignment->qual_len = 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write an alignment (line) to a SAM stream.
 *
 *      If field_mask is not BL_SAM_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate placeholder such as '.'
 *      rather than stored in alignment.  Possible mask values are:
 *
 *      BL_SAM_FIELD_ALL
 *      BL_SAM_FIELD_QNAME
 *      BL_SAM_FIELD_FLAG
 *      BL_SAM_FIELD_RNAME
 *      BL_SAM_FIELD_POS
 *      BL_SAM_FIELD_MAPQ
 *      BL_SAM_FIELD_CIGAR
 *      BL_SAM_FIELD_RNEXT
 *      BL_SAM_FIELD_PNEXT
 *      BL_SAM_FIELD_TLEN
 *      BL_SAM_FIELD_SEQ
 *      BL_SAM_FIELD_QUAL
 *
 *  Arguments:
 *      sam_stream  A FILE stream to which to write the line
 *      alignment   Pointer to a bl_sam_t structure
 *      field_mask  Bit mask indicating which fields to store in alignment
 *
 *  Returns:
 *      Number of items written (per fprintf() output)
 *
 *  See also:
 *      bl_sam_read(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-09  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_write(bl_sam_t *alignment, FILE *sam_stream,
			   sam_field_mask_t field_mask)

{
    int     count;
    
    // FIXME: Respect field_mask
    count = fprintf(sam_stream, "%s\t%u\t%s\t%" PRId64
		    "\t%u\t%s\t%s\t%" PRId64 "zu\t%zu\t%s\t%s\t%zu\t%zu\n",
		    alignment->qname,
		    alignment->flag,
		    alignment->rname,
		    alignment->pos,
		    alignment->mapq,
		    alignment->cigar,
		    alignment->rnext,
		    alignment->pnext,
		    alignment->tlen,
		    alignment->seq,
		    alignment->qual,
		    alignment->seq_len,
		    alignment->qual_len);
    return count;
}


/***************************************************************************
 *  Library:

 *      -lxtend
 *
 *  Description:
 *      Open a raw SAM file using fopen() or a gzipped, bzipped, or
 *      xzipped SAM file or BAM or CRAM file using popen().  If the
 *      file extension is .bam or .cram, or samtools_args is not
 *      NULL or "", data will be piped through "samtools view" with
 *      the given samtools_args as arguments.  The flag --with-header
 *      is always added for consistency with the case of reading a
 *      raw SAM file without piping through samtools.  Programs that
 *      don't want the header can filter it out by other means, such
 *      as bl_sam_skip_header(3).
 *
 *      bl_sam_fopen() must be used in conjunction with
 *      bl_sam_fclose() to ensure that fclose() or pclose() is called where
 *      appropriate.
 *
 *  Arguments:
 *      filename:       Name of the file to be opened
 *      mode:           "r" or "w", passed to fopen() or popen()
 *      samtools_args   Flags to pass to samtools view
 *
 *  Returns:
 *      A pointer to the FILE structure or NULL if open failed
 *
 *  See also:
 *      fopen(3), popen(3), gzip(1), bzip2(1), xz(1)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-05  Jason Bacon Derived from xt_fclose()
 ***************************************************************************/

FILE    *bl_sam_fopen(const char *filename, const char *mode,
		      char *samtools_args)

{
    char    *ext = strrchr(filename, '.'),
	    cmd[XT_CMD_MAX_CHARS + 1];
    struct stat sb;
    
    if ( samtools_args == NULL )
	samtools_args = "";
    
    if ( (strcmp(mode, "r") != 0 ) && (strcmp(mode, "w") != 0) )
    {
	fprintf(stderr, "bl_sam_fopen(): Only \"r\" and \"w\" modes supported.\n");
	return NULL;
    }
    
    if ( ext == NULL )
    {
	fprintf(stderr, "bl_sam_fopen(): No filename extension on %s.\n", filename);
	return NULL;
    }

    // popen() does not return NULL when the file does not exist
    if ( stat(filename, &sb) != 0 )
	return NULL;
    
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
	else if ( (strcmp(ext, ".bam") == 0) || (strcmp(ext, ".cram") == 0) 
		  || ! strblank(samtools_args) )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --with-header %s %s",
		    samtools_args, filename);
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
	else if ( strcmp(ext, ".bam") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --bam --with-header %s %s",
		    samtools_args, filename);
	    return popen(cmd, mode);
	}
	else if ( strcmp(ext, ".cram") == 0 )
	{
	    snprintf(cmd, XT_CMD_MAX_CHARS, "samtools view --cram --with-header %s %s",
		    samtools_args, filename);
	    return popen(cmd, mode);
	}
	else
	    return fopen(filename, mode);
    }
}


/***************************************************************************
 *  Library:

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
 *  2022-04-05  Jason Bacon Derived from xt_fclose()
 ***************************************************************************/

int     bl_sam_fclose(FILE *stream)

{
    return xt_fclose(stream);
}



/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <biolibc/gff.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Return the amount of overlap between a GFF feature and a SAM
 *      alignment.
 *  
 *  Arguments:
 *      alignment   Pointer to a bl_sam_t object
 *      feature     Pointer to a bl_gff_t object
 *
 *  Returns:
 *      The number of bases of overlap between the feature and alignment.
 *      A zero or negative return value indicates no overlap.
 *
 *  See also:
 *      bl_gff_sam_overlap(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-07  Jason Bacon Begin
 ***************************************************************************/

int64_t bl_sam_gff_overlap(bl_sam_t *alignment, bl_gff_t *feature)

{
    return bl_gff_sam_overlap(feature, alignment);
}


/***************************************************************************
 *  Use auto-c2man to generate a man page from this comment
 *
 *  Library:
 *      #include <sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Compare the positions of a SAM alignment and a GFF feature and
 *      return a status value much like strcmp().  0 is returned if the
 *      alignment and feature overlap.  A value < 0 is returned if the
 *      alignment is entirely "before" the feature, i.e. it is on an
 *      earlier chromosome according to bl_chrom_name_cmp(3), or on the
 *      same chromosome at a lower position.  A value > 0 is returned
 *      if the alignment is entirely "after" the feature, i.e. on a later
 *      chromosome or same chromosome and higher position.
 *
 *      This function is mainly intended for programs that sweep properly
 *      sorted GFF and SAM files locating overlaps in a single pass.
 *  
 *      A converse function, bl_gff_sam_cmp(3) is also provided so that
 *      the programmer can choose the more intuitive interface.
 *  
 *  Arguments:
 *      alignment   Pointer to a bl_sam_t object
 *      feature     Pointer to a bl_gff_t object
 *
 *  Returns:
 *      A value < 0 if the the alignment is entirely before the feature
 *      A value > 0 if the the alignment is entirely after the feature
 *      0 if the alignment and the feature overlap
 *
 *  See also:
 *      bl_gff_sam_cmp(3), bl_chrom_name_cmp(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-06  Jason Bacon Begin
 ***************************************************************************/

int     bl_sam_gff_cmp(bl_sam_t *alignment, bl_gff_t *feature)

{
    int     status = bl_chrom_name_cmp(BL_SAM_RNAME(alignment),
					    BL_GFF_SEQID(feature));
    
    if ( status != 0 )
	// Different chromosomes
	return status;
    else if ( BL_SAM_POS(alignment) + BL_SAM_SEQ_LEN(alignment) - 1
		< BL_GFF_START(feature) )
	// Alignment ends before the start of feature
	return -1;
    else if ( BL_SAM_POS(alignment) > BL_GFF_END(feature) )
	// Alignment starts after the end of feature
	return 1;
    else
	// Overlap
	return 0;
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
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of qname member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->qname[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the qname array
 *      new_qname_element The new value for qname[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char            new_qname_element;
 *
 *      if ( bl_sam_set_qname_ae(&bl_sam, c, new_qname_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QNAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qname_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_qname_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qname[c] = new_qname_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qname member in a bl_sam_t structure.
 *      Use this function to set qname in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_qname to bl_sam_ptr->qname.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qname       The new value for qname
 *      array_size      Size of the qname array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char            new_qname;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_qname_cpy(&bl_sam, new_qname, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QNAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qname_cpy(bl_sam_t *bl_sam_ptr, char new_qname[], size_t array_size)

{
    if ( new_qname == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->qname, new_qname, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for flag member in a bl_sam_t structure.
 *      Use this function to set flag in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      flag is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_flag        The new value for flag
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      unsigned        new_flag;
 *
 *      if ( bl_sam_set_flag(&bl_sam, new_flag)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_flag(bl_sam_t *bl_sam_ptr, unsigned new_flag)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->flag = new_flag;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of rname member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->rname[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the rname array
 *      new_rname_element The new value for rname[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char            new_rname_element;
 *
 *      if ( bl_sam_set_rname_ae(&bl_sam, c, new_rname_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNAME_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rname_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_rname_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->rname[c] = new_rname_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for rname member in a bl_sam_t structure.
 *      Use this function to set rname in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_rname to bl_sam_ptr->rname.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_rname       The new value for rname
 *      array_size      Size of the rname array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char            new_rname;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_rname_cpy(&bl_sam, new_rname, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNAME(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rname_cpy(bl_sam_t *bl_sam_ptr, char new_rname[], size_t array_size)

{
    if ( new_rname == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->rname, new_rname, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for pos member in a bl_sam_t structure.
 *      Use this function to set pos in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      pos is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_pos         The new value for pos
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      int64_t         new_pos;
 *
 *      if ( bl_sam_set_pos(&bl_sam, new_pos)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_pos(bl_sam_t *bl_sam_ptr, int64_t new_pos)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->pos = new_pos;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for mapq member in a bl_sam_t structure.
 *      Use this function to set mapq in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      mapq is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_mapq        The new value for mapq
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      unsigned char   new_mapq;
 *
 *      if ( bl_sam_set_mapq(&bl_sam, new_mapq)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_mapq(bl_sam_t *bl_sam_ptr, unsigned char new_mapq)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->mapq = new_mapq;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar member in a bl_sam_t structure.
 *      Use this function to set cigar in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      cigar is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar       The new value for cigar
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_cigar;
 *
 *      if ( bl_sam_set_cigar(&bl_sam, new_cigar)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar(bl_sam_t *bl_sam_ptr, char * new_cigar)

{
    if ( new_cigar == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar = new_cigar;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of cigar member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->cigar[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the cigar array
 *      new_cigar_element The new value for cigar[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char *          new_cigar_element;
 *
 *      if ( bl_sam_set_cigar_ae(&bl_sam, c, new_cigar_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_CIGAR_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_ae(bl_sam_t *bl_sam_ptr, size_t c, char  new_cigar_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar[c] = new_cigar_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar member in a bl_sam_t structure.
 *      Use this function to set cigar in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_cigar to bl_sam_ptr->cigar.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar       The new value for cigar
 *      array_size      Size of the cigar array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_cigar;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_cigar_cpy(&bl_sam, new_cigar, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_CIGAR(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_cpy(bl_sam_t *bl_sam_ptr, char * new_cigar, size_t array_size)

{
    if ( new_cigar == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->cigar, new_cigar, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of rnext member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->rnext[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the rnext array
 *      new_rnext_element The new value for rnext[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char            new_rnext_element;
 *
 *      if ( bl_sam_set_rnext_ae(&bl_sam, c, new_rnext_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNEXT_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rnext_ae(bl_sam_t *bl_sam_ptr, size_t c, char new_rnext_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->rnext[c] = new_rnext_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for rnext member in a bl_sam_t structure.
 *      Use this function to set rnext in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_rnext to bl_sam_ptr->rnext.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_rnext       The new value for rnext
 *      array_size      Size of the rnext array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char            new_rnext;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_rnext_cpy(&bl_sam, new_rnext, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_RNEXT(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_rnext_cpy(bl_sam_t *bl_sam_ptr, char new_rnext[], size_t array_size)

{
    if ( new_rnext == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->rnext, new_rnext, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for pnext member in a bl_sam_t structure.
 *      Use this function to set pnext in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      pnext is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_pnext       The new value for pnext
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      int64_t         new_pnext;
 *
 *      if ( bl_sam_set_pnext(&bl_sam, new_pnext)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_pnext(bl_sam_t *bl_sam_ptr, int64_t new_pnext)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->pnext = new_pnext;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for tlen member in a bl_sam_t structure.
 *      Use this function to set tlen in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      tlen is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_tlen        The new value for tlen
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      long            new_tlen;
 *
 *      if ( bl_sam_set_tlen(&bl_sam, new_tlen)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_tlen(bl_sam_t *bl_sam_ptr, long new_tlen)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->tlen = new_tlen;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_sam_t structure.
 *      Use this function to set seq in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq         The new value for seq
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_seq;
 *
 *      if ( bl_sam_set_seq(&bl_sam, new_seq)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq(bl_sam_t *bl_sam_ptr, char * new_seq)

{
    if ( new_seq == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq = new_seq;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of seq member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->seq[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the seq array
 *      new_seq_element The new value for seq[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char *          new_seq_element;
 *
 *      if ( bl_sam_set_seq_ae(&bl_sam, c, new_seq_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_SEQ_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_ae(bl_sam_t *bl_sam_ptr, size_t c, char  new_seq_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq[c] = new_seq_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq member in a bl_sam_t structure.
 *      Use this function to set seq in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_seq to bl_sam_ptr->seq.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq         The new value for seq
 *      array_size      Size of the seq array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_seq;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_seq_cpy(&bl_sam, new_seq, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_SEQ(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_cpy(bl_sam_t *bl_sam_ptr, char * new_seq, size_t array_size)

{
    if ( new_seq == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->seq, new_seq, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual member in a bl_sam_t structure.
 *      Use this function to set qual in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual        The new value for qual
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_qual;
 *
 *      if ( bl_sam_set_qual(&bl_sam, new_qual)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual(bl_sam_t *bl_sam_ptr, char * new_qual)

{
    if ( new_qual == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual = new_qual;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of qual member in a bl_sam_t
 *      structure. Use this function to set bl_sam_ptr->qual[c]
 *      in a bl_sam_t object from non-member functions.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      c               Subscript to the qual array
 *      new_qual_element The new value for qual[c]
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          c;
 *      char *          new_qual_element;
 *
 *      if ( bl_sam_set_qual_ae(&bl_sam, c, new_qual_element)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QUAL_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_ae(bl_sam_t *bl_sam_ptr, size_t c, char  new_qual_element)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual[c] = new_qual_element;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual member in a bl_sam_t structure.
 *      Use this function to set qual in a bl_sam_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_qual to bl_sam_ptr->qual.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual        The new value for qual
 *      array_size      Size of the qual array.
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      char *          new_qual;
 *      size_t          array_size;
 *
 *      if ( bl_sam_set_qual_cpy(&bl_sam, new_qual, array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_SAM_SET_QUAL(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_cpy(bl_sam_t *bl_sam_ptr, char * new_qual, size_t array_size)

{
    if ( new_qual == NULL )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_sam_ptr->qual, new_qual, array_size);
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar_array_size member in a bl_sam_t structure.
 *      Use this function to set cigar_array_size in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      cigar_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar_array_size The new value for cigar_array_size
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_cigar_array_size;
 *
 *      if ( bl_sam_set_cigar_array_size(&bl_sam, new_cigar_array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_array_size(bl_sam_t *bl_sam_ptr, size_t new_cigar_array_size)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar_array_size = new_cigar_array_size;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for cigar_len member in a bl_sam_t structure.
 *      Use this function to set cigar_len in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      cigar_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_cigar_len   The new value for cigar_len
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_cigar_len;
 *
 *      if ( bl_sam_set_cigar_len(&bl_sam, new_cigar_len)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_cigar_len(bl_sam_t *bl_sam_ptr, size_t new_cigar_len)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->cigar_len = new_cigar_len;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_array_size member in a bl_sam_t structure.
 *      Use this function to set seq_array_size in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq_array_size The new value for seq_array_size
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_seq_array_size;
 *
 *      if ( bl_sam_set_seq_array_size(&bl_sam, new_seq_array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_array_size(bl_sam_t *bl_sam_ptr, size_t new_seq_array_size)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq_array_size = new_seq_array_size;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for seq_len member in a bl_sam_t structure.
 *      Use this function to set seq_len in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      seq_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_seq_len     The new value for seq_len
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_seq_len;
 *
 *      if ( bl_sam_set_seq_len(&bl_sam, new_seq_len)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_seq_len(bl_sam_t *bl_sam_ptr, size_t new_seq_len)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->seq_len = new_seq_len;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual_array_size member in a bl_sam_t structure.
 *      Use this function to set qual_array_size in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual_array_size The new value for qual_array_size
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_qual_array_size;
 *
 *      if ( bl_sam_set_qual_array_size(&bl_sam, new_qual_array_size)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_array_size(bl_sam_t *bl_sam_ptr, size_t new_qual_array_size)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual_array_size = new_qual_array_size;
	return BL_SAM_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/sam.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual_len member in a bl_sam_t structure.
 *      Use this function to set qual_len in a bl_sam_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      qual_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_sam_ptr      Pointer to the structure to set
 *      new_qual_len    The new value for qual_len
 *
 *  Returns:
 *      BL_SAM_DATA_OK if the new value is acceptable and assigned
 *      BL_SAM_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_sam_t        bl_sam;
 *      size_t          new_qual_len;
 *
 *      if ( bl_sam_set_qual_len(&bl_sam, new_qual_len)
 *              == BL_SAM_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-30  gen-get-set Auto-generated from sam.h
 ***************************************************************************/

int     bl_sam_set_qual_len(bl_sam_t *bl_sam_ptr, size_t new_qual_len)

{
    if ( false )
	return BL_SAM_DATA_OUT_OF_RANGE;
    else
    {
	bl_sam_ptr->qual_len = new_qual_len;
	return BL_SAM_DATA_OK;
    }
}
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <sysexits.h>
#include <stdbool.h>




/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over meta-data lines in VCF input stream, leaving the FILE
 *      structure pointing to the first character in the first line of data
 *      or the first character of the header line starting with #CHROM if
 *      one is present.  The header line is typically read using
 *      bl_vcf_get_sample_ids(3). The skipped meta-data is copied to a
 *      temporary file whose FILE pointer is returned.
 *
 *  Arguments:
 *      vcf_stream  FILE pointer of VCF stream to be read
 *
 *  Returns:
 *      BL_READ_OK upon success, BL_READ_TRUNCATED if read fails
 *
 *  See also:
 *      bl_vcf_get_sample_ids(3), bl_vcf_read_static_fields(3), bl_vcf_read_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_vcf_skip_meta_data(FILE *vcf_stream)

{
    int     ch,
	    c,
	    count;
    char    start[6];
    FILE    *meta_stream;

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like vcf-split to replicate the
     *  header in output files.
     */

    meta_stream = tmpfile();
    
    while ( (ch = getc(vcf_stream)) == '#' )
    {
	count = fread(start, 1, 5, vcf_stream);
	
	// Put back "CHROM" or whatever was read regardless
	for (c = count - 1; c >= 0; --c)
	    ungetc(start[c], vcf_stream);
	// Something is seriously wrong if we don't find at least 5 chars
	if ( count != 5 )
	{
	    fclose(meta_stream);
	    return NULL;
	}
	
	if ( memcmp(start, "CHROM", 5) == 0 )
	{
	    // After return, read should start with #CHROM
	    ungetc(ch, vcf_stream);
	    rewind(meta_stream);
	    return meta_stream;
	}
	else
	{
	    // No #CHROM, transfer entire line to temp file
	    putc('#', meta_stream);
	    do
	    {
		ch = getc(vcf_stream);
		putc(ch, meta_stream);
	    }   while ( (ch != '\n') && (ch != EOF) );
	    if ( ch == EOF )
	    {
		fprintf(stderr,
		    "bl_vcf_skip_meta_data(): EOF reached reading meta-data.\n");
		fclose(meta_stream);
		return NULL;
	    }
	}
    }
    
    fprintf(stderr, "bl_vcf_skip_meta_data(): Warning: No #CHROM found in header.\n");
    rewind(meta_stream);
    return meta_stream;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Skip over meta-data lines and header line (beginning with #CHROM)
 *      in a VCF input stream, leaving the FILE structure pointing to the
 *      first character in the first line of data.
 *      The skipped meta-data and header are copied to a temporary
 *      file whose FILE pointer is returned.
 *
 *      Note that the header line (beginning with #CHROM and containing
 *      sample IDs) is typically read using bl_vcf_get_sample_ids(3).
 *      If you wish to do this, call bl_vcf_skip_meta_data() instead of
 *      bl_vcf_skip_header().
 *
 *  Arguments:
 *      vcf_stream  FILE pointer of VCF stream to be read
 *
 *  Returns:
 *      BL_READ_OK upon success, BL_READ_TRUNCATED if read fails
 *
 *  See also:
 *      bl_vcf_get_sample_ids(3), bl_vcf_read_static_fields(3), bl_vcf_read_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

FILE    *bl_vcf_skip_header(FILE *vcf_stream)

{
    int     ch;
    FILE    *meta_stream;

    /*
     *  Copy header to a nameless temp file and return the FILE *.
     *  This can be used by tools like vcf-split to replicate the
     *  header in output files.
     */

    meta_stream = bl_vcf_skip_meta_data(vcf_stream);
    if ( meta_stream != NULL )
    {
	if ( getc(vcf_stream) == '#' )  // #CHROM line?
	{
	    fseek(meta_stream, 0L, SEEK_END);   // Append header line
	    putc('#', meta_stream);
	    while ( ((ch = getc(vcf_stream)) != '\n') && (ch != EOF) )
		putc(ch, meta_stream);
	    putc(ch, meta_stream);
	    rewind(meta_stream);
	}
	else
	    ungetc('#', vcf_stream);
    }
    return meta_stream;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Extract sample IDs from a VCF input header line.  This is typically
 *      done following bl_vcf_skip_meta_data(3), which will leave the FILE
 *      pointer pointing to the beginning of the header line, if one is
 *      present.
 *
 *      The arguments first_col and last_col represent the first and
 *      last sample columns, both inclusive, from which sample IDs should
 *      be extracted.  A value of 1 represents the first sample column.
 *      This feature allows a VCF file with many columns to be processed
 *      in multiple stages.  For example, the vcf-split tool, based on
 *      biolibc, cannot efficiently process more than abou1 10,000 samples
 *      at once, since each sample requires an open output file.  A VCF
 *      with 150,000 samples can be processed in 15 separate passes.
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to the VCF input stream
 *      sample_ids  Array if character pointers to receive sample IDs
 *      first_col   First column from which a sample ID should be saved
 *      last_col    Last column from which a sample ID should be saved
 *
 *  See also:
 *      bl_vcf_skip_meta_data(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-06  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_get_sample_ids(FILE *vcf_stream, char *sample_ids[],
			   size_t first_col, size_t last_col)

{
    size_t  c,
	    len;
    char    temp_sample_id[BL_VCF_SAMPLE_ID_MAX_CHARS + 1];
    int     delimiter = 0;
    
    // Skip standard header tags to get to sample IDs
    for (c = 0; c < 9; ++c)
	tsv_skip_field(vcf_stream, &len);
    
    // Skip sample IDs before first_col
    for (c = 1; c < first_col; ++c)
	tsv_skip_field(vcf_stream, &len);
    
    for (; (c <= last_col) &&
	   (delimiter = tsv_read_field(vcf_stream, temp_sample_id,
				     BL_VCF_SAMPLE_ID_MAX_CHARS, &len)) != EOF; ++c)
    {
	sample_ids[c - first_col] = strdup(temp_sample_id);
	// fprintf(stderr, "'%s'\n", temp_sample_id);
    }
    
    if ( delimiter == 0 )
    {
	fprintf(stderr, "Reached last_col before reading any sample IDs.\n");
	fprintf(stderr, "Check your first_col and last_col values.\n");
	exit(EX_DATAERR);
    }
    
    // Skip any remaining fields after last_col
    if ( delimiter != '\n' )
	tsv_skip_rest_of_line(vcf_stream);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read static fields (columns 1 to 9) from one line of a VCF file.
 *      This function does not read any of the sample data in columns 10
 *      and on.  Samples can be read using a loop with tsv_read_field(3).
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in bed_feature.
 *      Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE stream for VCF input
 *      vcf_call    Pointer to bl_vcf_t structure to receive fields
 *      field_mask  Bit mask indicating which fields should be stored
 *
 *  Returns:
 *      BL_READ_OK upon success
 *      BL_READ_TRUNCATED if EOF is encountered while reading a call
 *      BL_READ_EOF if EOF is encountered between calls as it should be
 *
 *  Examples:
 *      FILE        *stream;
 *      bl_vcf_t  vcf_call;
 *      char        sample_data[MAX_CHARS + 1];
 *      size_t      len;
 *
 *      bl_vcf_read_static_fields(stream, &vcf_call, BL_VCF_FIELD_ALL);
 *      while ( tsv_read_field(stream, sample_data, MAX_CHARS, &len) != '\n' )
 *      {
 *          ...
 *      }
 *
 *  See also:
 *      bl_vcf_write_static_fields(3), bl_vcf_read_ss_call(3), bl_vcf_write_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-08  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_read_static_fields(bl_vcf_t *vcf_call, FILE *vcf_stream, 
	    vcf_field_mask_t field_mask)

{
    char    *end,
	    pos_str[BL_POSITION_MAX_DIGITS + 1];
    size_t  len;
    int     delim;
    
    vcf_call->ref_count = vcf_call->alt_count = vcf_call->other_count = 0;
    
    // Chromosome
    if ( field_mask & BL_VCF_FIELD_CHROM )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->chrom,
			&vcf_call->chrom_array_size, &vcf_call->chrom_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &len);
	vcf_call->chrom = strdup(".");
	vcf_call->chrom_array_size = 2;
	vcf_call->chrom_len = 1;
    }
    if ( delim == EOF )
    {
	// fputs("bl_vcf_read_static_fields(): Info: Got EOF reading CHROM, as expected.\n", stderr);
	return BL_READ_EOF;
    }
    
    // Call position
    if ( field_mask & BL_VCF_FIELD_POS )
	delim = tsv_read_field(vcf_stream, pos_str,
			BL_POSITION_MAX_DIGITS, &len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &len);
	strlcpy(pos_str, "0", 2);
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading POS: %s.\n",
		pos_str);
	return BL_READ_TRUNCATED;
    }
    else
    {
	vcf_call->pos = strtoul(pos_str, &end, 10);
	if ( *end != '\0' )
	{
	    fprintf(stderr,
		    "bl_vcf_read_static_fields(): Invalid call position: %s\n",
		    pos_str);
	    return BL_READ_TRUNCATED;
	}
    }
    
    // ID
    if ( field_mask & BL_VCF_FIELD_ID )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->id,
			&vcf_call->id_array_size, &vcf_call->id_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &len);
	vcf_call->id = strdup(".");
	vcf_call->id_array_size = 2;
	vcf_call->id_len = 1;
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading ID.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Ref
    if ( field_mask & BL_VCF_FIELD_REF )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->ref,
			&vcf_call->ref_array_size, &vcf_call->ref_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &len);
	vcf_call->ref = strdup(".");
	vcf_call->ref_array_size = 2;
	vcf_call->ref_len = 1;
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading REF.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Alt
    if ( field_mask & BL_VCF_FIELD_ALT )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->alt,
		   &vcf_call->alt_array_size, &vcf_call->alt_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &len);
	vcf_call->alt = strdup(".");
	vcf_call->alt_array_size = 2;
	vcf_call->alt_len = 1;
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading ALT.\n");
	return BL_READ_TRUNCATED;
    }

    // Qual
    if ( field_mask & BL_VCF_FIELD_QUAL )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->qual,
		   &vcf_call->qual_array_size, &vcf_call->qual_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &len);
	vcf_call->qual = strdup(".");
	vcf_call->qual_array_size = 2;
	vcf_call->qual_len = 1;
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading QUAL.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Filter
    if ( field_mask & BL_VCF_FIELD_FILTER )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->filter,
		   &vcf_call->filter_array_size, &vcf_call->filter_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &len);
	vcf_call->filter = strdup(".");
	vcf_call->filter_array_size = 2;
	vcf_call->filter_len = 1;
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading FILTER.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Info
    if ( field_mask & BL_VCF_FIELD_INFO )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->info,
		   &vcf_call->info_array_size, &vcf_call->info_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &vcf_call->info_len);
	vcf_call->info = strdup(".");
	vcf_call->info_array_size = 2;
	vcf_call->info_len = 1;
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading INFO.\n");
	return BL_READ_TRUNCATED;
    }
    
    // Format
    if ( field_mask & BL_VCF_FIELD_FORMAT )
	delim = tsv_read_field_malloc(vcf_stream, &vcf_call->format,
		   &vcf_call->format_array_size, &vcf_call->format_len);
    else
    {
	delim = tsv_skip_field(vcf_stream, &vcf_call->format_len);
	strlcpy(vcf_call->format, ".", 2);
	vcf_call->format_len = 1;
    }
    if ( delim == EOF )
    {
	fprintf(stderr, "bl_vcf_read_static_fields(): Got EOF reading FORMAT.\n");
	return BL_READ_TRUNCATED;
    }

    return BL_READ_OK;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Read a single-sample VCF call (static fields and one sample column).
 *      This should only be used with VCF inputs that have exactly one
 *      sample column.  For multisample VCFs, use bl_vcf_read_static_fields()
 *      followed by a loop to read the sample data.
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are discarded rather than stored in bed_feature.
 *      Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to VCF input stream
 *      vcf_call    bl_vcf_t structure to receive VCF data
 *      field_mask  Bit mask to indicate which fields to store
 *
 *  Returns:
 *      BL_READ_OK upon success
 *      BL_READ_TRUNCATED if EOF is encountered while reading a call
 *      BL_READ_EOF if EOF is encountered between calls as it should be
 *
 *  See also:
 *      bl_vcf_read_static_fields(3), tsv_read_field(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2019-12-11  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_read_ss_call(bl_vcf_t *vcf_call, FILE *vcf_stream,
	    vcf_field_mask_t field_mask)

{
    int     status;
    
    status = bl_vcf_read_static_fields(vcf_call, vcf_stream, field_mask);
    if ( status == BL_READ_OK )
    {
	if ( tsv_read_field_malloc(vcf_stream, &vcf_call->single_sample,
			&vcf_call->single_sample_array_size,
			&vcf_call->single_sample_len) != EOF )
	    return BL_READ_OK;
	else
	{
	    fprintf(stderr, "bl_vcf_read_ss_call(): Got EOF reading sample.\n");
	    return BL_READ_TRUNCATED;
	}
    }
    else
	return status;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write static fields from one line of a single-entry VCF file.
 *      Does not write sample data.
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate placeholder such as '.'
 *      rather than the actual data.  Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to the VCF output stream
 *      vcf_call    Pointer to the bl_vcf_t structure to output
 *      field_mask  Bit mask indicating which fields to output
 *
 *  Returns:
 *      The number of items output (as returned by fprintf())
 *
 *  See also:
 *      bl_vcf_read_static_fields(3), bl_vcf_write_ss_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_write_static_fields(bl_vcf_t *vcf_call, FILE *vcf_stream,
	    vcf_field_mask_t field_mask)

{
    char    *chrom = ".",
	    pos_str[BL_POSITION_MAX_DIGITS+1] = ".",
	    *id = ".",
	    *ref = ".",
	    *alt = ".",
	    *qual = ".",
	    *filter = ".",
	    *info = ".",
	    *format = ".";
    
    if ( field_mask & BL_VCF_FIELD_CHROM )
	chrom = vcf_call->chrom;
    if ( field_mask & BL_VCF_FIELD_POS )
	ltostrn(pos_str, vcf_call->pos, 10, BL_POSITION_MAX_DIGITS);
    if ( field_mask & BL_VCF_FIELD_ID )
	id = vcf_call->id;
    if ( field_mask & BL_VCF_FIELD_REF )
	ref = vcf_call->ref;
    if ( field_mask & BL_VCF_FIELD_ALT )
	alt = vcf_call->alt;
    if ( field_mask & BL_VCF_FIELD_QUAL )
	qual = vcf_call->qual;
    if ( field_mask & BL_VCF_FIELD_FILTER )
	filter = vcf_call->filter;
    if ( field_mask & BL_VCF_FIELD_INFO )
	info = vcf_call->info;
    if ( field_mask & BL_VCF_FIELD_FORMAT )
	format = vcf_call->format;
    
    return fprintf(vcf_stream,
	    "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t",
	    chrom, pos_str,
	    id, ref, alt, 
	    qual, filter, info,
	    format);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Write a single-sample VCF call to vcf_stream.
 *      This should only be used with VCF calls that have exactly one
 *      sample column.  For multisample VCFs, use bl_vcf_write_static_fields()
 *      followed by a loop to write the sample data.
 *
 *      If field_mask is not BL_VCF_FIELD_ALL, fields not indicated by a 1
 *      in the bit mask are written as an appropriate placeholder such as '.'
 *      rather than the actual data.  Possible mask values are:
 *
 *      BL_VCF_FIELD_ALL
 *      BL_VCF_FIELD_CHROM
 *      BL_VCF_FIELD_POS
 *      BL_VCF_FIELD_ID
 *      BL_VCF_FIELD_REF
 *      BL_VCF_FIELD_ALT
 *      BL_VCF_FIELD_QUAL
 *      BL_VCF_FIELD_FILTER
 *      BL_VCF_FIELD_INFO
 *      BL_VCF_FIELD_FORMAT
 *
 *  Arguments:
 *      vcf_stream  FILE pointer to the VCF output stream
 *      vcf_call    Pointer to the bl_vcf_t structure to output
 *      field_mask  Bit mask indicating which fields to output
 *
 *  Returns:
 *      The number of items output (as returned by fprintf())
 *
 *  See also:
 *      bl_vcf_read_ss_call(3), bl_vcf_write_static_fields(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

int     bl_vcf_write_ss_call(bl_vcf_t *vcf_call, FILE *vcf_stream,
	    vcf_field_mask_t field_mask)

{
    bl_vcf_write_static_fields(vcf_call, vcf_stream, field_mask);
    return fprintf(vcf_stream, "%s\n", vcf_call->single_sample);
}


// FIXME: Write a new function bl_vcf_read_multi-samples() that uses
// tsv_read_field_malloc() and extends the pointer array on-the-fly
#if 0
char    **bl_vcf_sample_alloc(bl_vcf_t *vcf_call, size_t samples)

{
    size_t  c;
    
    if ( (vcf_call->multi_samples =
	 (char **)xt_malloc(samples,
		    sizeof(*vcf_call->multi_samples))) != NULL )
    {
	for (c = 0; c < samples; ++c)
	{
	    if ( (vcf_call->multi_samples[c] =
		 (char *)xt_malloc(vcf_call->single_sample_array_size,
				sizeof(*vcf_call->multi_samples[c]))) == NULL )
		return NULL;
	}
    }
    return vcf_call->multi_samples;
}
#endif


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

#if 0
int     vcf_phred_add(bl_vcf_t *vcf_call, unsigned char score)

{
    if ( vcf_call->phreds == NULL )
    {
	// fprintf(stderr, "vcf_phred_add(): Allocating initial buffer.\n");
	if ( (vcf_call->phreds = xt_malloc(vcf_call->phred_buff_size,
				    sizeof(*vcf_call->phreds))) == NULL )
	{
	    fprintf(stderr, "vcf_phred_add(): Could not allocate phreds.\n");
	    exit(EX_UNAVAILABLE);
	}
    }
    
    // fprintf(stderr, "vcf_phred_add(): Adding '%c' at %zu\n", score, vcf_call->phred_count);
    vcf_call->phreds[vcf_call->phred_count++] = score;
    vcf_call->phreds[vcf_call->phred_count] = '\0';
    
    if ( vcf_call->phred_count == vcf_call->phred_buff_size )
    {
	vcf_call->phred_buff_size *= 2;
	if ( (vcf_call->phreds = xt_realloc(vcf_call->phreds,
		    vcf_call->phred_buff_size,
		    sizeof(*vcf_call->phreds))) == NULL )
	{
	    fprintf(stderr, "vcf_phred_add(): Could not reallocate phreds.\n");
	    exit(EX_UNAVAILABLE);
	}
    }
    return BL_READ_OK;
}


/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_phred_blank(bl_vcf_t *vcf_call)

{
    memcpy(vcf_call->phreds, "z", 2);
    vcf_call->phred_count = 0;
}

    
/***************************************************************************
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    vcf_phred_free(bl_vcf_t *vcf_call)

{
    if ( vcf_call->phreds != NULL )
    {
	free(vcf_call->phreds);
	vcf_call->phreds = NULL;
	vcf_call->phred_buff_size = BL_VCF_PHRED_BUFF_SIZE;
    }
    vcf_phred_blank(vcf_call);
}
#endif


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Free all memory associated with a VCF call.
 *
 *  Arguments:
 *      vcf_call    Pointer to the bl_vcf_t structure to free.
 *
 *  See also:
 *      bl_vcf_init(3), bl_vcf_sample_alloc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_free(bl_vcf_t *vcf_call)

{
    int     c;
    
    free(vcf_call->chrom);
    free(vcf_call->id);
    free(vcf_call->ref);
    free(vcf_call->alt);
    free(vcf_call->qual);
    free(vcf_call->filter);
    free(vcf_call->info);
    free(vcf_call->format);
    free(vcf_call->single_sample);
    if ( vcf_call->multi_samples != NULL )
    {
	for (c = 0; c < vcf_call->multi_sample_count; ++c)
	    free(vcf_call->multi_samples[c]);
	free(vcf_call->multi_sample_array_sizes);
	free(vcf_call->multi_sample_lens);
	free(vcf_call->multi_samples);
    }
    
    // Is this necessary?
    //bl_vcf_init(vcf_call);
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Initialize a bl_vcf_t structure, allocating default buffer
 *      sizes for some fields.
 *
 *  Arguments:
 *      vcf_call            Pointer to the bl_vcf_t structure to initialize
 *      info_array_size     Maximum size of INFO field in bytes
 *      format_array_size   Maximum size of FORMAT field in bytes
 *      single_sample_array_size   Maxixum size of SAMPLE field in bytes
 *
 *  See also:
 *      bl_vcf_free(3), vcf_read_call(3), bl_vcf_sample_alloc(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_init(bl_vcf_t *vcf_call)

{
    vcf_call->chrom_array_size = 0;
    vcf_call->chrom_len = 0;
    vcf_call->chrom = NULL;

    vcf_call->id_array_size = 0;
    vcf_call->id_len = 0;
    vcf_call->id = NULL;

    vcf_call->ref_array_size = 0;
    vcf_call->ref_len = 0;
    vcf_call->ref = NULL;

    vcf_call->alt_array_size = 0;
    vcf_call->alt_len = 0;
    vcf_call->alt = NULL;
    
    vcf_call->qual_array_size = 0;
    vcf_call->qual_len = 0;
    vcf_call->qual = NULL;
    
    vcf_call->filter_array_size = 0;
    vcf_call->filter_len = 0;
    vcf_call->filter = NULL;
    
    vcf_call->pos = 0;
    vcf_call->info_len = 0;
    vcf_call->ref_count = 0;
    vcf_call->alt_count = 0;
    vcf_call->other_count = 0;
    
    vcf_call->info_array_size = 0;
    vcf_call->info_len = 0;
    vcf_call->info = NULL;
    
    vcf_call->format_array_size = 0;
    vcf_call->format_len = 0;
    vcf_call->format = NULL;
    
    vcf_call->single_sample_array_size = 0;
    vcf_call->single_sample_len = 0;
    vcf_call->single_sample = NULL;
    
    vcf_call->multi_samples = NULL;
    vcf_call->multi_sample_count = 0;
    vcf_call->multi_sample_pointer_array_size = 0;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Convert a comma-separated list of VCF fields to a field bit mask
 *
 *  Arguments:
 *      spec    Character string containing comma-separated field list
 *
 *  Returns:
 *      A vcf_field_mask_t value with bits set for specified fields
 *
 *  See also:
 *      vcf_read_call(3), vcf_write_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-01-22  Jason Bacon Begin
 ***************************************************************************/

vcf_field_mask_t    bl_vcf_parse_field_spec(char *spec)

{
    vcf_field_mask_t    field_mask;
    char            *field_name;
    
    if ( strcmp(spec, "all") == 0 )
    {
	field_mask = BL_VCF_FIELD_ALL;
    }
    else
    {
	field_mask = 0x0;
	while ((field_name = strsep(&spec, ",")) != NULL)
	{
	    if ( strcmp(field_name, "chrom") == 0 )
		field_mask |= BL_VCF_FIELD_CHROM;
	    else if ( strcmp(field_name, "pos") == 0 )
		field_mask |= BL_VCF_FIELD_POS;
	    else if ( strcmp(field_name, "id") == 0 )
		field_mask |= BL_VCF_FIELD_ID;
	    else if ( strcmp(field_name, "ref") == 0 )
		field_mask |= BL_VCF_FIELD_REF;
	    else if ( strcmp(field_name, "alt") == 0 )
		field_mask |= BL_VCF_FIELD_ALT;
	    else if ( strcmp(field_name, "qual") == 0 )
		field_mask |= BL_VCF_FIELD_QUAL;
	    else if ( strcmp(field_name, "filter") == 0 )
		field_mask |= BL_VCF_FIELD_FILTER;
	    else if ( strcmp(field_name, "info") == 0 )
		field_mask |= BL_VCF_FIELD_INFO;
	    else if ( strcmp(field_name, "format") == 0 )
		field_mask |= BL_VCF_FIELD_FORMAT;
	    else
		field_mask = BL_VCF_FIELD_ERROR;
	}
    }
    return field_mask;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Determine if a VCF call is within a SAM alignment, i.e. on the
 *      same chrom and between the start and end positions of the
 *      alignment.
 *
 *  Arguments:
 *      vcf_call    Pointer to bl_vcf_t structure containing VCF call
 *      sam_alignment   Pointer to bl_sam_t structure containing alignment
 *
 *  Returns:
 *      true if the call is within the alignment
 *      false otherwise
 *
 *  See also:
 *      bl_vcf_call_downstream_of_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/


bool    bl_vcf_call_in_alignment(bl_vcf_t *vcf_call, bl_sam_t *sam_alignment)

{
    if ( (strcmp(BL_VCF_CHROM(vcf_call), BL_SAM_RNAME(sam_alignment)) == 0) &&
	 (BL_VCF_POS(vcf_call) >= BL_SAM_POS(sam_alignment)) &&
	 (BL_VCF_POS(vcf_call) <
	    BL_SAM_POS(sam_alignment) + BL_SAM_SEQ_LEN(sam_alignment)) )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Determine if a VCF call is downstream of a SAM alignment.
 *      For the purpose of this function, this could mean on the same
 *      chrom and higher position, or on a later chrom.
 *
 *  Arguments:
 *      vcf_call    Pointer to bl_vcf_t structure containing VCF call
 *      sam_alignment   Pointer to bl_sam_t structure containing alignment
 *
 *  Returns:
 *      true if the call is downstream of the alignment
 *      false otherwise
 *
 *  See also:
 *      bl_vcf_call_in_alignment(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-26  Jason Bacon Begin
 ***************************************************************************/

bool    bl_vcf_call_downstream_of_alignment(bl_vcf_t *vcf_call,
	    bl_sam_t *alignment)

{
    /*fprintf(stderr, "bl_vcf_call_downstream_of_alignment(): %s,%zu,%zu %s,%zu\n",
	    BL_SAM_RNAME(sam_alignment),BL_SAM_POS(sam_alignment),
	    BL_SAM_SEQ_LEN(sam_alignment),
	    BL_VCF_CHROM(vcf_call),BL_VCF_POS(vcf_call));*/
    if ( (BL_SAM_POS(alignment) + BL_SAM_SEQ_LEN(alignment) <= BL_VCF_POS(vcf_call)) &&
	  (strcmp(BL_SAM_RNAME(alignment), BL_VCF_CHROM(vcf_call)) == 0) )
	return true;
    else if ( bl_chrom_name_cmp(BL_SAM_RNAME(alignment), BL_VCF_CHROM(vcf_call)) < 0 )
	return true;
    else
	return false;
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Report VCF input sort error and terminate the process.
 *
 *  Arguments:
 *      vcf_call        Pointer to bl_vcf_t structure with latest call
 *      previous_chrom  Chromosome of previous VCF call
 *      previous_pos    Position of previous VCF call
 *
 *  Returns:
 *      Does not return
 *
 *  See also:
 *      vcf_read_call(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2020-05-27  Jason Bacon Begin
 ***************************************************************************/

void    bl_vcf_call_out_of_order(bl_vcf_t *vcf_call,
	    char *previous_chrom, int64_t previous_pos)

{
    fprintf(stderr, "ad2vcf: Error: VCF input must be sorted by chrom and then position.\n");
    fprintf(stderr, "Found %s,%" PRId64 " after %s,%" PRId64 ".\n",
	    BL_VCF_CHROM(vcf_call), BL_VCF_POS(vcf_call),
	    previous_chrom, previous_pos);
    exit(EX_DATAERR);
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
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of chrom member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->chrom[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the chrom array
 *      new_chrom_element The new value for chrom[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char            new_chrom_element;
 *
 *      if ( bl_vcf_set_chrom_ae(&bl_vcf, c, new_chrom_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_CHROM_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_chrom_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char new_chrom_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->chrom[c] = new_chrom_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for chrom member in a bl_vcf_t structure.
 *      Use this function to set chrom in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_chrom to bl_vcf_ptr->chrom.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_chrom       The new value for chrom
 *      array_size      Size of the chrom array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char            new_chrom;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_chrom_cpy(&bl_vcf, new_chrom, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_CHROM(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_chrom_cpy(bl_vcf_t *bl_vcf_ptr, char new_chrom[], size_t array_size)

{
    if ( new_chrom == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->chrom, new_chrom, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of id member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->id[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the id array
 *      new_id_element  The new value for id[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char            new_id_element;
 *
 *      if ( bl_vcf_set_id_ae(&bl_vcf, c, new_id_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_ID_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_id_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char new_id_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->id[c] = new_id_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for id member in a bl_vcf_t structure.
 *      Use this function to set id in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_id to bl_vcf_ptr->id.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_id          The new value for id
 *      array_size      Size of the id array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char            new_id;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_id_cpy(&bl_vcf, new_id, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_ID(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_id_cpy(bl_vcf_t *bl_vcf_ptr, char new_id[], size_t array_size)

{
    if ( new_id == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->id, new_id, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of ref member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->ref[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the ref array
 *      new_ref_element The new value for ref[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char            new_ref_element;
 *
 *      if ( bl_vcf_set_ref_ae(&bl_vcf, c, new_ref_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_REF_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_ref_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char new_ref_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->ref[c] = new_ref_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for ref member in a bl_vcf_t structure.
 *      Use this function to set ref in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_ref to bl_vcf_ptr->ref.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_ref         The new value for ref
 *      array_size      Size of the ref array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char            new_ref;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_ref_cpy(&bl_vcf, new_ref, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_REF(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_ref_cpy(bl_vcf_t *bl_vcf_ptr, char new_ref[], size_t array_size)

{
    if ( new_ref == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->ref, new_ref, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of alt member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->alt[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the alt array
 *      new_alt_element The new value for alt[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char            new_alt_element;
 *
 *      if ( bl_vcf_set_alt_ae(&bl_vcf, c, new_alt_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_ALT_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_alt_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char new_alt_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->alt[c] = new_alt_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for alt member in a bl_vcf_t structure.
 *      Use this function to set alt in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_alt to bl_vcf_ptr->alt.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_alt         The new value for alt
 *      array_size      Size of the alt array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char            new_alt;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_alt_cpy(&bl_vcf, new_alt, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_ALT(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_alt_cpy(bl_vcf_t *bl_vcf_ptr, char new_alt[], size_t array_size)

{
    if ( new_alt == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->alt, new_alt, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of qual member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->qual[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the qual array
 *      new_qual_element The new value for qual[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char            new_qual_element;
 *
 *      if ( bl_vcf_set_qual_ae(&bl_vcf, c, new_qual_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_QUAL_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_qual_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char new_qual_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->qual[c] = new_qual_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for qual member in a bl_vcf_t structure.
 *      Use this function to set qual in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_qual to bl_vcf_ptr->qual.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_qual        The new value for qual
 *      array_size      Size of the qual array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char            new_qual;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_qual_cpy(&bl_vcf, new_qual, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_QUAL(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_qual_cpy(bl_vcf_t *bl_vcf_ptr, char new_qual[], size_t array_size)

{
    if ( new_qual == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->qual, new_qual, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of filter member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->filter[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the filter array
 *      new_filter_element The new value for filter[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char            new_filter_element;
 *
 *      if ( bl_vcf_set_filter_ae(&bl_vcf, c, new_filter_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_FILTER_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_filter_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char new_filter_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->filter[c] = new_filter_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for filter member in a bl_vcf_t structure.
 *      Use this function to set filter in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_filter to bl_vcf_ptr->filter.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_filter      The new value for filter
 *      array_size      Size of the filter array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char            new_filter;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_filter_cpy(&bl_vcf, new_filter, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_FILTER(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_filter_cpy(bl_vcf_t *bl_vcf_ptr, char new_filter[], size_t array_size)

{
    if ( new_filter == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->filter, new_filter, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for info member in a bl_vcf_t structure.
 *      Use this function to set info in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      info is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_info        The new value for info
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char *          new_info;
 *
 *      if ( bl_vcf_set_info(&bl_vcf, new_info)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_info(bl_vcf_t *bl_vcf_ptr, char * new_info)

{
    if ( new_info == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->info = new_info;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of info member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->info[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the info array
 *      new_info_element The new value for info[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char *          new_info_element;
 *
 *      if ( bl_vcf_set_info_ae(&bl_vcf, c, new_info_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_INFO_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_info_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char  new_info_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->info[c] = new_info_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for info member in a bl_vcf_t structure.
 *      Use this function to set info in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_info to bl_vcf_ptr->info.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_info        The new value for info
 *      array_size      Size of the info array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char *          new_info;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_info_cpy(&bl_vcf, new_info, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_INFO(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_info_cpy(bl_vcf_t *bl_vcf_ptr, char * new_info, size_t array_size)

{
    if ( new_info == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->info, new_info, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for format member in a bl_vcf_t structure.
 *      Use this function to set format in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      format is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_format      The new value for format
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char *          new_format;
 *
 *      if ( bl_vcf_set_format(&bl_vcf, new_format)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_format(bl_vcf_t *bl_vcf_ptr, char * new_format)

{
    if ( new_format == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->format = new_format;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of format member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->format[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the format array
 *      new_format_element The new value for format[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char *          new_format_element;
 *
 *      if ( bl_vcf_set_format_ae(&bl_vcf, c, new_format_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_FORMAT_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_format_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char  new_format_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->format[c] = new_format_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for format member in a bl_vcf_t structure.
 *      Use this function to set format in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_format to bl_vcf_ptr->format.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_format      The new value for format
 *      array_size      Size of the format array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char *          new_format;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_format_cpy(&bl_vcf, new_format, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_FORMAT(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_format_cpy(bl_vcf_t *bl_vcf_ptr, char * new_format, size_t array_size)

{
    if ( new_format == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->format, new_format, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for single_sample member in a bl_vcf_t structure.
 *      Use this function to set single_sample in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      single_sample is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_single_sample The new value for single_sample
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char *          new_single_sample;
 *
 *      if ( bl_vcf_set_single_sample(&bl_vcf, new_single_sample)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_single_sample(bl_vcf_t *bl_vcf_ptr, char * new_single_sample)

{
    if ( new_single_sample == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->single_sample = new_single_sample;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of single_sample member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->single_sample[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the single_sample array
 *      new_single_sample_element The new value for single_sample[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char *          new_single_sample_element;
 *
 *      if ( bl_vcf_set_single_sample_ae(&bl_vcf, c, new_single_sample_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_SINGLE_SAMPLE_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_single_sample_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char  new_single_sample_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->single_sample[c] = new_single_sample_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for single_sample member in a bl_vcf_t structure.
 *      Use this function to set single_sample in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_single_sample to bl_vcf_ptr->single_sample.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_single_sample The new value for single_sample
 *      array_size      Size of the single_sample array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char *          new_single_sample;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_single_sample_cpy(&bl_vcf, new_single_sample, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_SINGLE_SAMPLE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_single_sample_cpy(bl_vcf_t *bl_vcf_ptr, char * new_single_sample, size_t array_size)

{
    if ( new_single_sample == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	// FIXME: Assuming char array is a null-terminated string
	strlcpy(bl_vcf_ptr->single_sample, new_single_sample, array_size);
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_samples member in a bl_vcf_t structure.
 *      Use this function to set multi_samples in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      multi_samples is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_samples The new value for multi_samples
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char **         new_multi_samples;
 *
 *      if ( bl_vcf_set_multi_samples(&bl_vcf, new_multi_samples)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_samples(bl_vcf_t *bl_vcf_ptr, char ** new_multi_samples)

{
    if ( new_multi_samples == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_samples = new_multi_samples;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of multi_samples member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->multi_samples[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the multi_samples array
 *      new_multi_samples_element The new value for multi_samples[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      char **         new_multi_samples_element;
 *
 *      if ( bl_vcf_set_multi_samples_ae(&bl_vcf, c, new_multi_samples_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_MULTI_SAMPLES_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_samples_ae(bl_vcf_t *bl_vcf_ptr, size_t c, char * new_multi_samples_element)

{
    if ( new_multi_samples_element == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_samples[c] = new_multi_samples_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_samples member in a bl_vcf_t structure.
 *      Use this function to set multi_samples in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_multi_samples to bl_vcf_ptr->multi_samples.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_samples The new value for multi_samples
 *      array_size      Size of the multi_samples array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      char **         new_multi_samples;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_multi_samples_cpy(&bl_vcf, new_multi_samples, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_MULTI_SAMPLES(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_samples_cpy(bl_vcf_t *bl_vcf_ptr, char ** new_multi_samples, size_t array_size)

{
    if ( new_multi_samples == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_vcf_ptr->multi_samples[c] = new_multi_samples[c];
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for pos member in a bl_vcf_t structure.
 *      Use this function to set pos in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      pos is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_pos         The new value for pos
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      int64_t         new_pos;
 *
 *      if ( bl_vcf_set_pos(&bl_vcf, new_pos)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_pos(bl_vcf_t *bl_vcf_ptr, int64_t new_pos)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->pos = new_pos;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for info_array_size member in a bl_vcf_t structure.
 *      Use this function to set info_array_size in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      info_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_info_array_size The new value for info_array_size
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_info_array_size;
 *
 *      if ( bl_vcf_set_info_array_size(&bl_vcf, new_info_array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_info_array_size(bl_vcf_t *bl_vcf_ptr, size_t new_info_array_size)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->info_array_size = new_info_array_size;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for info_len member in a bl_vcf_t structure.
 *      Use this function to set info_len in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      info_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_info_len    The new value for info_len
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_info_len;
 *
 *      if ( bl_vcf_set_info_len(&bl_vcf, new_info_len)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_info_len(bl_vcf_t *bl_vcf_ptr, size_t new_info_len)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->info_len = new_info_len;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for format_array_size member in a bl_vcf_t structure.
 *      Use this function to set format_array_size in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      format_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_format_array_size The new value for format_array_size
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_format_array_size;
 *
 *      if ( bl_vcf_set_format_array_size(&bl_vcf, new_format_array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_format_array_size(bl_vcf_t *bl_vcf_ptr, size_t new_format_array_size)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->format_array_size = new_format_array_size;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for format_len member in a bl_vcf_t structure.
 *      Use this function to set format_len in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      format_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_format_len  The new value for format_len
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_format_len;
 *
 *      if ( bl_vcf_set_format_len(&bl_vcf, new_format_len)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_format_len(bl_vcf_t *bl_vcf_ptr, size_t new_format_len)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->format_len = new_format_len;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for single_sample_array_size member in a bl_vcf_t structure.
 *      Use this function to set single_sample_array_size in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      single_sample_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_single_sample_array_size The new value for single_sample_array_size
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_single_sample_array_size;
 *
 *      if ( bl_vcf_set_single_sample_array_size(&bl_vcf, new_single_sample_array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_single_sample_array_size(bl_vcf_t *bl_vcf_ptr, size_t new_single_sample_array_size)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->single_sample_array_size = new_single_sample_array_size;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for single_sample_len member in a bl_vcf_t structure.
 *      Use this function to set single_sample_len in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      single_sample_len is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_single_sample_len The new value for single_sample_len
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_single_sample_len;
 *
 *      if ( bl_vcf_set_single_sample_len(&bl_vcf, new_single_sample_len)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_single_sample_len(bl_vcf_t *bl_vcf_ptr, size_t new_single_sample_len)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->single_sample_len = new_single_sample_len;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_sample_pointer_array_size member in a bl_vcf_t structure.
 *      Use this function to set multi_sample_pointer_array_size in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      multi_sample_pointer_array_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_sample_pointer_array_size The new value for multi_sample_pointer_array_size
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_multi_sample_pointer_array_size;
 *
 *      if ( bl_vcf_set_multi_sample_pointer_array_size(&bl_vcf, new_multi_sample_pointer_array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_pointer_array_size(bl_vcf_t *bl_vcf_ptr, size_t new_multi_sample_pointer_array_size)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_sample_pointer_array_size = new_multi_sample_pointer_array_size;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_sample_count member in a bl_vcf_t structure.
 *      Use this function to set multi_sample_count in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      multi_sample_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_sample_count The new value for multi_sample_count
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_multi_sample_count;
 *
 *      if ( bl_vcf_set_multi_sample_count(&bl_vcf, new_multi_sample_count)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_count(bl_vcf_t *bl_vcf_ptr, size_t new_multi_sample_count)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_sample_count = new_multi_sample_count;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_sample_array_sizes member in a bl_vcf_t structure.
 *      Use this function to set multi_sample_array_sizes in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      multi_sample_array_sizes is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_sample_array_sizes The new value for multi_sample_array_sizes
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t *        new_multi_sample_array_sizes;
 *
 *      if ( bl_vcf_set_multi_sample_array_sizes(&bl_vcf, new_multi_sample_array_sizes)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_array_sizes(bl_vcf_t *bl_vcf_ptr, size_t * new_multi_sample_array_sizes)

{
    if ( new_multi_sample_array_sizes == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_sample_array_sizes = new_multi_sample_array_sizes;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of multi_sample_array_sizes member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->multi_sample_array_sizes[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the multi_sample_array_sizes array
 *      new_multi_sample_array_sizes_element The new value for multi_sample_array_sizes[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      size_t *        new_multi_sample_array_sizes_element;
 *
 *      if ( bl_vcf_set_multi_sample_array_sizes_ae(&bl_vcf, c, new_multi_sample_array_sizes_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_MULTI_SAMPLE_ARRAY_SIZES_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_array_sizes_ae(bl_vcf_t *bl_vcf_ptr, size_t c, size_t  new_multi_sample_array_sizes_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_sample_array_sizes[c] = new_multi_sample_array_sizes_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_sample_array_sizes member in a bl_vcf_t structure.
 *      Use this function to set multi_sample_array_sizes in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_multi_sample_array_sizes to bl_vcf_ptr->multi_sample_array_sizes.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_sample_array_sizes The new value for multi_sample_array_sizes
 *      array_size      Size of the multi_sample_array_sizes array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t *        new_multi_sample_array_sizes;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_multi_sample_array_sizes_cpy(&bl_vcf, new_multi_sample_array_sizes, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_MULTI_SAMPLE_ARRAY_SIZES(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_array_sizes_cpy(bl_vcf_t *bl_vcf_ptr, size_t * new_multi_sample_array_sizes, size_t array_size)

{
    if ( new_multi_sample_array_sizes == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_vcf_ptr->multi_sample_array_sizes[c] = new_multi_sample_array_sizes[c];
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_sample_lens member in a bl_vcf_t structure.
 *      Use this function to set multi_sample_lens in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      multi_sample_lens is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_sample_lens The new value for multi_sample_lens
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t *        new_multi_sample_lens;
 *
 *      if ( bl_vcf_set_multi_sample_lens(&bl_vcf, new_multi_sample_lens)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_lens(bl_vcf_t *bl_vcf_ptr, size_t * new_multi_sample_lens)

{
    if ( new_multi_sample_lens == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_sample_lens = new_multi_sample_lens;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of multi_sample_lens member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->multi_sample_lens[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the multi_sample_lens array
 *      new_multi_sample_lens_element The new value for multi_sample_lens[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      size_t *        new_multi_sample_lens_element;
 *
 *      if ( bl_vcf_set_multi_sample_lens_ae(&bl_vcf, c, new_multi_sample_lens_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_MULTI_SAMPLE_LENS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_lens_ae(bl_vcf_t *bl_vcf_ptr, size_t c, size_t  new_multi_sample_lens_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->multi_sample_lens[c] = new_multi_sample_lens_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for multi_sample_lens member in a bl_vcf_t structure.
 *      Use this function to set multi_sample_lens in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_multi_sample_lens to bl_vcf_ptr->multi_sample_lens.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_multi_sample_lens The new value for multi_sample_lens
 *      array_size      Size of the multi_sample_lens array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t *        new_multi_sample_lens;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_multi_sample_lens_cpy(&bl_vcf, new_multi_sample_lens, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_MULTI_SAMPLE_LENS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_multi_sample_lens_cpy(bl_vcf_t *bl_vcf_ptr, size_t * new_multi_sample_lens, size_t array_size)

{
    if ( new_multi_sample_lens == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_vcf_ptr->multi_sample_lens[c] = new_multi_sample_lens[c];
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for ref_count member in a bl_vcf_t structure.
 *      Use this function to set ref_count in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      ref_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_ref_count   The new value for ref_count
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      unsigned        new_ref_count;
 *
 *      if ( bl_vcf_set_ref_count(&bl_vcf, new_ref_count)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_ref_count(bl_vcf_t *bl_vcf_ptr, unsigned new_ref_count)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->ref_count = new_ref_count;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for alt_count member in a bl_vcf_t structure.
 *      Use this function to set alt_count in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      alt_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_alt_count   The new value for alt_count
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      unsigned        new_alt_count;
 *
 *      if ( bl_vcf_set_alt_count(&bl_vcf, new_alt_count)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_alt_count(bl_vcf_t *bl_vcf_ptr, unsigned new_alt_count)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->alt_count = new_alt_count;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for other_count member in a bl_vcf_t structure.
 *      Use this function to set other_count in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      other_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_other_count The new value for other_count
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      unsigned        new_other_count;
 *
 *      if ( bl_vcf_set_other_count(&bl_vcf, new_other_count)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_other_count(bl_vcf_t *bl_vcf_ptr, unsigned new_other_count)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->other_count = new_other_count;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for phreds member in a bl_vcf_t structure.
 *      Use this function to set phreds in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      phreds is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_phreds      The new value for phreds
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      unsigned char * new_phreds;
 *
 *      if ( bl_vcf_set_phreds(&bl_vcf, new_phreds)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_phreds(bl_vcf_t *bl_vcf_ptr, unsigned char * new_phreds)

{
    if ( new_phreds == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->phreds = new_phreds;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for an array element of phreds member in a bl_vcf_t
 *      structure. Use this function to set bl_vcf_ptr->phreds[c]
 *      in a bl_vcf_t object from non-member functions.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      c               Subscript to the phreds array
 *      new_phreds_element The new value for phreds[c]
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          c;
 *      unsigned char * new_phreds_element;
 *
 *      if ( bl_vcf_set_phreds_ae(&bl_vcf, c, new_phreds_element)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_PHREDS_AE(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_phreds_ae(bl_vcf_t *bl_vcf_ptr, size_t c, unsigned char  new_phreds_element)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->phreds[c] = new_phreds_element;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for phreds member in a bl_vcf_t structure.
 *      Use this function to set phreds in a bl_vcf_t object
 *      from non-member functions.  This function copies the array pointed to
 *      by new_phreds to bl_vcf_ptr->phreds.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_phreds      The new value for phreds
 *      array_size      Size of the phreds array.
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      unsigned char * new_phreds;
 *      size_t          array_size;
 *
 *      if ( bl_vcf_set_phreds_cpy(&bl_vcf, new_phreds, array_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      BL_VCF_SET_PHREDS(3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_phreds_cpy(bl_vcf_t *bl_vcf_ptr, unsigned char * new_phreds, size_t array_size)

{
    if ( new_phreds == NULL )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	size_t  c;
	
	// FIXME: Assuming all elements should be copied
	for (c = 0; c < array_size; ++c)
	    bl_vcf_ptr->phreds[c] = new_phreds[c];
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for phred_count member in a bl_vcf_t structure.
 *      Use this function to set phred_count in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      phred_count is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_phred_count The new value for phred_count
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_phred_count;
 *
 *      if ( bl_vcf_set_phred_count(&bl_vcf, new_phred_count)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_phred_count(bl_vcf_t *bl_vcf_ptr, size_t new_phred_count)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->phred_count = new_phred_count;
	return BL_VCF_DATA_OK;
    }
}


/***************************************************************************
 *  Library:
 *      #include <biolibc/vcf.h>
 *      -lbiolibc -lxtend
 *
 *  Description:
 *      Mutator for phred_buff_size member in a bl_vcf_t structure.
 *      Use this function to set phred_buff_size in a bl_vcf_t object
 *      from non-member functions.  This function performs a direct
 *      assignment for scalar or pointer structure members.  If
 *      phred_buff_size is a pointer, data previously pointed to should
 *      be freed before calling this function to avoid memory
 *      leaks.
 *
 *  Arguments:
 *      bl_vcf_ptr      Pointer to the structure to set
 *      new_phred_buff_size The new value for phred_buff_size
 *
 *  Returns:
 *      BL_VCF_DATA_OK if the new value is acceptable and assigned
 *      BL_VCF_DATA_OUT_OF_RANGE otherwise
 *
 *  Examples:
 *      bl_vcf_t        bl_vcf;
 *      size_t          new_phred_buff_size;
 *
 *      if ( bl_vcf_set_phred_buff_size(&bl_vcf, new_phred_buff_size)
 *              == BL_VCF_DATA_OK )
 *      {
 *      }
 *
 *  See also:
 *      (3)
 *
 *  History: 
 *  Date        Name        Modification
 *  2022-04-24  gen-get-set Auto-generated from vcf.h
 ***************************************************************************/

int     bl_vcf_set_phred_buff_size(bl_vcf_t *bl_vcf_ptr, size_t new_phred_buff_size)

{
    if ( false )
	return BL_VCF_DATA_OUT_OF_RANGE;
    else
    {
	bl_vcf_ptr->phred_buff_size = new_phred_buff_size;
	return BL_VCF_DATA_OK;
    }
}
