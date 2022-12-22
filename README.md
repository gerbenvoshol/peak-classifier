# peak-classifier

## Description

peak-classifier classify ChIP/ATAC-Seq peaks based on features provided in
a GFF file.

Peaks are provided in a BED file sorted by chromosome and position.  Typically
these are output from a peak caller such as MACS2, or the differential
analysis that follows.  The GFF must also be sorted by chromosome, position,
and subfeature, which is the default for common data sources.

Peak-classifier generates features that are not explicitly identified in the
GFF, such as introns and potential promoter regions, and outputs the augmented
feature list to a BED file.  It then identifies overlapping features by
running bedtools intersect on the augmented feature list and peak list,
outputting an annotated BED-like TSV file with additional columns to describe
the feature.  If a peak overlaps multiple features, a separate line is output
for each.

Alternative approaches to this problem include R scripting with a tool such
as ChIPpeakAnno or multistage processing of the GFF using awk and bedtools.

In contrast, peak-classifier is a simple Unix command that takes a BED file
and a GFF file as inputs and reports all peak classifications in a matter of
seconds.

Admittedly, an optimal C program isn't really necessary to solve this problem,
since the crappiest implementation I can imagine would not take more than
hours to run for a typical ATAC-Seq peak set.  However:

    * It's an opportunity to develop and test biolibc code that will be
      useful for other problems and bigger data
    * It's more about making peak classification convenient than fast
    * It never hurts to hone your C skills
    * There's no such thing as a program that's too fast

## Design and Implementation

The code is organized following basic object-oriented design principals, but
implemented in C to minimize overhead and keep the source code accessible to
scientists who don't have time to master the complexities of C++.

Structures are treated as classes, with accessor macros and mutator functions
provided, so dependent applications and libraries need not access
structure members directly.  Since the C language cannot enforce this, it's
up to application programmers to exercise self-discipline.

## Building and installing

### Building peak-classifier locally

1. Clone the repository
2. Run "make" to buils
3. Run "make install"

The default install prefix is /usr/local.  

View the Makefile for full details.
