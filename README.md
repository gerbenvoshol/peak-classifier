# peak-classifier

## Description

peak-classifier classify ChIP/ATAC-Seq peaks based on features provided in
a GFF file. Based on [peak-classifier by Paul Auer](https://github.com/auerlab/peak-classifier)

An excerpt from Pauls README:

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

## Differences compared to Pauls version

* Amalgamated libxtend and biolibc for easier build
* Add feature name and IDs to the feature type ouput
* Add optional argument to bedtools location

## Building and installing

### Building peak-classifier locally

1. Clone the repository
2. Run "make" to buils
3. Run "make install"

The default install prefix is /usr/local.  

View the Makefile for full details.
