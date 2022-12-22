#############################################################################
#   Description:
#       Split a BED file containing feature info into major features
#       such as genes and pseudogenes.  Major features are separated
#       by lines containing ###.  The BED file is typically the augmented
#       conversion of a GFF generated by peak-classifier.
#
#   History: 
#   Date        Name        Modification
#   2021-06-13  Jason Bacon Begin
#############################################################################

BEGIN {
    OFS=FS
}
{
    # Skip leading header lines
    while ( getline && ($0 ~ "^#") ) {
    }
    
    if ( $4 ~ "gene" ) {
	filename=dir "/" $1 "-" $2 "-" $3 ".bed";
	print $0 > filename;
	while ( getline && ($0 != "###") ) {
	    print $0 >> filename;
	}
	close(filename);
    }
    else {
	# Discard non-gene features
	while ( getline && ($0 != "###") ) {
	}
    }
}