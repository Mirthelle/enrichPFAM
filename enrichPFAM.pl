#!/usr/bin/perl -w 

use strict; 
use LWP::Simple;

######################################################################
sub getFASTA {                                                      ##
######################################################################
## Gets FASTA from PFAM url using PFAM accession ID as input.       ##
######################################################################
    my($acc) = @_;
    
    # Getting fasta entry
    my $loc = "http://pfam.xfam.org/family/"; 
    my $data = "/alignment/seed/format?format=fasta&alnType=seed&order=a&case=l&gaps=none&download=0";
    my $all = $loc . $acc . $data;
    my $entry = get($all);
    my @fasta = split("\n", $entry);    
}

######################################################################
sub humanGenesNames {                                               ##
######################################################################
## Using a fasta file as input writes a text file with human        ##
## gene names.                                                      ## 
######################################################################
    my(@fasta) = @_;
     
    # Selecting headers correspondig to human genes and storing them in @headers.
    my @headers = ();
    foreach my $line (@fasta) {
        if (($line =~  m/^>/)  and ($line =~ m/_HUMAN/)) {
            my $header = $line;
            push @headers, $header;
        }
    }
    
    # Selecting only gene names and storing them in @gene_names.
    my @gene_names = ();
    foreach my $name (@headers) {
        $name =~ s/>//;     # Remove first >
        my @parts = split("_", $name);  # Split header in "_"
        push @gene_names, $parts[0];    # Select only first element which is gene name
    }
    
    return @gene_names;
}

##############################################
### 	   M A I N 	  P R O G R A M         ##
##############################################
{
	my $PFAM_ID = $ARGV[0];

	my @fasta = getFASTA($PFAM_ID);
	my @gene_names = humanGenesNames(@fasta);
    print join(",", @gene_names);
}
