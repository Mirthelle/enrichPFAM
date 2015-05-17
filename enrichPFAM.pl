#!/usr/bin/perl -w 

use strict; 
use LWP::Simple;

use IO::File;

use GO::TermFinder;
use GO::AnnotationProvider::AnnotationParser;
use GO::OntologyProvider::OboParser;

use GO::Utils::File    qw (GenesFromFile);
use GO::Utils::General qw (CategorizeGenes);

######################################################################
sub getFASTA {                                                      ##
######################################################################
## Gets FASTA from PFAM url using PFAM accession ID as input.       ##
######################################################################
    my($acc) = @_;
    
    # Getting fasta entry
    my $loc="http://pfam.xfam.org/family/"; 
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

    # Writing text file with one gene name per line.
    my $filename = "names.txt";
    open(my $fh, '>', $filename) or die "Could not open file '$filename' $!";
    print $fh join("\n", @gene_names);
    close $fh;
    return $filename;
}

######################################################################
sub GOterms {                                                       ##
######################################################################
## Takes 3 file names as inputs: gene_list, annotation_file         ##
## (annotated genome) and obo_file (ontology_file).                 ##
## Returns hash with key being aspects (P, C, F) and values being   ##
## arrays with annotated GO terms in that aspect.                   ##
######################################################################
    my($annotation_file, $genes_list) = @_;
    
    my $annotation = GO::AnnotationProvider::AnnotationParser->new(annotationFile=>$annotation_file);
    
    ## Gene categorization (and error control)
        my @genes = GenesFromFile($genes_list);
        unlink $genes_list;
        
        my @standard_genes = ();
        foreach my $gene (@genes) {
            my $st_name = $annotation->standardNameByName($gene);
            push @standard_genes, $st_name;
        }
    return @standard_genes;
}

##############################################
### 	   M A I N 	  P R O G R A M         ##
##############################################
{
	my $PFAM_ID = $ARGV[0];
	my $anno_human_file = "files/gene_association.goa_human";

	my @fasta = getFASTA($PFAM_ID);
	my $gene_names = humanGenesNames(@fasta);
	
	my @gene_st_names = GOterms($anno_human_file, $gene_names);
    print join(",", @gene_st_names);

}