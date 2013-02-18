#!/usr/bin/perl 

## 03/12/2010 12:20:35 PM CET
## 02/18/2013 04:25:17 PM
## TODO: print strict PHYLIP output

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
Getopt::Long::Configure("bundling_override");


#---------------------------------------------------------------------------
#  Global variables
#---------------------------------------------------------------------------
my %HoH              = ();   # 
my %seqid_HoH        = ();   # 
my %nseq_hash        = ();   # key:infile, val:nseq
my %seqid_count_hash = ();   # key:seqid, val:count
my $term             = $/;   # input record separator
my @hash_ref_array   = ();   # array with hash references
my $nfiles           = 0;    # count number of files
my $space            = "\t"; # spacer for aligned print
my $nchar            = 0;    # nchar for phyml header.
my $nseq;                    # nseq for phyml header. Do not initiate!
my $first_name       = q{};  # First name in matrix.
my $fasta            = 0;    # Print phyml format by default
my $man              = 0;    # Manual
my $help             = 0;    # Help
my $verbose          = 0;    # Verbose
my $noprint          = 0;    # Do not print the concatenation
my $sequential       = 0;    # Print sequential with line breaks in in sequence (default is interleaved) 
my $strict_phylip    = 0;    # Print strict phylip format (http://evolution.genetics.washington.edu/phylip/doc/sequence.html)
my $lwidth           = 60;   # default line width for fasta


#---------------------------------------------------------------------------
#  Handle arguments
#---------------------------------------------------------------------------
if (@ARGV < 1) {
    die "No arguments. Try:\n\n $0 -man\n\n";
}
else {
    GetOptions('h|help|?'         => sub { pod2usage(1) },
               'm|man'            => sub { pod2usage(-exitstatus => 0, -verbose => 2) },
               'f|fasta'          => \$fasta,
               'v|verbose'        => \$verbose,
               'n|noprint'        => \$noprint,
               's|sequential'     => \$sequential,
               'p|phylip'         => \$strict_phylip,
              );
}


#---------------------------------------------------------------------------
#  Read all infiles to count sequences
#---------------------------------------------------------------------------
print STDERR "\nCecking sequences in infiles...\n\n" if ($verbose);

foreach my $arg (@ARGV) {

    my $infile  = $arg;
    my %seq_hash = parse_fasta($infile); # key: seqid, value:sequence

    print STDERR "  File $infile: " if ($verbose);

    $nfiles++;

    ## Save sequences in array with hash references. Does this work for really large number of fasta files?
    my $hash_ref     = \%seq_hash;
    push(@hash_ref_array, $hash_ref);

    ## Add nseqs to global nseq_hash:
    $nseq_hash{$infile} = scalar(keys(%seq_hash));

    ## Get length of sequence for all tax labels. Put in hashes.
    foreach my $tax_key (keys %seq_hash) {
        $seqid_count_hash{$tax_key}++;
        $HoH{$infile}{$tax_key} = length($seq_hash{$tax_key});
        $seqid_HoH{$infile}{$tax_key}++;
    }

    ## Check all seqs are same length
    my $length;
    my $lname;
    foreach my $name (keys %seq_hash) {
        my $l = length $seq_hash{$name};
        if (defined $length) {
            if ($length != $l) {
                print STDERR "Error!\nSequences in $infile not all same length ($lname is $length, $name is $l)\n";
                exit(1);
            }
        }
        else {
            $length = length $seq_hash{$name};
            $lname  = $name;
        }
    }
    print STDERR " ntax=$nseq_hash{$infile} nchar=$length\n" if ($verbose);

} # Done with file


#---------------------------------------------------------------------------
#  Check if the same number of sequences
#---------------------------------------------------------------------------
my $lname;
foreach my $file (keys %nseq_hash) {
    my $l = $nseq_hash{$file}; # val is a length
    if (defined $nseq) {
        if ($nseq != $l) {
            print STDERR "Error!\nNumber of sequences in files differ ($lname has $nseq, $file has $l)\n";
            exit(1);
        }
    }
    else {
        $nseq = $nseq_hash{$file};
        $lname  = $file;
    }
}


#---------------------------------------------------------------------------
#  Check sequence id's
#---------------------------------------------------------------------------
if (scalar((keys %seqid_count_hash)) != $nseq) { # number of unique seqid's not eq to nseqs
    foreach my $key (sort { $seqid_count_hash{$b} <=> $seqid_count_hash{$a} } (keys %seqid_count_hash)) {
        print STDERR "$key --> $seqid_count_hash{$key}\n";
    }
    print STDERR "\nError!\nSome sequence labels does not occur in all files.\n";
    print STDERR "That is, sequence id's needs to be identical for concatenation.\n\n";
    exit(1);
}
else {
    ## Find the longest taxon name for aligned printing
    my @sorted_names = sort { length($b) <=> length($a) } keys %seqid_count_hash;
    $space = length( shift(@sorted_names) ) + 2;
    $first_name = $sorted_names[0];
}


#---------------------------------------------------------------------------
#  Get nchar
#---------------------------------------------------------------------------
foreach my $h_ref (@hash_ref_array) {
    $nchar = $nchar + length($h_ref->{$first_name});
} 


#---------------------------------------------------------------------------
#  Print everything to STDOUT
#---------------------------------------------------------------------------
if ($verbose) {
    print STDERR "\nChecked $nfiles files -- sequence labels and lengths seems OK.\n";
    print STDERR "Concatenated $nseq sequences, length $nchar.\n";
    print STDERR "Printing concatenation to STDOUT.\n\n" unless $noprint;
}

## Noprint?
if ($noprint) {
    print STDERR "\nEnd of script.\n\n" if ($verbose);
    exit(1);
}

## Phyml header?
if ($strict_phylip) {
    print STDOUT "   $nseq    $nchar\n";
}
else {
    print STDOUT "$nseq $nchar\n" unless $fasta;
}

## Print the array with hash references (does this work with really large number of files (hashes))?
if ($fasta or $sequential) {
    ## First, concatenate all sequences from hashes
    my %print_hash = (); # key:label, value:sequence
    foreach my $h_ref (@hash_ref_array) {
        foreach my $seqid (sort keys %$h_ref) {
            $print_hash{$seqid} .= $h_ref->{$seqid};
        }
    }
    ## Then print, and add line breaks in sequences
    foreach my $label (sort keys  %print_hash) {
        if ($fasta) {
            print STDOUT ">$label\n";
        }
        elsif ($strict_phylip) {
            my $phylip_label = phylip_label($label);
            print STDOUT "$phylip_label\n";
        }
        else {
            print STDOUT "$label\n";
        }
        ## Print sequence
        ## TODO: phylip strict printing of sequence in blocks of 10
        $print_hash{$label} =~ s/\S{$lwidth}/$&\n/gs; ## replace word of size $lwidth with itself and "\n"
        print STDOUT $print_hash{$label}, "\n";
    }
}
else { # default: phyml interleaved
    my $did_first = 0;
    foreach my $h_ref (@hash_ref_array) {
        foreach my $seqid (sort keys %$h_ref) {
            if ($strict_phylip) {
                my $phylip_seqid = phylip_label($seqid);
                print STDOUT $phylip_seqid unless $did_first;
            }
            else {
                printf STDOUT "%-${space}s ", $seqid unless $did_first;
            }
            ## Print sequence
            ## TODO: phylip strict printing of sequence in blocks of 10
            ## TODO: print length of 60
            print STDOUT "$h_ref->{$seqid}\n";
        }
        print "\n";
        $did_first = 1;
    } 
}

print STDERR "End of script.\n\n" if ($verbose);


#===  FUNCTION  ================================================================
#         NAME:  parse_fasta
#      VERSION:  11/07/2011 05:25:35 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  filename
#      RETURNS:  hash
#         TODO:  ???
#===============================================================================
sub parse_fasta {

    my ($infile) = @_;

    my $term     = $/; # input record separator;
    my %seq_hash = (); # key:seqid, val:seq

    open my $INFILE, "<", $infile or die "could not open infile '$infile' : $! \n"; 
    $/ = ">";
    while(<$INFILE>) {
        chomp;
        next if($_ eq '');
        my ($id, @sequencelines) = split /\n/;
        foreach my $line (@sequencelines) {
            $seq_hash{$id} .= $line;
        }
    }
    $/ = $term;

    return(%seq_hash);

} # end of parse_fasta


#===  FUNCTION  ================================================================
#         NAME:  phylip_label
#      VERSION:  02/18/2013 04:43:00 PM
#  DESCRIPTION:  manipulates input string to be of length 10, possibly padded with
#                white space.
#   PARAMETERS:  string
#      RETURNS:  string with the length of 10
#         TODO:  ???
#===============================================================================
sub phylip_label {

    my ($string) = @_;

    if (length($string) > 10) {
        $string = substr($string, 0, 10); # truncate label!
    }
    else {
        my $string_length = length($string);
        my $pad = ' ' x ((10 - $string_length)); # pad end with white space
        $string = $string . $pad;
    }

    return($string);

} # end of phylip_label





#===  POD DOCUMENTATION  =======================================================
#      VERSION:  02/18/2013 04:44:16 PM
#  DESCRIPTION:  Documentation
#         TODO:  ?
#===============================================================================
=pod

=head1 NAME

catfasta2phyml.pl -- Concatenate FASTA alignments to PHYML or FASTA format


=head1 SYNOPSIS

catfasta2phyml.pl [options] [files]


=head1 OPTIONS

=over 8

=item B<-h, -?, --help>

Print a brief help message and exits.


=item B<-m, --man>

Prints the manual page and exits.


=item B<-f, --fasta>

Print output in FASTA format. Default is PHYML format.


=item B<-p, --phylip>

[Not fully implemented]
Print output in a strict PHYLIP format.
See http://evolution.genetics.washington.edu/phylip/doc/sequence.html.



=item B<-v, --verbose>

Be verbose by showing some useful output. See the combination with B<-n>.


=item B<-n, --noprint>

Do not print the concatenation, just check if all files have the same sequence lables and lengths.
Program returns 1 on exit. See also the combination with B<-v>.


=back



=head1 DESCRIPTION

B<catfasta2phyml.pl> will concatenate FASTA alignments to one
file (interleaved PHYML or FASTA format) after checking that
all sequence labels are present in all files, and that sequences
are aligned (of same length).

Prints to STDOUT.


=head1 USAGE

To concatenate fasta files to a phyml readable format:

    catfasta2phyml.pl file1.fas file2.fas > out.phy
    catfasta2phyml.pl *.fas > out.phy
    catfasta2phyml.pl --sequential *.fas > out.phy
    catfasta2phyml.pl --verbose *.fas > out.phy

To concatenate fasta files to fasta format:

    catfasta2phyml.pl -f file1.fas file2.fas > out.fasta
    catfasta2phyml.pl -f *.fas > out.fasta

To check fasta alignments

    catfasta2phyml.pl --noprint --verbose *.fas
    catfasta2phyml.pl -nv *.fas
    catfasta2phyml.pl -n *.fas


=head1 AUTHOR

Written by Johan A. A. Nylander



=head1 DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage



=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010, 2011, 2012, 2013 Johan Nylander. All rights reserved.

This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details. 
http://www.gnu.org/copyleft/gpl.html 


=cut


__END__

