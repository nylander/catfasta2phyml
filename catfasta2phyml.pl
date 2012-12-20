#!/usr/bin/perl 
#===============================================================================
#
#         FILE:  catfasta2phyml.pl
#
#        USAGE:  ./catfasta2phyml.pl *.fasta > out
#
#  DESCRIPTION:  Concatenates fasta alignments to a phyml file
#
#      OPTIONS:  ---
# REQUIREMENTS:  ---
#         BUGS:  ---
#        NOTES:  ---
#       AUTHOR:  Johan A. A. Nylander (JN), <jnylander @ users.sourceforge.net>
#      COMPANY:  SU
#      VERSION:  1.0
#      CREATED:  03/12/2010 12:20:35 PM CET
#     REVISION:  12/20/2012 10:38:42 AM
#===============================================================================

use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;


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
my $dontprint        = 0;    # Do not print the concatenation
my $relaxed_phylip   = 0;    # Print relaxed phylip format
my $lwidth           = 60;   # default line width for fasta


#---------------------------------------------------------------------------
#  Handle arguments
#---------------------------------------------------------------------------
if (@ARGV < 1) {
    die "No arguments. Try:\n\n $0 -man\n\n";
}
else {
    GetOptions('help|?'         => sub { pod2usage(1) },
               'man'            => sub { pod2usage(-exitstatus => 0, -verbose => 2) },
               'fasta'          => \$fasta,
               'verbose!'       => \$verbose,
               'dont-print'     => \$dontprint,
               'relaxed-phylip' => \$relaxed_phylip,
              );
}


#---------------------------------------------------------------------------
#  Read all infiles to count sequences
#---------------------------------------------------------------------------
print STDERR "\nCecking sequences in infiles...\n" if ($verbose);

foreach my $arg (@ARGV) {
    my $infile       = $arg;
    my %seq_hash     = parse_fasta($infile); # key: seqid, value:sequence

    $nfiles++;

    ## Save sequences in array with hash references.
    ## Does this work for really large number of fasta files?
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
    print STDERR "Printing concatenation to STDOUT.\n\n" unless $dontprint;
}

## Noprint?
if ($dontprint) {
    print STDERR "\nEnd of script!\n\n" if ($verbose);
    exit(1);
}

## Phyml header?
print STDOUT "$nseq $nchar\n" unless $fasta;

## Print the array with hash references
## (Does this work with really large number of files (hashes))?
if ($fasta) {
    ## First, concatenate all sequences from hashes
    my %print_hash = (); # key:label, value:sequence
    foreach my $h_ref (@hash_ref_array) {
        foreach my $seqid (sort keys %$h_ref) {
            $print_hash{$seqid} .= $h_ref->{$seqid};
        }
    }
    ## Then print, first add line breaks in sequences
    foreach my $label (sort keys  %print_hash) {
        $print_hash{$label} =~ s/\S{$lwidth}/$&\n/gs; ## replace word of size $lwidth with itself and "\n"
        print STDOUT ">$label\n";
        print STDOUT $print_hash{$label}, "\n";
    }
}
else {
    my $did_first = 0;
    foreach my $h_ref (@hash_ref_array) {
        foreach my $seqid (sort keys %$h_ref){
            if ($relaxed_phylip) {
                printf STDOUT "%-${space}s ", $seqid unless $did_first;
            }
            else {
                printf STDOUT "%-${space}s ", $seqid;
            }
            print  STDOUT "$h_ref->{$seqid}\n";
        }
        print "\n\n";
        $did_first = 1;
    } 
}

print STDERR "End of script!\n\n" if ($verbose);


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


#===  POD DOCUMENTATION  =======================================================
#      VERSION:  11/07/2011 06:25:23 PM
#  DESCRIPTION:  Documentation
#         TODO:  ?
#===============================================================================
=pod

=head1 NAME

catfasta2phyml.pl  -- Concatenate FASTA alignments to PHYML or FASTA format



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

=item B<-r, --relaxed-phylip>

Print output in relaxed PHYLIP format.
That is, sequence labels are only printed once, and, the "relaxed", labels can have more than 8 characters.

=item B<-v, --verbose>

Be verbose.

=item B<-d, --dont-print>

Do not print the concatenation. Program returns 1 on exit.


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
    catfasta2phyml.pl --verbose *.fas > out.phy

To concatenate fasta files to fasta format:

    catfasta2phyml.pl -f file1 file2 > out.fasta
    catfasta2phyml.pl -f *.fas > out.fasta

To check fasta alignments

    catfasta2phyml.pl --dont-print --verbose *.fasta
    catfasta2phyml.pl -d *.fasta


=head1 AUTHOR

Written by Johan A. A. Nylander



=head1 DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage



=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010, 2011, 2012 Johan Nylander. All rights reserved.

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

