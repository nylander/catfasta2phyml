#!/usr/bin/env perl 


## Fri 28 Sep 2018 01:24:05 PM CEST
## TODO:
## - Print strict interleaved PHYLIP output:
##   ./catfasta2phyml.pl -v -c -p testing2/* > outfile
## - Avoid tmp print hash.
##


use strict;
use warnings;
use Data::Dumper;
use Pod::Usage;
use Getopt::Long;
Getopt::Long::Configure("bundling_override");


#---------------------------------------------------------------------------
#  Global variables
#---------------------------------------------------------------------------
# my %HoH = (
#     'file' => {
#                 'aligned' => 1,
#                 'nseqs'   => 10,
#                 'nchars'  => 4,
#                 'seqs'    => {
#                               'apa' => 'ACGT',
#                               'bpa' => 'ACGT',
#                           },
# 
#               },
# 
# );
# ${HoH}{'file'}{'seqs'}{'apa'};

my %HoH              = ();
my %seqids           = ();   # 
my %nseq_hash        = ();   # key:infile, val:nseq
my %nchar_hash       = ();   # key:infile, val:nchar (for aligned data)
my %seqid_count_hash = ();   # key:seqid, val:count
my $term             = $/;   # Input record separator
my @hash_ref_array   = ();   # Array with hash references
my @all_labels_found = ();   # All unique labels in all files
my @intersect_labels = ();   # Labels occuring in all files (intersection)
my @infiles          = ();   # Infiles on ARGV
my $nfiles           = 0;    # Count number of files
my $space            = "\t"; # Spacer for aligned print
my $nchar            = 0;    # nchar for phyml header
my $nseq             = 0;    # nseq for phyml header
my $fasta            = 0;    # Print phyml format by default
my $concatenate      = 0;    # Concatenate after adding missing taxa
my $intersect        = 0;    # Concatenate only intersecting taxa
my $man              = 0;    # Manual
my $help             = 0;    # Help
my $verbose          = 0;    # Verbose
my $noprint          = 0;    # Do not print the concatenation
my $sequential       = 0;    # Print sequential with line breaks in in sequence (default is interleaved) 
my $strict_phylip    = 0;    # Print strict phylip format (http://evolution.genetics.washington.edu/phylip/doc/sequence.html)
my $lwidth           = 60;   # Default line width for fasta
my $nt_counter       = 1;    # Counter for partitions
my $end_count        = 0;    # Counter for partitions


#---------------------------------------------------------------------------
#  Handle arguments
#---------------------------------------------------------------------------
if (@ARGV < 1) {
    die "No arguments. Try:\n\n $0 -man\n\n";
}
else {
    GetOptions('h|help|?'         => sub { pod2usage(1) },
               'm|man'            => sub { pod2usage(-exitstatus => 0, -verbose => 2) },
               'c|concatenate'    => \$concatenate,
               'f|fasta'          => \$fasta,
               'i|intersect'      => \$intersect,
               'n|noprint'        => \$noprint,
               'p|phylip'         => \$strict_phylip,
               's|sequential'     => \$sequential,
               'v|verbose'        => \$verbose,
              );
}


#---------------------------------------------------------------------------
#  Read all infiles, count labels and get sequence lengths
#---------------------------------------------------------------------------
print STDERR "\nChecking sequences in infiles...\n\n" if ($verbose);
foreach my $infile (@ARGV) {

    my $seq_hash_ref = parse_fasta($infile); # key:seqid, value:sequence
    print STDERR "  File $infile: " if ($verbose);
    push (@infiles, $infile);
    $nfiles++;

    ## Are sequences aligned?
    my (@ret) = aligned($seq_hash_ref);
    if (scalar(@ret) > 1) {
        print STDERR "\n\nError! Expecting aligned input sequences.\n";
        print STDERR "Sequences in $infile are not all of the same length:\n$ret[1] is $ret[0], $ret[2] is $ret[3]\n";
        exit(0);
    }
    elsif (scalar(@ret) == 1) {
        $HoH{$infile}{'aligned'} = 1;
        $HoH{$infile}{'nchars'} = $ret[0];
    }

    ## Save sequences from file(s) in one hash. TODO: Mem limit?
    foreach my $key (keys %$seq_hash_ref) {
        $HoH{$infile}{'seqs'}{$key} = ${$seq_hash_ref}{$key};
    }

    ## Get nseqs for file
    $HoH{$infile}{'nseqs'} = scalar(keys %{${HoH}{$infile}{'seqs'}});
    $nseq_hash{$infile} = $HoH{$infile}{'nseqs'};

    ## Get length of sequence for all tax labels.
    foreach my $tax_key (keys %$seq_hash_ref) {
        $seqid_count_hash{$tax_key}++;
    }

    print STDERR " ntax=$HoH{$infile}{'nseqs'} nchar=$HoH{$infile}{'nchars'}\n" if ($verbose);

} # Done with file


#---------------------------------------------------------------------------
# Collect taxon labels and find the length of longest name for printing
#---------------------------------------------------------------------------
(@all_labels_found) = sort { length($b) <=> length($a) } keys %seqid_count_hash;
$space = length($all_labels_found[0]) + 2;


#---------------------------------------------------------------------------
# Get nseq. First check if nseqs are equal among input files
#---------------------------------------------------------------------------
my %string = map { $_, 1 } values %nseq_hash;
if (keys %string == 1) { # all values equal
    my (@nseqs) = values %nseq_hash;
    $nseq = shift(@nseqs);
}
elsif ($concatenate) {
    foreach my $file (keys %HoH) {
        my %second = map {$_=>1} (keys %{${HoH}{$file}{'seqs'}}); # Get all seq labels in file,
        my @only_in_all_labels_found = grep { !$second{$_} } @all_labels_found; # and compare with all labels found in all files
        if (@only_in_all_labels_found) {
            print STDERR "\n  Need to fill missing sequences for file $file:\n" if ($verbose);
            my $allgapseq = '-' x $HoH{$file}{'nchars'};
            foreach my $seqid (@only_in_all_labels_found) {
                print STDERR "    Adding all gaps for seqid $seqid\n" if ($verbose);
                $HoH{$file}{'seqs'}{$seqid} = $allgapseq;
            }
        }
        $nseq = scalar(keys %{$HoH{$file}{'seqs'}});
        $HoH{$file}{'nseqs'} = $nseq;
    }
}
elsif ($intersect) {
    my ($max_nfiles_for_label) = sort { $b <=> $a } values %seqid_count_hash;
    (@intersect_labels) = grep{ $seqid_count_hash{$_} eq $max_nfiles_for_label } keys %seqid_count_hash;
    $nseq = scalar(@intersect_labels);
}
else {
    print STDERR "\n";
    foreach my $key (sort { $seqid_count_hash{$a} <=> $seqid_count_hash{$b} } (keys %seqid_count_hash)) {
        printf STDERR "%-${space}s --> %d\n", $key, $seqid_count_hash{$key};
    }
    print STDERR "\n\nError! Some sequence labels does not occur in all files.\n";
    print STDERR "That is, sequence labels needs to be identical for safe concatenation.\n";
    print STDERR "Use the --concatenate (or -c) to concatenate anyway.\n";
    print STDERR "Empty (all gap) sequences will be added where needed.\n\n";
    print STDERR "Alternatively, you may use --intersect (or -i) to only\n";
    print STDERR "concatenate sequences with labels present in all files.\n\n";
    exit(1);
}


#---------------------------------------------------------------------------
# Check if names can be abbreviated
#---------------------------------------------------------------------------
if ($strict_phylip) {
    my %test = ();
    foreach my $label (@all_labels_found) {
        my $phylabel = phylip_label($label);
        if ($test{$phylabel}++) {
            print STDERR "\nWarning! Strict Phylip format results in duplicate labels for these data (e.g., $phylabel)!\n";
            exit(0);
        }
    }
}


#---------------------------------------------------------------------------
#  Get nchar
#---------------------------------------------------------------------------
foreach my $file (keys %HoH) {
    $nchar = $nchar + ($HoH{$file}{'nchars'});
}


#---------------------------------------------------------------------------
#  Print everything to STDOUT
#---------------------------------------------------------------------------

print STDERR "\nChecked $nfiles files -- sequence labels and lengths seems OK.\n\n" if ($verbose);

if ($noprint) {
    print STDERR "\n\nEnd of script.\n\n" if ($verbose);
    exit(1);
}
elsif($verbose) {
    print STDERR "Printing concatenation to STDOUT, and partition information to STDERR.\n\n";
    print STDERR "Concatenating sequences for $nseq sequence labels, total length $nchar.\n\n";
}

if ($strict_phylip) {
    print STDOUT "   $nseq    $nchar\n";
}
else {
    print STDOUT "$nseq $nchar\n" unless $fasta;
}


#---------------------------------------------------------------------------
# Print the hash via intermediate hash
# TODO: Try to circumvent the intermediate hash (mem limit?!)
#---------------------------------------------------------------------------
if ($fasta or $sequential) {
    ## First, concatenate all sequences from hashes
    my %print_hash = (); # key:label, value:sequence
    foreach my $file (@infiles) { # Keep input order
        die "\n\nError: $file not in HoH\n" unless exists(${HoH}{$file});
        $end_count = $nt_counter + ${HoH}{$file}{'nchars'} - 1;
        print STDERR "$file = $nt_counter-$end_count\n";
        $nt_counter = $nt_counter + ${HoH}{$file}{'nchars'};
        my @seq_ids = ();
        if ($intersect) {
            @seq_ids = @intersect_labels;
        }
        else {
            @seq_ids = keys %{$HoH{$file}{'seqs'}};
        }
        foreach my $seqid (@seq_ids) {
            $print_hash{$seqid} .= $HoH{$file}{'seqs'}{$seqid}; # Concatenate seqs from files
        }
    }
    ## Then print, and add line breaks in sequences
    foreach my $label (sort keys %print_hash) {
        if ($fasta) {
            print STDOUT ">$label\n";
            $print_hash{$label} =~ s/\S{$lwidth}/$&\n/gs; # replace word of size $lwidth with itself and "\n"
            print STDOUT $print_hash{$label}, "\n";
        }
        elsif ($strict_phylip) {
            my $phylip_label = phylip_label($label);
            printf STDOUT "%-10s ", $phylip_label;
            my $s = phylip_blocks($print_hash{$label});
            print $s, "\n";
        }
        else {
            printf STDOUT "%-${space}s ", $label;
            print STDOUT $print_hash{$label}, "\n";
        }
    }
}
else { # default: phyml interleaved, file by file
    my $did_first = 0;
    foreach my $file (@infiles) { # Keep input order
        die "\n\nError: $file not in HoH\n" unless exists(${HoH}{$file});
        $end_count = $nt_counter + ${HoH}{$file}{'nchars'} - 1;
        print STDERR "$file = $nt_counter-$end_count\n";
        $nt_counter = $nt_counter + ${HoH}{$file}{'nchars'};
        my @seq_ids = ();
        if ($intersect) {
            @seq_ids = @intersect_labels;
        }
        else {
            @seq_ids = sort keys %{$HoH{$file}{'seqs'}};
        }
        foreach my $seqid (@seq_ids) {
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
            print STDOUT "$HoH{$file}{'seqs'}{$seqid}\n";
        }
        print "\n";
        $did_first = 1;
    } 
}

print STDERR "\nEnd of script.\n\n" if ($verbose);


#===  FUNCTION  ================================================================
#         NAME:  aligned
#      VERSION:  08/31/2015 07:05:54 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  ref to hash with seqs
#      RETURNS:  0 if aligned, array with names and lengths of the first encountered
#                seqs of unequal length otherwise.
#                "Sequences in $infile are not all of the same length ($lname is $length, $name is $l)"
#         TODO:  ???
#===============================================================================
sub aligned {

    my ($h_ref) = shift(@_);

    my $length;
    my $lname;
    my @aligned = ();

    foreach my $name (keys %$h_ref) {
        my $l = length($h_ref->{$name});
        if (defined $length) {
            if ($length != $l) {
                @aligned = ($length,$lname,$name,$l);
                last;
            }
        }
        else {
            $length = length($h_ref->{$name});
            $lname  = $name;
            @aligned = ($length);
        }
    }

    return @aligned;

} # end of aligned


#===  FUNCTION  ================================================================
#         NAME:  parse_fasta
#      VERSION:  09/27/2018 05:17:10 PM
#  DESCRIPTION:  ???
#   PARAMETERS:  filename
#      RETURNS:  hash ref
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
    close($INFILE);
    $/ = $term;

    return(\%seq_hash);

} # end of parse_fasta


#===  FUNCTION  ================================================================
#         NAME:  phylip_label
#      VERSION:  02/18/2013 04:43:00 PM
#  DESCRIPTION:  manipulates input string to be of length 10, possibly padded
#                with white space.
#   PARAMETERS:  string
#      RETURNS:  string with the length of 10
#         TODO:  ???
#===============================================================================
sub phylip_label {

    my ($string) = @_;

    if (length($string) > 10) {
        $string = substr($string, 0, 10); # Truncate label!
    }
    else {
        my $string_length = length($string);
        my $pad = ' ' x ((10 - $string_length)); # Pad end with white space
        $string = $string . $pad;
    }

    return($string);

} # end of phylip_label


#===  FUNCTION  ================================================================
#         NAME:  phylip_blocks
#      VERSION:  09/03/2015 10:35:12 PM
#  DESCRIPTION:  return string in blocks of ten characters separated by spaces.
#                No more than 6 blocks wide, and 5 for the first row (providing
#                space to sequence label).
#   PARAMETERS:  string
#      RETURNS:  string
#         TODO:  ???
#===============================================================================
sub phylip_blocks {

    my ($string) = @_;

    my $ret_seq = '';
    my $first   = 1;
    my $i       = 0;

    my @foo = unpack("(A10)*", $string);

    foreach my $p (@foo) {
        $i++;
        $ret_seq .= $p;
        if ($i == 5) {
            if ($first) {
                $ret_seq .= "\n";
                $first = 0;
                $i = 0;
            }
            else {
               $ret_seq .= ' ';
            }
        }
        elsif ($i == 6) {
            $ret_seq .= "\n" unless $first;
            $i = 0;
        }
        else {
            $ret_seq .= ' ';
        }
    }

    return $ret_seq;

} # end of phylip_blocks


#===  POD DOCUMENTATION  =======================================================
#      VERSION:  Wed 26 Sep 2018 05:53:29 PM CEST
#  DESCRIPTION:  Documentation
#         TODO:  ?
#===============================================================================
=pod

=head1 NAME

catfasta2phyml.pl -- Concatenate FASTA alignments to PHYML, PHYLIP, or FASTA format


=head1 SYNOPSIS

catfasta2phyml.pl [options] [files]


=head1 OPTIONS

=over 8

=item B<-h, -?, --help>

Print a brief help message and exits.


=item B<-m, --man>

Prints the manual page and exits.


=item B<-c, --concatenate>

Concatenate files even when number of taxa differ among alignments.
Missing data will be filled with all gap (-) sequences.


=item B<-i, --intersect>

Concatenate sequences for sequence labels occuring in all input files
(intersection).


=item B<-f, --fasta>

Print output in FASTA format (default is PHYML format).


=item B<-p, --phylip>

Print output in a strict PHYLIP format.
See L<http://evolution.genetics.washington.edu/phylip/doc/sequence.html>.

B<Note:> The current output is not entirely strict for the
interleaved format. Left to do is to efficiently print sequences
in blocks of 10 characters. The sequential PHYLIP format works,
on the other hand (use B<-s> in combination with B<-p>).


=item B<-s, --sequential>

Print output in sequential format (default is interleaved).


=item B<-v, --verbose>

Be verbose by showing some useful output. See the combination with B<-n>.


=item B<-n, --noprint>

Do not print the concatenation, just check if all files have the same
sequence lables and lengths. Program returns 1 on exit.
See also the combination with B<-v>.


=back

=head1 DESCRIPTION

B<catfasta2phyml.pl> will concatenate FASTA alignments to one file
(interleaved PHYML or FASTA format) after checking that all sequences
are aligned (of same length). If there are sequence labels that are not
present in all files, a warning will be issued. Sequenced can, however,
still be concatenated (and missing sequences be filled with missing data
(gaps)) if the argument B<--concatenate> is used.

In addition, only sequences with sequence labels present in all files
(the intersection) can be printed using the B<--intersect> argument.

The program prints the concatenated data to B<STDOUT>. A table with
information about partitions is printed to B<STDERR>. Example: 

    file1.fas = 1-625
    file2.fas = 626-1019
    file3.fas = 1020-2061
    file4.fas = 2062-3364
    file5.fas = 3365-3796




=head1 USAGE

To concatenate fasta files to a phyml readable format:

    catfasta2phyml.pl file1.fas file2.fas > out.phy
    catfasta2phyml.pl *.fas > out.phy 2> partitions.txt
    catfasta2phyml.pl --sequential *.fas > out.phy
    catfasta2phyml.pl --verbose *.fas > out.phy

To concatenate fasta files to fasta format:

    catfasta2phyml.pl -f file1.fas file2.fas > out.fasta
    catfasta2phyml.pl -f *.fas > out.fasta

To check fasta alignments:

    catfasta2phyml.pl --noprint --verbose *.fas
    catfasta2phyml.pl -nv *.fas
    catfasta2phyml.pl -n *.fas

To concatenate fasta files, while filling in missing taxa:

    catfasta2phyml.pl --concatenate --verbose *.fas

To concatenate sequences for sequence labels occuring in all files:

    catfasta2phyml.pl --intersect *.fas


=head1 AUTHOR

Written by Johan A. A. Nylander


=head1 DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage


=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010-2020 Johan Nylander

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.


=head1 DOWNLOAD

https://github.com/nylander/catfasta2phyml


=cut


__END__

