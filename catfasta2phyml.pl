#!/usr/bin/env perl

use strict;
use warnings;
use Pod::Usage;
use Getopt::Long;
use File::Basename qw(basename);

Getopt::Long::Configure("bundling_override", "no_ignore_case");

my $VERSION = '2.0.0';

# Input/output settings.
my $fasta         = 0;
my $concatenate   = 0;
my $intersect     = 0;
my $noprint       = 0;
my $strict_phylip = 0;
my $sequential    = 0;
my $verbose       = 0;
my $basename_suffix;

my $line_width = 60;
my $chunk_size = 64 * 1024;

#---------------------------------------------------------------------------
# Arguments
#---------------------------------------------------------------------------
die "No arguments. Try:\n\n $0 -man\n\n" unless @ARGV;

# A bare basename option must not consume the next input filename.
# Normalize it to an explicitly empty optional argument.
# Do not modify arguments following the option terminator "--".
for my $arg (@ARGV) {
    last if $arg eq '--';
    if ($arg eq '-b' || $arg eq '--basename') {
        $arg = '--basename=';
    }
}

GetOptions(
    'h|help|?'      => sub { pod2usage(1) },
    'm|man'         => sub {
        pod2usage(-exitstatus => 0, -verbose => 2);
    },
    'b|basename:s'  => \$basename_suffix,
    'c|concatenate' => \$concatenate,
    'f|fasta'       => \$fasta,
    'i|intersect'   => \$intersect,
    'n|noprint'     => \$noprint,
    'p|phylip'      => \$strict_phylip,
    's|sequential'  => \$sequential,
    'v|verbose'     => \$verbose,
    'V|version'     => sub {
        print STDOUT "$0 v.$VERSION\n";
        exit(0);
    },
) or pod2usage(2);

die "Error: No input files given\n" unless @ARGV;

# Each entry contains metadata only:
#
# {
#     path    => filename,
#     nseq    => number of labels,
#     nchar   => alignment length,
#     records => {
#         label => {
#             start  => byte offset immediately after the header,
#             end    => byte offset before the next header, or EOF,
#             length => number of sequence characters,
#         },
#     },
# }
#
# No sequence strings are stored here.
my @inputs;
my %label_count;

#---------------------------------------------------------------------------
# Pass 1: index labels, byte ranges, and sequence lengths
#---------------------------------------------------------------------------
print STDERR "\nChecking sequences in infiles...\n\n" if $verbose;

for my $path (@ARGV) {
    my $input = scan_fasta($path);
    push @inputs, $input;

    # Records are keyed by unique labels, so each label contributes
    # exactly once per input alignment.
    $label_count{$_}++ for keys %{ $input->{records} };

    printf STDERR "  File %s: ntax=%d nchar=%d\n",
        $path, $input->{nseq}, $input->{nchar}
        if $verbose;
}

my $nfiles = scalar @inputs;
my @all_labels = sort keys %label_count;
my $space = max_label_length(\@all_labels) + 2;

#---------------------------------------------------------------------------
# Select output labels
#---------------------------------------------------------------------------
my @labels;
my $nseq;
my $selection;

if ($intersect) {
    # Do this BEFORE checking whether taxon counts are equal.
    @labels = grep { $label_count{$_} == $nfiles } @all_labels;

    unless (@labels) {
        print STDERR "Warning: no labels occur in all files\n";
        exit(1);
    }

    $selection = 'intersection';
    $nseq = scalar @labels;
}
elsif ($concatenate) {
    @labels = @all_labels;
    $selection = 'union';
    $nseq = scalar @labels;

    if ($verbose) {
        for my $input (@inputs) {
            my @missing = grep {
                !exists $input->{records}{$_}
            } @labels;

            next unless @missing;

            print STDERR
                "\n  Need to fill missing sequences for file ",
                "$input->{path}:\n";

            print STDERR "    Adding all gaps for seqid $_\n"
                for @missing;
        }
    }
}
else {
    # Equal numbers of taxa are sufficient, even if labels differ.
    my %taxon_counts = map { $_->{nseq} => 1 } @inputs;

    if (scalar(keys %taxon_counts) != 1) {
        print STDERR "\n";

        for my $label (
            sort {
                $label_count{$a} <=> $label_count{$b}
                    || $a cmp $b
            } @all_labels
        ) {
            printf STDERR "%-${space}s --> %d\n",
                $label, $label_count{$label};
        }

        print STDERR
            "\nError: Input files contain different numbers of taxa.\n",
            "Use --concatenate (-c) to fill missing sequences with gaps,\n",
            "or --intersect (-i) to retain labels present in all files.\n";

        exit(1);
    }

    $selection = 'legacy';
    @labels = @all_labels;
    $nseq = $inputs[0]->{nseq};
}

# FASTA takes precedence over PHYLIP formatting.
if ($strict_phylip && !$fasta) {
    my %seen;

    for my $label (@labels) {
        my $short = phylip_label($label);

        if (exists $seen{$short}) {
            die
                "Error: Strict PHYLIP labels collide: ",
                "'$seen{$short}' and '$label' both become '$short'\n";
        }

        $seen{$short} = $label;
    }
}

my $nchar = 0;
$nchar += $_->{nchar} for @inputs;

print STDERR
    "\nChecked $nfiles files -- input alignments and selection checked.\n\n"
    if $verbose;

if ($noprint) {
    print STDERR "End of script.\n" if $verbose;
    exit(0);
}

if ($verbose) {
    print STDERR
        "Printing concatenation to STDOUT, ",
        "and partition information to STDERR.\n",
        "Header dimensions: $nseq sequence labels, $nchar characters.\n\n";

    if ($selection eq 'legacy' && @labels != $nseq) {
        print STDERR
            "Warning: Equal taxon counts were accepted, but labels differ.\n",
            "Default-mode output may have inconsistent labels or lengths.\n",
            "Use --concatenate to pad missing segments, or --intersect ",
            "to retain shared labels.\n\n";
    }
}

#---------------------------------------------------------------------------
# Output header and partition table
#---------------------------------------------------------------------------
unless ($fasta) {
    if ($strict_phylip) {
        print STDOUT "   $nseq    $nchar\n";
    }
    else {
        print STDOUT "$nseq $nchar\n";
    }
}

print_partitions(\@inputs, $basename_suffix);

#---------------------------------------------------------------------------
# Pass 2: retrieve and print sequence data on demand
#---------------------------------------------------------------------------
if ($fasta || $sequential) {
    # One output label at a time, across all input alignments.
    #
    # Open only one input handle at a time. This avoids descriptor limits
    # for large numbers of files, at the cost of repeated opens/seeks.
    my $mode = $fasta
        ? 'fasta'
        : $strict_phylip
            ? 'phylip'
            : 'plain';

    for my $label (@labels) {
        if ($fasta) {
            print STDOUT ">$label\n";
        }
        elsif ($strict_phylip) {
            print STDOUT phylip_label($label), ' ';
        }
        else {
            printf STDOUT "%-${space}s ", $label;
        }

        # Keep the writer alive across partitions so wrapping is based
        # on the concatenated sequence, not on individual input files.
        my $writer = make_writer($mode, $line_width);

        for my $input (@inputs) {
            if (exists $input->{records}{$label}) {
                my $fh = open_input($input->{path});

                stream_sequence(
                    $fh,
                    $input,
                    $label,
                    $writer->{put},
                    $chunk_size,
                );

                close_input($fh, $input->{path});
            }
            elsif ($selection eq 'union') {
                stream_gaps(
                    $input->{nchar},
                    $writer->{put},
                    $chunk_size,
                );
            }
            elsif ($selection eq 'intersection') {
                # This would indicate an internal metadata error.
                die
                    "Error: Selected intersection label '$label' ",
                    "is missing from '$input->{path}'\n";
            }

            # In legacy default mode, missing segments are skipped,
        }

        $writer->{finish}->();
    }
}
else {
    # Interleaved: each input alignment is one output block.
    my $first_block = 1;

    for my $input (@inputs) {
        my $fh = open_input($input->{path});

        # Preserve the original default behavior when labels differ:
        # each block uses that input file's own sorted labels.
        my @block_labels = $selection eq 'legacy'
            ? sort keys %{ $input->{records} }
            : @labels;

        for my $label (@block_labels) {
            if ($first_block) {
                if ($strict_phylip) {
                    print STDOUT phylip_label($label);
                }
                else {
                    printf STDOUT "%-${space}s ", $label;
                }
            }

            my $put = sub { print STDOUT $_[0]; };

            if (exists $input->{records}{$label}) {
                stream_sequence(
                    $fh,
                    $input,
                    $label,
                    $put,
                    $chunk_size,
                );
            }
            elsif ($selection eq 'union') {
                stream_gaps($input->{nchar}, $put, $chunk_size);
            }
            else {
                die
                    "Error: Selected label '$label' ",
                    "is missing from '$input->{path}'\n";
            }

            print STDOUT "\n";
        }

        close_input($fh, $input->{path});
        print STDOUT "\n";
        $first_block = 0;
    }
}

# Flush and report output errors, including a failure detected on close.
close STDOUT or die "Error: Could not close STDOUT: $!\n";

print STDERR "\nEnd of script.\n\n" if $verbose;
exit(0);

#---------------------------------------------------------------------------
# Open inputs in raw mode so tell/seek/read all use byte offsets.
# Inputs must be regular, seekable files.
#---------------------------------------------------------------------------
sub open_input {
    my ($path) = @_;

    open my $fh, '<:raw', $path
        or die "Error: Could not open infile '$path': $!\n";

    unless (-f $fh) {
        close $fh;
        die "Error: Input '$path' must be a regular, seekable file\n";
    }

    return $fh;
}

sub close_input {
    my ($fh, $path) = @_;

    close $fh
        or die "Error: Could not close infile '$path': $!\n";
}

sub input_position {
    my ($fh, $path) = @_;

    my $position = tell($fh);
    die "Error: Could not determine position in '$path': $!\n"
        if $position < 0;

    return $position;
}

#---------------------------------------------------------------------------
# Pass 1: scan FASTA without retaining sequence data.
#
# Full header text is used as the label.
# Sequence whitespace is removed consistently in both passes.
#---------------------------------------------------------------------------
sub scan_fasta {
    my ($path) = @_;
    my $fh = open_input($path);

    my %records;
    my $label;

    while (1) {
        my $line_start = input_position($fh, $path);

        # Clear errno so EOF can be distinguished from a read error.
        $! = 0;
        my $line = <$fh>;

        unless (defined $line) {
            die "Error: Could not read infile '$path': $!\n" if $!;
            last;
        }

        $line =~ s/\r?\n\z//;

        if ($line =~ /^>(.*)\z/) {
            my $new_label = $1;

            #die "Error: Empty FASTA header in '$path'\n"
            #    unless $new_label =~ /\S/;
            if ($new_label !~ /\S/) {
                print STDERR "Error: Empty FASTA header in $path \n" if ($verbose);
                exit 1;
            }

            #die "Error: Duplicate FASTA header '$new_label' in '$path'\n"
            #    if exists $records{$new_label};
            if (exists $records{$new_label}) {
                print STDERR "Error: Duplicate FASTA header '$new_label' in '$path'\n" if ($verbose);
                exit 1;
            }

            # Close the preceding record before starting the next.
            $records{$label}{end} = $line_start
                if defined $label;

            $label = $new_label;
            $records{$label} = {
                start  => input_position($fh, $path),
                length => 0,
            };
        }
        elsif (defined $label) {
            $line =~ s/\s+//g;
            $records{$label}{length} += length($line);
        }
        elsif ($line =~ /\S/) {
            #die
            #    "Error: Sequence data before the first FASTA header ",
            #    "in '$path'\n";
            print STDERR "Error: Sequence data before the first FASTA header in '$path'\n" if ($verbose);
            exit 1;
        }
    }

    $records{$label}{end} = input_position($fh, $path)
        if defined $label;

    close_input($fh, $path);

    die "Error: Could not read FASTA sequences in '$path'\n"
        unless keys %records;

    my @labels = sort keys %records;
    my $reference_label = $labels[0];
    my $nchar = $records{$reference_label}{length};

    for my $name (@labels) {
        my $length = $records{$name}{length};

        #die "Error: No sequence for header '$name' in '$path'\n"
        #    unless $length;
        if (! $length) {
            print STDERR "Error: No sequence for header '$name' in '$path'\n" if ($verbose);
        }

        if ($length != $nchar) {
            #die
            #    "Error: Expecting aligned input sequences.\n",
            #    "Sequences in '$path' are not all of the same length:\n",
            #    "$reference_label is $nchar, $name is $length\n";
           print STDERR "Error: Expecting aligned input sequences.\n" if ($verbose);
           print STDERR "Sequences in '$path' are not all of the same length:\n" if ($verbose);
           print STDERR "$reference_label is $nchar, $name is $length\n" if ($verbose);
           exit 1;
        }
    }

    return {
        path    => $path,
        nseq    => scalar(@labels),
        nchar   => $nchar,
        records => \%records,
    };
}

#---------------------------------------------------------------------------
# Pass 2: read only the indexed byte range, in bounded chunks.
#
# A sequence can span many lines, or be written on one very long line.
# Neither case requires storing the full sequence during output.
#---------------------------------------------------------------------------
sub stream_sequence {
    my ($fh, $input, $label, $put, $size) = @_;

    my $record = $input->{records}{$label};
    my $path = $input->{path};

    seek($fh, $record->{start}, 0)
        or die "Error: Could not seek in '$path': $!\n";

    my $remaining = $record->{end} - $record->{start};
    my $observed_length = 0;

    while ($remaining > 0) {
        my $wanted = $remaining > $size ? $size : $remaining;
        my $chunk = '';

        my $read = read($fh, $chunk, $wanted);

        die "Error: Could not read infile '$path': $!\n"
            unless defined $read;

        die "Error: Unexpected end of file in '$path'; input may have changed\n"
            if $read == 0;

        $remaining -= $read;

        $chunk =~ s/\s+//g;
        $observed_length += length($chunk);

        $put->($chunk) if length($chunk);
    }

    if ($observed_length != $record->{length}) {
        die
            "Error: Sequence length changed for '$label' in '$path'; ",
            "do not modify inputs while the script is running\n";
    }
}

#---------------------------------------------------------------------------
# Stream missing data without constructing a full gap sequence.
#---------------------------------------------------------------------------
sub stream_gaps {
    my ($length, $put, $size) = @_;

    my $gaps = '-' x ($length > $size ? $size : $length);

    while ($length > 0) {
        my $take = $length > $size ? $size : $length;

        if ($take == length($gaps)) {
            $put->($gaps);
        }
        else {
            $put->(substr($gaps, 0, $take));
        }

        $length -= $take;
    }
}

#---------------------------------------------------------------------------
# Incremental sequence formatter.
#
# plain:  print chunks unchanged, then a newline
# fasta:  wrap at the configured width
# phylip: groups of 10, five groups on the first row, six thereafter
#
# The formatter buffers at most one FASTA line or PHYLIP group.
#---------------------------------------------------------------------------
sub make_writer {
    my ($mode, $width) = @_;

    if ($mode eq 'plain') {
        return {
            put    => sub { print STDOUT $_[0]; },
            finish => sub { print STDOUT "\n"; },
        };
    }

    my $buffer = '';
    my $unit = $mode eq 'fasta' ? $width : 10;

    my $row = 0;
    my $blocks_on_row = 0;

    my $emit = sub {
        my ($text) = @_;

        if ($mode eq 'fasta') {
            print STDOUT $text, "\n";
            return;
        }

        my $row_limit = $row == 0 ? 5 : 6;

        if ($blocks_on_row == $row_limit) {
            print STDOUT "\n";
            $row++;
            $blocks_on_row = 0;
        }
        elsif ($blocks_on_row > 0) {
            print STDOUT ' ';
        }

        print STDOUT $text;
        $blocks_on_row++;
    };

    return {
        put => sub {
            # Consume the chunk in slices. Do not append the entire
            # chunk to the formatting buffer.
            my $offset = 0;
            my $length = length($_[0]);

            while ($offset < $length) {
                my $take = $unit - length($buffer);
                my $available = $length - $offset;
                $take = $available if $available < $take;

                $buffer .= substr($_[0], $offset, $take);
                $offset += $take;

                if (length($buffer) == $unit) {
                    $emit->($buffer);
                    $buffer = '';
                }
            }
        },
        finish => sub {
            $emit->($buffer) if length($buffer);
            $buffer = '';

            # FASTA's emit already terminates each output line.
            print STDOUT "\n" if $mode eq 'phylip';
        },
    };
}

#---------------------------------------------------------------------------
# Partition coordinates retain input argument order.
#---------------------------------------------------------------------------
sub print_partitions {
    my ($inputs, $suffix) = @_;

    my $start = 1;

    for my $input (@$inputs) {
        my $end = $start + $input->{nchar} - 1;
        my $name = $input->{path};

        if (defined $suffix) {
            $name = length($suffix)
                ? basename($name, $suffix)
                : basename($name);
        }

        print STDERR "$name = $start-$end\n";
        $start = $end + 1;
    }
}

sub phylip_label {
    my ($label) = @_;
    return sprintf("%-10s", substr($label, 0, 10));
}

sub max_label_length {
    my ($labels) = @_;

    my $maximum = 0;

    for my $label (@$labels) {
        my $length = length($label);
        $maximum = $length if $length > $maximum;
    }

    return $maximum;
}

__END__

=pod

=head1 NAME

catfasta2phyml.pl -- Concatenate FASTA alignments to PHYML, PHYLIP, or FASTA format

=head1 SYNOPSIS

catfasta2phyml.pl [options] [files]

=head1 OPTIONS

=over 8

=item B<-h, -?, --help>

Print a brief help message and exit.

=item B<-m, --man>

Print the manual page and exit.

=item B<-c, --concatenate>

Concatenate the union of labels across input files. Missing sequences are
filled with gap (-) characters of the corresponding alignment length.
This applies even when input files contain equal numbers of taxa.

=item B<-i, --intersect>

Concatenate only sequences whose labels occur in every input file.

If no labels occur in all files, print a warning to STDERR and exit with
status 1 without printing alignment data.

This option takes precedence over B<--concatenate> if both are supplied.

=item B<-f, --fasta>

Print FASTA output, wrapped at 60 characters per line. This option takes
precedence over PHYLIP formatting.

=item B<-p, --phylip>

Use ten-character PHYLIP labels, truncated or padded as necessary.
Labels that collide after truncation cause an error.

Sequential output uses groups of ten characters, with five groups on
the first row and six on subsequent rows.

Note: Interleaved output is not fully strict PHYLIP (see
L<https://phylipweb.github.io/phylip/doc/sequence.html>).

Use B<-s -p> for sequential PHYLIP output.

=item B<-s, --sequential>

Print sequential output. The default is interleaved output.

=item B<-b, --basename=suffix>

Use file basenames in partition definitions. Remove the supplied suffix
(optional) when it matches the end of the basename.

=item B<-v, --verbose>

Print progress and selection information to STDERR.

=item B<-n, --noprint>

Validate input alignments and the requested label selection without
printing alignment data or partition definitions.

Return status 0 on success, or a nonzero status on failure.

Without B<-c> or B<-i>, validation preserves the original equal-taxon-count
rule; it does not require identical label sets.

=item B<-V, --version>

Print the version number and exit.

=back

=head1 DESCRIPTION

Each input file must contain aligned FASTA sequences: all sequences within
that file must have the same nonzero length.

The first pass records labels, sequence lengths, and byte ranges. Sequence
data are not retained. During output, the script seeks to these ranges and
streams sequence data in bounded chunks.

FASTA headers must start with C<E<gt>>. The complete header text following
that character is used as the sequence label. Duplicate labels within an
input file and empty headers or sequences are rejected.

Whitespace in sequence data, including CRLF line endings, is removed.

Inputs must be regular, seekable, uncompressed files. Standard input,
pipes, and process substitutions are not supported. Input files must not
change while the script runs.

B<--intersect> retains only labels present in all input files.
B<--concatenate> retains all labels and pads missing segments with gaps.

Without either option files with equal numbers of taxa are accepted even when
their labels differ.  Interleaved blocks use each file's own sorted labels.
Sequential and FASTA output omit missing segments rather than padding them.
Consequently, default-mode output for differing label sets can have
inconsistent labels or sequence lengths. Use B<--concatenate> for gap-padded
output.

Alignment data are printed to STDOUT. Partition definitions are printed
to STDERR in input argument order:

    file1.fas = 1-625
    file2.fas = 626-1019
    file3.fas = 1020-2061

=head1 MEMORY AND I/O

Stored metadata grows with the number of input records and their labels,
not with the total number of sequence characters.

The metadata scan reads one input line at a time, so its temporary memory
depends on the longest input line. The output pass reads chunks of at most
64 KiB; formatting buffers hold at most 60 characters.

Sequential and FASTA output reopen input files as needed for each label.
Only one input filehandle is open at a time. Interleaved output opens each
file once during the output pass.

=head1 USAGE

Default interleaved output:

    catfasta2phyml.pl file1.fas file2.fas > out.phy 2> partitions.txt

Sequential output:

    catfasta2phyml.pl -s *.fas > out.phy

Sequential PHYLIP output:

    catfasta2phyml.pl -sp *.fas > out.phy

FASTA output with missing data padded:

    catfasta2phyml.pl -cf *.fas > out.fasta

Labels present in every input file:

    catfasta2phyml.pl -i *.fas > out.phy
    catfasta2phyml.pl -if *.fas > out.fasta

Validation:

    catfasta2phyml.pl -nv *.fas

Basename and suffix removal:

    catfasta2phyml.pl -b dat/file1.fas dat/file2.fas > out.phy
    catfasta2phyml.pl -b'.fas' dat/file1.fas dat/file2.fas > out.phy

=head1 AUTHOR

Written by Johan A. A. Nylander

=head1 DEPENDENCIES

Uses Perl modules Getopt::Long, Pod::Usage, and File::Basename.

=head1 LICENSE AND COPYRIGHT

Copyright (c) 2010-2026 Johan Nylander

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
