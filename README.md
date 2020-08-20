# catfasta2phyml

### NAME

`catfasta2phyml.pl` -- Concatenate FASTA alignments to PHYML, PHYLIP, or FASTA format

### SYNOPSIS

    catfasta2phyml.pl [options] [files]

### OPTIONS

- **-h, -?, --help**

    Print a brief help message and exits.

- **-m, --man**

    Prints the manual page and exits.

- **-c, --concatenate**

    Concatenate files even when number of taxa differ among alignments.
    Missing data will be filled with all gap (-) sequences.

- **-i, --intersect**

    Concatenate sequences for sequence labels occuring in all input files
    (intersection).

- **-f, --fasta**

    Print output in FASTA format (default is PHYML format).

- **-p, --phylip**

    Print output in a strict PHYLIP format.
    See [http://evolution.genetics.washington.edu/phylip/doc/sequence.html](http://evolution.genetics.washington.edu/phylip/doc/sequence.html).

    **Note:** The current output is not entirely strict for the
    interleaved format. Left to do is to efficiently print sequences
    in blocks of 10 characters. The sequential PHYLIP format works,
    on the other hand (use **-s** in combination with **-p**).

- **-s, --sequential**

    Print output in sequential format (default is interleaved).

- **-v, --verbose**

    Be verbose by showing some useful output. See the combination with **-n**.

- **-n, --noprint**

    Do not print the concatenation, just check if all files have the same
    sequence lables and lengths. Program returns 1 on exit.
    See also the combination with **-v**.

### DESCRIPTION

**catfasta2phyml.pl** will concatenate FASTA alignments to one file
(interleaved PHYML or FASTA format) after checking that all sequences
are aligned (of same length). If there are sequence labels that are not
present in all files, a warning will be issued. Sequenced can, however,
still be concatenated (and missing sequences be filled with missing data
(gaps)) if the argument **--concatenate** is used.

In addition, only sequences with sequence labels present in all files
(the intersection) can be printed using the **--intersect** argument.

The program prints the concatenated data to **STDOUT**. A table with
information about partitions is printed to **STDERR**. Example: 

    file1.fas = 1-625
    file2.fas = 626-1019
    file3.fas = 1020-2061
    file4.fas = 2062-3364
    file5.fas = 3365-3796

### USAGE

To concatenate fasta files to a phyml readable format:

    $ catfasta2phyml.pl file1.fas file2.fas > out.phy
    $ catfasta2phyml.pl *.fas > out.phy 2> partitions.txt
    $ catfasta2phyml.pl --sequential *.fas > out.phy
    $ catfasta2phyml.pl --verbose *.fas > out.phy

To concatenate fasta files to fasta format:

    $ catfasta2phyml.pl -f file1.fas file2.fas > out.fasta
    $ catfasta2phyml.pl -f *.fas > out.fasta

To check fasta alignments:

    $ catfasta2phyml.pl --noprint --verbose *.fas
    $ catfasta2phyml.pl -nv *.fas
    $ catfasta2phyml.pl -n *.fas

To concatenate fasta files, while filling in missing taxa:

    $ catfasta2phyml.pl --concatenate --verbose *.fas

To concatenate sequences for sequence labels occuring in all files:

    $ catfasta2phyml.pl --intersect *.fas

### TIPS

**1. "Argument list too long" error?**

If we run into the issue of "Argument list too long" (where we have a command
line longer than allowed on our system (`getconf ARG_MAX`) - which may happen
if we try to concatenate many files), we can still do it, but in steps. For
example (here with some help of [GNU
parallel](https://www.gnu.org/software/parallel/)):

    $ catfasta2phyml.pl -c $(find . -type f -name '*.ali') > concatenated.phy 2>/dev/null
    -bash: catfasta2phyml.pl: Argument list too long

Instead, start by concatenating to intermediate files using GNU parallel

    $ find . -type f -name '*.ali' | \
          parallel -N1000 'catfasta2phyml.pl -c -f '"{}"' > tmp.'"{#}"'.conc'

Then concatenate the intermediate files to one

    $ catfasta2phyml.pl -c tmp.*.conc > concatenated.phy 2>/dev/null
    $ rm tmp.*.conc


**2. But I want to split, not concatenate!**

Facing the "opposite" situation (having a large concatenated fasta file that
you want to split into individual alignments)? If you have a corresponding
partitions file, you may give FastEAR a try
([https://github.com/nylander/FastEAR](https://github.com/nylander/FastEAR))!



### AUTHOR

Written by Johan A. A. Nylander

### DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage

### LICENSE AND COPYRIGHT

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

### DOWNLOAD

<https://github.com/nylander/catfasta2phyml>
