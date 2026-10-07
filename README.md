# catfasta2phyml

### NAME

`catfasta2phyml.pl` -- Concatenate FASTA alignments to PHYML, PHYLIP, or FASTA
format

### SYNOPSIS

    catfasta2phyml.pl [options] [files]

### OPTIONS

- **-h, -?, --help**

Print a brief help message and exits.

- **-m, --man**

Prints the manual page and exits.

- **-c, --concatenate**

Concatenate files even when number of taxa differ among alignments. Missing
data will be filled with all gap (-) sequences.

- **-i, --intersect**

Concatenate sequences for sequence labels occuring in all input files
(intersection). This option takes precedence over --concatenate if both are
supplied.

- **-f, --fasta**

Print output in FASTA format (default is PHYML format).

- **-p, --phylip**

Print output in a strict PHYLIP format.

Note: The current output is not entirely strict for the interleaved format.
Left to do is to efficiently print sequences in blocks of 10 characters. The
sequential PHYLIP format works, on the other hand (use **-s** in combination
with **-p**). See section "Data file format" on
[https://phylipweb.github.io/phylip/doc/main.html#inputfiles](https://phylipweb.github.io/phylip/doc/main.html#inputfiles)

- **-s, --sequential**

Print output in sequential format. The default is interleaved output.

- **-b, --basename=suffix**

Use file basenames in partition definitions. Remove the supplied suffix
(optional) when it matches the end of the basename.

- **-v, --verbose**

Be verbose by showing some useful output. See the combination with **-n**.

- **-n, --noprint**

Validate input alignments and the requested label selection without printing
alignment data or partition definitions.

Return status 0 on success, or a nonzero status on failure.

See also the combination with **-v**.

- **-V, --version**

    Print version number and exit.

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

See below for how this table can be used in other software (e.g., IQ-Tree,
RAxML-ng).

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
    $ for f in dat/*.fas ; do
        if ./catfasta2phyml.pl -n "$f" 2> /dev/null ; then
          echo "$f OK"
        else
          echo "$f FAIL"
        fi
      done

To concatenate fasta files, while filling in missing taxa:

    $ catfasta2phyml.pl --concatenate --verbose *.fas

To concatenate sequences for sequence labels occuring in all files:

    $ catfasta2phyml.pl --intersect *.fas

To ensure basename as name and suffix removal in partition definition:

    $ catfasta2phyml.pl -b dat/file1.fas dat/file2.fas > out.phy
    $ catfasta2phyml.pl -b'.fas' dat/file1.fas dat/file2.fas > out.phy


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


**2. Prepare a RAxML-style partitions file**

Catfasta2phyml does not check what data type (DNA, PROTEIN, etc) that is being
concatenated. It only checks the sequence labels and sequence lengths. When
running catfasta2phyml, a list of partition names and relative positions are
written to standard error.  A partition file (for, e.g.,
[IQ-Tree](http://www.iqtree.org/) and
[RAxML-ng]((https://github.com/amkozlov/raxml-ng)) does require, however, a
data type to be given in front of the partition specification.

Assuming that we are concatenating the same kind of data type, the preparation
of a partitions file is straightforward.  Below is an example using `sed` (GNU
Linux). Let us assume that that the data type is "DNA":

    $ ./catfasta2phyml.pl -b'.fas' -c dat/*.fas > out.phy 2> partitions.txt
    $ cat partitions.txt
    file1 = 1-625
    file2 = 626-1019
    file3 = 1020-2061
    file4 = 2062-3364
    file5 = 3365-3796

We can now add `DNA, ` on each line:

    $ sed -i -e 's/^/DNA, /' partitions.txt
    $ cat partitions.txt
    DNA, file1 = 1-625
    DNA, file2 = 626-1019
    DNA, file3 = 1020-2061
    DNA, file4 = 2062-3364
    DNA, file5 = 3365-3796


**3. But I want to split, not concatenate!**

Facing the "opposite" situation (having a large concatenated fasta file that
you want to split into individual alignments)? If you have a corresponding
partitions file, you may give FastEAR a try
([https://github.com/nylander/FastEAR](https://github.com/nylander/FastEAR))!


### AUTHOR

Written by Johan A. A. Nylander

### DEPENDENCIES

Uses Perl modules Getopt::Long and Pod::Usage

### LICENSE AND COPYRIGHT

Copyright (c) 2010-2024 Johan Nylander, [MIT license](LICENSE)

### DOWNLOAD

<https://github.com/nylander/catfasta2phyml>


### INSTALL WITH CONDA

    $ conda install -c bioconda catfasta2phyml


