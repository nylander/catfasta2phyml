catfasta2phyml
==============

NAME

    catfasta2phyml.pl -- Concatenate FASTA alignments to PHYML, PHYLIP, or FASTA format


SYNOPSIS

    catfasta2phyml.pl [options] [files]


OPTIONS

    -h, -?, --help
            Print a brief help message and exits.

    -m, --man
            Prints the manual page and exits.

    -f, --fasta
            Print output in FASTA format. Default is PHYML format.

    -p, --phylip
            [Working, but not yet entirely strict...] Print output in a
            strict PHYLIP format. See
            http://evolution.genetics.washington.edu/phylip/doc/sequence.htm
            l.

    -s, --sequential
            Print output in sequential format. Default is interleaved.

    -v, --verbose
            Be verbose by showing some useful output. See the combination
            with -n.

    -n, --noprint
            Do not print the concatenation, just check if all files have the
            same sequence lables and lengths. Program returns 1 on exit. See
            also the combination with -v.


DESCRIPTION

    catfasta2phyml.pl will concatenate FASTA alignments to one file
    (interleaved PHYML or FASTA format) after checking that all sequence
    labels are present in all files, and that sequences are aligned (of same
    length).

    Prints to STDOUT.


USAGE

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


AUTHOR

    Written by Johan A. A. Nylander


DEPENDENCIES

    Uses Perl modules Getopt::Long and Pod::Usage


LICENSE AND COPYRIGHT

    Copyright (c) 2010, 2011, 2012, 2013 Johan Nylander. All rights
    reserved.

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details. http://www.gnu.org/copyleft/gpl.html

