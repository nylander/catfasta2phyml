catfasta2phyml
==============

Concatenates FASTA formatted files to one "phyml" (PHYLIP) formatted file


NAME
    catfasta2phyml.pl -- Concatenate FASTA alignments to PHYML or FASTA
    format

SYNOPSIS
    catfasta2phyml.pl [options] [files]

OPTIONS
    -h, -?, --help
            Print a brief help message and exits.

    -m, --man
            Prints the manual page and exits.

    -f, --fasta
            Print output in FASTA format. Default is PHYML format.

    -r, --relaxed-phylip
            Print output in relaxed PHYLIP format. That is, sequence labels
            are only printed once, and, the "relaxed", labels can have more
            than 8 characters.

    -v, --verbose
            Be verbose.

    -d, --dont-print
            Do not print the concatenation. Program returns 1 on exit.

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
        catfasta2phyml.pl --verbose *.fas > out.phy

    To concatenate fasta files to fasta format:

        catfasta2phyml.pl -f file1 file2 > out.fasta
        catfasta2phyml.pl -f *.fas > out.fasta

    To check fasta alignments

        catfasta2phyml.pl --dont-print --verbose *.fasta
        catfasta2phyml.pl -d *.fasta

AUTHOR
    Written by Johan A. A. Nylander

DEPENDENCIES
    Uses Perl modules Getopt::Long and Pod::Usage

LICENSE AND COPYRIGHT
    Copyright (c) 2010, 2011, 2012 Johan Nylander. All rights reserved.

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details. http://www.gnu.org/copyleft/gpl.html

