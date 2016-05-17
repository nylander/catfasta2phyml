catfasta2phyml
==============

NAME

    catfasta2phyml.pl -- Concatenate FASTA alignments to PHYML, PHYLIP, or
    FASTA format

SYNOPSIS

    catfasta2phyml.pl [options] [files]

OPTIONS

    -h, -?, --help
            Print a brief help message and exits.

    -m, --man
            Prints the manual page and exits.

    -c, --concatenate
            Concatenate files even when number of taxa differ among
            alignments. Missing data will be filled with all gap (-)
            sequences.

    -f, --fasta
            Print output in FASTA format. Default is PHYML format.

    -p, --phylip
            Print output in a strict PHYLIP format. See
            http://evolution.genetics.washington.edu/phylip/doc/sequence.htm
            l.
            
            Note: the current output format is not entirely strict. Left
            to do is to efficiently print sequences in blocks of 10 characters.
            Check output if using this option.

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
    (interleaved PHYML or FASTA format) after checking that all sequences
    are aligned (of same length). If there are sequence labels that are not
    present in all files, a warning will be issued. Sequenced can, however,
    still be concatenated (and missing sequences be filled with missing data
    (gaps)) if the argument --concatenate is used.

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

    To check fasta alignments:

        catfasta2phyml.pl --noprint --verbose *.fas
        catfasta2phyml.pl -nv *.fas
        catfasta2phyml.pl -n *.fas

    To concatenate fasta files, while filling in missing taxa:

        catfasta2phyml.pl --concatenate --verbose *.fas

AUTHOR

    Written by Johan A. A. Nylander

DEPENDENCIES

    Uses Perl modules Getopt::Long and Pod::Usage

LICENSE AND COPYRIGHT

    Copyright (c) 2010, 2011, 2012, 2013, 2014, 2015, 2016 Johan Nylander. All
    rights reserved.

    This program is free software; you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation; either version 2 of the License, or (at your
    option) any later version.

    This program is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General
    Public License for more details. http://www.gnu.org/copyleft/gpl.html

DOWNLOAD

    https://github.com/nylander/catfasta2phyml

