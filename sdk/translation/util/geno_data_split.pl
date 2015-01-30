#!/usr/bin/env perl
###############################################################################
#
#  geno_data_split.pl
#
#  Revision Control
#  ----------------
#  [$Id: geno_data_split.pl,v 1.5 2009-01-13 18:00:01 mspald Exp $]
#
#  Synopsis:
#  ----------------
#
#  Split a Genotype Short report experiment data file into its various
#  constituencies for testing purposes. The purpose is to iterate
#  the calling of the allele-translation-dmet3 program over all
#  experiment genes split by this program.
#  This is required because invalid data within the Genotype short report
#  experiment input files can cause the allele-translation-dmet3 program
#  to exit. By iterating each experiment gene then we are guaranteed that
#  all experiment genes can be tested independently, especially from those
#  invalid data experiments.
#
#
###############################################################################
#------------------------------------------------------------------------------
# INCLUDES
////////////////////////////////////////////////////////////////
//
// Copyright (C) 2008 Affymetrix, Inc.
//
// This program is free software; you can redistribute it and/or modify 
// it under the terms of the GNU General Public License (version 2) as 
// published by the Free Software Foundation.
// 
// This program is distributed in the hope that it will be useful, 
// but WITHOUT ANY WARRANTY; without even the implied warranty of 
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
// General Public License for more details.
// 
// You should have received a copy of the GNU General Public License 
// along with this program;if not, write to the 
// 
// Free Software Foundation, Inc., 
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////
use English qw( -no_match_vars );
use Getopt::Long;
use Pod::Usage;
use IO::File;
use Data::Dumper;
use Cwd ();
use warnings;
use strict;

#------------------------------------------------------------------------------
# CONSTANTS

#------------------------------------------------------------------------------
# GLOBAL VARIABLES

my ( $use_output, $use_verbose );

my $outDir      = ".";
my $genoFile    = "";
my $ttFile      = "";
my %ttFileGenes = ();

#------------------------------------------------------------------------------
# ROUTINES

#-----------------------------------------------------------------------------
#
#  initialize
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Sets up the environment for running.
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Parameters:
#
#  Returns:
#    exit value
#      - 0 = OK, 1 = user input error.
#
#  Environment:
#    croaks
#      - system errors
#
#------------------------------------------------------------------------------
sub initialize {

    if ( @ARGV != 2 ) {
        pod2usage();
    }

    if ( !-f $ARGV[0] ) {
        croak "File not found: $ARGV[0]\n";
    }

    $genoFile = $ARGV[0];

    if ( !-r $ARGV[0] ) {
        croak "Permission denied: $genoFile\n";
    }

    if ( !-f $ARGV[1] ) {
        croak "File not found: $ARGV[1]\n";
    }

    $ttFile = $ARGV[1];

    if ( !-r $ARGV[1] ) {
        croak "Permission denied: $ttFile\n";
    }

    if ($use_output) {
        if ( !-d $use_output ) {
            croak "Directory not found: $use_output\n";
        }
        $outDir = $use_output;
        if ( ( !-w $outDir ) || ( !-r $outDir ) ) {
            croak "Permission denied: $outDir\n";
        }
    }

    # Suck the translation table file into memory.

    my $tt_fh = IO::File->new( $ttFile, "<" );

    if ( !$tt_fh ) {
        croak "Can't open $ttFile, $OS_ERROR\n";
    }

TTFILE_LINE:
    foreach my $line (<$tt_fh>) {
        if ( $line =~ m{ ^(\w+)\t }xms ) {
            my $gene = $1;

            if ( !( $gene =~ m{Gene}xms ) ) {
                push @{ $ttFileGenes{$gene} }, $line;
            }
        }
    }

    close $tt_fh;

    return;

}

# end intialize
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  geno_data_split
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  For each line of the input file, look for qualifying lines and
#  when found output the records to an experiment file.
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Parameters:
#
#  Returns:
#    exit value
#      - 0 = OK, 1 = user input error.
#
#  Environment:
#    croaks
#      - system errors
#
#------------------------------------------------------------------------------
sub geno_data_split() {

    # Experiment records may not be in sorted order so we set up a hash
    # of all records.
    my %experiment_gene_lines    = ();
    my %experiment_gene_dir_name = ();
    my %experiment_short_name    = ();
    my %experiment_names         = ();
    my %skipped_genes            = ();

    my $geno_fh = IO::File->new( $genoFile, "<" );

    # READ records
    if ( !$geno_fh ) {
        croak "Can't open file $genoFile; $OS_ERROR\n";
    }

    my $header_line = <$geno_fh>;

    my $total_lines          = 1;
    my $skipped_lines        = 1;
    my $num_experiment_genes = 0;
READ_LINE:
    foreach my $line (<$geno_fh>) {

        $total_lines++;

        if ( $line =~ m{ ^([^\t]+\t)([^\t]+\t)([^\t]+\t)([^\t]+\t){3} }xms ) {

            my $experiment_name = $2;
            chop $experiment_name;

            if ( !defined( $experiment_names{$experiment_name} ) ) {
                $experiment_names{$experiment_name} = 1;
            }

            my $gene = $3;
            chop $gene;

            if ( !defined( $ttFileGenes{$gene} ) ) {
                if ( !defined( $skipped_genes{$gene} ) ) {
                    $skipped_genes{$gene} = 1;
                }
                else {
                    $skipped_genes{$gene}++;
                }
                $skipped_lines++;
                next READ_LINE;
            }

            my $experiment_gene = "$experiment_name:$gene";

            if ( !defined( $experiment_gene_lines{$experiment_gene} ) ) {

                if ( !defined( $experiment_short_name{$experiment_name} ) ) {
                    my $num_chars = 7;
                    my $experiment_short_name = substr $experiment_name, 0,
                        $num_chars;

                ADD_CHAR_TO_SHORT_NAME:
                    while (
                        defined $experiment_short_name{$experiment_short_name}
                        )
                    {
                        my $delete_name = $experiment_short_name;
                        $num_chars++;
                        $experiment_short_name = substr $experiment_name, 0,
                            $num_chars;
                        $experiment_short_name
                            =~ s{ [\(\)\!\.\#\?\/\\] }{+}gxms;
                        my $swap = $experiment_short_name{$delete_name};
                        $experiment_short_name{$swap} = substr $swap, 0,
                            $num_chars;
                        $experiment_short_name{$swap}
                            =~ s{ [\(\)\!\.\#\?\/\\] }{+}gxms;
                        delete $experiment_short_name{$delete_name};
                        $experiment_short_name{ $experiment_short_name{$swap}
                            } = $swap;
                    }
                    $experiment_short_name{$experiment_name}
                        = $experiment_short_name;
                    $experiment_short_name{$experiment_short_name}
                        = $experiment_name;
                }
                $experiment_gene_dir_name{$experiment_gene}
                    = $experiment_short_name{$experiment_name};
                push @{ $experiment_gene_lines{$experiment_gene} }, $line;
                $num_experiment_genes++;
            }
            else {
                push @{ $experiment_gene_lines{$experiment_gene} }, $line;
            }
        }
        else {
            $skipped_lines++;
        }

    }    # foreach line

    close $geno_fh;

EXPERIMENT_SHORT_NAME_PRINT:
    foreach my $experiment ( sort( keys(%experiment_names) ) ) {
        print "$experiment -> " . $experiment_short_name{$experiment} . "\n";
    }

    # WRITE records

    my $total_experiment_lines = 0;

EXPERIMENT_GENE_WRITE:
    foreach my $experiment_gene ( sort keys(%experiment_gene_lines) ) {

        my $num_lines = @{ $experiment_gene_lines{$experiment_gene} };

        $total_experiment_lines += $num_lines;

        print "$experiment_gene: $num_lines\n" if ($use_verbose);
        my $egDir = "$outDir/$experiment_gene_dir_name{$experiment_gene}";
        if ( !-d $egDir ) {
            if ( system("mkdir $egDir ") != 0 ) {
                croak "$egDir unable to create: $OS_ERROR\n";
            }
        }

        my $gene;
        if ( $experiment_gene =~ m{ :([^:]+)$ }xms ) {
            $gene = $1;
        }
        else {
            croak "Programmer error!\n";
        }

        my $gene_fh = IO::File->new( "$egDir/$gene", ">" );

        if ( !$gene_fh ) {
            croak "$gene can't create file: $OS_ERROR\n";
        }
        print {$gene_fh} $header_line;

    EXPERIMENT_LINE_WRITE:
        foreach my $line ( @{ $experiment_gene_lines{$experiment_gene} } ) {
            print {$gene_fh} $line;
        }
        close $gene_fh;

    }

    print
        "$total_lines lines : $total_experiment_lines experiment lines: $skipped_lines lines skipped.\n";
    print "$num_experiment_genes experiment genes with files created.\n";

    my $num_skipped_genes = keys(%skipped_genes);

    print "$num_skipped_genes: $num_skipped_genes\n";

SKIPPD_GENE:
    foreach my $skipped_gene ( sort keys(%skipped_genes) ) {
        print "$skipped_gene: $skipped_genes{$skipped_gene} skipped lines.\n";
    }
    return;
}

# end geno_data_split
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  main
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Simple control loop. If its not simple it don't belong hair.
#
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Parameters:
#
#  Returns:
#
#  Environment:
#
#------------------------------------------------------------------------------
# sub main

my ( $time_stamp, $exit_value );

GetOptions(
    'output=s' => \$use_output,
    'help'     => sub {
        pod2usage(
            -verbose => 2,
            -section => "SYNOPSIS|DESCRIPTION|CONFIGURATION",
            -exitval => 2
        );

    },
    'verbose' => \$use_verbose,
) or pod2usage(2);

initialize();

geno_data_split();

exit(0);

#------------------------------------------------------------------------------
#                               EOF
#------------------------------------------------------------------------------
1;
##############################################################################

=head1 NAME

geno_data_split.pl - split Genotype short report TSV file into various single experiment gene files.

=cut

=head1 SYNOPSIS

geno_data_split.pl [options] <Geno type short report TSV file> <Translation Table file>

Options:

    --output, -o     Output directory where all files will be written too.
                     NOTE: Each experiment will have its own subdirectory
                     created under the output directory. 
    --help, -h       Show this message

Example:

    geno_data_split.pl -o output 20080109_31set_DMET3_Genotypes_Short_qc_16various.txt DMET3_TTable_v20080110_EarlyAccess.txt

=head1 DESCRIPTION

C<geno_data_split.pl> is a command-line interface designed to facilitate
testing allele-translation-dmet3 by splitting the input Genotype Short Report
files into one file/experiment gene. In this way individual experiment genes
can be tested in addition to isolating experiments from each other. Any invalid
data in a Geno type short report input file will cause the
allele-translation-dmet3 program to exit. Call the allele-translation-dmet3
program on each individual experiment gene input file inside some shell script
"for loop" to avoid one bad experiment gene data set from invalidating all
experiment gene analysis. 

=head1 CONFIGURATION

None


=cut

##############################################################################
