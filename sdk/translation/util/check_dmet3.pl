#!/usr/bin/env perl
###############################################################################
#
#  check_dmet3.pl
#
#  Revision Control
#  ----------------
#  [$Id: check_dmet3.pl,v 1.5 2009-01-13 18:00:00 mspald Exp $]
#
#  Synopsis:
#  ----------------
#
#  At Affy, the "make check" target runs a bunch of unit test cases.
#  This script here runs regression, end-to-end testing.
#
#
###############################################################################
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
#------------------------------------------------------------------------------
# INCLUDES
use Carp;
use English qw( -no_match_vars );
use Env;
use Getopt::Long;
use Pod::Usage;
use IO::File;
use Data::Dumper;
use Cwd ();
use warnings;
use strict;

#------------------------------------------------------------------------------
# CONSTANTS

my $A3BIN = "bin/allele-translation-dmet3foo";

my $OUTDIR = "/tmp";

my $DMET2_MARKER_REPORT_EXT    = "dmet2_marker.reg";
my $DMET2_HAPLOTYPE_REPORT_EXT = "dmet2_haplotype.reg";
my $DMET3_MARKER_REPORT_EXT    = "dmet3_marker.reg";
my $DMET3_HAPLOTYPE_REPORT_EXT = "dmet3_haplotype.reg";

my %HEADERS = (
    "Experiment" => 0,
    "Gene"       => 1,
    "AssayId"    => 5,
    "Call"       => 6,
    "Ref"        => 7,
    "Var"        => 8,
    "Allele"     => 9,
);

my @HEADER_ORDER
    = ( "Experiment", "Gene", "AssayId", "Call", "Ref", "Var", "Allele" );

#------------------------------------------------------------------------------
# GLOBAL VARIABLES

my ( $use_verbose, $use_first, $use_marker_only, $use_haplotype_only );
my ($use_binary, $use_data_dir);

my $genotype_root_dir = "";

my @genotype_files = ();

my $a3_binary = "";
my $a3_data_dir = "";

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

    if ( $use_marker_only && $use_haplotype_only ) {
        croak "Can't request marker only with haplotype only.\n";
    }

    if ( @ARGV != 1 ) {
        pod2usage();
    }

    if ( !-d $ARGV[0] ) {
        croak "Directory not found: $ARGV[0]\n";
    }

    if ( !-r $ARGV[0] ) {
        croak "Permission denied: $ARGV[0]\n";
    }

    $genotype_root_dir = $ARGV[0];

    if ( -f $A3BIN ) {
        $a3_binary = $A3BIN;
    }
    else {
        if ( !$use_binary ) {
            print STDERR
                "DMET3 binary \"allele-translation-dmet3\" not specified.\n\n";
            pod2usage();
        }
        elsif ( !-f $use_binary ) {
            croak "$use_binary: file not found.\n\n";
        }
        elsif ( !-e $use_binary ) {
            croak "$use_binary: permission denied.\n\n";
        }
        $a3_binary = $use_binary;
    }

    if ( !$use_data_dir ) {
        print {\*STDERR} "Data directory required and missing for Translation Table and Copy Number files.\n";
        pod2usage();
    }
    elsif ( ! -d $use_data_dir ) {
        croak "$use_data_dir: directory not found.\n";
    }
    $a3_data_dir = $use_data_dir;

    my $pipe_fh = IO::File->new("find $genotype_root_dir -type f|");
    if ( !$pipe_fh ) {
        croak "Can't open pipe: $OS_ERROR.\n";
    }

            
READFILE_LABEL:
    foreach my $geno_file (<$pipe_fh>) {

        chomp $geno_file;

        if ( $geno_file =~ m{ \.(rpt) }xms ) {
            next READFILE_LABEL;
        }

        my $line = `head -1 $geno_file`;

        #print "$geno_file: $line";
        chomp $line if ($line);

        if ( $line && ( $line =~ m{^Sample\sName}xms ) ) {

            if (   ( -f "${geno_file}.$DMET2_MARKER_REPORT_EXT" )
                && ( -s "${geno_file}.$DMET2_MARKER_REPORT_EXT" ) )
            {
                push @genotype_files, $geno_file;
            }
            else {
                print "$geno_file: no marker report file found\n";
            }
        }

    }

    close($pipe_fh);

    if (@genotype_files) {
        print "Found " . ( @genotype_files + 0 ) . " genotype data files.\n";
    }
    else {
        croak "No genotype files found!\n";
    }

    return;

}

# end intialize
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  run_a3_on_genotype_files
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Given the list of known genotype files, call the a3 program and
#  create the dmet2 comparison output log.
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Parameters:
#
#  @genotype_files: initialized earlier.
#
#  Returns:
#    Nothing
#
#  Environment:
#    croaks
#      - system errors
#
#------------------------------------------------------------------------------
sub run_a3_on_genotype_files( ) {

    print "run_a3_on_genotype_files: BEGIN\n";

INPUT_GENOTYPE_FILE_LABEL:
    foreach my $genotype_file (@genotype_files) {

        my $genotype_dir = $genotype_file;
        if ( $genotype_dir =~ m{ \/ }xms ) {
            $genotype_dir =~ s{/[^/]+$}{}xms;
        }

        my $cmd = "$a3_binary -r 1 -i $a3_data_dir -g $genotype_file";

        $cmd = "bash --norc --noprofile -f -c \'$cmd\' >>/dev/null 2>&1";

        #print "$cmd\n";

        if ( system($cmd) == 0 ) {

            my $marker_regression_log_file    = `ls *marker.reg`;
            my $haplotype_regression_log_file = `ls *haplotype.reg`;
            chomp $marker_regression_log_file;
            chomp $haplotype_regression_log_file;

            $marker_regression_log_file    =~ s{[^\/]+\/}{}gxms;
            $marker_regression_log_file    =~ s{^\/}{}xms;
            $haplotype_regression_log_file =~ s{[^\/]+\/}{}gxms;
            $haplotype_regression_log_file =~ s{^\/}{}xms;

            if ( !$use_haplotype_only ) {
                if ( !$marker_regression_log_file ) {
                    croak
                        "Can't find regression marker report for file: ${genotype_file}\n";
                }
                if ( !-f $marker_regression_log_file ) {
                    croak "Can't find file $marker_regression_log_file.\n";
                }
                $cmd
                    = "/bin/cp -f $marker_regression_log_file ${genotype_file}.$DMET3_MARKER_REPORT_EXT";
                if ( system($cmd) != 0 ) {
                    croak "$cmd, $OS_ERROR\n";
                }
            }
            if ( !$use_marker_only ) {
                if ( !$haplotype_regression_log_file ) {
                    croak
                        "Can't find regression haplotype report for file: ${genotype_file}\n";
                }

                if ( !-f $haplotype_regression_log_file ) {
                    croak "Can't find file $haplotype_regression_log_file.\n";
                }
                $cmd
                    = "sort $haplotype_regression_log_file > ${genotype_file}.$DMET3_HAPLOTYPE_REPORT_EXT";
                if ( system($cmd) != 0 ) {
                    croak "$cmd, $OS_ERROR\n";
                }
            }

            #print "$cmd\n";

#$cmd = "sort $regression_log_file >& ${genotype_file}.$DMET3_MARKER_REPORT_EXT";
        }
        elsif ($use_first) {
            croak "Exiting on first failure: $genotype_file\n";
        }

        print ".";

    }
    print "\n";

    print "run_a3_on_genotype_files: END\n";
}

# end run_a3_on_genotype_files
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  analyze_dmet2_with_dmet3_results
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Run diff on the two outputs and report the differences.
#
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
sub analyze_dmet2_with_dmet3_results() {

    my %results = ();

    my $total_tests             = 0;
    my $total_marker_passed     = 0;
    my $total_marker_failed     = 0;
    my $total_marker_aborted    = 0;
    my $total_haplotype_passed  = 0;
    my $total_haplotype_failed  = 0;
    my $total_haplotype_aborted = 0;

ANALYZE_FILES_LABEL:
    foreach my $genotype_file (@genotype_files) {

        $total_tests++;

        if ( !$use_haplotype_only ) {

            if ( !-f "${genotype_file}.$DMET3_MARKER_REPORT_EXT" ) {
                $total_marker_aborted++;
                $total_haplotype_aborted++;
                if ($use_first) {
                    croak "Exiting on first failure: $genotype_file.\n";
                }
                next ANALYZE_FILES_LABEL;
            }

            my $cmd
                = "diff  ${genotype_file}.$DMET2_MARKER_REPORT_EXT ${genotype_file}.$DMET2_MARKER_REPORT_EXT >& ${genotype_file}.marker.diff";
            my $temp;
            system($cmd);

            if ( !-f "${genotype_file}.marker.diff" ) {
                croak "FAILED: $cmd\n";
            }

            if ( -z "${genotype_file}.marker.diff" ) {
                $results{$genotype_file} = 1;
                $total_marker_passed++;
            }
            else {
                $total_marker_failed++;
                $results{$genotype_file} = 0;
                if ($use_first) {
                    croak "Existing on first failure: $genotype_file.\n";
                }
            }

        }

        if ( !$use_marker_only ) {

            if ( !-f "${genotype_file}.$DMET3_HAPLOTYPE_REPORT_EXT" ) {
                $total_haplotype_aborted++;
                if ($use_first) {
                    croak "Exiting on first failure: $genotype_file.\n";
                }
            }

            my $cmd
                = "diff  ${genotype_file}.$DMET2_HAPLOTYPE_REPORT_EXT ${genotype_file}.$DMET3_HAPLOTYPE_REPORT_EXT >& ${genotype_file}.haplotype.diff";
            my $temp;
            system($cmd);

            if ( !-f "${genotype_file}.haplotype.diff" ) {
                croak "FAILED: $cmd\n";
            }

            if ( -z "${genotype_file}.haplotype.diff" ) {
                $results{$genotype_file} = 1;
                $total_haplotype_passed++;
            }
            else {
                $total_haplotype_failed++;
                $results{$genotype_file} = 0;
                if ($use_first) {
                    croak "Existing on first failure: $genotype_file.\n";
                }
            }

        }

    }    # foreach test file.

    if ( !$use_haplotype_only ) {
        print
            "MARKER Summary:    $total_tests tests | $total_marker_passed passed | $total_marker_failed failed | $total_marker_aborted aborted\n";
    }
    if ( !$use_marker_only ) {
        print
            "HAPLOTYPE Summary: $total_tests tests | $total_haplotype_passed passed | $total_haplotype_failed failed | $total_haplotype_aborted aborted\n";
    }

    return;

}

# end analyze_dmet2_with_dmet3_results
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
    'binary=s'  => \$use_binary,
    'data=s'    => \$use_data_dir,
    'first'     => \$use_first,
    'haplotype' => \$use_haplotype_only,
    'help'      => sub {
        pod2usage(
            -verbose => 2,
            -section => "SYNOPSIS|DESCRIPTION|CONFIGURATION",
            -exitval => 2
        );

    },
    'marker' => \$use_marker_only,

    'verbose' => \$use_verbose,
) or pod2usage(2);

initialize();

run_a3_on_genotype_files();

analyze_dmet2_with_dmet3_results();

exit(0);

#------------------------------------------------------------------------------
#                               EOF
#------------------------------------------------------------------------------
1;
##############################################################################

=head1 NAME

check_dmet3.pl - DMET3 regression testing for allele translation markers 

=cut

=head1 SYNOPSIS

check_dmet3.pl [options] <directory of Genotype Short Report files> 

Options:
    --binary, -b     REQUIRED: allele-translation-dmet3 binary
    --data,  -d      REQUIRED: DMET3 data directory with translation table, copy number report file. 
    --first, -f      First to fail. Stop execution on the first test to fail.
    --haplotype      Haplotype reports only. 
    --help           Show this message
    --marker, -m     Marker reports only. 
    --verbose, -v    Output the results to console as well as to file. 
Example:

    check_marker_dmet3.pl -b bin/allele-translation-dmet3 testdata_dmet2b/

=head1 DESCRIPTION

C<check_marker_dmet3.pl> is a command-line interface designed to 
compare DMET3 results with known DMET2 results as a form of regression. 

Every Genotype Short Report file found will have DMET2 output files
of the same name but with an extension ".dmet2_XXXX.reg". In addition,
the DMET3 report will be the same name with an extension ".reg".

Regression is comprised of running "diff" on two files and reporting
the results. 


=head1 CONFIGURATION

None


=cut

##############################################################################
