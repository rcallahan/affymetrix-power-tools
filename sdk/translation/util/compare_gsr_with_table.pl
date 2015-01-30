#!/usr/bin/env perl
###############################################################################
#
#  compare_gsr_with_table.pl
#
#  Revision Control
#  ----------------
#  [$Id: compare_gsr_with_table.pl,v 1.2 2009-01-13 18:00:00 mspald Exp $]
#
#  Synopsis:
#  ----------------
#
#  Compare a Genotype Short Report with a translation table.
#  To do this requires two operations.
#
#  1.) Filter data not in the translation table out of the GSR.
#  2.) Order the AssayIDs in translation table order.
#
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

#------------------------------------------------------------------------------
# GLOBAL VARIABLES

my $tt_file  = "";
my $gsr_file = "";

my %assayId_order = ();

my @gsr_lines       = ();
my $gsr_header_line = "";

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

    if ( !-r $ARGV[0] ) {
        croak "Permission denied: $ARGV[0]\n";
    }

    if ( !-f $ARGV[1] ) {
        croak "File not found: $ARGV[1]\n";
    }

    if ( !-r $ARGV[1] ) {
        croak "Permission denied: $ARGV[1]\n";
    }

    $tt_file  = $ARGV[0];
    $gsr_file = $ARGV[1];

    my $tt_fh = IO::File->new( $tt_file, "<" );
    if ( !$tt_fh ) {
        croak "$tt_file: can't open, $OS_ERROR\n";
    }

    my $header_line = <$tt_fh>;

    my @splits = split( m{\t}xms, $header_line );

    if ( @splits < 15 ) {
        croak "$tt_file: invalid Translation Table file.\n";
    }
    if ( $splits[2] ne "AssayID" ) {
        croak "$tt_file: invalid Translation Table file.\n";
    }

    my $order = 1;

TT_LINE:
    foreach my $line (<$tt_fh>) {

        @splits = split( m{\t}xms, $line );
        if ( @splits < 15 ) {
            croak "$tt_file: invalid Translation Table file.\n";
        }

        if ( !( $splits[2] =~ m{^\s*\#}xms ) ) {
            $assayId_order{ $splits[2] } = $order;
            $order++;
        }

    }
    close $tt_fh;

    my $gsr_fh = IO::File->new( $gsr_file, "<" );
    if ( !$gsr_fh ) {
        croak "$gsr_file: can't open file, $OS_ERROR\n";
    }

    $gsr_header_line = <$gsr_fh>;
    @splits = split( m{\t}xms, $gsr_header_line );

    if ( @splits != 7 ) {
        croak "$gsr_file: invalid Genotype Short Report file.\n";
    }

    if (   ( $splits[0] ne "Sample Name" )
        || ( $splits[1] ne "Experiment Name" )
        || ( $splits[2] ne "Gene" )
        || ( $splits[4] ne "Assay Id" ) )
    {
        print "$gsr_header_line";
        croak "$gsr_file: invalid Genotype Short Report file.\n;";
    }

    my $count = 1;
GSR_LINE:
    foreach my $line (<$gsr_fh>) {
        @splits = split m{\t}xms, $line;
        my %gsr_record = ();

        $gsr_record{'EXPERIMENT'} = $splits[1];
        $gsr_record{'GENE'}       = $splits[2];
        $gsr_record{'ASSAYID'}    = $splits[4];
        $gsr_record{'LINE'}       = $line;

        push @gsr_lines, {%gsr_record};
    }

    close $gsr_fh;

    return;

}

# end intialize
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  compare_gsr_with_table
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Sort the input %gsr_lines by experiment name, gene, and assay Id where
#  the assay Id is ordered from the input translation table.
#
#  If the assay Id is not found in the translation table then the line is
#  filtered.
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Returns:
#    Nothing
#
#  Environment:
#    Outputs the comparison to console.
#
#
#------------------------------------------------------------------------------
sub compare_gsr_with_table( ) {

    # FILTER for only those Assay IDs found int the translation table.
    my @gsr_filtered_lines = ();

GSR_RECORD:
    foreach my $gsr_record_ref (@gsr_lines) {
        if ( defined( $assayId_order{ $gsr_record_ref->{'ASSAYID'} } ) ) {
            push @gsr_filtered_lines, $gsr_record_ref;
        }
    }

    # SORT

    my @gsr_sorted_lines = sort {
        my $a_experiment    = $a->{'EXPERIMENT'};
        my $a_gene          = $a->{'GENE'};
        my $a_assayId_order = $assayId_order{ $a->{'ASSAYID'} };

        my $b_experiment    = $b->{'EXPERIMENT'};
        my $b_gene          = $b->{'GENE'};
        my $b_assayId_order = $assayId_order{ $b->{'ASSAYID'} };

        if ( $a_experiment ne $b_experiment ) {
            return ( $a_experiment cmp $b_experiment );
        }

        if ( $a_gene ne $b_gene ) {
            return ( $a_gene cmp $b_gene );
        }

        return ( $a_assayId_order <=> $b_assayId_order );

    } @gsr_filtered_lines;

    my $gsr_compare_file = $gsr_file . ".compare_gsr_with_table";
    my $gsr_compare_fh = IO::File->new( $gsr_compare_file, ">");

    if ( ! $gsr_compare_fh ) {
        croak "$gsr_compare_file: can't open file for write, $OS_ERROR\n";
    }

    
    print $gsr_header_line;
    print {$gsr_compare_fh} $gsr_header_line;
    
    foreach my $gsr_sorted_lines (@gsr_sorted_lines) {

        print $gsr_sorted_lines->{'LINE'};
        print {$gsr_compare_fh} $gsr_sorted_lines->{'LINE'};
    }

    close $gsr_compare_fh;
    
    return;
    

}

# end compare_gsr_with_trable.
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
    'help' => sub {
        pod2usage(
            -verbose => 2,
            -section => "SYNOPSIS|DESCRIPTION|CONFIGURATION",
            -exitval => 2
        );

    },
) or pod2usage(2);

initialize();

compare_gsr_with_table();

exit(0);

#------------------------------------------------------------------------------
#                               EOF
#------------------------------------------------------------------------------
1;
##############################################################################

=head1 NAME

compare_gsr_with_table.pl - DMET3 regression testing for allele translation.

=cut

=head1 SYNOPSIS

check_dmet3.pl [options] <translation table file> <Genotype Short Report file> 

Options:
    --help           Show this message
Example:

    compare_gsr_with_table.pl DMET3_TTable_v20080110_EarlyAccess.txt  CYP2D6

=head1 DESCRIPTION

C<compare_gsr_with_table.pl> is a command-line interface designed to 
compare DMET3 Genotype Short Report TSV files with the allele translation
table file. 

This tool is NOT used by regression itself, but instead is designed
to aid humans evaluate data by giving a visually easy alignment of
Assay ID comparison for determining if a call in a report is correct. 

=head1 CONFIGURATION

None


=cut

##############################################################################
