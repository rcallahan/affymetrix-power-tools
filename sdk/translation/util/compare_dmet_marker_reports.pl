#!/usr/bin/env perl
###############################################################################
#
#  create_dmet_marker_reports.pl
#
#  Revision Control
#  ----------------
#  [$Id: compare_dmet_marker_reports.pl,v 1.5 2009-01-13 18:00:00 mspald Exp $]
#
#  Synopsis:
#  ----------------
#
#  HACK: DMET2 marker reports are not to included as part of the regression.
#  Instead, filtered, normalized versions are to be kept. This script
#  is intended soley for the purposes for QA of the initial release of
#  DMET3.
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

my $OUTDIR = "/tmp";

my %MARKER_HEADERS = (
    "Experiment",            => 0,
    "Gene"                   => 1,
    "Sample"                 => 2,
    "Functional Change"      => 3,
    "External ID"            => 4,
    "Assay ID"               => 5,
    "Basecall"               => 6,
    "Reference Base"         => 7,
    "Variant Base"           => 8,
    "Call"                   => 9,
    "Haplotype Marker"       => 10,
    "Change for Variant"     => 11,
    "Variant cDNA Change"    => 12,
    "Variant DNA Change"     => 13,
    "dbSNP ID"               => 14,
    "Validated"              => 15,
    "Allele Defining Marker" => 16,
    "Relevant Alleles"       => 17,
);

#Experiment	Gene	Sample	Functional_Change	External_ID	Assay_ID	Basecall	Reference_Base	Variant_Base	Call	Haplotype_Marker	Change_for_Variant	Variant_cDNA_Change	Variant_DNA_Change	dbSNP_ID	Validated	Allele_Defining_Marker	Relevant_Alleles

#------------------------------------------------------------------------------
# GLOBAL VARIABLES

my ( $use_known_discrepancies, $use_count, $use_summary, $use_lines );

my $dmet2_file = "";
my $dmet3_file = "";

my %dmet2              = ();
my @dmet2_marker_lines = ();

my %dmet3              = ();
my @dmet3_marker_lines = ();

my %missing_experiments          = ();
my %missing_experiment_assay_ids = ();
my %diffs                        = ();
my %diff_stats                   = ();

my %known_discrepancies = ( 'QUOTES' => 0, );

my %known_discrepancy_description
    = ( 'QUOTES' => "Variant DNA Change quoted DMET2",  );


my %dmet2_stats = ( 'TOTAL_MARKERS' => 0,
                    'INVALID_LINES' => 0, );
my %dmet3_stats = ( 'TOTAL_MARKERS' => 0,
                    'INVALID_LINES' => 0, );
my %dmet3_duplicate_markers = ();

#------------------------------------------------------------------------------
# ROUTINES

#-----------------------------------------------------------------------------
#
#  resolve_known_discrepancies
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Some differences are trivial, such as DMET2 quoting things occasionally
#  that are not quoted in the original data set (Translation Table).
#
#  This routine calculates such things and records them.
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
sub resolve_known_discrepancies {
    my ( $col, $dmet2_ref, $dmet3_ref ) = @ARG;

    return
        if ( $dmet2_ref->[ $MARKER_HEADERS{$col} ] eq
        $dmet3_ref->[ $MARKER_HEADERS{$col} ] );

    if ( $col eq "Variant DNA Change" ) {
        my $quote_test = '"' . $dmet3_ref->[ $MARKER_HEADERS{$col} ] . '"';
        if ( $dmet2_ref->[ $MARKER_HEADERS{$col} ] eq $quote_test ) {
            $known_discrepancies{'QUOTES'}++;
            $dmet3_ref->[ $MARKER_HEADERS{$col} ] = $quote_test;
        }
    }

}

# end resolve_known_discrepancies
#-----------------------------------------------------------------------------
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
        croak "$ARGV[0]: file not found.\n";
    }

    if ( !-f $ARGV[1] ) {
        croak "$ARGV[1]: file not found.\n";
    }

    if ( !( $ARGV[0] =~ m{ \.marker\.rpt }xms ) ) {
        croak "$ARGV[0]: DMET2 files end with '.marker.rpt'\n";
    }

    if ( !( $ARGV[1] =~ m{ dmet3_3\.0\.marker\.rpt$ }xms ) ) {
        croak "$ARGV[0]: DMET3 files end with 'dmet_3.0..marker.rpt'\n";
    }

    $dmet2_file = $ARGV[0];
    $dmet3_file = $ARGV[1];

    foreach my $col ( keys %MARKER_HEADERS ) {
        $diff_stats{$col} = 0;
    }

    return;

}

# end intialize
#-----------------------------------------------------------------------------

#-----------------------------------------------------------------------------
#
#  read_marker_files
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Read the marker files completely into memory has hashes keyed by
#  experiment : assay id.
#
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#
#------------------------------------------------------------------------------
sub read_marker_files() {

    my $num_fields = scalar keys(%MARKER_HEADERS);

    my $dmet2_fh = IO::File->new( $dmet2_file, "<" );

    if ( !$dmet2_fh ) {
        croak "${dmet2_file}: $OS_ERROR\n";
    }

    my $line_num = 0;

DMET2_LINE:
    foreach my $line (<$dmet2_fh>) {
        push @dmet2_marker_lines, $line;
        $line_num++;
        if ( $line =~ m{^\s*\#}xms ) {
            next DMET2_LINE;
        }
        $dmet2_stats{'TOTAL_MARKERS'}++;
        while ( chomp $line ) { }
        $line =~ s{\t$}{}xms;
        my @splits = split m{\t}xms, $line;

        if ( $num_fields != ( scalar @splits ) ) {
            if ( ( $num_fields - 1 ) == ( scalar @splits ) ) {
                push @splits, "";
            }
            else {
                print "$line\n";
                print "Only had "
                    . scalar @splits
                    . " of $num_fields fields\n";
                $dmet2_stats{'INVALID_LINES'}++;
                next DMET2_LINE;
            }
        }

        $dmet2{ $splits[ $MARKER_HEADERS{"Experiment"} ] }
            { $splits[ $MARKER_HEADERS{"Assay ID"} ] }{"cols"} = [@splits];
        $dmet2{ $splits[ $MARKER_HEADERS{"Experiment"} ] }
            { $splits[ $MARKER_HEADERS{"Assay ID"} ] }{"line"} = $line_num;

    }

    close $dmet2_fh;
    $dmet2_fh = 0;

    my $dmet3_fh = IO::File->new( $dmet3_file, "<" );

    if ( !$dmet3_fh ) {
        croak "${dmet3_file}: $OS_ERROR\n";
    }

    $line_num = 0;
DMET3_LINE:
    foreach my $line (<$dmet3_fh>) {
        push @dmet3_marker_lines, $line;
        $line_num++;
        if ( $line =~ m{^\s*\#}xms ) {
            next DMET3_LINE;
        }
        $dmet3_stats{'TOTAL_MARKERS'}++;
        while ( chomp $line ) { }

        my @splits = split m{\t}xms, $line;

        if ( $num_fields != ( scalar @splits ) ) {
            if ( ( $num_fields - 1 ) == ( scalar @splits ) ) {
                push @splits, "";
            }
            else {
                print "$line\n";
                print "Only had "
                    . scalar @splits
                    . " of $num_fields fields\n";
                $dmet3_stats{'INVALID_LINES'}++;
                next DMET3_LINE;
            }
        }

        if ( defined $dmet3{ $splits[ $MARKER_HEADERS{"Experiment"} ] }{$splits[ $MARKER_HEADERS{"Assay ID"} ] }) {
            if ( !defined $dmet3_duplicate_markers{$splits[ $MARKER_HEADERS{"Assay ID"}]} ) {
                $dmet3_duplicate_markers{$splits[ $MARKER_HEADERS{"Assay ID"}]} = 1;
            }
            else {
                $dmet3_duplicate_markers{$splits[ $MARKER_HEADERS{"Assay ID"}]}++;
            }
        }

        $dmet3{ $splits[ $MARKER_HEADERS{"Experiment"} ] }
            { $splits[ $MARKER_HEADERS{"Assay ID"} ] }{"cols"} = [@splits];
        $dmet3{ $splits[ $MARKER_HEADERS{"Experiment"} ] }
            { $splits[ $MARKER_HEADERS{"Assay ID"} ] }{"line"} = $line_num;

    }
    close $dmet3_fh;
    $dmet3_fh = 0;

    return;
}

# end read_marker_files
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  compare_files
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#
#
#
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#
#------------------------------------------------------------------------------
sub compare_files {

    # DMET2 is the authority. Missing data is missing relative to DMET2
    # as well as extra being extra.

DMET2_EXPERIMENT:
    foreach my $dmet2_experiment ( sort( keys %dmet2 ) ) {

        if ( !defined $dmet3{$dmet2_experiment} ) {
            $missing_experiments{$dmet2_experiment} = 1;
            next DMET2_EXPERIMENT;
        }

    DMET2_EXPERIMENT_ASSAYID:
        foreach
            my $dmet2_assay_id ( sort( keys %{ $dmet2{$dmet2_experiment} } ) )
        {

            if ( !defined $dmet3{$dmet2_experiment}{$dmet2_assay_id} ) {
                $missing_experiment_assay_ids{$dmet2_experiment}
                    {$dmet2_assay_id}
                    = $dmet2{$dmet2_experiment}{$dmet2_assay_id}{"line"};
                next DMET2_EXPERIMENT_ASSAYID;
            }

            my @dmet2
                = @{ $dmet2{$dmet2_experiment}{$dmet2_assay_id}{"cols"} };
            my $dmet2_line
                = $dmet2{$dmet2_experiment}{$dmet2_assay_id}{"line"};
            my @dmet3
                = @{ $dmet3{$dmet2_experiment}{$dmet2_assay_id}{"cols"} };
            my $dmet3_line
                = $dmet3{$dmet2_experiment}{$dmet2_assay_id}{"line"};

            my $okRecord = 1;

        COLUMN_TO_MATCH:
            foreach my $col (
                sort { return $MARKER_HEADERS{$a} <=> $MARKER_HEADERS{$b}; }
                keys(%MARKER_HEADERS) )
            {

                if ( ( $col eq "Basecall" ) || ( $col eq "Call" ) ) {
                    my $dmet2 = $dmet2[ $MARKER_HEADERS{$col} ];
                    my @splits = split( m{/}xms, $dmet2 );
                    if ( @splits > 1 ) {
                        @splits = sort @splits;
                        $dmet2  = $splits[0] . "/" . $splits[1];
                        $dmet2[ $MARKER_HEADERS{$col} ] = $dmet2;
                    }
                    my $dmet3 = $dmet3[ $MARKER_HEADERS{$col} ];
                    @splits = split( m{/}xms, $dmet3 );
                    if ( @splits > 1 ) {
                        @splits = sort @splits;
                        $dmet3  = $splits[0] . "/" . $splits[1];
                        $dmet3[ $MARKER_HEADERS{$col} ] = $dmet3;
                    }

                }
                elsif ( ( $col eq "Variant DNA Change" ) ) {

                    if ($use_known_discrepancies) {
                        resolve_known_discrepancies( $col, \@dmet2, \@dmet3 );
                    }
                    my $dmet2 = $dmet2[ $MARKER_HEADERS{$col} ];
                    my @splits = split( m{,}xms, $dmet2 );
                    if ( @splits > 1 ) {
                        @splits = sort @splits;
                        $dmet2  = $splits[0] . "," . $splits[1];
                        $dmet2[ $MARKER_HEADERS{$col} ] = $dmet2;
                    }
                    my $dmet3 = $dmet3[ $MARKER_HEADERS{$col} ];
                    @splits = split( m{,}xms, $dmet3 );
                    if ( @splits > 1 ) {
                        @splits = sort @splits;
                        $dmet3  = $splits[0] . "," . $splits[1];
                        $dmet3[ $MARKER_HEADERS{$col} ] = $dmet3;
                    }
                }
                elsif (($col eq "Variant Base" ) ||
                       ($col eq "Change for Variant") ||
                       ($col eq "Functional Change") ||
                       ($col eq "Variant DNA Change") ||
                       ($col eq "Variant cDNA Change") ||
                       ($col eq "Allele Defining Marker")){
                    my $dmet2_value = $dmet2[ $MARKER_HEADERS{$col} ];
                    my @splits = split( m{,}xms, $dmet2_value );
                    if ( @splits > 1 ) {
                        @splits = sort @splits;
                        $dmet2_value  = $splits[0] . "," . $splits[1];
                        $dmet2[ $MARKER_HEADERS{$col} ] = $dmet2_value;
                    }
                    $dmet2_value =~ s{,$}{}xms;
                    $dmet2[ $MARKER_HEADERS{$col} ] = $dmet2_value;
                    my $dmet3_value = $dmet3[ $MARKER_HEADERS{$col} ];
                    @splits = split( m{,}xms, $dmet3_value );
                    if ( @splits > 1 ) {
                        @splits = sort @splits;
                        $dmet3_value  = $splits[0] . "," . $splits[1];
                        $dmet3[ $MARKER_HEADERS{$col} ] = $dmet3_value;
                    }

                }
                elsif (( $col eq "Relevant Alleles" )
                    && ( defined $dmet2[ $MARKER_HEADERS{$col} ] ) )
                {
                    my @splits
                        = split( m{,}xms, $dmet2[ $MARKER_HEADERS{$col} ] );
                    if ( @splits > 1 ) {
                        my %sort_order = ();
                        foreach my $split (@splits) {
                            $sort_order{$split} = 1;
                            if ( ( $split =~ m{(\d+)}xms ) ) {
                                $sort_order{$split} = $1;
                            }
                        }
                        @splits = sort( {
                                if ( ( $sort_order{$a} != $sort_order{$b} ) )
                                {
                                    return $sort_order{$a} <=>
                                        $sort_order{$b};
                                }
                                return $a cmp $b;
                        } @splits );
                        my $newRelevantAlleles = "";
                        foreach my $split (@splits) {
                            $newRelevantAlleles .= $split . ",";
                        }
                        chop $newRelevantAlleles;
                        $dmet2[ $MARKER_HEADERS{$col} ] = $newRelevantAlleles;
                    }
                    @splits
                        = split( m{,}xms, $dmet3[ $MARKER_HEADERS{$col} ] );
                    if ( @splits > 1 ) {
                        my %sort_order = ();
                        foreach my $split (@splits) {
                            $sort_order{$split} = 1;
                            if ( ( $split =~ m{(\d+)}xms ) ) {
                                $sort_order{$split} = $1;
                            }
                        }
                        @splits = sort( {
                                if ( ( $sort_order{$a} != $sort_order{$b} ) )
                                {
                                    return $sort_order{$a} <=>
                                        $sort_order{$b};
                                }
                                return $a cmp $b;
                        } @splits );
                        my $newRelevantAlleles = "";
                        foreach my $split (@splits) {
                            $newRelevantAlleles .= $split . ",";
                        }
                        chop $newRelevantAlleles;
                        $dmet3[ $MARKER_HEADERS{$col} ] = $newRelevantAlleles;
                    }
                }

                if ( $dmet2[ $MARKER_HEADERS{$col} ] ne
                    $dmet3[ $MARKER_HEADERS{$col} ] )
                {
                    $okRecord = 0;
                    if ( !defined( $diffs{$dmet3_line} ) ) {
                        $diffs{$dmet3_line}
                            = "### DMET3 [${dmet3_line}]: DMET2 [${dmet2_line}]: ";
                    }
                    $diff_stats{$col} += 1;
                    $diffs{$dmet3_line} .= " $col 2=("
                        . $dmet2[ $MARKER_HEADERS{$col} ] . ") 3=("
                        . $dmet3[ $MARKER_HEADERS{$col} ] . ") |";

                }
            }
        }
    }

}

# end compare_files
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  report
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#
#
#
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#
#------------------------------------------------------------------------------
sub report {

    my $dmet2_missing_marker_count      = 0;
    my $dmet3_extra_marker_count        = 0;
    my %dmet3_extra_experiments = ();
    
    if ( !$use_count ) {
        print "----------------------------------------------------\n";
        print "$dmet2_file <=> $dmet3_file\n";
        print "----------------------------------------------------\n";
    }

DMET2_EXPERIMENT:
    foreach my $dmet2_experiment ( sort( keys %dmet2 ) ) {

        if ( $dmet2_experiment eq "Experiment" ) {
            next DMET2_EXPERIMENT;
        }

        if ( $missing_experiments{$dmet2_experiment} ) {
            if ( !$use_count && !$use_summary ) {
                print "${dmet2_experiment}: DMET3 IS MISSING EXPERIMENT \n";
            }
            next DMET2_EXPERIMENT;
        }

        if ( defined $missing_experiment_assay_ids{$dmet2_experiment} ) {
        MISSING_EXPERIMENT_ASSAYID:
            foreach my $missing_assay_id (
                sort(
                    keys %{ $missing_experiment_assay_ids{$dmet2_experiment} }
                ) )
            {
                $dmet2_missing_marker_count++;
                my $line_num
                    = $missing_experiment_assay_ids{$dmet2_experiment}
                    {$missing_assay_id};
                if ( !$use_count && !$use_summary ) {

                    print
                        "${dmet2_experiment}: $line_num: ${missing_assay_id}: DMET3 IS MISSING ASSAY ID\n";
                }
            }
        }
    }

DMET3_EXPERIMENT:
    foreach my $dmet3_experiment ( sort( keys %dmet3 ) ) {

        if ( $dmet3_experiment eq "Experiment" ) {
            next DMET3_EXPERIMENT;
        }
        if ( !$use_count && !$use_summary ) {
            print "DMET3: $dmet3_experiment\n";
        }

        if ( !defined $dmet2{$dmet3_experiment} ) {
            $dmet3_extra_experiments{$dmet3_experiment} =
                scalar keys(  %{ $dmet3{$dmet3_experiment}} );
            if ( !$use_count && !$use_summary ) {
                print "${dmet3_experiment}: DMET3 EXPERIMENT NOT IN DMET2\n";
            }
            next DMET3_EXPERIMENT;
        }

    DMET3_EXPERIMENT_ASSAYID:
        foreach
            my $dmet3_assay_id ( sort( keys %{ $dmet3{$dmet3_experiment} } ) )
        {

            my $dmet3_line
                = $dmet3{$dmet3_experiment}{$dmet3_assay_id}{"line"};
            if ( !defined $dmet2{$dmet3_experiment}{$dmet3_assay_id} ) {
                my $gene = $dmet3{$dmet3_experiment}{$dmet3_assay_id}{"cols"}
                    [ $MARKER_HEADERS{"Gene"} ];

                if ( !$use_count && !$use_summary ) {
                    print
                        "${dmet3_experiment}: $dmet3_line: $gene: ${dmet3_assay_id}: DMET3 ASSAY ID NOT IN DMET2\n";
                    $dmet3_extra_marker_count++;
                }
                next DMET3_EXPERIMENT_ASSAYID;
            }

            if ( defined $diffs{$dmet3_line} ) {
                if ( !$use_count && !$use_summary ) {
                    if ($use_lines) {
                        my $dmet2_line
                            = $dmet2{$dmet3_experiment}{$dmet3_assay_id}
                            {"line"};
                        print "\n";
                        print "DMET2 [$dmet2_line]:"
                            . $dmet2_marker_lines[ $dmet2_line - 1 ];
                        print "DMET3 [$dmet3_line]:"
                            . $dmet3_marker_lines[ $dmet3_line - 1 ];
                    }

                    print $diffs{$dmet3_line} . "\n";
                }
            }
        }

    }

    if ($use_count) {
        my $count = ( scalar keys(%diffs) ) + $dmet2_missing_marker_count
            + $dmet3_extra_marker_count;
        print "$count\n";
        exit(0);
    }
    my $dmet3_extra_experiment_markers = 0;
    foreach my $experiment ( keys( %dmet3_extra_experiments ) ) {
        $dmet3_extra_experiment_markers += $dmet3_extra_experiments{$experiment};
    }
    my $dmet3_duplicate_markers_total = scalar keys( %dmet3_duplicate_markers) ;
    print "----------------------------------------------------\n";
    print "SUMMARY\n";
    print "----------------------------------------------------\n";
    print "DMET2 missing experiment total: ";
    print scalar keys ( %missing_experiments ) . "\n";
    foreach my $experiment ( sort keys( %missing_experiments ) ) {
        print "$experiment\n";
    }
    print "DMET2 marker             total: "
        . $dmet2_stats{'TOTAL_MARKERS'} . "\n";
    print "DMET2 missing marker     total: $dmet2_missing_marker_count\n";
    print "DMET2 invalid line       total: " . $dmet2_stats{'INVALID_LINES'} . "\n";
    print "DMET3 extra   experiment total: ";
    print scalar keys( %dmet3_extra_experiments) . "\n";
    foreach my $experiment ( sort keys (%dmet3_extra_experiments) ) {
        print "DMET3 $experiment\n";
    }
    print "DMET3 marker             total: "
        . $dmet3_stats{'TOTAL_MARKERS'} . "\n";
    print "DMET3 invalid line       total: " . $dmet3_stats{'INVALID_LINES'} . "\n";
    print "DMET3 extra marker       total: $dmet3_extra_marker_count\n";
    print "DMET3 extra expr. marker total: $dmet3_extra_experiment_markers\n";
    
    print "DMET3 duplicate markers  total: $dmet3_duplicate_markers_total\n";
    print "DMET2 compared  markers  total: ";
    my $dmet2_compared_markers = ($dmet2_stats{'TOTAL_MARKERS'} - $dmet2_missing_marker_count) - $dmet2_stats{'INVALID_LINES'};
    print  "$dmet2_compared_markers ";
    print  "($dmet2_stats{'TOTAL_MARKERS'} - $dmet2_missing_marker_count ";
    print "-  " . $dmet2_stats{'INVALID_LINES'} . ")\n";
    print "DMET3 compared markers   total: ";
    my $dmet3_compared_markers = ($dmet3_stats{'TOTAL_MARKERS'} - $dmet3_extra_marker_count - $dmet3_extra_experiment_markers - $dmet3_duplicate_markers_total ) - $dmet3_stats{'INVALID_LINES'};
    print "$dmet3_compared_markers ";
    print "($dmet3_stats{'TOTAL_MARKERS'} - $dmet3_extra_marker_count - $dmet3_extra_experiment_markers - $dmet3_duplicate_markers_total ";
    print "- " . $dmet3_stats{'INVALID_LINES'} . ")\n";
    foreach my $assayId ( sort keys( %dmet3_duplicate_markers ) ) {
        print "$assayId duplicate         count: " . $dmet3_duplicate_markers{$assayId} . "\n";
    }
    print "Diff                     total: " . scalar keys(%diffs) . "\n";

    foreach my $col (
        sort( { return ( $MARKER_HEADERS{$a} <=> $MARKER_HEADERS{$b} );
            } keys(%MARKER_HEADERS) )
        )
    {

        if ( $diff_stats{$col} ) {
            printf "%-24s diffs: %.0d\n", $col, $diff_stats{$col};
        }
    }
    print "Discrepancies            total: ";
    if ($use_known_discrepancies) {
        my $total = 0;
        foreach my $key ( keys(%known_discrepancies) ) {
            $total += $known_discrepancies{$key};
        }
        print "$total\n";
        foreach my $key ( keys(%known_discrepancies) ) {
            if ( $known_discrepancies{$key} ) {
                printf "%-30s", $known_discrepancy_description{$key};
                print " : " . $known_discrepancies{$key} . "\n";
            }
        }
    }
    else {
        print "N/A (not specified, use -k)\n";
    }

    print "----------------------------------------------------\n";
    if ( $dmet2_compared_markers != $dmet3_compared_markers ) {
        print "!!! WARNING !!! Possible reporting bug!\n";
        print "!!! WARNING !!! Number of DMET2 and DMET3 markers did not match!\n";
        print "!!! WARNING !!! DMET2 $dmet2_compared_markers != DMET3 $dmet3_compared_markers!\n";
    }

}

# end report
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
    'count' => \$use_count,
    'help'  => sub {
        pod2usage(
            -verbose => 2,
            -section => "SYNOPSIS|DESCRIPTION|CONFIGURATION",
            -exitval => 2
        );

    },
    'known'   => \$use_known_discrepancies,
    'line'    => \$use_lines,
    'summary' => \$use_summary,
) or pod2usage();

initialize();
read_marker_files();
compare_files();
report();

exit(0);

#------------------------------------------------------------------------------
#                               EOF
#------------------------------------------------------------------------------
1;
##############################################################################

=head1 NAME

compare_dmet_marker_reports.pl - compare DMET2 & DMET3 Marker report files.

=cut    

=head1 SYNOPSIS

compare_dmet_marker_reports.pl [options] <DMET2 file> <DMET3 file>


Options:
    --count, -c      Just print one line with a count of the diffs. 
    --summary, -s    Don't print detail data, just the summary report. 
    --known, -k      known discrepancies, compare and report on known
                     but otherwise acceptable discrepancies.
    --line, -l       Print the lines of marker records with diffs. 
    --help           Show this message

Example:

    compare_dmet_marker_reports.pl ABCB1.marker.rpt ABCB1.dmet3_3.0.marker.rpt

=head1 DESCRIPTION

C<compare_dmet_marker_reports.pl> is a command-line interface designed to 
ease testing DMET2 output with DMET3 output. This script normalizes the
output between the two versions and reports discrepancies per line. 


=head1 CONFIGURATION

None

=cut

##############################################################################
