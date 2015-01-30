#!/usr/bin/env perl
###############################################################################
#
#  query_tsv_file.pl
#
#  Revision Control
#  ----------------
#  [$Id: probeset_tsv_file_compare.pl,v 1.6 2009-03-23 22:13:55 mspald Exp $]
#
#  Synopsis:
#  ----------------
#
#
###############################################################################
#////////////////////////////////////////////////////////////////
#//
#// Copyright (C) 2008 Affymetrix, Inc.
#//
#// This program is free software; you can redistribute it and/or modify 
#// it under the terms of the GNU General Public License (version 2) as 
#// published by the Free Software Foundation.
#// 
#// This program is distributed in the hope that it will be useful, 
#// but WITHOUT ANY WARRANTY; without even the implied warranty of 
#// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
#// General Public License for more details.
#// 
#// You should have received a copy of the GNU General Public License 
#// along with this program;if not, write to the 
#// 
#// Free Software Foundation, Inc., 
#// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#//
#////////////////////////////////////////////////////////////////
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

my $tsv_file_a = "";
my $tsv_file_a_fh = undef;
my $tsv_file_b = "";
my $tsv_file_b_fh = undef;

my %gene_probe_set = ();

my %exper_probe_set
    = ( 'A' => {},
        'B' => {},
    );

my %exper_probe_linenum
    = ( 'A' => {},
        'B' => {},
    );

my %exper_probe_show
    = ( );

my @exper_probe_sets;

my @headers = ();
my $experiment_col = -1;
my $gene_col = -1;
my $probeset_col = -1;
my $line_num = 0;

my ( $use_verbose, $use_number, $use_show );

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
        warn "Two TSV files of identical type must be supplied\n\n";
        pod2usage();
    }

    if ( !-f $ARGV[0] ) {
        croak "$ARGV[0]: file not found\n";
    }
    if ( !-f $ARGV[1] ) {
        croak "$ARGV[1]: file not found\n";
    }

    $tsv_file_a = $ARGV[0];
    $tsv_file_b = $ARGV[1];

    return;

}

# end intialize
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  initialize
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Read the header lines of both TSV files. Insure the number of columns
#  and names are identical.
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
sub initialize_headers {

    $tsv_file_a_fh = IO::File->new( $tsv_file_a, "<" );

    if ( $tsv_file_a_fh == 0 ) {
        croak "$tsv_file_a: failed open $OS_ERROR\n";
    }

    $tsv_file_b_fh = IO::File->new( $tsv_file_b, "<" );

    if ( $tsv_file_b_fh == 0 ) {
        croak "$tsv_file_b: failed open $OS_ERROR\n";
    }

    # skip the comment lines.

    my $a_header = <$tsv_file_a_fh>;
    my $b_header = <$tsv_file_b_fh>;

    $line_num++;
    while ( $a_header =~ m{^\#}xms ) {
        $a_header = <$tsv_file_a_fh>;
        $b_header = <$tsv_file_b_fh>;
        $line_num++;
    }

    if ( $a_header ne $b_header ) {
        warn "$a_header";
        warn "$b_header";
        #croak "File headers do not match.\n\n";
        warn "File headers do not match.\n\n";
    }

    # DOS or UNIX files can be compared. 
    while ( chomp $a_header ) { }

    @headers = split m{\t}xms, $a_header;
    my $count = 0;
  HEADER_COLUMN:
    foreach my $column ( @headers ) {

        if ( $column =~ m{ ^(?:assay\s*id|probe\s*set\s*id) }ixms ) {
            $probeset_col = $count;
        }
        elsif ( $column =~ m{ ^(?:chp\s*filename|experiment\s*name) }ixms ) {
            $experiment_col = $count;
        }
        elsif ( $column =~ m{ ^Gene$ }ixms ) {
            $gene_col = $count;
        }
        $count++;

        last HEADER_COLUMN 
            if ( (($probeset_col >= 0 ) && ($experiment_col >= 0)) &&($gene_col >= 0) );
    }

    if ( $probeset_col < 0 ) {
        croak "Failed to find a header named \"Assay\" or \"Probe Set\"\n\n";
    }
    if ( $experiment_col < 0 ) {
        croak "Failed to find a header named \"CHP filename\" or \"Experiment Name\"\n\n";
    }

    if ( $use_verbose ) {
        if ( $use_number ) {
            print "$line_num: ";
        }
        print $a_header . "\n";
        print "Experiment Column: $headers[$experiment_col]\n";
        print "Probe Set Column: $headers[$probeset_col]\n";
    }
    return;

}

# end initialize_headers 
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  probeset_tsv_file_read
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Returns:
#    Nothing
#
#  Environment:
#    Outputs the query results to console.
#
#
#------------------------------------------------------------------------------
sub probeset_tsv_file_read {

    my $okToCompare = 1;

    # READ file a;
    my $a_line_num = $line_num;
    
  READ_A_LINE:
    foreach my $line ( <$tsv_file_a_fh> ) {

        $a_line_num++;
        
        # DOS or UNIX
        while ( chomp $line ) {}

        if ( $line =~ m{^\s*$}xms ) {
            next READ_A_LINE;
        }

            
        
        my @cols = split m{\t}xms, $line;

        if ( @cols > @headers ) {
            my @temp = ();
            my $count = 0;
            foreach my $c ( @cols ) {
                push @temp,$c;
                $count++;
                if ( $count == scalar @headers ) {
                    last;
                }
            }
            @cols = @temp;
            
            #$okToCompare = 0;
            #warn "Invalid row:\n";
            #warn "$a_line_num: $line";
           # warn "$tsv_file_a expecting " . scalar @headers . " columns and found " . scalar @cols  . "\n";
            #last READ_A_LINE;
        }
        elsif ( @cols < @headers ) {
            foreach my $i ( @cols .. @headers ) {
                push @cols, "";
            }
        }

        my $probeset_id = $cols[$probeset_col];
        my $experiment_id = $cols[$experiment_col];
        my $gene = $cols[$gene_col];
        
        my $expr_probe_key = "${experiment_id}::$probeset_id";

        if ( defined $exper_probe_set{'A'}{$expr_probe_key} ) {
            croak "$a_line_num: $tsv_file_a $expr_probe_key: duplicate experiment, probeset detected.\n";
        }
        push @exper_probe_sets, $expr_probe_key;
        if ( $use_show ) {
            if ( $line =~ m{ $use_show }xms ) {
                $exper_probe_show{$expr_probe_key} = $line;
            }
        }
        
        $exper_probe_set{'A'}{$expr_probe_key} = [@cols];
        $exper_probe_linenum{'A'}{$expr_probe_key} = $a_line_num;
        $gene_probe_set{$expr_probe_key} = $gene;
    }
    

    close( $tsv_file_a_fh);

    if ( ! $okToCompare ) {
        close( $tsv_file_b_fh);
        return $okToCompare;
    }

    # READ file b;
    my $b_line_num = $line_num;
  READ_B_LINE:
    foreach my $line ( <$tsv_file_b_fh> ) {

        $b_line_num++;
        
        if ( $line =~ m{^\s*$}xms ) {
            next READ_A_LINE;
        }

        # DOS or UNIX
        while ( chomp $line ) {}


        my @cols = split m{\t}xms, $line;

        if ( @cols > @headers ) {
            my @temp = ();
            my $count = 0;
            foreach my $c ( @cols ) {
                push @temp,$c;
                $count++;
                if ( $count == scalar @headers ) {
                    last;
                }
            }
            @cols = @temp;
            
            #$okToCompare = 0;
            #warn "Invalid row:\n";
            #warn "$b_line_num: $line";
           # warn "$tsv_file_a expecting " . scalar @headers . " columns and found " . scalar @cols  . "\n";
            #last READ_A_LINE;
        }
        elsif ( @cols < @headers ) {
            foreach my $i ( @cols .. @headers ) {
                push @cols, "";
            }
        }

        my $probeset_id = $cols[$probeset_col];
        my $experiment_id = $cols[$experiment_col];

        my $expr_probe_key = "${experiment_id}::$probeset_id";

        if ( ! defined $exper_probe_set{'A'}{$expr_probe_key} ) {
            push  @exper_probe_sets, $expr_probe_key;
        }
        if ( defined $exper_probe_set{'B'}{$expr_probe_key} ) {
            warn "Invalid row:\n";
            warn "$b_line_num: ";
            croak "$tsv_file_b $expr_probe_key: duplicate experiment, probeset id detected.\n";
        }
        if ( $use_show ) {
            if ( $line =~ m{ $use_show }xms ) {
                $exper_probe_show{$expr_probe_key} = $line;
            }
        }
        $exper_probe_set{'B'}{$expr_probe_key} = [@cols];
        $exper_probe_linenum{'B'}{$expr_probe_key} = $b_line_num;

    }
    close( $tsv_file_b_fh);

    return $okToCompare;

}

# end probeset_tsv_file_read
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  probeset_tsv_file_compare
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Returns:
#    Nothing
#
#  Environment:
#    Outputs the query results to console.
#
#
#------------------------------------------------------------------------------
sub probeset_tsv_file_compare {

    @exper_probe_sets = sort @exper_probe_sets;

  PROBESET_TO_COMPARE:
    foreach my $exper_probeset_id ( @exper_probe_sets ) {

        if ( ! defined $exper_probe_set{'A'}{$exper_probeset_id} ) {
            if ( $use_number ) {
                print $exper_probe_linenum{'B'}{$exper_probeset_id} . ":";
            }
            print "> " . @{$exper_probe_set{'B'}{$exper_probeset_id}} . "\n";
            next PROBESET_TO_COMPARE;
        }
        elsif ( !defined $exper_probe_set{'B'}{$exper_probeset_id} ) {
            if ( $use_number ) {
                print $exper_probe_linenum{'A'}{$exper_probeset_id} . ":";
            }
            print "< " . @{$exper_probe_set{'A'}{$exper_probeset_id}} . "\n";
            next PROBESET_TO_COMPARE;
        }

        my @a_cols = @{$exper_probe_set{'A'}{$exper_probeset_id}};
        my @b_cols = @{$exper_probe_set{'B'}{$exper_probeset_id}};

        my @diff_cols = ();

        foreach my $i (0..(@headers -1)) {
            if ( $a_cols[$i] ne $b_cols[$i] ) {
                push @diff_cols, "[ " . $headers[$i] . "]" .  " < " . $a_cols[$i] . " > " . $b_cols[$i] . "\t";
            }
        }

        if ( @diff_cols ) {
            if ( $use_number ) {
                print $exper_probe_linenum{'A'}{$exper_probeset_id} . " " .
                    $exper_probe_linenum{'B'}{$exper_probeset_id} . ": ";
            }
            print "$exper_probeset_id:\t" . $gene_probe_set{$exper_probeset_id} . "\t";
            print @diff_cols;
            print "\n";
        }
        elsif ( $use_show && defined $exper_probe_show{$exper_probeset_id} ) {
            print $exper_probe_show{$exper_probeset_id} . "\n";
        }
    }

}


# end probeset_tsv_file_compare
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
    'number'  => \$use_number,
    'show=s'    => \$use_show, 
    'verbose' => \$use_verbose,
) or pod2usage(2);

initialize();

initialize_headers();

if ( probeset_tsv_file_read() ) {
    probeset_tsv_file_compare();
}
else {
    exit(1);
}

exit(0);

#------------------------------------------------------------------------------
#                               EOF
#------------------------------------------------------------------------------
1;
##############################################################################

=head1 NAME

probeset_tsv_file_compare.pl - Compare two TSV files of identical format, different data

=cut

=head1 SYNOPSIS

probeset_tsv_file_compare.pl [options]  <TSV file1, TSV file2> 

Options:
    --number, -n  Line number. Output the line number as "n:";
    --help        Show this message
    --show, -s    Show any record which matches the a given Perl regular expression.
    --verbose     Output header information.

Example:

    probeset_tsv_file_comparte.pl a_comprehensive.rpt b_comprehensive.rpt

=head1 DESCRIPTION

C<probeset_tsv_file_compare.pl> is a command-line interface designed to
comparte two TSV files of the same type ( same columns, etc.). The script
inspects the files for the experiment and  probe set header
column names. Experiment probe sets must be unique for this to work.
The script builds an internal hash keyed by experiment, probeset of both
files and outputs and column diffs per experiment probeset. 

=head1 CONFIGURATION

None


=cut

##############################################################################
