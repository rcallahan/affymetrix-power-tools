#!/usr/bin/env perl
################################################################
##
## Copyright (C) 2008 Affymetrix, Inc.
##
## This program is free software; you can redistribute it and/or modify 
## it under the terms of the GNU General Public License (version 2) as 
## published by the Free Software Foundation.
## 
## This program is distributed in the hope that it will be useful, 
## but WITHOUT ANY WARRANTY; without even the implied warranty of 
## MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU 
## General Public License for more details.
## 
## You should have received a copy of the GNU General Public License 
## along with this program;if not, write to the 
## 
## Free Software Foundation, Inc., 
## 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
##
################################################################
###############################################################################
#
#  extract_comprehensive_for_override.pl
#
#  Revision Control
#  ----------------
#  [$Id: extract_comprehensive_for_override.pl,v 1.1 2009-01-13 18:00:00 mspald Exp $]
#
#  Synopsis:
#  ----------------
#  The override file only contains those genes with no calls. This
#  script will extract any experiment gene from a comprehensive report
#  and output a suitable override file or append to an existing one. 
#
#  This script simply scrapes the comprehensive report and does not
#  read the annotation file. Therefore the data produced will not work
#  for aliases. 
#
###############################################################################
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

my $comprehensive_file = "";
my $comprehensive_fh = undef;
my $override_file = "";
my $override_fh = undef;
my $extract_gene = "";
my $extract_experiment = "";

my %exper_probe_set = ();
my %exper_gene_records = ();




my @headers = ();
my %header_index = ();

my $experiment_col = -1;
my $probeset_col = -1;
my $gene_col = -1;
my $line_num = 0;

my ( $use_verbose, $use_number, $use_show );

my $OVERRIDE_HEADER_LINE = "CHP Filename	Gene	Common Name	Probe Set ID	Basecall	Override Comment	Reference Allele	Variant Allele";
    
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

    if ( @ARGV < 3 || @ARGV > 4 ) {
        pod2usage();
    }

    if ( !-f $ARGV[0] ) {
        croak "$ARGV[0]: file not found\n";
    }

    $comprehensive_file = $ARGV[0];
    $override_file = $ARGV[1];

    $extract_experiment = $ARGV[2];
    if ( @ARGV > 3 ) {
        $extract_gene = $ARGV[3];
    }
    

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

    my $override_output_headers = 0;
    $comprehensive_fh = IO::File->new( $comprehensive_file, "<" );

    if ( $comprehensive_fh == 0 ) {
        croak "$comprehensive_file: failed open $OS_ERROR\n";
    }

    $override_fh = 0;
    
    if ( -e $override_file ) {
        $override_fh = IO::File->new( $override_file, ">>" );
    }
    else {
        $override_fh = IO::File->new( $override_file, ">" );
        $override_output_headers = 1;
    }
        

    if ( $override_fh == 0 ) {
        croak "$override_file: failed open $OS_ERROR\n";
    }

    # skip the comment lines in the comprehensive report.


    my $num_override_headers = 8;
    my $override_header_count = 0;
    
  HEADER_LINE:
    while ( 1 ) {
        my $comprehensive_header = <$comprehensive_fh>;
        print $comprehensive_header;
        if ( $comprehensive_header =~ m{^\#}xms ) {
            if ( $override_output_headers &&
                     ($override_header_count < $num_override_headers) ) {
                $override_header_count++;
                print {$override_fh} $comprehensive_header;
            }
            next HEADER_LINE;
        }
        if ( $override_output_headers ) {
            print {$override_fh} "$OVERRIDE_HEADER_LINE\n";
        }

        while ( chomp $comprehensive_header ) {}

        @headers = split(m{ \t }xms, $comprehensive_header);
        my $count = 0;
      HEADER_NAME:
        foreach my $header_name ( @headers ) {
            $header_index{$header_name} = $count;
            $count++;
        }
        last HEADER_LINE;
        
    }

}

# end initialize_headers 
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  comprehensive_report_read
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
#
#------------------------------------------------------------------------------
sub comprehensive_report_read {

    my $line_num = 0;
    
  READ_COMPREHENSIVE_LINE:
    foreach my $line (<$comprehensive_fh>) {

        $line_num++;

        
        # DOS or UNIX
        while ( chomp $line ) {}

        if ( $line =~ m{^\s*$}xms ) {
            next READ_COMPREHENSIVE_LINE;
        }
        
        my @cols = split m{\t}xms, $line;

        if ( @cols > @headers ) {
            warn "Invalid row:\n";
            warn "$line_num: $line";
            warn "$comprehensive_file expecting " . scalar @headers . " columns and found " . scalar @cols  . "\n";
            last READ_A_LINE;
        }
        elsif ( @cols < @headers ) {
            foreach my $i ( @cols .. @headers ) {
                push @cols, "";
            }
        }

        my $probeset_id = $cols[$header_index{'Probe Set ID'}];
        my $experiment_id = $cols[$header_index{'CHP Filename'}];
        $experiment_id =~ s{ \.[Cc][hH][pP]$ }{}xms;

        my $gene = $cols[$header_index{'Gene'}];
        my $expr_probe_key = "${experiment_id}::$probeset_id";

        if ( defined $exper_probe_set{$expr_probe_key} ) {
            croak "$line_num: $comprehensive_file $expr_probe_key: duplicate experiment, probeset detected.\n";
        }

        $exper_probe_set{$expr_probe_key} = [@cols];
        push( @{$exper_gene_records{$experiment_id}}, [@cols]);
    }

    close( $comprehensive_fh);


}
# end comprehensive_report_read
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  override_report_write
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
#
#
#------------------------------------------------------------------------------
sub override_report_write {

    if ( not defined($exper_gene_records{$extract_experiment}) ) {
        print "$extract_experiment: not found in $comprehensive_file\n";
        print "Possible experiments: " . keys( %exper_gene_records) . "\n";
        exit(1);
    }

    my @experiment_records = @{$exper_gene_records{$extract_experiment}};

    my @override_cols = split( m{ \t }xms, $OVERRIDE_HEADER_LINE);
    
  EXPERIMENT_RECORD:
    foreach my $record_ref ( @experiment_records ) {
        
        my @record = @{$record_ref};
        if (!$extract_gene ||
                ($record[$header_index{'Gene'}] eq $extract_gene)) {

            my $count = 0;
          COLUMN:
            foreach my $col ( @override_cols ) {
                if ( $col eq "Reference Allele" ) {
                    $col = "Reference Base";
                }
                elsif ( $col eq "Variant Allele" ) {
                    $col = "Variant Base";
                }
                
                if ( defined( $header_index{$col} ) ) {
                    print {$override_fh} $record[$header_index{$col}];
                }
                if ( $count < (@record -1)) {
                    print {$override_fh} "\t";
                }
                $count++;
            }
            print {$override_fh}  "\n";
        }
    }

    close( $override_fh);
    
}
# end override_report_write
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
    'verbose' => \$use_verbose,
) or pod2usage(2);

initialize();

initialize_headers();

comprehensive_report_read();
override_report_write();


exit(0);

#------------------------------------------------------------------------------
#                               EOF
#------------------------------------------------------------------------------
1;
##############################################################################

=head1 NAME

extract_comprehensive_for_override.pl - Extract records from a comprehensive


=cut

=head1 SYNOPSIS

extrac_comprehensive_for_override.pl [options]  <input_comprehesive.rpt > <output_override.rpt> [<Experiment>] [<Gene>]

Options:
    --help        Show this message
    --verbose     Output header information.

Example:

    extract_comprehensive_for_override.pl test_comprehensive.rpt 

=head1 DESCRIPTION

C<extract_comprehensive_for_override.pl> is a command-line script
to convert comprehensive records to override records. The script
will create the output file or append to the file if it already exists.

Since the script only reads the comprehensive report then alleles will
aliases will fail since aliases cannot be input into the translation script.

=head1 CONFIGURATION

None


=cut

##############################################################################
