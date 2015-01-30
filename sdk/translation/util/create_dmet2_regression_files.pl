#!/usr/bin/env perl
###############################################################################
#
#  create_dmet2_regression_files.pl
#
#  Revision Control
#  ----------------
#  [$Id: create_dmet2_regression_files.pl,v 1.7 2009-01-13 18:00:00 mspald Exp $]
#
#  Synopsis:
#  ----------------
#
#  Therefore this script takes the output of DMET2 reports and
#  strips all but the data required for regression QA with the new
#  DMET3 program.
#
#  The objective is to create a normalized output that can be "diff'd"
#  against the output of DMET3 and the where no "diffs" exist then the QA
#  is considered passed.
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

my $A2ROOT = "$ENV{HOME}/devel/affymetrix/dmet2bad";
my $A3ROOT = "$ENV{HOME}/devel/affymetrix/dmet3";
my $A2BIN  = "$A2ROOT/bin/allele_translation_dmet2.pl";
my $OUTDIR = "/tmp";

my $TRANSLATION_TABLE_FILE
    = "$A3ROOT/testdata_dmet2b/DMET3_TTable_v20080110_EarlyAccess.db";
my $COPYNUMBER_FILE
    = "$A3ROOT/testdata_dmet2b/20080109_31set_DMET3_cn_16various.txt";

my $A2_MARKER_REPORT_FILE    = "$A2ROOT/test/20080109_31set_DMET2.marker.rpt";
my $A2_HAPLOTYPE_REPORT_FILE = "$A2ROOT/test/20080109_31set_DMET2.hap.rpt";

my $A2OPTIONS
    = " -use $TRANSLATION_TABLE_FILE -cn $COPYNUMBER_FILE  --outdir $OUTDIR -gsh ";

my %MARKER_HEADERS = (
    "Experiment" => 0,
    "Gene"       => 1,
    "AssayId"    => 5,
    "Call"       => 6,
    "Ref"        => 7,
    "Var"        => 8,
    "Allele"     => 9,
);

my %HAPLOTYPE_HEADERS = (
    "Experiment"  => 0,
    "Gene"        => 1,
    "Call"        => 3,
    "Call_Count"  => 4,
    "Known_Count" => 5,
);

#Experiment	Gene	Sample	Functional_Change	External_ID	Assay_ID	Basecall	Reference_Base	Variant_Base	Call	Haplotype_Marker	Change_for_Variant	Variant_cDNA_Change	Variant_DNA_Change	dbSNP_ID	Validated	Allele_Defining_Marker	Relevant_Alleles

my @MARKER_HEADER_ORDER
    = ( "Experiment", "Gene", "AssayId", "Call", "Ref", "Var", "Allele" );

my $MARKER_HEADER = "Experiment\tGene\tAssayId\tA1\tA2\tRef\tVar\tCall\n";
my $HAPLOTYPE_HEADER = "Experiment\tGene\tCall\tCall_Count\tKnown_Count\n";

my $DMET2_MARKER_REPORT_EXT    = "dmet2_marker.reg";
my $DMET2_HAPLOTYPE_REPORT_EXT = "dmet2_haplotype.reg";

#Experiment	Gene	Sample	Call	Call_Count	Known_Count	UNK_Exists	Basecall_Rate	Basecall_Count	NoCall_Count	PossibleRareAllele_Count	NotAvailable_Count

my @HAPLOTYPE_HEADER_ORDER
    = ( "Experiment", "Gene", "Call", "Call_Count", "Known_Count" );

#------------------------------------------------------------------------------
# GLOBAL VARIABLES

my ( $use_verbose, $use_command, $use_force, $use_marker_file,
    $use_haplotype_file );

my $genotype_root_dir = "";

my $input_haplotype_file = "";
my $input_marker_file    = "";

my @genotype_files = ();

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

    if ( -f $A2_HAPLOTYPE_REPORT_FILE ) {
        $input_haplotype_file = $A2_HAPLOTYPE_REPORT_FILE;
    }
    else {
        if ( !$use_haplotype_file ) {
            print STDERR "No haplotype report file (hapt.rpt) specified\n";
            pod2usage();
        }
        if ( !-f $use_haplotype_file ) {
            croak "$use_haplotype_file: file not found\n";
        }
        if ( !-r $use_haplotype_file ) {
            croak "$use_haplotype_file: persmission denied\n";
        }
        $input_haplotype_file = $use_haplotype_file;
    }

    if ( -f $A2_MARKER_REPORT_FILE ) {
        $input_marker_file = $A2_MARKER_REPORT_FILE;
    }
    else {
        if ( !$use_marker_file ) {
            print STDERR "No marker report file (hapt.rpt) specified\n";
            pod2usage();
        }
        if ( !-f $use_marker_file ) {
            croak "$use_marker_file: file not found\n";
        }
        if ( !-r $use_marker_file ) {
            croak "$use_marker_file: persmission denied\n";
        }
        $input_marker_file = $use_marker_file;
    }

    my $pipe_fh = IO::File->new("find $genotype_root_dir -type f|");

    if ( !$pipe_fh ) {
        croak "Can't open pipe: $OS_ERROR.\n";
    }

READFILE_LABEL:
    foreach my $geno_file (<$pipe_fh>) {

        chomp $geno_file;

        my $line = `head -1 $geno_file`;

        #print "$geno_file: $line";
        chomp $line if ($line);

        if ( $line && ( $line =~ m{^Sample\sName}xms ) ) {
            push @genotype_files, "$geno_file";

            #print "Found file: $geno_file\n" if ( $use_verbose );
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
#  run_a2_on_genotype_files
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Given the list of known genotype files, call the a2 script and
#  create the marker.rpt file in the same directory as the file.
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
sub run_a2_on_genotype_files( ) {

    print "run_a2_on_genotype_files: BEGIN\n";

INPUT_GENOTYPE_FILE_LABEL:
    foreach my $genotype_file (@genotype_files) {

        my $genotype_dir = $genotype_file;
        if ( $genotype_dir =~ m{/}xms ) {
            $genotype_dir =~ s{/[^/]+$}{}xms;
        }
        my $marker_report_file    = "${genotype_file}.marker.rpt";
        my $haplotype_report_file = "${genotype_file}.hap.rpt";
        my $marker_log_file       = "${genotype_file}.log";

        if ( !$use_force && ( -f $marker_report_file ) ) {
            next INPUT_GENOTYPE_FILE_LABEL;
        }

        $marker_report_file    =~ s{[^/]+/}{}gxms;
        $marker_log_file       =~ s{[^/]+/}{}gxms;
        $haplotype_report_file =~ s{[^/]+/}{}gxms;

        #my $cmd = "$A2BIN $A2OPTIONS $genotype_file >>/dev/null 2>&1 ";
        #if ( system($cmd) != 0 ) {
        #    croak "$cmd failed: $OS_ERROR\n";
        #}
        my %exprs = ();

        my $genotype_fh = IO::File->new($genotype_file, "<");
        if ( ! $genotype_fh ) {
            croak "$genotype_file: can't open file, $OS_ERROR\n";
        }
      GENOTYPE_LINE:
        foreach my $line ( <$genotype_fh> ) {
            while ( chomp $line ) {};
            next GENOTYPE_LINE if ( !$line ) ;
            next GENOTYPE_LINE if ( $line =~ m{^\s*\#}xms ) ;
            
            my @splits = split( m{\t}xms, $line);
            next GENOTYPE_LINE if ( (scalar @splits) < 6 ) ;

            my $experiment = $splits[1];
            my $gene       = $splits[2];
            $experiment =~ s{\s}{\.}gxms;
            $experiment =~ s{[\(\)\]\[!\$\+\*\&\^\%\#\@\:\;\'\"\>\<\?]}{\.}gxms;
            my $expr = "$experiment\\t$gene\\t";
            $exprs{$expr} = 1;
            
        }
        close $genotype_fh;
        $genotype_fh = 0;

        my $marker_tmp_fh = IO::File->new("$OUTDIR/$marker_report_file", ">");
        if ( ! $marker_tmp_fh ) {
            croak "$OUTDIR/$marker_report_file, can't open for write, $OS_ERROR\n";
        }
        my $marker_input_fh = IO::File->new($input_marker_file, "<");
        if ( !$marker_input_fh ) {
            croak "$input_marker_file: can't open, $OS_ERROR\n";
        }
      MARKER_LINE:
        foreach my $line ( <$marker_input_fh> ) {
            foreach my $expr ( keys(%exprs) ) {
                if ( $line =~ m{$expr}xms ) {
                    print {$marker_tmp_fh} $line;
                    next MARKER_LINE;
                }
            }
        }
        close $marker_input_fh;
        $marker_input_fh = 0;
        close $marker_tmp_fh;
        $marker_tmp_fh = 0;

        my $haplotype_tmp_fh = IO::File->new("$OUTDIR/$haplotype_report_file", ">");
        if ( ! $haplotype_tmp_fh ) {
            croak "$OUTDIR/$marker_report_file, can't open for write, $OS_ERROR\n";
        }

        my $haplotype_input_fh = IO::File->new($input_haplotype_file, "<");
        if ( !$haplotype_input_fh ) {
            croak "$input_haplotype_file: can't open, $OS_ERROR\n";
        }

      HAPLOTYPE_LINE:
        foreach my $line ( <$haplotype_input_fh> ) {
            foreach my $expr ( keys(%exprs) ) {
                if ( $line =~ m{$expr}xms ) {
                    print {$haplotype_tmp_fh} $line;
                    next HAPLOTYPE_LINE;
                }
            }
        }
        close $haplotype_input_fh;
        $haplotype_input_fh = 0;
        close $haplotype_tmp_fh;
        $haplotype_tmp_fh = 0;

        my $cmd = "mv -f $OUTDIR/$marker_report_file $genotype_dir";
        if ( system($cmd) != 0 ) {
            croak "$cmd, $OS_ERROR\n";
        }

        $cmd = "mv -f $OUTDIR/$haplotype_report_file $genotype_dir/$haplotype_report_file";
        if ( system($cmd) != 0 ) {
            croak "$cmd, $OS_ERROR\n";
        }

        #$cmd = "mv -f $OUTDIR/$marker_log_file $genotype_dir";
        #if ( system($cmd) != 0 ) {
        #    croak "$cmd, $OS_ERROR\n";
        #}

        if ( !-f "$genotype_dir/$marker_report_file" ) {
            croak "Can't find file $genotype_dir/$marker_report_file.\n";
        }
        my $temp = $genotype_file;
        $temp =~ s{$genotype_dir}{}xms;

        $cmd = "rm -f $OUTDIR/$temp.*";

        if ( system($cmd) != 0 ) {
            croak "$cmd, $OS_ERROR\n";
        }

        print ".";

    }
    print "\n";

    print "run_a2_on_genotype_files: END\n";
}

# end run_a2_on_genotype_files
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  convert_a2_marker_report_files
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Filter the input file and create an output file of the same name
#  with a new extension.
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
sub convert_a2_marker_report_files() {

INPUT_GENOTYPE_FILE_LABEL:
    foreach my $genotype_file (@genotype_files) {

        #print "processing...$genotype_file\n" if ( $use_verbose );

        my $marker_report_file = "${genotype_file}.marker.rpt";
        my $marker_qa_file     = "${genotype_file}.$DMET2_MARKER_REPORT_EXT";

        my $marker_in_fh  = IO::File->new( $marker_report_file, "<" );
        my $marker_out_fh = IO::File->new( $marker_qa_file,     ">" );

        if ( !$marker_in_fh ) {
            croak "Can't open file $marker_report_file; $OS_ERROR\n";
        }

        if ( !$marker_out_fh ) {
            croak "Can't open file $marker_qa_file; $OS_ERROR\n";
        }

        my @qa_lines = ();

    READ_LINE:
        foreach my $line (<$marker_in_fh>) {

            next READ_LINE if ( $line =~ m{ ^\s*\# }xms );

            next READ_LINE if ( $line =~ m{ ^Experiment }xms );

            chomp $line;

            my @columns = split( /\t/, $line );

            next if ( @columns < 11 );

            next if ( $columns[ $MARKER_HEADERS{"Call"} ] eq "NotAvailable" );

            my $qa_line = "";

            foreach my $header (@MARKER_HEADER_ORDER) {
                if ( $header eq "Ref" ) {
                    $columns[ $MARKER_HEADERS{$header} ]
                        = "R|" . $columns[ $MARKER_HEADERS{$header} ];
                }
                elsif ( $header eq "Var" ) {
                    my $variant = $columns[ $MARKER_HEADERS{$header} ];
                    if ( $variant =~ m{,}xms ) {
                        my @splits = split( m{,}xms, $variant );
                        @splits  = sort @splits;
                        $variant = "";
                        foreach my $split (@splits) {
                            $variant .= $split;
                            if ( $split ne $splits[$#splits] ) {
                                $variant .= ",";
                            }
                        }
                    }
                    $columns[ $MARKER_HEADERS{$header} ] = "V|" . $variant;
                }
                elsif ( $header eq "Call" ) {
                    if ( $columns[ $MARKER_HEADERS{$header} ] ) {
                        my @splits = split( /\//,
                            $columns[ $MARKER_HEADERS{$header} ] );
                        if ( @splits == 2 ) {
                            $columns[ $MARKER_HEADERS{$header} ]
                                = "A1|$splits[0]\tA2|$splits[1]";
                        }
                        if ( $columns[ $MARKER_HEADERS{$header} ]
                            =~ m{PossibleRareAllele|NoCall}xms )
                        {
                            $columns[ $MARKER_HEADERS{$header} ] .= "\t";
                        }
                    }
                }
                elsif ( $header eq "Allele" ) {
                    my $call = $columns[ $MARKER_HEADERS{$header} ];
                    if ( $call =~ m{ ([^/]+)/([^/]+) }xms ) {
                        my $first       = $1;
                        my $second      = $2;
                        my $sort_first  = $first;
                        my $sort_second = $second;
                        if ( $first =~ m{ (\d+) }xms ) {
                            $sort_first = $1;
                        }
                        if ( $second =~ m{ (\d+) }xms ) {
                            $sort_second = $1;
                        }
                        if ( $first eq "UNK" ) {
                            $call = "$second/$first";
                        }
                        elsif ( $second eq "UNK" ) {
                            $call = "$first/$second";
                        }
                        elsif ( ( $sort_first cmp $sort_second ) == -1 ) {
                            $call = "$first/$second";
                        }
                        else {
                            $call = "$second/$first";
                        }
                    }
                    elsif ($call) {
                        croak "$call: Bad Allele column!\n";
                    }
                    $columns[ $MARKER_HEADERS{$header} ] = $call;
                }

                $qa_line .= $columns[ $MARKER_HEADERS{$header} ] . "\t";
            }

            chop $qa_line;
            $qa_line .= "\n";

            push @qa_lines, $qa_line;

        }    # foreach line

        print {$marker_out_fh} $MARKER_HEADER;
        if (@qa_lines) {
        SORT_LINES:
            @qa_lines = sort @qa_lines;
            print {$marker_out_fh} @qa_lines;
            print $MARKER_HEADER if ($use_verbose);
            print @qa_lines      if ($use_verbose);
        }

        close $marker_in_fh;
        close $marker_out_fh;

    }    # foreach $geno_file, process the output.

    return;
}

# end convert_a2_marker_report_files
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  convert_a2_haplotype_report_files
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Filter the input file and create an output file of the same name
#  with a new extension.
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
sub convert_a2_haplotype_report_files() {

INPUT_GENOTYPE_FILE_LABEL:
    foreach my $genotype_file (@genotype_files) {

        #print "processing...$genotype_file\n" if ( $use_verbose );

        my $haplotype_report_file = "${genotype_file}.hap.rpt";
        my $haplotype_qa_file
            = "${genotype_file}.$DMET2_HAPLOTYPE_REPORT_EXT";

        print "$haplotype_qa_file, $haplotype_report_file\n";
        
        my $haplotype_in_fh  = IO::File->new( $haplotype_report_file, "<" );
        my $haplotype_out_fh = IO::File->new( $haplotype_qa_file,     ">" );

        if ( !$haplotype_in_fh ) {
            croak "Can't open file $haplotype_report_file; $OS_ERROR\n";
        }

        if ( !$haplotype_out_fh ) {
            croak "Can't open file $haplotype_qa_file; $OS_ERROR\n";
        }

        my @qa_lines = ();

    READ_LINE:
        foreach my $line (<$haplotype_in_fh>) {
            next READ_LINE if ( $line =~ m{ ^\s*\# }xms );

            next READ_LINE if ( $line =~ m{ ^Experiment }xms );

            chomp $line;

            my @columns = split( /\t/, $line );

            next if ( @columns < 11 );

            #next if ( $columns[ $HEADERS{"Call"} ] eq "NotAvailable" );

            my $qa_line = "";

            foreach my $header (@HAPLOTYPE_HEADER_ORDER) {

                if ( $header eq "Call" ) {

                    # Order the two calls by any number in the system.
                    my $call = $columns[ $HAPLOTYPE_HEADERS{$header} ];

                    if ( $call
                        =~ m{ [^\d/]*(\d+)[^/]*/[^\d/]*(\d+)[^\d]* }xms )
                    {
                        my $first  = $1;
                        my $second = $2;
                        if ( $first > $second ) {
                            my @splits = split( m{/}xms, $call, 2 );
                            $call = "$splits[1]/$splits[0]";
                            $columns[ $HAPLOTYPE_HEADERS{$header} ] = $call;
                        }
                    }
                }

                $qa_line .= $columns[ $HAPLOTYPE_HEADERS{$header} ] . "\t";
            }

            chop $qa_line;
            $qa_line .= "\n";

            push @qa_lines, $qa_line;

        }    # read each input file line.

        print {$haplotype_out_fh} $HAPLOTYPE_HEADER;
        
        if (@qa_lines) {
        SORT_LINES:
            
            
            #@qa_lines = sort {
            #    my @a_splits = split( m{\t}xms, $a );
            #    my @b_splits = split( m{\t}xms, $b );

             #   if ( @a_splits < 3 ) {
              #croak "$a";
             #   }
             #   if ( @b_splits < 3 ) {
             #       croak "$b";
             #   }
             #   my $a_call = $a_splits[2];
             #
             #   my $b_call     = $b_splits[2];
             #   my $a_call_num = 0;
             #   my $b_call_num = 0;
             #   if ( $a_call =~ m{(\d+)}xms ) {
             #       $a_call_num = $1;
             #   }
             #   if ( $b_call =~ m{(\d+)}xms ) {
             #       $b_call_num = $1;
             #   }
             #   if (   ( $a_call_num && $b_call_num )
             #       && ( $a_call_num != $b_call_num ) )
             #   {
             #
             #      return ( $a_call_num <=> $b_call_num );
             #   }
             #   return ( $a cmp $b );
             #} @qa_lines;

            @qa_lines = sort @qa_lines;

            print {$haplotype_out_fh} @qa_lines;
            print @qa_lines if ($use_verbose);
        }

        close $haplotype_in_fh;
        close $haplotype_out_fh;

        my $cmd = "sort $haplotype_qa_file >  ${haplotype_qa_file}.sort";
        if ( system( $cmd ) != 0 ) {
            croak "${cmd}: FAILED, $OS_ERROR\n";
        }
        $cmd = "/bin/mv -f ${haplotype_qa_file}.sort $haplotype_qa_file";
        if ( system( $cmd ) != 0 ) {
            croak "${cmd}: FAILED, $OS_ERROR\n";
        }

    }

    # foreach $geno_file, process the output.

    return;
}

# end convert_a2_haplotype_report_files
#-----------------------------------------------------------------------------
#-----------------------------------------------------------------------------
#
#  cleanup_a2_marker_files
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Synopsis:
#
#  Remove 0 length files and report files with only headers.
#  . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . . .
#
#  Parameters:
#
#
#
#  Returns:
#    Nothing
#
#  Environment:
#    croaks
#      - system errors
#
#------------------------------------------------------------------------------
sub cleanup_a2_files( ) {

INPUT_GENOTYPE_FILE_LABEL:
    foreach my $genotype_file (@genotype_files) {

        my $marker_report_file = "${genotype_file}.marker.rpt";
        my $marker_qa_file     = "${genotype_file}.$DMET2_MARKER_REPORT_EXT";
        my $haplotype_qa_file
            = "${genotype_file}.$DMET2_HAPLOTYPE_REPORT_EXT";

        if ( -z "$marker_qa_file" ) {
            unlink($marker_qa_file)
                || croak "Can't remove file: $marker_qa_file\n";
        }
        if ( -z "$haplotype_qa_file" ) {
            unlink($haplotype_qa_file)
                || croak "Can't remove file: $haplotype_qa_file\n";
        }

        #my $line = `tail -1 $marker_report_file`;

        #if ( $line =~ m{Sample\s+Functional_Change}xms ) {
        #    unlink( $marker_report_file ) ||
        #   croak "Can't remove file: $marker_report_file\n";
        #}

    }

}

# end cleanup_a2_files
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
    'command'     => \$use_command,
    'force'       => \$use_force,
    'haplotype=s' => \$use_haplotype_file,
    'help'        => sub {
        pod2usage(
            -verbose => 2,
            -section => "SYNOPSIS|DESCRIPTION|CONFIGURATION",
            -exitval => 2
        );

    },
    'marker=s' => \$use_marker_file,
    'verbose'  => \$use_verbose,
) or pod2usage();

if ($use_command) {
    print "$A2BIN $A2OPTIONS\n";
    exit(0);
}

initialize();

run_a2_on_genotype_files();

convert_a2_marker_report_files();

convert_a2_haplotype_report_files();

#cleanup_a2_files();

exit(0);

#------------------------------------------------------------------------------
#                               EOF
#------------------------------------------------------------------------------
1;
##############################################################################

=head1 NAME

create_dmet2_regression_files.pl - normalize DMET2 reports for QA with DMET3. 

=cut    

=head1 SYNOPSIS

create_dmet2_regression_files.pl [options] [directory of Genotype Short Report files]


Options:
    --command, -c    Output the DMET2 command used.
    --force, -f      Overwrite existing qa files.
    --haplotype      REQUIRED: Haplotype file (hap.rpt) to use. 
    --help           Show this message
    --marker,-m      REQUIRED: Marker file (marker.rpt) to use. 
    --verbose, -v    Output the results to console as well as to file. 

Example:

    create_dmet2_regression_files.pl -o 20080109_31set --marker 20080109_31set_DMET2.marker.rpt --haplotype 20080109_31set_DMET2.hap.rpt testdata_dmet2b 

=head1 DESCRIPTION

C<create_dmet2_regression_files.pl> is a command-line interface designed to 
ease testing DMET2 output with DMET3 output. The output of this script
should be directly testable against DMET3 output.

Every Genotype Short Report file found will have additional DMET2 output files
of the same name but with extensions "dmet2_marker.reg" or "dmet2_haplotype.reg".



=head1 CONFIGURATION

None

=cut

##############################################################################
