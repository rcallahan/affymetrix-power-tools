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
/**
 * @file   translation/regression/README.txt
 * @author Mybrid Spalding
 * @date   Fri Jun  6 11:19:52 PDT 2008
 * @brief  Regression notes. 
 */

###################################################################################

This document describes the differences between this regression test
suite and the typical affy/sdk regression tests using the
"RegressionTest.h" object.

In general, running "make regression" in the top level "translation"
directory will run the tests and give a pass/fail message as with the
typical regression run.

==================================================================================

1.) Data directory needs to be moved.

The data directory currently reside under the regression source code
"input" directory ("translation/regression/input").  Root privileges
not available to me are required to move these directories to the
typical regression data storage path.

See: affy/sdk/regression-data/data/README.txt

     To prevent accidental updates to this data, it is owned by
     root. To add data to this directory, you will need to login
     to a "root-permitted" (no_squash_root) host such as
     "oratino" to modify the data.


Once the input is moved to the final regression place then the INPUT_DIR constant
in the "ATDRegession.h" header file will need to be change.


   adams:{mybrid}-bash-[ ~/devel/affymetrix/affy/sdk/translation/regression ]

   2008-06-06 11:30:33 % grep INPUT_DIR ATDRegression.h

   const std::string INPUT_DIR           = "input";

   adams:{mybrid}-bash-[ ~/devel/affymetrix/affy/sdk/translation/regression ]

==================================================================================

2.) DMET2

DMET3 C++ was first ported from DMET2 Perl code. As such the DMET3 C++
code was necessarily tested to meet the output for the DMET2 Perl
code. It was decided to make this a permanent test going forward as
the DMET3 code does differ with respect to the call logic. Therefore
the DMET2 regression ONLY tests the minimal call data elements and no
other reporting data.

In the event that the DMET3 calling logic (like with customer
overrides and copy number) vary such that it is untenable to maintain
backwards DMET2 regression testing then the DMET2 regression should be
removed. To facilitate this the DMET2 regression code was broken out
into its own set of files.

Finally, the data set tested can be added to any time by simply
dropping files with the appropriate names in the input directory.

What's required is the input Genotype Short report file with
corresponding ".dmet2_marker.reg" and ".dmet2_haplotype.reg" files of
the DMET2 results. These files can be created using a Perl script found
in the "util" subdirectory given the raw ".marker.rpt" and ".hap.rpt"
files. 

The report files must have the complete Genotype Short report file
name before the report extensions. In addition. the Genotype Short
Report must have "DMET2" in the name to signify the file is to be used
for DMET2 regression.


==================================================================================

3.) DMET3

At the time of this writing there is no DMET3 only regression testing. Presumably
this should be added as time moves on.

==================================================================================

3.) Business logic required over line-by-line comparison

The regression testing for allele translation is non-numeric. In
addition when entire experiments are missing then the line-by-line
regression report is not helpful. Therefore a custom report with
sophisticated business logic is used. In this way when the failed test
is the result of a missing an entire experiment, or an extra
experiment, only one exception of this erroneous experiment is
reported and not all the lines associated.



###################################################################################


