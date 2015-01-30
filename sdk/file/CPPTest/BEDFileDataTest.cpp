////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License 
// (version 2.1) as published by the Free Software Foundation.
// 
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
// 
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA 
//
////////////////////////////////////////////////////////////////
#include "file/CPPTest/BEDFileDataTest.h"
//
#include "file/BEDFileData.h"
//

using namespace affxbed;

CPPUNIT_TEST_SUITE_REGISTRATION( BEDFileDataTest );


void BEDFileDataTest::setUp()
{
}

void BEDFileDataTest::tearDown()
{
}

void BEDFileDataTest::testIntervalGroupClear()
{
	TagValuePairType params;
	IntervalEntry intervals;
	IntervalGroup group;
	group.intervals.push_back(intervals);
	group.intervals.push_back(intervals);
	CPPUNIT_ASSERT(group.intervals.size() == 2);
	group.parameters.push_back(params);
	group.parameters.push_back(params);
	CPPUNIT_ASSERT(group.parameters.size() == 2);
	group.Clear();
	CPPUNIT_ASSERT(group.intervals.size() == 0);
	CPPUNIT_ASSERT(group.parameters.size() == 0);
}

void BEDFileDataTest::testFormatTrack()
{
	BEDFileData bed;
	bed.FormatTrack("AFFX", "desc");

	std::string t = bed.Track();
	std::string ref = "track name=\"AFFX\" description=\"desc\"";
	CPPUNIT_ASSERT(t == ref);
}

void BEDFileDataTest::testFormatBrowser()
{
	BEDFileData bed;
	bed.FormatBrowser("chr1", 1536283, 224672335);
	CPPUNIT_ASSERT(bed.Browser() == "browser position chr1:1536283-224672335");
}

void BEDFileDataTest::testRead_no_file()
{
	const char *file = "./data/no_file";
	
	BEDFileData bed;
	bed.FileName() = file;
	CPPUNIT_ASSERT(bed.Read() == false);
}

void BEDFileDataTest::testClear()
{
	const char *file = "./data/test1.bed";
	
	BEDFileData bed;
	bed.FileName() = file;
	CPPUNIT_ASSERT(bed.Read() == true);
	bed.Clear();
	CPPUNIT_ASSERT(bed.Browser() == "");
	CPPUNIT_ASSERT(bed.Track() == "");
	CPPUNIT_ASSERT(bed.IntervalGroups().size() == 0);
}

void BEDFileDataTest::testRead_with_no_params()
{
	const char *file = "./data/test1.bed";
	
	BEDFileData bed;
	bed.FileName() = file;
	CPPUNIT_ASSERT(bed.Read() == true);

	CPPUNIT_ASSERT(bed.Browser() == "browser position chr1:1536283-224672335");
	CPPUNIT_ASSERT(bed.Track() == "track name=\"AFFX track\" description=\"Thr=1.300000 Gap=4 Run=6 Less=false\"");

	CPPUNIT_ASSERT(bed.IntervalGroups().size() == 1);

	IntervalGroup &group = *bed.IntervalGroups().begin();

	CPPUNIT_ASSERT(group.parameters.size() == 0);
	CPPUNIT_ASSERT(group.intervals.size() == 27);

	IntervalEntry &an = *group.intervals.begin();
	CPPUNIT_ASSERT(an.seq == "chr1");
	CPPUNIT_ASSERT(an.start == 1536283);
	CPPUNIT_ASSERT(an.stop == 1536288);
	CPPUNIT_ASSERT(an.overlap == 0);
	CPPUNIT_ASSERT(an.probeSetName == "");
	CPPUNIT_ASSERT(an.strand == ' ');
}

void BEDFileDataTest::testRead_with_params()
{
	const char *file = "./data/test2.bed";
	
	BEDFileData bed;
	bed.FileName() = file;
	CPPUNIT_ASSERT(bed.Read() == true);

	CPPUNIT_ASSERT(bed.Browser() == "browser position chr1:1536283-224672335");
	CPPUNIT_ASSERT(bed.Track() == "track name=\"AFFX track\" description=\"Thr=1.300000 Gap=4 Run=6 Less=false\"");

	CPPUNIT_ASSERT(bed.IntervalGroups().size() == 3);

	IntervalGroupListIt git;
	TagValuePairTypeList::iterator pit;
	IntervalEntryList::iterator iit;
	IntervalGroup *group;
	IntervalEntry *interval;
	TagValuePairType *param;

	git = bed.IntervalGroups().begin();
	{
		group = &(*git);

		CPPUNIT_ASSERT(group->parameters.size() == 2);
		pit = group->parameters.begin();
		param = &(*pit);
		CPPUNIT_ASSERT(param->Tag == "param1");
		CPPUNIT_ASSERT(param->Value == "value1");
		++pit;
		param = &(*pit);
		CPPUNIT_ASSERT(param->Tag == "param2");
		CPPUNIT_ASSERT(param->Value == "value2");

		CPPUNIT_ASSERT(group->intervals.size() == 6);
		iit = group->intervals.begin();
		interval = &(*iit);
		CPPUNIT_ASSERT(interval->seq == "chr1");
		CPPUNIT_ASSERT(interval->start == 1536283);
		CPPUNIT_ASSERT(interval->stop == 1536288);
		CPPUNIT_ASSERT(interval->overlap == 0);
		CPPUNIT_ASSERT(interval->probeSetName == "");
		CPPUNIT_ASSERT(interval->strand == ' ');

		iit = group->intervals.end();
		--iit;
		interval = &(*iit);
		CPPUNIT_ASSERT(interval->seq == "chr1");
		CPPUNIT_ASSERT(interval->start == 224672329);
		CPPUNIT_ASSERT(interval->stop == 224672335);
		CPPUNIT_ASSERT(interval->overlap == 0);
		CPPUNIT_ASSERT(interval->probeSetName == "");
		CPPUNIT_ASSERT(interval->strand == ' ');
	}

	++git;
	{
		group = &(*git);

		CPPUNIT_ASSERT(group->parameters.size() == 1);
		pit = group->parameters.begin();
		param = &(*pit);
		CPPUNIT_ASSERT(param->Tag == "param3");
		CPPUNIT_ASSERT(param->Value == "value3");


		CPPUNIT_ASSERT(group->intervals.size() == 1);
		iit = group->intervals.begin();
		interval = &(*iit);
		CPPUNIT_ASSERT(interval->seq == "chr10");
		CPPUNIT_ASSERT(interval->start == 86451955);
		CPPUNIT_ASSERT(interval->stop == 86451960);
		CPPUNIT_ASSERT(interval->overlap == 0);
		CPPUNIT_ASSERT(interval->probeSetName == "");
		CPPUNIT_ASSERT(interval->strand == ' ');
	}

	++git;
	{
		group = &(*git);

		CPPUNIT_ASSERT(group->parameters.size() == 2);
		pit = group->parameters.begin();
		param = &(*pit);
		CPPUNIT_ASSERT(param->Tag == "param4");
		CPPUNIT_ASSERT(param->Value == "value4");
		++pit;
		param = &(*pit);
		CPPUNIT_ASSERT(param->Tag == "param5");
		CPPUNIT_ASSERT(param->Value == "value5");


		CPPUNIT_ASSERT(group->intervals.size() == 20);
		iit = group->intervals.begin();
		interval = &(*iit);
		CPPUNIT_ASSERT(interval->seq == "chr11");
		CPPUNIT_ASSERT(interval->start == 955420);
		CPPUNIT_ASSERT(interval->stop == 955425);
		CPPUNIT_ASSERT(interval->overlap == 0);
		CPPUNIT_ASSERT(interval->probeSetName == "");
		CPPUNIT_ASSERT(interval->strand == ' ');

		iit = group->intervals.end();
		--iit;
		interval = &(*iit);
		CPPUNIT_ASSERT(interval->seq == "chr11");
		CPPUNIT_ASSERT(interval->start == 1394265);
		CPPUNIT_ASSERT(interval->stop == 1394270);
		CPPUNIT_ASSERT(interval->overlap == 0);
		CPPUNIT_ASSERT(interval->probeSetName == "");
		CPPUNIT_ASSERT(interval->strand == ' ');
	}
}

void BEDFileDataTest::testRead()
{
	const char *file = "./data/test.promoter.annotations.txt";
	
	BEDFileData bed;
	bed.FileName() = file;
	CPPUNIT_ASSERT(bed.Read() == true);

	CPPUNIT_ASSERT(bed.IntervalGroups().size() == 1);
	IntervalEntryList &an = bed.IntervalGroups().begin()->intervals;
	CPPUNIT_ASSERT(an.size() == 14);

	IntervalEntryListIt it = an.begin();

	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "1553567_s_at");
	CPPUNIT_ASSERT(it->start == 1);
	CPPUNIT_ASSERT(it->stop == 14787);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);

	++it;
	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "242131_at");
	CPPUNIT_ASSERT(it->start == 1);
	CPPUNIT_ASSERT(it->stop == 14787);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 1, 0.00001);

	++it;
	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "238199_x_at");
	CPPUNIT_ASSERT(it->start == 1);
	CPPUNIT_ASSERT(it->stop == 14787);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 2, 0.00001);
	++it;
	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "224372_at");
	CPPUNIT_ASSERT(it->start == 1);
	CPPUNIT_ASSERT(it->stop == 14787);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);

	++it;
	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "224373_s_at");
	CPPUNIT_ASSERT(it->start == 1);
	CPPUNIT_ASSERT(it->stop == 14787);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);
	++it;

	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "1555653_at");
	CPPUNIT_ASSERT(it->start == 1);
	CPPUNIT_ASSERT(it->stop == 14787);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);

	++it;
	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "211600_at");
	CPPUNIT_ASSERT(it->start == 1);
	CPPUNIT_ASSERT(it->stop == 14787);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);

	++it;
	CPPUNIT_ASSERT(it->seq == "chrM");
	CPPUNIT_ASSERT(it->probeSetName == "224375_at");
	CPPUNIT_ASSERT(it->start == 6759);
	CPPUNIT_ASSERT(it->stop == 16571);
	CPPUNIT_ASSERT(it->strand == '-');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);

	++it;
	CPPUNIT_ASSERT(it->seq == "chr20");
	CPPUNIT_ASSERT(it->probeSetName == "1552998_at");
	CPPUNIT_ASSERT(it->start == 55851);
	CPPUNIT_ASSERT(it->stop == 65800);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);

	it = an.end();
	--it;
	CPPUNIT_ASSERT(it->seq == "chr20");
	CPPUNIT_ASSERT(it->probeSetName == "204432_at");
	CPPUNIT_ASSERT(it->start == 293715);
	CPPUNIT_ASSERT(it->stop == 303688);
	CPPUNIT_ASSERT(it->strand == '+');
	CPPUNIT_ASSERT_DOUBLES_EQUAL(it->overlap, 0, 0.00001);
}
