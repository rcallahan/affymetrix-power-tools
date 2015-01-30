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
#include "file/CPPTest/BEDFileWriterTest.h"
//
#include "file/BEDFileWriter.h"
//

using namespace affxbed;

CPPUNIT_TEST_SUITE_REGISTRATION( BEDFileWriterTest );


void BEDFileWriterTest::setUp()
{
}

void BEDFileWriterTest::tearDown()
{
}

void BEDFileWriterTest::testWrite()
{
	const char *file = "./data/new.bed";
	BEDFileWriter bed;
	IntervalGroup group;
	IntervalEntry interval;
	TagValuePairType param;

	bed.FormatTrack("AFFX", "desc");
	bed.FormatBrowser("chr1", 1, 2);
	bed.FileName() = file;

	group.Clear();
	param.Tag = "param1";
	param.Value = "value1";
	group.parameters.push_back(param);
	param.Tag = "param2";
	param.Value = "value2";
	group.parameters.push_back(param);

	interval.seq = "chr1";
	interval.start = 1;
	interval.stop = 2;
	interval.overlap = 0;
	interval.probeSetName = "";
	interval.strand = ' ';
	group.intervals.push_back(interval);

	interval.seq = "chr1";
	interval.start = 3;
	interval.stop = 4;
	interval.overlap = 0;
	interval.probeSetName = "";
	interval.strand = ' ';
	group.intervals.push_back(interval);

	bed.IntervalGroups().push_back(group);

	group.Clear();
	param.Tag = "param3";
	param.Value = "value3";
	group.parameters.push_back(param);

	interval.seq = "chr2";
	interval.start = 11;
	interval.stop = 12;
	interval.overlap = 0;
	interval.probeSetName = "";
	interval.strand = ' ';
	group.intervals.push_back(interval);

	bed.IntervalGroups().push_back(group);

	CPPUNIT_ASSERT(bed.Write() == true);


	bed.Clear();

	bed.FileName() = file;
	CPPUNIT_ASSERT(bed.Read() == true);
	CPPUNIT_ASSERT(bed.Browser() == "browser position chr1:1-4");

	std::string t = bed.Track();
	std::string ref = "track name=\"AFFX\" description=\"desc\"";
	CPPUNIT_ASSERT(t == ref);


	CPPUNIT_ASSERT(bed.IntervalGroups().size() == 2);

	IntervalGroupListIt git;
	TagValuePairTypeList::iterator pit;
	IntervalEntryList::iterator iit;

	git = bed.IntervalGroups().begin();
	{
		group = *git;

		CPPUNIT_ASSERT(group.parameters.size() == 2);
		pit = group.parameters.begin();
		param = *pit;
		CPPUNIT_ASSERT(param.Tag == "param1");
		CPPUNIT_ASSERT(param.Value == "value1");
		++pit;
		param = *pit;
		CPPUNIT_ASSERT(param.Tag == "param2");
		CPPUNIT_ASSERT(param.Value == "value2");

		CPPUNIT_ASSERT(group.intervals.size() == 2);
		iit = group.intervals.begin();
		interval = *iit;
		CPPUNIT_ASSERT(interval.seq == "chr1");
		CPPUNIT_ASSERT(interval.start == 1);
		CPPUNIT_ASSERT(interval.stop == 2);
		CPPUNIT_ASSERT(interval.overlap == 0);
		CPPUNIT_ASSERT(interval.probeSetName == "");
		CPPUNIT_ASSERT(interval.strand == ' ');

		++iit;
		interval = *iit;
		CPPUNIT_ASSERT(interval.seq == "chr1");
		CPPUNIT_ASSERT(interval.start == 3);
		CPPUNIT_ASSERT(interval.stop == 4);
		CPPUNIT_ASSERT(interval.overlap == 0);
		CPPUNIT_ASSERT(interval.probeSetName == "");
		CPPUNIT_ASSERT(interval.strand == ' ');
	}

	++git;
	{
		group = *git;

		CPPUNIT_ASSERT(group.parameters.size() == 1);
		pit = group.parameters.begin();
		param = *pit;
		CPPUNIT_ASSERT(param.Tag == "param3");
		CPPUNIT_ASSERT(param.Value == "value3");


		CPPUNIT_ASSERT(group.intervals.size() == 1);
		iit = group.intervals.begin();
		interval = *iit;
		CPPUNIT_ASSERT(interval.seq == "chr2");
		CPPUNIT_ASSERT(interval.start == 11);
		CPPUNIT_ASSERT(interval.stop == 12);
		CPPUNIT_ASSERT(interval.overlap == 0);
		CPPUNIT_ASSERT(interval.probeSetName == "");
		CPPUNIT_ASSERT(interval.strand == ' ');
	}
}
