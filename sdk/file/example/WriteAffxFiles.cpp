////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

//
#include "file/BARFileWriter.h"
#include "file/CDFFileData.h"
#include "file/CELFileWriter.h"
#include "file/CHPFileWriter.h"
//
#include <cstdio>
#include <cstring>
#include <iostream>
#include <string.h>
#include <string>
//

using namespace std;

void write_bar(const char *file)
{
	affxbarwriter::CBARFileWriter bar;
	bar.SetFileName(file);
	bar.CreateNewFile();

	TagValuePairType param;
	param.Tag="Param1";
	param.Value="NoVal";
	bar.AddAlgorithmParameter(param.Tag.c_str(), param.Value.c_str());

	// Affy BAR files currently expect 1 int for the position and 1 or more floats.
	bar.AddColumn(affxbar::BAR_DATA_INTEGER);
	bar.AddColumn(affxbar::BAR_DATA_FLOAT);

	bar.SetNumberSequences(1);
	affxbar::BarSequenceResultData data;
	srand(999999);
	for (int k=0; k<bar.GetNumberSequences(); k++)
	{
		affxbar::CGDACSequenceResultItem *p = bar.GetResultsPtr(k);
		p->SetNumberDataPoints(5);
		p->SetName("Seq1");
		p->SetGroupName("Group1");
		p->SetVersion("1.0");
		p->AddParameter("SeqParam1", "NoVal");
		for (int l=0; l<p->GetNumberDataPoints(); l++)
		{
			data.iValue = l;
			p->SetDataPoint(l, 0, data);
			data.fValue = (rand() % 100) + (0.2f * l);
			p->SetDataPoint(l, 1, data);
		}
	}

	if (bar.Save() == false)
	{
		printf("failed to save file\n");
		return;
	}
	bar.Close();
}

void write_cel(const char *file)
{
	affxcel::CCELFileWriter cel;
	cel.SetFileName(file);

	//cel.CreateNewFile();


	// Due to a limitation with GCOS, the name of the algorithm must be "Percentile"
	cel.SetAlgorithmName("Percentile");
	cel.AddAlgorithmParameter("Percentile", "50");
	cel.AddAlgorithmParameter("CellMargin", "2");
	cel.AddAlgorithmParameter("OutlierHigh", "1.500");
	cel.AddAlgorithmParameter("OutlierLow", "1.004");


	cel.SetMargin(2);
	cel.SetDimensions(126, 126);
	cel.SetChipType("Test3");
	GridCoordinatesType grid;
	grid.upperleft.x = 0;
	grid.upperleft.y = 0;

	grid.upperright.x = 1000;
	grid.upperright.y = 0;

	grid.lowerleft.x = 0;
	grid.lowerleft.y = 1000;

	grid.lowerright.x = 1000;
	grid.lowerright.y = 1000;

	cel.SetGridCorners(grid);

	cel.AllocateEntries();
	affxcel::CELFileEntryType entry;
	srand(999999);
	for (int ir=0; ir<cel.GetRows(); ir++)
	{
		for (int ic=0; ic<cel.GetCols(); ic++)
		{
			entry.Intensity = (float)(rand() % 1000);
			entry.Pixels = 10;
			entry.Stdv = 1.3f;
			cel.SetCellEntry(ic, ir, &entry);
		}
	}

	cel.SetMask(0, 0, true);
	cel.SetMask(0, 1, true);
	cel.SetOutlier(0, 2, true);

	if (cel.WriteXDABCel() == false)
	{
		printf("failed to save file\n");
		return;
	}
	cel.Close();
}

void write_chp(const char *file, const char *cdfFile)
{
	affxcdf::CCDFFileData cdf;
	cdf.SetFileName(cdfFile);
	cdf.Read();

	affxchpwriter::CCHPFileWriter chp;
	chp.SetFileName(file);
	if (chp.CreateNewFile() == false)
	{
		printf("failed to create a chp file\n");
		return;
	}
	chp.InitializeForWriting(cdf);
	chp.SetParentCelFileName("no_cel.cel");

	// Due to a limitation with GCOS, the name of the alg must be "Stat" for expression arrays.
	chp.SetAlgorithmName("Stat");
	chp.SetAlgorithmVersion("1.0.02");
	chp.AddAlgorithmParameter("p1", "1.0");
	chp.AddAlgorithmParameter("p2", "1.1");
	chp.AddAlgorithmParameter("p3", "1.2");
	chp.AddChipSummaryParameter("s1", "2.0");
	chp.AddChipSummaryParameter("s2", "2.1");

	// Due to a limitation with GCOS, the name of the progid must be "GeneChip.CallGEBaseCall.1" for expression arrays.
	chp.SetProgID("GeneChip.CallGEBaseCall.1");

	affxchp::CExpressionProbeSetResults entry;
	affxcdf::CCDFProbeSetInformation info;
	for (int i=0; i<chp.GetHeader().GetNumProbeSets(); i++)
	{
		cdf.GetProbeSetInformation(i, info);
		entry.m_HasCompResults = false;
		entry.Detection = i % ABS_NO_CALL;
		entry.Signal = (float)i;
		entry.DetectionPValue = (float) i / chp.GetHeader().GetNumProbeSets();
		entry.NumPairs = info.GetNumLists();
		entry.NumUsedPairs = info.GetNumLists();
		chp.SetExpressionEntry(i, &entry);
	}

	if (chp.Save() == false)
	{
		printf("failed to save the chp file\n");
		return;
	}
	chp.Clear();
}

void test_file_writers(std::string file, std::string cdfFile)
{
	if (strstr(file.c_str(), ".cel") != NULL || strstr(file.c_str(), ".CEL") != NULL)
	{
		write_cel(file.c_str());
	}

	else if (strstr(file.c_str(), ".chp") != NULL || strstr(file.c_str(), ".CHP") != NULL)
	{
		write_chp(file.c_str(), cdfFile.c_str());
	}

	else if (strstr(file.c_str(), ".bar") != NULL || strstr(file.c_str(), ".bar") != NULL)
	{
		write_bar(file.c_str());
	}
}


