////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#include "mas5-stat/src/ExpressionAlgorithmImplementation.h"
//
#include <cstring>
#include <string.h>
//

std::string get_file(int argc, char *argv[], int &i)
{
	std::string file = "";

	// Get the file name.
	int j=i+1;
	while (j<argc && argv[j][0] != '-')
	{
		file += argv[j];
		if (j<argc-1 && argv[j+1][0] != '-')
			file += " ";
		++j;
	}
	i = j-1;

	return file;
}

int main(int argc, char* argv[])
{
	CExpressionAlgorithmImplementation exp;
	CExpStatAlgSettings &params = exp.GetParameters();

	// Get the inputs from the command line.
	std::string cdfFile;
	std::string celFile;
	std::string baselineFile;
	bool writeCell = false;
	std::string outputDir;
	int i=1;
	while(i<argc)
	{
		if (strcmp(argv[i], "-cdf") == 0)
			cdfFile = get_file(argc, argv, i);

		else if (strcmp(argv[i], "-i") == 0)
			celFile = get_file(argc, argv, i);

		else if (strcmp(argv[i], "-b") == 0)
			baselineFile = get_file(argc, argv, i);

		else if (strcmp(argv[i], "-Alpha1") == 0)
			params.Alpha1 = atof(argv[i+1]);

		else if (strcmp(argv[i], "-Alpha2") == 0)
			params.Alpha2 = atof(argv[i+1]);

		else if (strcmp(argv[i], "-Tau") == 0)
			params.Tau = atof(argv[i+1]);

		else if (strcmp(argv[i], "-TGT") == 0)
			params.TGT = atof(argv[i+1]);

		else if (strcmp(argv[i], "-Gamma1H") == 0)
			params.Gamma1H = atof(argv[i+1]);

		else if (strcmp(argv[i], "-Gamma1L") == 0)
			params.Gamma1L = atof(argv[i+1]);

		else if (strcmp(argv[i], "-Gamma2H") == 0)
			params.Gamma2H = atof(argv[i+1]);

		else if (strcmp(argv[i], "-Gamma2L") == 0)
			params.Gamma2L = atof(argv[i+1]);

		else if (strcmp(argv[i], "-Perturbation") == 0)
			params.Perturbation = atof(argv[i+1]);

		else if (strcmp(argv[i], "-SFMethod") == 0)
			params.SFMethod = (CExpStatAlgSettings::ScalingOptionsEnum) atoi(argv[i+1]);

		else if (strcmp(argv[i], "-NormMethod") == 0)
			params.NormMethod = (CExpStatAlgSettings::NormalizationOptionsEnum) atoi(argv[i+1]);

		else if (strcmp(argv[i], "-NormFactor") == 0)
			params.NormFactor = atof(argv[i+1]);

		else if (strcmp(argv[i], "-ScaleFactor") == 0)
			params.ScaleFactor = atof(argv[i+1]);

		else if (strcmp(argv[i], "-ScaleProbeSet") == 0)
			params.ScaleGenes.push_back(get_file(argc, argv, i));

		else if (strcmp(argv[i], "-NormProbeSet") == 0)
			params.NormGenes.push_back(get_file(argc, argv, i));

		else if (strcmp(argv[i], "-write-cell") == 0)
			writeCell = true;

		else if (strcmp(argv[i], "-out-dir") == 0)
			outputDir = get_file(argc, argv, i);
	
		++i;
	}



	////////////////////////////////////////////////////////////////////////
	////
	//// Run the expression algorithm.
	////
	////////////////////////////////////////////////////////////////////////
    if (exp.RunStat(celFile.c_str(), baselineFile.c_str(), cdfFile.c_str()) == false)
	{
		std::cerr << exp.GetError() << endl;
		return 0;
	}

	// Output the background zone information
	std::cout << "Number of Zones=" << exp.GetBackgroundZoneInfo().number_zones << endl;
	std::cout << "Smooth factor=" << exp.GetBackgroundZoneInfo().smooth_factor << endl;
	if (exp.GetBackgroundZoneInfo().pZones != NULL)
	{
		for (int iZone=0; iZone<exp.GetBackgroundZoneInfo().number_zones; ++iZone)
		{
			std::cout << "Zone" << iZone << "X=" << exp.GetBackgroundZoneInfo().pZones[iZone].center.x << endl;
			std::cout << "Zone" << iZone << "Y=" << exp.GetBackgroundZoneInfo().pZones[iZone].center.y << endl;
			std::cout << "Zone" << iZone << "B=" << exp.GetBackgroundZoneInfo().pZones[iZone].background << endl;
		}
		std::cout << endl;
	}

	// Output the header
	std::cout << "Name" << "\t";
	std::cout << "Signal" << "\t";
	std::cout << "Detection" << "\t";
	std::cout << "p-value" << "\t";
	std::cout << "pairs" << "\t";
	std::cout << "pairs used";
	if (exp.DoesCompDataExists() == true)
	{
		std::cout << "\t";
		std::cout << "Change" << "\t";
		std::cout << "Change p-value" << "\t";
		std::cout << "SLR" << "\t";
		std::cout << "SLR-Low" << "\t";
		std::cout << "SLR-High" << "\t";
		std::cout << "Common";
	}
	std::cout << endl;


	// Output the results.
	int n=exp.GetNumResults();
	for (int i=0; i<n; i++)
	{
		AbsStatExpressionProbeSetResultType *abs_result = exp.GetAbsStatResult(i);
		std::cout << exp.GetProbeSetName(i) << "\t";
		std::cout << abs_result->Signal << "\t";
		std::cout << ExpressionAbsCallString[abs_result->Detection] << "\t";
		std::cout << abs_result->DetectionPValue << "\t";
		std::cout << abs_result->NumPairs << "\t";
		std::cout << abs_result->NumUsedPairs;

		if (exp.DoesCompDataExists() == true)
		{
			CompStatExpressionProbeSetResultType *comp_result = exp.GetCompStatResult(i);
			std::cout << "\t";
			std::cout << ExpressionCompCallString[comp_result->Change] << "\t";
			std::cout << comp_result->ChangePValue << "\t";
			std::cout << comp_result->SignalLogRatio << "\t";
			std::cout << comp_result->SignalLogRatioLow << "\t";
			std::cout << comp_result->SignalLogRatioHigh << "\t";
			std::cout << comp_result->NumCommonPairs;
		}
		std::cout << endl;
	}
	std::cout << endl;

	// Output the rawq information
	std::cout << "RawQ\t" << exp.GetRawQ() << endl;

	// Output the bg stats
	AvgStdvMinMaxType bg = exp.GetBgStats();
	std::cout << "BG Avg\t" << bg.avg << endl;
	std::cout << "BG Stdv\t" << bg.stdv << endl;
	std::cout << "BG Min\t" << bg.min << endl;
	std::cout << "BG Max\t" << bg.max << endl;

	// Output the noise stats
	AvgStdvMinMaxType noise = exp.GetNoiseStats();
	std::cout << "Noise Avg\t" << noise.avg << endl;
	std::cout << "Noise Stdv\t" << noise.stdv << endl;
	std::cout << "Noise Min\t" << noise.min << endl;
	std::cout << "Noise Max\t" << noise.max << endl;

	// The control information.
	ControlInformationList &controls = exp.GetControlInfo();
	ControlInformationList::iterator iter;
	for (iter=controls.begin(); iter!=controls.end(); iter++)
	{
		ControlInformationType info = *iter;
		switch (info.qcType)
		{
		case affxcdf::CheckerboardPositiveQCProbeSetType:
			std::cout << "Corner+\t";
			break;

		case affxcdf::CheckerboardNegativeQCProbeSetType:
			std::cout << "Corner-\t";
			break;

		case affxcdf::CentralCrossPositiveQCProbeSetType:
			std::cout << "Central+\t";
			break;

		case affxcdf::CentralCrossNegativeQCProbeSetType:
			std::cout << "Central-\t";
			break;

		default:
			break;
		}
		std::cout << info.avg << "\t" << info.count << endl;
	}


	return 0;
}
