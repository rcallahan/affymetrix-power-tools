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
#include "copynumber/CNReferenceMethodPDNN.h"
//
#include "copynumber/MersenneTwister.h"
//
#include "file/TsvFile/TsvFile.h"
#include "file5/File5.h"
#include "file5/File5_Tsv.h"
#include "util/AffxMultiDimensionalArray.h"
#include "util/Verbose.h"
//
#include <algorithm>
#include <limits>
//


	
/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNReferenceMethodPDNN::explainSelf() 
	{
	CNReferenceMethodPDNN obj;
	SelfDoc doc;
	doc.setDocName(obj.getType());
	doc.setDocDescription(obj.getDescription());
	doc.setDocOptions(obj.getDefaultDocOptions());
	return doc;
	}

/**
 * @brief Default Getter method for parameters and their documentation.
 * @return map of parameters and their descriptions.
 */
std::vector<SelfDoc::Opt> CNReferenceMethodPDNN::getDefaultDocOptions() 
	{
	std::vector<SelfDoc::Opt> opts;

	// SelfDoc::Opt(name, type, value, default, min, max, description)

	return opts;
	}


/**
 * @brief This static function should be overridden by child classes
 * to return an object of the correct type initialized correctly
 * with the parameters in the string, string map. All objects
 * created this way should be deleted when finished using.
 *
 * @param param - Map of key/value pairs to initialize the object.
 *
 * @return Pointer toCreate object, this should be sub casted as necessary.
 */
SelfCreate * CNReferenceMethodPDNN::newObject(
	std::map<std::string,std::string>& params) 
	{
	SelfDoc doc = explainSelf();
	std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
	CNReferenceMethodPDNN* pMethod = new CNReferenceMethodPDNN();
	std::string strPrefix = getPrefix();
	
	return pMethod;
	}

/**
 * @brief Constructor
 */
CNReferenceMethodPDNN::CNReferenceMethodPDNN()
{
}

/**
 * @brief Destructor
 */
CNReferenceMethodPDNN::~CNReferenceMethodPDNN()
{
}

/**
 * @brief Run the calculations.
 */
void CNReferenceMethodPDNN::run()
{
	AffxString strReferenceFileName = m_pEngine->getOpt("reference-output");
	AffxString strProbeFileName = m_pEngine->getOpt("probe-file");
	Verbose::out(1, "CNReferenceMethodPDNN::run(...) start");
	AffxArray<PDNNProbe> vPDNNProbes;
	loadProbes(strReferenceFileName, strProbeFileName, vPDNNProbes);

	for (unsigned int uiProbeLength = 2; (uiProbeLength <= 100); uiProbeLength++)
	{
		int iPopulationSize = 0;
		for (int iIndex = 0; (iIndex < vPDNNProbes.getCount()); iIndex++)
		{
			PDNNProbe* p = vPDNNProbes.getAt(iIndex);
			if (p->Sequence.length() == uiProbeLength) {iPopulationSize++;}
		}
		if (iPopulationSize > (uiProbeLength - 1))
		{
			// Load up the population.
			AffxArray<PDNNProbe> vProbes;
			vProbes.resize(iPopulationSize);
			int i = 0;
			for (int iIndex = 0; (iIndex < vPDNNProbes.getCount()); iIndex++)
			{
				PDNNProbe* p = vPDNNProbes.getAt(iIndex);
				if (p->Sequence.length() == uiProbeLength) {vProbes[i] = p; i++;}
			}
			int iSampleSize = Min(200000, iPopulationSize);
			std::vector<int> vProbeIndexes;
			getRandomSample(iPopulationSize, iSampleSize, vProbeIndexes);
			ColumnVector vPopulationIntensities(iPopulationSize);
			AffxMultiDimensionalArray<char> mxPopulationDinucleotideTypes(iPopulationSize, (uiProbeLength - 1));
			for (int iIndex = 0; (iIndex < iPopulationSize); iIndex++)
			{
				PDNNProbe* p = vProbes.getAt(iIndex);
				vPopulationIntensities.element(iIndex) = log(p->MedianIntensity);
				for (int iBaseIndex = 0; (iBaseIndex < (uiProbeLength - 1)); iBaseIndex++)
				{
					mxPopulationDinucleotideTypes.set(iIndex, iBaseIndex, getDinucleotideType(p->Sequence.charAt(iBaseIndex), p->Sequence.charAt(iBaseIndex + 1)));
				}
			}

			// Load up the sampled data
			AffxMultiDimensionalArray<char> mxSampledDinucleotideTypes(iSampleSize, (uiProbeLength - 1));
			ColumnVector vSampledIntensities(iSampleSize);
//			Verbose::out(1, "*");
			for (int iIndex = 0; (iIndex < vProbeIndexes.size()); iIndex++)
			{
				int iProbeIndex = vProbeIndexes[iIndex];
				PDNNProbe* p = vProbes.getAt(iProbeIndex);
//				Verbose::out(1, "PDNN Random Sample ProbeID = " + ToStr(p->ProbeID));
				vSampledIntensities(iIndex + 1) = log(p->MedianIntensity);
				for (int iBaseIndex = 0; (iBaseIndex < (uiProbeLength - 1)); iBaseIndex++)
				{
					mxSampledDinucleotideTypes.set(iIndex, iBaseIndex, getDinucleotideType(p->Sequence.charAt(iBaseIndex), p->Sequence.charAt(iBaseIndex + 1)));
				}
			}
//			Verbose::out(1, "*");

			// Initialize the run
			ColumnVector vInitialEnergy(10);
			for (int iIndex = 0; (iIndex < 10); iIndex++) {vInitialEnergy[iIndex] = 0.1;}
			ColumnVector vEnergy = vInitialEnergy;
			ColumnVector vPrevEnergy = vInitialEnergy;

			ColumnVector vInitialWeights(uiProbeLength - 1);
			for (int iIndex = 0; (iIndex < (uiProbeLength - 1)); iIndex++) {vInitialWeights[iIndex] = (1.0 / (double)(uiProbeLength - 1));}
			ColumnVector vWeights = vInitialWeights;
			ColumnVector vPrevWeights = vInitialWeights;
			
			double dPrevFit = 1;
			double dTOL = 1e-2;

			Matrix X, U, V; DiagonalMatrix S;

			// Run to convergence
			for (int iIteration = 1; (iIteration <= 10); iIteration++)
			{
				X.ReSize(iSampleSize, 10);
				for (int iSampleIndex = 0; (iSampleIndex < iSampleSize); iSampleIndex++)
				{
					for (int iDinucleotideType = 1; (iDinucleotideType <= 10); iDinucleotideType++)
					{
						double dSum = 0;
						for (int iBaseIndex = 0; (iBaseIndex < (uiProbeLength - 1)); iBaseIndex++)
						{
							if (mxSampledDinucleotideTypes.get(iSampleIndex, iBaseIndex) == iDinucleotideType)
							{
								dSum += vWeights[iBaseIndex];
							}
						}
						X((iSampleIndex + 1), iDinucleotideType) = dSum;
					}
				}
				svd(X, S, U, V, vSampledIntensities, vEnergy);
				X.ReSize(iSampleSize, (uiProbeLength - 1));
				for (int iSampleIndex = 0; (iSampleIndex < iSampleSize); iSampleIndex++)
				{
					for (int iBaseIndex = 0; (iBaseIndex < (uiProbeLength - 1)); iBaseIndex++)
					{
						X(iSampleIndex + 1, iBaseIndex + 1) = vEnergy(mxSampledDinucleotideTypes.get(iSampleIndex, iBaseIndex));
					}
				}
				svd(X, S, U, V, vSampledIntensities, vWeights);
				CNAnalysisMethod::memory("CNReferenceMethodPDNN::run(...)");
				flipSignIfNeeded(vWeights, vEnergy, X);
				if (hasConverged(uiProbeLength, iIteration, vSampledIntensities, X, vEnergy, vPrevEnergy, vWeights, vPrevWeights, dPrevFit, dTOL)) {break;}
			}

			// Calculate Predicted Intensities
			ColumnVector vPredictedIntensities(iPopulationSize);
			for (int k = 1; (k <= iPopulationSize); k++)
			{
				vPredictedIntensities(k) = 0;
				for (int j = 1; (j <= (uiProbeLength - 1)); j++)
				{
					vPredictedIntensities(k) += (vEnergy(mxPopulationDinucleotideTypes.get((k - 1), (j - 1))) * vWeights(j));
					vProbes.getAt(k - 1)->PredictedIntensity = exp(vPredictedIntensities(k));
				}
			}

			// Calculate statistics
			ColumnVector v = (vPopulationIntensities - vPredictedIntensities);
			double dSSE = 0;
			for (int i = 1; (i <= v.Nrows()); i++)
			{
				dSSE += v(i) * v(i);
			}
			double dSSTO = var(vPopulationIntensities) * (iPopulationSize - 1);
			double dR2 = 1.0 - dSSE / dSSTO;
			double dOverallFit = corr(vPopulationIntensities, vPredictedIntensities);
			Verbose::out(1, "PDNN Model Fit\tProbeLength = " + ToStr(uiProbeLength) + "\tR2 = " + ToStr(dR2) + "\tOverallFit = " + ToStr(dOverallFit));
			vProbes.nullAll();
		}
	}
	updateReference(strReferenceFileName, vPDNNProbes);

        m_pEngine->setOpt("pdnn-reference-values-exist","true");

	vPDNNProbes.deleteAll();
	Verbose::out(1, "CNReferenceMethodPDNN::run(...) end");
}

void CNReferenceMethodPDNN::loadProbes(	const std::string& strReferenceFileName, 
                                        const std::string& strProbeFileName, 
                                        AffxArray<PDNNProbe>& vPDNNProbes)
{
  affx::File5_File file5;
  affx::File5_Group* group5 = NULL;
  affx::File5_Tsv* tsv5 = NULL;

  int iCount=0;
  if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
    {
      iCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "Cyto2", "ProbeEffects");
      file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
      group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
      tsv5 = group5->openTsv("ProbeEffects", affx::FILE5_OPEN);

    }
  else
    {
      iCount = affx::File5_Tsv::getFileTsvLineCount(strReferenceFileName, "CN5", "CN5.plier.feature-response");
      file5.open(strReferenceFileName, affx::FILE5_OPEN_RO);
      group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
      tsv5 = group5->openTsv("CN5.plier.feature-response", affx::FILE5_OPEN);
    }
  if (iCount <= 0) {throw(Except("Cannot find any probes in the CN reference file to process."));}

  vPDNNProbes.reserve(iCount);
  int i = 0;
  double d = 0;
  float f = 0;
  while (tsv5->nextLine() == affx::FILE5_OK)
    {
      PDNNProbe* p = new PDNNProbe;
      tsv5->get(0, 0, &i); p->ProbeID = i;
      tsv5->get(0, 1, &d); p->ProbeEffect = d;
      tsv5->get(0, 2, &i); p->ProbeSetIndex = i;
      tsv5->get(0, 3, &i); p->Allele = (char)i;
      tsv5->get(0, 4, &f); p->MedianIntensity = f;
      vPDNNProbes.add(p);
    }
  tsv5->close();
  delete tsv5;

  group5->close();
  delete group5;
  file5.close();

  vPDNNProbes.quickSort(0);

  affx::TsvFile tsv;
  tsv.m_optAutoTrim = true;
        
  AffxString strSequence;
  PDNNProbe objSearch;

  int iXCount = m_pEngine->getOptInt("num-cols");
  int iYCount = m_pEngine->getOptInt("num-rows");
  m_pEngine->setOpt("array-size", ToStr(iYCount));

  std::string strX, strY;
  tsv.open(strProbeFileName);
  tsv.bind(0,1,&strX, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,2,&strY, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,4,&strSequence, affx::TSV_BIND_REQUIRED);
        
  while (tsv.nextLevel(0) == affx::TSV_OK){
    int x = AffxByteArray(strX).parseInt();
    int y = AffxByteArray(strY).parseInt();
    if (x >= m_pEngine->getOptInt("num-cols")) {continue;}
    if (y >= m_pEngine->getOptInt("num-rows")) {continue;}
    strSequence = strSequence.reverse();
    switch(strSequence.length()) {
    case 23: strSequence = "AT" + strSequence; break;
    case 50: strSequence = strSequence.substring(1, strSequence.length()); break;
    case 48: strSequence = "A" + strSequence; break;
    case 47: strSequence = "AT" + strSequence; break;
    case 46: strSequence = "ATG" + strSequence; break;
    case 44: strSequence = "A" + strSequence; break;
    case 43: strSequence = "AT" + strSequence; break;
    case 42: strSequence = "ATG" + strSequence; break;
    case 41: strSequence = "ATGC" + strSequence; break;
    case 39: strSequence = "A" + strSequence; break;
    case 38: strSequence = "AT" + strSequence; break;
    }
    objSearch.ProbeID = y * iYCount + x + 1; // Using iYCount intentionally, because X does not cover the whole array.
    int iSearchIndex = vPDNNProbes.binarySearch(objSearch, 0);
    if (iSearchIndex != -1)   {
      PDNNProbe* p = vPDNNProbes.getAt(iSearchIndex);
      p->Sequence = strSequence;
    }

  }
  tsv.clear();
}

void CNReferenceMethodPDNN::getRandomSample(int iPopulationSize, int iSampleSize, std::vector<int>& vRandomSample)
{
	vRandomSample.clear();
	MTRand mt(5489);
	for (int iIndex = 0; (iIndex < iSampleSize); iIndex++)
	{
		int i = mt.randInt(iPopulationSize - 1);
		bool bAlreadyUsed = false;
		for (int j = 0; (j < vRandomSample.size()); j++)
		{
			if (i == vRandomSample[j]) {bAlreadyUsed = true; break;}
		}
		if (bAlreadyUsed) {iIndex--; continue;}
		vRandomSample.push_back(i);
	}
	std::sort(vRandomSample.begin(), vRandomSample.end());
}

char CNReferenceMethodPDNN::getDinucleotideType(char c1, char c2)
{
	if ((c1 == 'A') && (c2 == 'A')) {return 1;}
	else if ((c1 == 'T') && (c2 == 'T')) {return 1;}
	else if ((c1 == 'A') && (c2 == 'T')) {return 2;}
	else if ((c1 == 'A') && (c2 == 'G')) {return 3;}
	else if ((c1 == 'C') && (c2 == 'T')) {return 3;}
	else if ((c1 == 'A') && (c2 == 'C')) {return 4;}
	else if ((c1 == 'G') && (c2 == 'T')) {return 4;}
	else if ((c1 == 'T') && (c2 == 'A')) {return 5;}
	else if ((c1 == 'T') && (c2 == 'G')) {return 6;}
	else if ((c1 == 'C') && (c2 == 'A')) {return 6;}
	else if ((c1 == 'T') && (c2 == 'C')) {return 7;}
	else if ((c1 == 'G') && (c2 == 'A')) {return 7;}
	else if ((c1 == 'G') && (c2 == 'G')) {return 8;}
	else if ((c1 == 'C') && (c2 == 'C')) {return 8;}
	else if ((c1 == 'G') && (c2 == 'C')) {return 9;}
	else if ((c1 == 'C') && (c2 == 'G')) {return 10;}
	else return 0;
}

void CNReferenceMethodPDNN::svd(        Matrix& mxX, 
                                        DiagonalMatrix& mxDiagonal, 
                                        Matrix& mxU, 
                                        Matrix& mxV, 
                                        Matrix& mxData, 
                                        ColumnVector& vResult)
{
	SVD(mxX,mxDiagonal,mxU,mxV);
	invertDiagonal(mxDiagonal);
	vResult = mxV * mxDiagonal * mxU.t() * mxData;
}

void CNReferenceMethodPDNN::invertDiagonal(DiagonalMatrix& mxDiagonal)
{
	DiagonalMatrix mxTemp = mxDiagonal;
	for (int x = 0; ( x < mxTemp.Ncols()); x++)
	{
		for (int y = 0; (y < mxTemp.Nrows()); y++)
		{
			if (x == y) 
			{
				double d = 1.0 / mxTemp.element(y, x);
				if (fabs(d) > 1e6) {d = 0;}
				mxDiagonal(y + 1, x + 1) = d;
			}
		}
	}
}

void CNReferenceMethodPDNN::flipSignIfNeeded(ColumnVector& vWeights, ColumnVector& vEnergy, Matrix& X)
{
	int iSum = 0;
	for (int iElementIndex = 0; (iElementIndex < vWeights.Nrows()); iElementIndex++)
	{
		if (vWeights.element(iElementIndex) < 0) {iSum++;}
	}
	if (iSum > (vWeights.Nrows() / 2.0))
	{
		for (int iElementIndex = 0; (iElementIndex < vWeights.Nrows()); iElementIndex++)
		{
			vWeights.element(iElementIndex) *= -1;
		}
		for (int iElementIndex = 0; (iElementIndex < vEnergy.Nrows()); iElementIndex++)
		{
			vEnergy.element(iElementIndex) *= -1;
		}
		for (int y = 0; (y < X.Nrows()); y++)
		{
			for (int x = 0; (x < X.Ncols()); x++)
			{
				X.element(y, x) *= -1;
			}
		}
	}
}

bool CNReferenceMethodPDNN::hasConverged(       unsigned int uiProbeLength, 
                                                int iIteration, 
                                                ColumnVector& vSampledIntensities, 
                                                Matrix& X, 
                                                ColumnVector& vEnergy, 
                                                ColumnVector& vPrevEnergy, 
                                                ColumnVector& vWeights, 
                                                ColumnVector& vPrevWeights, 
                                                double& dPrevFit, 
                                                double dTOL)
{
	bool bConverged = false;
	Matrix mx = X * vWeights;
	double dFit = corr(vSampledIntensities, mx);
	mx = vEnergy - vPrevEnergy;
	double delta1 = norm(mx);
	mx = vWeights - vPrevWeights;
	double delta2 = norm(mx);
	double delta3 = fabs(dPrevFit - dFit);

	//Verbose::out(1, "ProbeLength = " + ToStr(uiProbeLength) + "\tInteration = " + ToStr(iIteration) + "\tEnergyDelta = " + ::getDouble(delta1, 10) + "\tWeightsDelta = " + ::getDouble(delta2, 10) + "\tfit = " + ToStr(dFit));
	if ((delta3 < dTOL) || ((delta1 < dTOL) && (delta2 < dTOL)))
	{
		//Verbose::out(1, "ProbeLength = " + ToStr(uiProbeLength) + "\tMethod has converged. TOL = " + ToStr(dTOL));
		bConverged = true;
	}
	else
	{
		vPrevWeights = vWeights;
		vPrevEnergy = vEnergy;
		dPrevFit = dFit;
	}
	return bConverged;
}

void CNReferenceMethodPDNN::updateReference(const std::string& strReferenceFileName, AffxArray<PDNNProbe>& vPDNNProbes)
{

	affx::File5_File file5;
	affx::File5_Group* group5 = NULL;
	affx::File5_Tsv* tsv5 = NULL;
	file5.open(strReferenceFileName, affx::FILE5_OPEN);

        if ((m_pEngine->isOptDefined("cyto2")) && (m_pEngine->getOptBool("cyto2")))
        {
                group5 = file5.openGroup("Cyto2", affx::FILE5_OPEN);
	        tsv5 = group5->openTsv("ProbeEffects", affx::FILE5_REPLACE);

        }
        else
        {
                group5 = file5.openGroup("CN5", affx::FILE5_OPEN);
	        tsv5 = group5->openTsv("CN5.plier.feature-response", affx::FILE5_REPLACE);
        }
	tsv5->defineColumn(0, 0, "probe_id", affx::FILE5_DTYPE_INT);
	tsv5->defineColumn(0, 1, "feature_response", affx::FILE5_DTYPE_DOUBLE);
	tsv5->defineColumn(0, 2, "ProbeSetIndex", affx::FILE5_DTYPE_INT);
	tsv5->defineColumn(0, 3, "AlleleCode", affx::FILE5_DTYPE_INT);
	tsv5->defineColumn(0, 4, "MedianNormalizedIntensity", affx::FILE5_DTYPE_FLOAT);
	tsv5->defineColumn(0, 5, "PredictedIntensity", affx::FILE5_DTYPE_FLOAT);
        //tsv5->gotoLine(0);
	for (int iIndex = 0; (iIndex < vPDNNProbes.getCount()); iIndex++)
	{
		PDNNProbe* p = vPDNNProbes.getAt(iIndex);
		tsv5->set_i(0, 0, p->ProbeID);
		tsv5->set_d(0, 1, p->ProbeEffect);
		tsv5->set_i(0, 2, p->ProbeSetIndex);
		tsv5->set_i(0, 3, (int)p->Allele);
		tsv5->set_f(0, 4, p->MedianIntensity);

		tsv5->set_f(0, 5, p->PredictedIntensity);
		tsv5->writeLevel(0);
	}

	tsv5->close();
	delete tsv5;
	group5->close();
	delete group5;
	file5.close();
}

