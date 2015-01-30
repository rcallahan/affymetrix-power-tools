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

//

#include "algorithm/em/PrimeEM.h"
//
#include <cppunit/CompilerOutputter.h>
#include <cppunit/extensions/HelperMacros.h>
#include <cppunit/extensions/TestFactoryRegistry.h>
#include <cppunit/ui/text/TestRunner.h>
//
#include <cmath>
#include <cstring>
#include <string>
//
#define SAMPLES 89
#define GROUP    3
// This is a double.
#define K      4.0

using namespace std;

// test data, real data with 1, 2 and 3 clusters respectively
double gTestData[GROUP][2][SAMPLES] = 
{
	{{242.61956,270.64284,269.49567,299.96376,261.05227,221.69487,285.47846,245.4567,237.04691,246.65369,265.34452,213.1842,227.63572,258.95084,243.43004,236.62447,254.31584,220.51687,235.60744,287.53137,247.25758,229.39466,247.93611,269.73922,227.08859,237.65537,276.88458,221.48915,288.2868,260.69428,266.17449,295.44566,257.15763,234.30802,297.78891,303.62558,310.17827,272.44682,268.51199,268.88533,253.42991,248.84577,222.60804,265.84817,289.57178,281.36162,284.75008,243.79828,283.43191,260.02168,219.39858,320.15289,259.2002,246.89022,287.78357,296.78309,230.50912,286.53247,256.49754,239.34509,270.82639,258.54595,250.26837,245.57994,328.72683,261.96118,278.92814,270.78793,274.39088,212.48815,228.65898,253.73635,283.85794,245.54647,272.59564,226.71105,299.38768,250.73921,252.25402,239.90874,232.3826,252.86261,328.9986,257.37016,251.69094,256.22824,234.43024,266.18346,293.23704},
	{1067.11268,1146.56683,1330.81423,1202.15579,1210.11218,1287.00868,1218.36179,1226.94554,1189.6987,1317.7677,1213.11281,1042.07794,1115.80626,1191.64105,1191.38541,1028.58037,1042.07727,1169.94908,1314.16842,1202.19091,1250.81307,1286.53378,1091.5206,1105.93964,1185.22099,1159.71319,1218.04992,917.12593,1026.85217,1164.1149,1078.10207,1001.8706,1195.16741,1113.46569,1263.48274,946.99375,929.22137,975.61752,1401.66689,1193.05439,1172.51379,1062.04988,1251.62095,1049.2366,1177.33782,1076.92501,1211.0713,1181.7589,1097.94243,1157.0031,1094.84841,1235.23178,1160.93264,1282.8976,1206.90408,945.56217,1047.83142,1285.03406,1201.40977,1042.3606,1249.87819,1164.93862,1170.09888,1155.28786,1311.64892,941.61324,1049.19215,1114.46877,1047.25301,1161.34585,1123.74644,1005.23794,1121.68841,1249.10349,1114.08839,1182.87323,1271.01918,1263.7495,1127.21591,1158.99916,1106.50432,1135.02975,1014.80834,1046.92886,1091.89957,1033.40892,1167.18985,1169.50034,850.03539}},
	{{223.52361,246.20328,237.04925,174.89228,200.49079,201.35273,226.48263,215.10963,237.36351,271.55366,276.32857,220.27987,224.04764,279.16638,201.86962,259.46803,253.22125,245.20754,194.35652,247.04461,252.05059,239.83847,220.22843,188.34149,228.42472,221.95042,212.03998,218.19286,268.04149,267.08575,248.9427,209.93155,687.65916,187.97899,210.34984,229.9774,229.52198,719.88569,232.6964,666.81055,166.4902,224.73468,231.27519,228.26813,209.45474,217.39949,237.02103,234.72981,166.19081,232.35986,257.45438,224.36186,188.9343,257.92521,217.80739,194.27333,188.57137,732.84487,176.31323,277.98629,685.78873,199.52979,193.49481,180.2313,174.86265,216.54933,217.25495,200.57251,214.15635,203.10093,250.79349,255.42111,213.28143,203.07288,231.10296,213.96174,215.68057,250.21006,242.67768,219.43702,231.01168,742.68303,263.59973,197.65844,557.33989,179.72383,722.49598,221.50704,193.11055},
	{993.24887,1103.92584,1161.89875,894.85955,1193.72044,1087.0304,1035.83509,1037.1036,1203.37657,1268.66071,1254.99834,1150.21808,1171.61531,1241.07015,1124.26841,1113.77964,1283.34058,1114.79056,1172.12473,1181.33745,1069.04691,1070.47158,1108.97993,1098.59535,1243.73865,1022.49776,1090.68659,1050.11268,1207.24005,1202.27728,1086.80342,1134.08837,711.63727,1044.77538,937.53947,1098.50999,1178.19675,669.44965,1089.03115,626.07226,1029.29672,1165.56032,1030.25694,1074.82029,948.2767,936.7401,1176.6176,1074.67108,930.69969,1045.41928,1311.10632,1228.07939,941.24979,1044.82656,1061.38917,1042.85413,1178.39355,702.55957,1028.38118,1045.16477,679.36273,1121.71534,1100.3766,1079.92249,1086.5865,1024.97727,1057.6307,1035.49286,929.52209,1016.63895,1021.63797,1189.98315,1166.82517,1276.12584,1124.40956,1134.58018,1016.64207,1232.13141,1156.98495,1063.35047,966.28309,645.77349,1160.0103,1004.34978,602.19311,1107.54563,698.45594,1118.50616,1039.69136}},
	{{531.98945,211.84991,529.85686,377.27812,468.36419,514.20691,231.30706,193.58304,504.66956,391.62333,357.2313,227.90794,434.31717,400.13258,503.45134,541.83654,494.89046,367.12623,185.48975,488.72675,480.53695,203.12775,392.70953,387.60362,208.60965,469.84432,216.34452,500.40041,375.87951,350.0387,228.42918,521.86889,504.04632,555.05553,461.62188,512.08418,555.44041,628.45268,550.33599,518.19244,511.01404,565.70197,469.80938,517.72573,495.12265,540.37437,452.97515,531.54148,479.46525,380.80784,534.08068,477.77029,490.51037,549.69219,537.4535,494.45929,552.9155,549.56007,559.92306,540.15235,476.3994,491.00794,373.42828,532.47597,390.00076,498.06107,514.20856,503.10566,470.74608,526.6147,524.48882,514.88508,482.79793,540.81715,453.10322,491.29187,400.86921,381.55002,219.43469,512.51438,404.15902,401.75395,534.75055,565.76322,528.34141,520.04727,373.29286,340.9595,486.68118},
	{198.51966,840.50556,184.71642,595.32747,167.51963,168.25056,902.35832,712.934,193.53201,557.63004,602.88193,966.23677,527.28457,645.58648,153.24265,177.23392,197.24908,608.74735,858.06396,189.01277,184.05243,764.79235,497.60216,581.95255,963.83963,190.78116,861.38058,175.2046,539.83057,473.90274,914.01921,176.65033,190.43275,200.07561,196.36426,184.53034,149.90521,198.37792,183.77075,186.09255,165.64512,191.3765,182.48727,191.59873,152.92068,187.54954,164.54044,186.04612,203.99143,519.88091,202.87835,181.61338,223.2761,206.13735,162.26493,203.70828,162.37148,170.91806,195.28529,182.61692,155.18046,210.26889,176.07442,161.17725,588.2448,174.03579,196.98436,211.73367,201.13841,214.4885,196.73094,215.35079,195.03719,144.49524,236.08333,170.7634,583.59349,557.11244,850.16736,195.29095,614.39178,589.66765,200.0847,168.48383,199.8152,186.21094,470.93076,538.48705,171.38638}}
};

// test data, real data with 1, 2 and 3 clusters respectively
int gRefCall[GROUP+2][SAMPLES] = 
{
	{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
	{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,1,2,1,2,2},
	{0,2,0,1,0,0,2,2,0,1,1,2,1,1,0,0,0,1,2,0,0,2,1,1,2,0,2,0,1,1,2,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,1,1,2,0,1,1,0,0,0,0,1,1,0},
	{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2},
	{2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,2,2,1,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,1,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,1,2,2,1,2,1,2,2}
};

// test data, real data with 1, 2 and 3 clusters respectively
double gRefConf[GROUP+2][SAMPLES] = 
{
	{0.000029,0.000149,0.003091,0.001567,0.000837,0.013483,0.000076,0.003663,0.003862,0.007458,0.000523,0.002645,0.002763,0.000663,0.002697,0,0.00081,0.007012,0.01043,0.000344,0.004264,0.010832,0.000034,0.000793,0.006006,0.002573,0.00003,0.000542,0.010629,0.000152,0.001172,0.017384,0.000904,0.00157,0.000132,0.032268,0.036939,0.010013,0.006016,0.000089,0.000789,0.000075,0.011019,0.002248,0.001045,0.004027,0.000107,0.002299,0.003263,0.000115,0.003581,0.003508,0.000187,0.005745,0.000292,0.029007,0.000413,0.000204,0.001118,0.000001,0.000729,0.000268,0.001066,0.001247,0.001749,0.009589,0.005299,0.000691,0.004222,0.008982,0.002877,0.002067,0.002189,0.004557,0.000885,0.005987,0.000126,0.004074,0.00016,0.002167,0.001637,0.000215,0.033536,0.00103,0.000001,0.001323,0.003472,0.000024,0.040892},
	{0.002393,0.001978,0.000007,0.000261,0.005808,0.001483,0.001214,0.000103,0.000149,0.000611,0.001461,0.000613,0.000644,0.00237,0.002592,0.004401,0.000145,0.001423,0.006538,0.000193,0.00526,0.002183,0.000086,0.004709,0.00176,0.000984,0.000339,0.00012,0.001788,0.00181,0.003335,0.001504,0.003938,0.002512,0.002249,0.000208,0.000308,0.001598,0.000571,0.00095,0.008044,0.00048,0.002274,0.000443,0.00158,0.004148,0.000009,0.001179,0.002818,0.001832,0.0002,0.001941,0.000021,0.009393,0.000028,0.00131,0.00873,0.000118,0.004707,0.01929,0.000429,0.002979,0.003483,0.006165,0.008367,0.000347,0.000033,0.000399,0.003683,0.000044,0.008815,0.000678,0.001923,0.009095,0.000036,0.000972,0.000422,0,0.000233,0.000061,0.006366,0.010584,0.002889,0.000173,0.010781,0.007842,0.000008,0.000109,0.001399},
	{0.000122,0.000274,0.000647,0.00357,0.000145,0.003248,0.00058,0.003181,0.000699,0.000212,0.010735,0.000778,0.021371,0.005575,0.008155,0.003266,0.002408,0.008547,0.00487,0.000982,0.000663,0.00188,0.011922,0.00062,0.00479,0.003638,0.000234,0.000539,0.000064,0.002939,0.000191,0.001631,0.000323,0.000062,0.008002,0.000065,0.019593,0.005459,0.00222,0.000098,0.003779,0.001655,0.001146,0.000041,0.007023,0.000769,0.000014,0.000547,0.008023,0.002257,0.000443,0.000459,0.018106,0.00019,0.008815,0.004795,0.0112,0.006508,0.000636,0.001681,0.003498,0.008793,0.025372,0.008604,0.000795,0.000588,0.000668,0.006826,0.008522,0.003868,0.000193,0.006193,0.00327,0.020688,0.04114,0.000728,0.000005,0.00002,0.000761,0.00052,0.001142,0.000069,0.000157,0.00997,0.000344,0.000131,0.012753,0.003641,0.000408},
	{2.91467E-05,0.00014919,0.003090859,0.001566768,0.000837266,0.013482928,7.56383E-05,0.003663003,0.003861547,0.007457614,0.000522733,0.002644658,0.002763093,0.000662684,0.002696753,0,0.000809848,0.007012188,0.010429561,0.000344455,0.004264176,0.010831654,3.39746E-05,0.00079298,0.006005585,0.002572596,3.01003E-05,0.000542521,0.010629356,0.000151753,0.001172245,0.017383814,0.000904202,0.001570046,0.000131965,0.032267869,0.036939263,0.010012627,0.006016254,8.86917E-05,0.000789285,7.54595E-05,0.011018634,0.00224793,0.001045466,0.004027307,0.000106931,0.002299368,0.003263176,0.000114739,0.003581464,0.003508449,0.000187218,0.005745232,0.000291944,0.0290066,0.00041306,0.000203788,0.001118064,6.55651E-07,0.000729203,0.000268102,0.001066089,0.00124675,0.001748681,0.009589374,0.005298734,0.000690579,0.004222333,0.008981466,0.002876639,0.00206697,0.002188623,0.004557431,0.000884712,0.005986631,0.000125587,0.004074395,0.000159681,0.002166867,0.00163728,0.000215232,0.033535779,0.001029551,9.53674E-07,0.001323283,0.003471971,2.40207E-05,0.040892243},
	{0.004845679,0.004006863,1.45435E-05,0.000528693,0.011743009,0.003005207,0.002460182,0.000208378,0.000301063,0.001239061,0.002959788,0.001243055,0.001305699,0.004801273,0.005250335,0.008905053,0.000293911,0.002883971,0.013215482,0.000391483,0.01063764,0.004422665,0.00017339,0.009527326,0.003565431,0.001993775,0.000687122,0.000242829,0.003621638,0.003667533,0.006752253,0.00304687,0.005517244,0.005087137,0.004554451,0.00042057,0.000623524,0.002239406,0.00115782,0.001331568,0.016246319,0.000973284,0.004605115,0.000898361,0.003201842,0.008393228,1.80006E-05,0.002390027,0.005705774,0.003711283,0.000404418,0.003933311,4.19617E-05,0.018955706,5.68628E-05,0.002655089,0.017625868,0.000165284,0.00952214,0.034523878,0.000601232,0.00603199,0.007051587,0.012463391,0.016896605,0.000702858,6.70552E-05,0.000808537,0.007455289,8.92282E-05,0.017794965,0.001373649,0.003896594,0.018359303,7.31349E-05,0.001970172,0.000856102,7.15256E-07,0.000472963,0.000124395,0.012867927,0.01480633,0.00584954,0.00035131,0.015081406,0.015840948,1.11461E-05,0.000221431,0.002835035}
};

// posteriors
double gPosterior[7][3] =
{
	{-0.77082276,-0.32579938,0.65786839},	 //mean
	{0.029999997,0.070000000,0.046110325},	 //std
	{0.11233905,0.20226787,0.68539310},	 //weight
	{-1.2500000,-0.58037776,0.25302526},	 //mean minimum
	{-0.59037775,0.24302526,1.2500000},	 //mean maximum
	{0.029999997,0.070000000,0.029999997},	 //std minimum
	{0.038696874,0.079571441,0.099403568}	 //std maximum
};

// other stats being tested, FLD, SEP, likelihood etc.
double gOthers[4] =
{
	// likelihood
	39.304608, // @todo: check this is correct. (Copied from the output after converting to doubles.)
	// FLD1
	5.8434429,
	// FLD2
	11.735169,
	// SEP
	8.2809296
};

// values to test  calcProbability()
double gProbInputs[6]=
{
	-0.78956616,	//x
	-0.78707337,	//mu
	0.032660954,	//sigma
	1.000000,	//weight
	2.500000,	//zCut
	0.000000	//useWeight
};
double gProbOutput = 30.528552;	// @todo: correct?

// values for testing calcGaussianWithExpTail()
double gGassianWithTail[] =
{
	2.332768, // z score
	5.000000, //
	0.250000 // weighting factor
};
double gGassianReturn = 0.506502;

// prior information for testing seed generation
double gPrior[4][3] = 
{
	{-0.66000003,0.00000000,0.66000003},	//m_mu
	{0.12500000,0.17500000,0.12500000},	//m_muSD
	{0.064999998,0.085000001,0.064999998},	//m_sigma
	{0.035000000,0.015000000,0.035000000}	//m_sigmaSD
};
// seed to be generated based on above prior and given data	in gTestData[3]
double gSeed[7][3]=
{
	{-0.77082235,-0.32575530,0.65786850},	//m_mu
	{0.016832022,0.068721704,0.046110392},	//m_sigma
	{0.11235955,0.20224719,0.68539327},	//m_weight
	{-1.2500000,-0.58037776,0.25302526},	//m_minMu
	{0.029999997,0.070000000,0.029999997},	//m_minSigma
	{-0.59037775,0.24302526,1.2500000},	//m_maxMu
	{0.099999994,0.10000000,0.099999994}	//m_maxSigma
};

double gProb[89][3] = {
	{1.6251690e-005,0.011946543,7.9777722},
	{5.5685539,0.0037552123,5.1084735e-006},
	{1.2811068e-005,0.0094173579,7.9466205},
	{0.71134406,0.17024176,0.00023159153},
	{1.3982775e-005,0.010278675,7.9916506},
	{1.0426107e-005,0.0076641766,7.6913648},
	{5.7206225,0.0039097439,5.3186836e-006},
	{6.2425394,0.0045120837,6.1380892e-006},
	{1.7970890e-005,0.013210320,7.9044242},
	{0.24496661,0.35935348,0.00048348811},
	{1.4423846,0.11005451,0.00014971466},
	{4.9843683,0.0032259957,4.3885398e-006},
	{0.066035569,2.8269153,0.0017935552},
	{0.91681290,0.14679764,0.00019969865},
	{8.3828791e-006,0.0061622104,7.2132883},
	{1.0415542e-005,0.0076564257,7.6896152},
	{2.0830897e-005,0.015312705,7.7069721},
	{1.2186280,0.12299483,0.00016731817},
	{4.2717328,0.0026830356,3.6499146e-006},
	{1.8552844e-005,0.013638112,7.8706818},
	{1.7888882e-005,0.013150036,7.9088564},
	{6.0425558,0.0042663408,5.8037940e-006},
	{0.095684744,1.7633659,0.0012377980},
	{0.37270331,0.24051082,0.00032718311},
	{4.2811627,0.0026896517,3.6589131e-006},
	{2.2419759e-005,0.016480664,7.5710506},
	{5.5367651,0.0037238966,5.0658682e-006},
	{1.2999578e-005,0.0095559359,7.9568300},
	{0.26138082,0.33309069,0.00045312598},
	{0.16536319,0.74849093,0.00071623223},
	{5.4922848,0.0036806243,5.0070016e-006},
	{1.1620597e-005,0.0085422480,7.8512626},
	{1.7003304e-005,0.012499046,7.9510751},
	{1.4365871e-005,0.010560282,7.9976816},
	{2.7133137e-005,0.019945428,7.1135569},
	{1.4350712e-005,0.010549144,7.9975162},
	{6.0434572e-006,0.0044425158,6.1876431},
	{9.3352837e-006,0.0068623247,7.4737020},
	{1.1120899e-005,0.0081749223,7.7931399},
	{1.4180343e-005,0.010423901,7.9952469},
	{1.0125279e-005,0.0074430499,7.6389031},
	{1.1598652e-005,0.0085261129,7.8489532},
	{1.8860615e-005,0.013864366,7.8513198},
	{1.5771098e-005,0.011593251,7.9897728},
	{8.7478275e-006,0.0064304913,7.3217998},
	{1.2621787e-005,0.0092782238,7.9351158},
	{1.4758730e-005,0.010849070,7.9999752},
	{1.2984994e-005,0.0095452107,7.9560833},
	{2.7153847e-005,0.019960660,7.1114497},
	{0.17672819,0.66573793,0.00067017257},
	{1.7347596e-005,0.012752127,7.9359436},
	{1.7392023e-005,0.012784803,7.9338670},
	{3.6514553e-005,0.026841670,6.1635323},
	{1.6545395e-005,0.012162445,7.9684319},
	{8.1879716e-006,0.0060189436,7.1505203},
	{2.3768223e-005,0.017471896,7.4463792},
	{7.5700509e-006,0.0055647152,6.9273090},
	{8.9295554e-006,0.0065640789,7.3716717},
	{1.2830349e-005,0.0094315363,7.9477215},
	{1.1574674e-005,0.0085084904,7.8464065},
	{1.0280375e-005,0.0075570606,7.6666231},
	{2.7913366e-005,0.020518970,7.0339847},
	{4.3018299e-005,0.031622570,5.5586052},
	{8.2489878e-006,0.0060637882,7.1705427},
	{0.39752990,0.23276787,0.00031664985},
	{1.2911743e-005,0.0094913635,7.9522247},
	{1.7900862e-005,0.013158836,7.9082150},
	{2.5944380e-005,0.019071596,7.2335467},
	{2.7647980e-005,0.020323886,7.0611043},
	{2.2696111e-005,0.016683824,7.5460920},
	{1.6559223e-005,0.012172598,7.9679570},
	{2.5285270e-005,0.018587070,7.2991457},
	{2.1965961e-005,0.016147071,7.6112638},
	{5.8910259e-006,0.0043304665,6.0963922},
	{7.1189403e-005,0.052331030,3.6337996},
	{1.2683689e-005,0.0093237273,7.9390197},
	{0.28924543,0.30100223,0.00040947375},
	{0.29563257,0.29449910,0.00040062700},
	{5.7829652,0.0039755027,5.4081420e-006},
	{1.7548402e-005,0.012899750,7.9263473},
	{0.44273165,0.22021416,0.00029957228},
	{0.30714694,0.28345886,0.00038560820},
	{1.6410595e-005,0.012063349,7.9728999},
	{7.8735293e-006,0.0057877959,7.0417042},
	{1.7067265e-005,0.012546057,7.9483943},
	{1.4036458e-005,0.010318132,7.9927325},
	{0.092172250,1.8567178,0.0012849686},
	{0.71872962,0.16924454,0.00023023479},
	{1.3256313e-005,0.0097446600,7.9688158}
};
double gNormClusterProb[] = {19.080069,9.0199823,60.899944};

/**
 * @class EMTest
 * @brief cppunit class for testing conversion functions.
 */
class EMTest : public CppUnit::TestFixture
{
	CPPUNIT_TEST_SUITE(EMTest);
	CPPUNIT_TEST(testDefaultEMwithSeedGeneration);
	CPPUNIT_TEST(test1ClusterEM);
	CPPUNIT_TEST(test2ClusterEM);  
	CPPUNIT_TEST(testEMPosteriorWith3Clusters);
	CPPUNIT_TEST(testcalcProbability);
	CPPUNIT_TEST(testcalcGaussianWithExpTail);
	CPPUNIT_TEST(testgenerateEmpiricalEMSeeds);
	CPPUNIT_TEST(testdoEstep);
	CPPUNIT_TEST(testupdateMu);
	CPPUNIT_TEST(testadjustMu);
	CPPUNIT_TEST(testupdateWeight);
	CPPUNIT_TEST(testupdateSigma);
	CPPUNIT_TEST(testadjustSigma);
	CPPUNIT_TEST(testcomputeOptimalEstimates);
	CPPUNIT_TEST(testcomputeFLD);	   
	CPPUNIT_TEST_SUITE_END();
public:
	// three pretend data sets
	vector<vector<double> > data;

	// reference data
	vector<vector<int> > refCall;
	vector<vector<double> > refConf;

	// prior
	CEMPrior prior;
	// seed
	vector<CEMSeed> outseeds;

	void setUp();
	void testDefaultEMwithSeedGeneration();
	void test1ClusterEM();
	void test2ClusterEM();
	void testEMPosteriorWith3Clusters();
	void testcalcProbability();
	void testcalcGaussianWithExpTail();
	void testgenerateEmpiricalEMSeeds();
	void testdoEstep();
	void testupdateMu();
	void testadjustMu();
	void testupdateWeight();
	void testupdateSigma();
	void testadjustSigma();
	void testcomputeOptimalEstimates();
	void testcomputeFLD();
	void check(CEMEst* pResults, const int g);
	void checkSeeds();
	static bool closeEnough(double d1, double d2, int digits=6);
};

bool EMTest::closeEnough(double d1, double d2, int digits) 
{
	double diff = fabs(d1 - d2);
	if(diff < 1 / (pow((double)10, digits)))
		return true;
	return false;
}

// Registers the fixture into the 'registry'
CPPUNIT_TEST_SUITE_REGISTRATION(EMTest);

void EMTest::setUp()
{
	data.resize(GROUP);
	refCall.resize(GROUP+2);
	refConf.resize(GROUP+2);

	for(int s = 0; s < SAMPLES; s ++)
	{
		// number of clusters with data k + 1
		for(int g = 0; g < GROUP; g ++)
		{
			double r = K*(gTestData[g][0][s]-gTestData[g][1][s])/(gTestData[g][0][s]+gTestData[g][1][s]);
			double r1 = log(r+sqrt(r*r+1));
			data[g].push_back(r1/log(K+sqrt(K*K+1)));
			refCall[g].push_back(gRefCall[g][s]);
			refConf[g].push_back(gRefConf[g][s]);
		}
		for(int g = GROUP; g < GROUP+2; g++)
		{
			refCall[g].push_back(gRefCall[g][s]);
			refConf[g].push_back(gRefConf[g][s]);
		}
	}
	prior.set();
	vector<double> v(3);
	copy(&gPrior[0][0],&gPrior[0][3],v.begin());
	prior.setMu(v);
	copy(&gPrior[1][0],&gPrior[1][3],v.begin());
	prior.setMuSD(v);
	copy(&gPrior[2][0],&gPrior[2][3],v.begin());
	prior.setSigma(v);
	copy(&gPrior[3][0],&gPrior[3][3],v.begin());
	prior.setSigmaSD(v);
	prior.setStatus(true);
}

void EMTest::check(CEMEst* pResults, const int g)
{
	for(int s = 0; s < SAMPLES; s ++)
	{
		CPPUNIT_ASSERT(closeEnough(pResults->m_class[s], refCall[g][s]));
		CPPUNIT_ASSERT(closeEnough(pResults->m_confidence[s], refConf[g][s]));
	}
}

void EMTest::checkSeeds()
{
	for(int i = 0; i < outseeds.size(); i ++)
	{
		for(int g = 0; g < 3; g ++)
		{
			CPPUNIT_ASSERT(closeEnough(outseeds[i].getMu()[g], gSeed[0][g]));
			CPPUNIT_ASSERT(closeEnough(outseeds[i].getSigma()[g], gSeed[1][g]));
			CPPUNIT_ASSERT(closeEnough(outseeds[i].getWeight()[g], gSeed[2][g]));
			CPPUNIT_ASSERT(closeEnough(outseeds[i].getMinMu()[g], gSeed[3][g]));
			CPPUNIT_ASSERT(closeEnough(outseeds[i].getMinSigma()[g], gSeed[4][g]));
			CPPUNIT_ASSERT(closeEnough(outseeds[i].getMaxMu()[g], gSeed[5][g]));
			CPPUNIT_ASSERT(closeEnough(outseeds[i].getMaxSigma()[g], gSeed[6][g]));
		}
	}
}

void EMTest::testDefaultEMwithSeedGeneration() 
{ 
	CPrimeEM myEM;
	CEMEst* pResults;
	for(int g = 0; g < GROUP; g ++)
	{
		myEM.setData(data[g]);
		myEM.goPrimeEM();
		pResults = myEM.getEMEstimates();
		check(pResults,g);
	}
}

void EMTest::test1ClusterEM() 
{ 
	CPrimeEM myEM;
	CEMEst* pResults;
	CEMParam param;
	param.m_nClusters = 1;
	CEMSeed seed;
	seed.set(1);

	vector<double> v;
	v.push_back(-0.66f);
	seed.setMu(v);
	v.clear();
	v.push_back(-1.0f);
	seed.setMinMu(v);
	v.clear();
	v.push_back(-0.2f);
	seed.setMaxMu(v);
	v.clear();
	v.push_back(0.125f);
	seed.setSigma(v);
	v.clear();
	v.push_back(0.5f);
	seed.setMaxSigma(v);
	v.clear();
	v.push_back(0.05f);
	seed.setMinSigma(v);
	v.clear();
	v.push_back(1.0f);
	seed.setWeight(v);

	// set up
	myEM.setParam(param);
	myEM.setData(data[0]);
	myEM.EMEstimate(seed);
	pResults = myEM.getEMEstimates();
	check(pResults,3);
}

void EMTest::test2ClusterEM() 
{
	CPrimeEM myEM;
	CEMEst* pResults;
	CEMParam param;
	param.m_nClusters = 2;
	CEMSeed seed;
	seed.set(2);

	vector<double> v;
	v.push_back(-0.66f);
	v.push_back(0.0);
	seed.setMu(v);
	v.clear();
	v.push_back(-1.0f);
	v.push_back(-0.3f);
	seed.setMinMu(v);
	v.clear();
	v.push_back(-0.33f);
	v.push_back(0.3f);
	seed.setMaxMu(v);
	v.clear();
	v.push_back(0.125f);
	v.push_back(0.125f);
	seed.setSigma(v);
	v.clear();
	v.push_back(0.5f);
	v.push_back(0.5f);
	seed.setMaxSigma(v);
	v.clear();
	v.push_back(0.05f);
	v.push_back(0.05f);
	seed.setMinSigma(v);
	v.clear();
	v.push_back(0.5f);
	v.push_back(0.5f);
	seed.setWeight(v);
	//set up
	myEM.setParam(param);
	myEM.setData(data[1]);
	myEM.EMEstimate(seed);
	pResults = myEM.getEMEstimates();
	check(pResults,4);
}

void EMTest::testEMPosteriorWith3Clusters()
{
	CPrimeEM myEM;
	CEMEst* pResults;

	myEM.setData(data[2]);
	myEM.goPrimeEM();
	pResults = myEM.getEMEstimates();

	// check posterior model
	for(int c = 0; c < 3; c ++)
	{
		CPPUNIT_ASSERT(closeEnough(pResults->m_mu[c], gPosterior[0][c]));
		CPPUNIT_ASSERT(closeEnough(pResults->m_sigma[c], gPosterior[1][c]));
		CPPUNIT_ASSERT(closeEnough(pResults->m_weight[c], gPosterior[2][c]));
		CPPUNIT_ASSERT(closeEnough(pResults->m_minMu[c], gPosterior[3][c]));
		CPPUNIT_ASSERT(closeEnough(pResults->m_maxMu[c], gPosterior[4][c]));
		CPPUNIT_ASSERT(closeEnough(pResults->m_minSigma[c], gPosterior[5][c]));
		CPPUNIT_ASSERT(closeEnough(pResults->m_maxSigma[c], gPosterior[6][c]));
	}
	// check other stats
  // @todo check why this is...
  // on vercelli pResults->m_ML = 39.304615020751953 with -O1 or above 
  // with "-O0" vercilli gets the same results as ada6. 
  // gOthers[0]= 39.304625999999999
  // ada5      = 39.3046188
  // vercelli  = 39.304615020751953
  // the orginal code
	// CPPUNIT_ASSERT(closeEnough(pResults->m_ML, gOthers[0],5));
  // printf("pResults->m_ML = %f  (ref=%f)\n",pResults->m_ML,gOthers[0]);
  CPPUNIT_ASSERT(closeEnough(pResults->m_ML,   gOthers[0],6));
	CPPUNIT_ASSERT(closeEnough(pResults->m_FLD1, gOthers[1],5));
	CPPUNIT_ASSERT(closeEnough(pResults->m_FLD2, gOthers[2],5));
	CPPUNIT_ASSERT(closeEnough(pResults->m_sep,  gOthers[3],5));
}

void EMTest::testcalcProbability()
{
	double x = gProbInputs[0];
	double mu = gProbInputs[1];
	double sigma = gProbInputs[2];
	double weight = gProbInputs[3];
	double zcut = gProbInputs[4];
	bool useWeight = (gProbInputs[5] == 0);
  // @todo: again, this was copied from the output.
  // converting to doubles changed the answer from 30.528553 to 30.528552
  double calcp=calcProbability(x,mu,sigma,weight,zcut,useWeight);
  // printf("calcp=%f  (%f)\n",calcp,gProbOutput);
	CPPUNIT_ASSERT(closeEnough(calcp,gProbOutput,6));
}

void EMTest::testcalcGaussianWithExpTail()
{
	CPPUNIT_ASSERT(closeEnough(calcGaussianWithExpTail(gGassianWithTail[0],gGassianWithTail[1],gGassianWithTail[2]),gGassianReturn));
}

void EMTest::testgenerateEmpiricalEMSeeds()
{
	CPrimeEM myEM;
	myEM.setData(data[2]);
	myEM.setPrior(prior);
	myEM.generateEmpiricalEMSeeds();
	outseeds = myEM.getSeeds();
	checkSeeds();
}
	
void EMTest::testdoEstep()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);

	CEMParam param;
	param.m_nClusters = 3;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);
	myEM.doEstep();

	// test the probabilities of each data point and normalized prob of each cluster
	for(int g = 0; g < GROUP; g ++)
	{
		for(int i = 0; i< SAMPLES; i ++)
		{
			CPPUNIT_ASSERT(closeEnough(myEM.m_p[i][g],gProb[i][g],4));
		}
		CPPUNIT_ASSERT(closeEnough(myEM.m_clusterNP[g],gNormClusterProb[g],4));
	}
}

void EMTest::testupdateMu()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);

	CEMParam param;
	param.m_nClusters = 3;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);

	// calculate initial probs
	myEM.doEstep();

	// starts from {-0.6600000,0.00000000,0.6600000}
	myEM.updateMu();
	// to 	
	double gmu[]={-0.57764363,-0.27485964,0.65778154};
	for(int i = 0; i< GROUP; i ++)
	{
		CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_mu[i],gmu[i],4));
	}
}

void EMTest::testadjustMu()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);

	CEMParam param;
	param.m_nClusters = 3;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);

	// calculate initial probs
	myEM.doEstep();

	// starts from {-0.6600000,0.00000000,0.6600000}
	myEM.updateMu();

	// starts from {-0.57764363,-0.27485964,0.65778154};
	myEM.adjustMu();
	// to
	double gmu[]={-0.57764363,-0.25000000,0.65778154};
	for(int i = 0; i< GROUP; i ++)
	{
		CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_mu[i],gmu[i],5));
	}
}

void EMTest::testupdateWeight()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);

	CEMParam param;
	param.m_nClusters = 3;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);

	// calculate initial probs
	myEM.doEstep();

	// starts from {-0.6600000,0.00000000,0.6600000}
	myEM.updateMu();

	// starts from {-0.57764363,-0.27485964,0.65778154};
	myEM.adjustMu();
	// starts from 
	// {0.450000,0.350000,0.200000};
	myEM.updateWeight();
	// to
	double gweight[] = {0.21438280,0.10134812,0.68426901};
	for(int i = 0; i< GROUP; i ++)
	{
		CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_weight[i],gweight[i],5));
	}
}

void EMTest::testupdateSigma()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);

	CEMParam param;
	param.m_nClusters = 3;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);

	// calculate initial probs
	myEM.doEstep();

	// starts from {-0.6600000,0.00000000,0.6600000}
	myEM.updateMu();

	// starts from {-0.57764363,-0.27485964,0.65778154};
	myEM.adjustMu();
	// starts from 
	// {0.450000,0.350000,0.200000};
	myEM.updateWeight();
	// starts from m_sigma=[3](0.12500000,0.12500000,0.12500000)
	myEM.updateSigma();
	// to 
	double gsigma[] = {0.20530988,0.12232611,0.047600299};
	for(int i = 0; i< GROUP; i ++)
	{
		CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_sigma[i],gsigma[i],5));
	}
}

void EMTest::testadjustSigma()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);

	CEMParam param;
	param.m_nClusters = 3;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);
	myEM.m_pSeed = &seed;

	// calculate initial probs
	myEM.doEstep();

	// starts from {-0.6600000,0.00000000,0.6600000}
	myEM.updateMu();

	// starts from {-0.57764363,-0.27485964,0.65778154};
	myEM.adjustMu();
	// starts from 
	// {0.450000,0.350000,0.200000};
	myEM.updateWeight();
	// starts from m_sigma=[3](0.12500000,0.12500000,0.12500000)
	myEM.updateSigma();
	// starts from {0.20530988,0.12232611,0.047600299}
	myEM.adjustSigma();
	// to
	double gsigma[] = {0.20530988,0.12232611,0.05000000};
	for(int i = 0; i< GROUP; i ++)
	{
		CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_sigma[i],gsigma[i],5));
	}
}

void EMTest::testcomputeOptimalEstimates()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);
	CEMParam param;
	param.m_nClusters = 3;

	// don't ask EM to coverage
	param.m_nMaxIter = 0;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);
	myEM.m_pSeed = &seed;

	// one step EM
	myEM.doEstep(true);

	// get optimal estimates
	myEM.computeOptimalEstimates();

	int gcalls[] = {0,2,0,-1,0,0,2,2,0,-1,-1,2,1,-1,0,0,0,-1,2,0,0,2,-1,-1,2,0,2,0,-1,-1,2,0,0,0,0,0,0,0,0,0,0,0,
			0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1,0,0,0,0,0,0,0,0,0,0,0,-1,-1,2,0,-1,-1,0,0,0,0,-1,-1,0};
	double gconfs[] = {0.001523036,0.004287965,0.001250545,0.21240127,0.001294957,0.001388397,0.004028716,0.003197967,
			0.001788416,0.42352676,0.08667279,0.005363928,0.032938771,0.15649047,0.001887273,0.001389944,0.002355253,
			0.10860845,0.006878324,0.001892464,0.0017743,0.003505756,0.065705597,0.41056785,0.006856641,0.002721782,
			0.004343187,0.0012536,0.45721707,0.20012775,0.004421104,0.001274263,0.001630746,0.00132161,0.003966388,
			0.001320418,0.003281183,0.001597019,0.001309597,0.001308016,0.001434711,0.001275511,0.001949993,0.001461726,
			0.001762209,0.001249218,0.001354337,0.001253301,0.003972336,0.22917998,0.001684596,0.001691744,0.00692931,
			0.001563489,0.001962016,0.00305604,0.002240159,0.001706588,0.001250785,0.00127682,0.001409886,0.004190936,
			0.009270554,0.001938015,0.38795012,0.001251984,0.001776384,0.003633594,0.004114087,0.002788686,0.001565395,
			0.003454086,0.002613921,0.003421645,0.0219457,0.001249462,0.50649649,0.51529223,0.003924767,0.001717179,
			0.35122657,0.49662569,0.001544603,0.002095308,0.001640576,0.001298297,0.061109096,0.20986737,0.001260352};
	for(int i = 0; i< SAMPLES; i ++)
	{
		CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_class[i],gcalls[i]));
		CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_confidence[i],gconfs[i]));
	}
}

void EMTest::testcomputeFLD()
{
	CEMSeed seed;
	double mu[] = {-0.660000,0.000000,0.660000};
	seed.setMu(mu);
	double minMu[] = {-1.000000,-0.250000,0.300000};
	seed.setMinMu(minMu);
	double maxMu[] = {-0.300000,0.250000,1.000000};
	seed.setMaxMu(maxMu);
	double sigma[] = {0.125000,0.125000,0.125000};
	seed.setSigma(sigma);
	double minSigma[] = {0.0500000,0.0500000,0.0500000};
	seed.setMinSigma(minSigma);
	double maxSigma[] = {0.500000,0.500000,0.500000};
	seed.setMaxSigma(maxSigma);
	double weight[] = {0.450000,0.350000,0.200000};
	seed.setWeight(weight);

	CEMParam param;
	param.m_nClusters = 3;
	// don't ask EM to coverage
	param.m_nMaxIter = 0;

	CPrimeEM myEM;
	myEM.setParam(param);
	myEM.setData(data[2]);
	myEM.validate(seed);

	// validate seed
	myEM.validate(seed);
	myEM.m_pSeed = &seed;

	// one step EM
	myEM.doEstep(true);

	// get optimal estimates
	myEM.computeOptimalEstimates();

	// compute FLD
	myEM.computeFLD();

	// check diff of FLDs, since no empirical information is used
	// all FLDs should be the same
	double gFLD = 3.7335241;
	CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_FLD1,gFLD));
	CPPUNIT_ASSERT(closeEnough(myEM.m_est.m_FLD2,gFLD));
}

