////////////////////////////////////////////////////////////////
//
// Copyright (C) 2006 Affymetrix, Inc.
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

#include "mas5-stat/src/mathLib.h"
//
#include "mas5-stat/src/pTable.h"
//

float median(const std::vector<float> & v) {
	std::vector <float> u = v;
	long len = (long) u.size();
	if (len < 1) {
		return 0.0f;
	}
	std::sort(u.begin(), u.end());
	long half = len / 2;
	if (len % 2 == 1) {
		return u[half];
	}
	else {
		return (u[half - 1] + u[half]) / 2.0f;
	}
}

DI oneSidedSignRank3(const std::vector<float> & x, const double alpha) 
{
//	int len = (int) x.size(); // pk unused
	std::vector <float> newdiff = x;
	int n = 0;
	int i;
	for (i = 0; i < (int) x.size(); ++i) {
		if (x[i] != 0.0) {
			newdiff[n] = x[i];
			n++;
		}
	}
	newdiff.resize(n);
	std::vector <float> ranks(n);
	for (i=0; i<n; ++i) {
		ranks[i] = (float)(i + 1);
	}

	if (n == 0)
		return DI(0.5, 0);
	else {
		std::vector <DI> absdiff(n);
		for (i = 0; i < n; ++i) {
			absdiff[i].d = fabs(newdiff[i]);
			absdiff[i].i = i;
		}
		std::sort(absdiff.begin(), absdiff.end());
		int nTies = 0;
		std::vector <int> ties(n-1);
		for (i = 0; i < n - 1; ++i) {
			if (absdiff[i].d == absdiff[i+1].d) {
				ties[nTies] = i;
				nTies ++;
			}
		}
		ties.resize(nTies);
		int tieGroup = 0;
		double doubleVarMod = 0;
		if (nTies) {
			i = 0;
			while ( i < n - 1) {
				double initElement = absdiff[i].d;
				int tieGroupSize = 1;
				for (int j = i + 1; j < n; ++j) {
					if (absdiff[j].d == initElement) {
						tieGroupSize ++;
						if (j == n - 1) {
							i = j;
							tieGroup ++;
							for (int m = j - tieGroupSize + 1; m <= j; ++m) {
								ranks[m] = (2*j - tieGroupSize + 3) / 2.0f;
							}
							doubleVarMod += tieGroupSize *
								((double)tieGroupSize * tieGroupSize - 1);
							break;
						}
					}
					else {
						i = j;
						if (tieGroupSize > 1) {
							tieGroup ++;
							for (int m = j - tieGroupSize; m <= j - 1; ++m) {
								ranks[m] = (2 * j - tieGroupSize + 1) / 2.0f;
							}
							doubleVarMod += tieGroupSize *
								((double)tieGroupSize * tieGroupSize - 1);
						}
						break;
					}
				}
			}
		}
		std::vector <float> invr(n);
		for (i = 0; i < n; ++i) {
			invr[absdiff[i].i] = ranks[i];
		}

		double w = 0;
		for (i = 0; i < n; ++i) {
			if (newdiff[i] > 0)
				w += invr[i];
		}
		DI ans;
		if (n > 11) {
			double dw = w - ((double) n) * (n + 1) / 4.0;
			double denom2 = (((double)n)*(n+1)*(2*n+1) - 0.5*doubleVarMod)/24.0;
			if (denom2 <=0) {
				std::cerr << "denom2=" << denom2 << " dw=" << dw << std::endl;
				exit(1);
			}
			double z = dw / sqrt(denom2);
			ans.d = 1 - normalCDF(z);
		}
		else {
			if (nTies == 0) {
				int myCode = 0;
				for (i = 0; i < n; ++i) {
					if (newdiff[i] > 0) {
						myCode += 1 << ((int) invr[i] - 1);
					}
				}
				ans.d = fGetPValue(n-1, myCode);
			}
			else {
				int twoToN = 1 << n;
				std::vector<int> mask(n);
				for (i = 0; i < n; ++i) {
					mask[i] = 1 << i;
				}
				std::vector <double> posRanks(twoToN);
				for (i = 0; i < twoToN; ++i) {
					double sum = 0;
					for (int j = 0; j < n; ++j) {
						if (i & mask[j])
							sum += ranks[j];
					}
					posRanks[i] = sum;
				}
				double tail = 0;
				for (i = 0; i < twoToN; ++i) {
					if (posRanks[i] > w) {
						tail ++;
					}
					else if (posRanks[i] == w) {
						tail += 0.5;
					}
				}
				ans.d = tail / (double) twoToN;
			}
		}
		ans.i = (ans.d < alpha) ? 1 : 0;
		return ans;
	}
}

double normalCDF(const double x) {
	const double sqrt2 = 1.414213562373095048801689;
	double unAdjusted = 0.5 - 0.5 * erf( - x / sqrt2);
	if (unAdjusted > 1) {
		return 1;
	}
	else if (unAdjusted < 0) {
		return 0;
	}
	else {
		return unAdjusted;
	}
}

double erf(const double x) 
{
	static const int erfCoeffSize[] = {3, 7, 4};
	static const double erfA1[] = {-3.5609843701815385e-2,
		6.9963834886191355,    2.1979261618294152e1,
		2.4266795523053175e2};
    static const double erfB1[] = {1.5082797630407787e1,
		9.1164905404514901e1,  2.1505887586986120e2};
	static const double erfA2[] = {-1.368648573827167067e-7,
		5.641955174789739711e-1, 7.211758250883093659,
		4.316222722205673530e1,	1.529892850469404039e2,
		3.393208167343436870e2,	4.519189537118729422e2,
		3.004592610201616005e2};
    static const double erfB2[] = {1.278272731962942351e1,
		7.700015293522947295e1,  2.775854447439876434e2,
		6.389802644656311665e2,  9.313540948506096211e2,
		7.909509253278980272e2,  3.004592609569832933e2};
	static const double erfA3[] ={2.23192459734184686e-2,
		2.78661308609647788e-1, 2.26956593539686930e-1,
		4.94730910623250734e-2, 2.99610707703542174e-3};
    static const double erfB3[] = {1.98733201817135256,
		1.05167510706793207,  1.91308926107829841e-1,
		1.06209230528467918e-2};
	static const double erfXBounds[] = {0.46875, 4.0};
	static const double invSqrtPi = 0.56418958354775627928; // 1 / sqrt(pi)

	double absX = fabs(x);
	double xSquared = absX * absX;
	double temp = 0;


	if (absX <= erfXBounds[0]) {
		double num = erfA1[0] * xSquared;
		double den = xSquared;
		int last = erfCoeffSize[0];
		for (int i = 1; i < last; i++) {
			num = (num + erfA1[i]) * xSquared;
			den = (den + erfB1[i - 1]) * xSquared;
		}
		return x * (num + erfA1[last]) / (den + erfB1[last - 1]);
	}
	else {
		if (absX <= erfXBounds[1]) {
			double num = erfA2[0] * absX;
			double den = absX;
			int last = erfCoeffSize[1];
			for (int i = 1; i < last; i++) {
				num = (num + erfA2[i]) * absX;
				den = (den + erfB2[i - 1]) * absX;
			}
			temp = (num + erfA2[last]) / (den + erfB2[last - 1]);
		}
		else {
			double xInvSquared = 1.0 / xSquared;
			double num = erfA3[0] * xInvSquared;
			double den = xInvSquared;
			int last = erfCoeffSize[2];
			for (int i = 1; i < last; i++) {
				num = (num + erfA3[i]) * xInvSquared;
				den = (den + erfB3[i - 1]) * xInvSquared;
			}
			temp = xInvSquared * (num + erfA3[last]) / (den + erfB3[last - 1]);
			temp = (invSqrtPi - temp) / absX;
		}
		temp = 1 - exp(-xSquared) * temp;
		if (x > 0)      // in fact, we may use if (x > erfXBounds[0])
			return temp;
		else
			return -temp;
	}
}


double cFraction(double x, double a, double b) {
	double ai = 1, bi = 1, y = 1;
	double aPlus1 = a + 1, aPlusB = a + b;
//	double = aMinus1 = a - 1; // pk unused
	double z0 = 1 - x * aPlusB / aPlus1;
	double error = 100;
	for (int i = 1; i < MaxIterations; ++i) {
		double aPlus2I = a + 2 * i;
		double xModified = x / aPlus2I;
		double c = xModified * i * (b - i) / (aPlus2I - 1);
		double d = - xModified * (a + i) * (aPlusB + i) / (aPlus2I + 1);
		double y1 = y + ai * c;
		double y2 = y1 + y * d;
		double z1 = 0, z2 = 0;
		if (i == 1) {
			z1 = bi * c + z0;
			z2 = z1 + d * z0;
		}
		else {
			z1 = bi * c + 1;
			z2 = z1 + d;
		}
		ai = y1 / z2;
		bi = z1 / z2;
		double yold = y;
		y = y2 / z2;
		error = fabs(y - yold);
		if (error < RelativeBound * fabs(y)) {
			return y;
		}
	}
	return y;
}

double incompleteBeta(const double x, const double a, const double b)  {
	if (x == 0) {
		return 0;
	}
	else if (x == 1) {
		return 1;
	}
	else {
		double coeff = pow(x, a) * pow(1 - x, b) /
			exp(logGamma(a) + logGamma(b) - logGamma(a + b));
		if (x < (1.0 + a) / (2.0 + a + b)) { 
			return coeff * cFraction(x, a, b) / a;
		}
		else {
			return 1.0 - coeff * cFraction(1.0 - x, b, a) / b;
		}
	}
}
double tCDF(double t, double df) {
	double p;
	if (df == 1) {
		p = 0.5 + atan(t) / PI;
	}
	else {
		if (t == 0) {
			p = 0.5;
		}
		else {
			double x = df / (t*t + df);
			double y = 1 - incompleteBeta(x, df/2.0, 0.5);
			p = (1 - y) / 2.0;
			if (t > 0) {
				p = 1 - p;
			}
		}
	}
	if (p > 1)
		return 1;
	else if (p < 0)
		return 0;
	else
		return p;
}

double tCDFinversed(double p, double df) {
	double t = 0;
	if (p == 0) {
		t = -9999.0;
	}
	else if (p == 1) {
		t = 9999.0;
	}
	else if (p == 0.5) {
		t = 0;
	}
	else if (df == 1) {
		t = tan(PI * (p - 0.5));
	}	
	else if (p < 0.5) {
		t = - tCDFinversed(1 - p, df);
	}
	else {
		int maxIter = 120;
		double rB = 0.00001;
		double aB = 0.00005;
//		double coeff = exp(logGamma((df + 1) / 2.0) - logGamma(df / 2.0))
//				/ sqrt(df * PI); // pk unused
//		bool success = false; // pk unused
		double t0 = 1;
		double t1 = 2;
		double p0 = tCDF(t0, df);
		if (p0 == p) {
			return t0;
		}
		double p1 = tCDF(t1, df);
		if (p1 == p) {
			return t1;
		}
		for (int i = 0; i < maxIter; ++i) {
			if (p1 == p) {
				return t1;
			}
			else if (p0 == p) {
				return t0;
			}
			else if (p0 > p) {
				t1 = t0;
				p1 = p0;
				t0 = t0 / 2.0;
				p0 = tCDF(t0, df);
			}
			else if (p1 < p) {
				t0 = t1;
				p0 = p1;
				t1 = 2 * t1;
				p1 = tCDF(t1, df);
			}
			else {
				double midT = (t0 + t1) / 2.0;
				double midP = tCDF(midT, df);
				double dT = fabs(t1 - midT);
				if (midP == p || dT < aB || dT < rB * t1) {
					return midT;
				}
				else if (midP < p) {
					t0 = midT;
					p0 = tCDF(t0, df);
				}
				else {
					t1 = midT;
					p1 = tCDF(t1, df);
				}
			}
		}
	}
	return t;
}

double logGamma(const double x) {
	double a[8] = {5.7083835261e-03, -1.910444077728e-03,
		8.4171387781295e-04, -5.952379913043012e-04,
		7.93650793500350248e-04, -2.777777777777681622553e-03,
		8.333333333333333331554247e-02, 0.9189385332046727417803297};

	double c1[9] = {4.945235359296727046734888e0,
		2.018112620856775083915565e2, 2.290838373831346393026739e3,
		1.131967205903380828685045e4, 2.855724635671635335736389e4,
		3.848496228443793359990269e4, 2.637748787624195437963534e4,
		7.225813979700288197698961e3, -5.772156649015328605195174e-1};
	
	double d1[8] = {6.748212550303777196073036e1,
		1.113332393857199323513008e3, 7.738757056935398733233834e3,
		2.763987074403340708898585e4, 5.499310206226157329794414e4,
		6.161122180066002127833352e4, 3.635127591501940507276287e4,
		8.785536302431013170870835e3};
     
	double c2[9] = {4.974607845568932035012064e0,
		5.424138599891070494101986e2, 1.550693864978364947665077e4,
		1.847932904445632425417223e5, 1.088204769468828767498470e6,
		3.338152967987029735917223e6, 5.106661678927352456275255e6,
		3.074109054850539556250927e6, 4.227843350984671393993777e-1};

	double d2[8] = {1.830328399370592604055942e2,
		7.765049321445005871323047e3, 1.331903827966074194402448e5,
		1.136705821321969608938755e6, 5.267964117437946917577538e6,
		1.346701454311101692290052e7, 1.782736530353274213975932e7,
		9.533095591844353613395747e6};
		
	double c4[9] = {-1.474502166059939948905062e4,
		-2.426813369486704502836312e6, -1.214755574045093227939592e8,
		-2.663432449630976949898078e9, -2.940378956634553899906876e10,
		-1.702665737765398868392998e11, -4.926125793377430887588120e11,
		-5.606251856223951465078242e11, 1.791759469228055000094023e0};

	double d4[8] = {-2.690530175870899333379843e3,
		-6.393885654300092398984238e5, -4.135599930241388052042842e7,
		-1.120872109616147941376570e9, -1.488613728678813811542398e10,
		-1.016803586272438228077304e11, -3.417476345507377132798597e11,
		-4.463158187419713286462081e11};

   	double numerator = 0;
	double denominator = 1;
	int n = 8;

	if (x <= Eps) {
		return - log(x);
	}
	else if (x <= 0.5) {
		for (int i = 0; i < n; ++i) {
			numerator = x * numerator + c1[i];
			denominator = x * denominator + d1[i];
		}
		return - log(x) + x * (x * numerator / denominator + c1[n]);
	}
	else if (x <= 0.6796875) {
		double z = x - 1;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c2[i];
			denominator = z * denominator + d2[i];
		}
		return -log(x) + z * (z * numerator / denominator + c2[n]);
	}
	else if (x <= 1.5) {
		double z = x - 1;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c1[i];
			denominator = z * denominator + d1[i];
		}
		return z * (z * numerator / denominator + c1[n]);
	}
	else if (x <= 4) {
		double z = x - 2;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c2[i];
			denominator = z * denominator + d2[i];
		}
		return z * (z * numerator / denominator + c2[n]);
	}
	else if (x <= 12) {
		double z = x - 4;
		for (int i = 0; i < n; ++i) {
			numerator = z * numerator + c4[i];
			denominator = z * denominator + d4[i];
		}
		return z * numerator / denominator + c4[n];
	}
	else {
		double z = x * x;
		numerator = a[0];
		for (int i = 1; i <= 6; ++i) {
			numerator = numerator / z + a[i];
		}
		numerator /= x;
		return numerator + log(x) * (x - 0.5) - x + a[7];
	}	
}

double t_distribution_lookup(const int df) {
	
	static const double dfFrom1To30[] = {
		6.314, 2.920, 2.353, 2.132, 2.015,
		1.943, 1.895, 1.860, 1.833, 1.812,
		1.796, 1.782, 1.771, 1.761, 1.753,
		1.746, 1.740, 1.734, 1.729, 1.725,
		1.721, 1.717, 1.714, 1.711,	1.708,
		1.706, 1.703, 1.701, 1.699, 1.697
	};

	static const double df40 = 1.684;
	static const double df50 = 1.676;
	static const double df75 = 1.665;
	static const double df100 = 1.660;
	static const double df200 = 1.653;
	static const double df1000 = 1.646;

	if (df >=1 && df <= 30)
		return dfFrom1To30[df - 1];
	else if (df > 30 && df <= 40)
		return df40;
	else if (df > 40 && df <= 50)
		return df50;
	else if (df > 50 && df <= 75)
		return df75;
	else if (df > 75 && df <= 100)
		return df100;
	else if (df > 100 && df <= 200)
		return df200;
	else
		return df1000;
}


