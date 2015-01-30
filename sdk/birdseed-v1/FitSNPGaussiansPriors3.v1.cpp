////////////////////////////////////////////////////////////////
//
// Copyright (C) 2007 The Broad Institute and Affymetrix, Inc.
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
#include "birdseed-v1/FitSNPGaussiansPriors3.h"
//
#include "birdseed-dev/Matrix.h"
#include "broadutil/APTUtil.h"
#include "broadutil/BroadUtil.h"
#include "chipstream/SelfCreate.h"
#include "chipstream/SelfDoc.h"
#include "util/PgOptions.h"
//
#include <algorithm>
#include <cassert>
#include <iostream>
#include <limits>
#include <memory>
#include <vector>

using namespace std;
using namespace birdseed::v1;
using namespace birdseed::dev;

// See below for documentation of these parameters.
static double std_slope;
static double epsilon;
static double eps;
static double var_start;
static double cluster_distance_ratio_cutoff;
// TODO: Consider removing this, and removing tempmeandist, and lambda3
static double mean_dist2;
static double lambda3;
static double merged_cluster_threshold;
static double small_cluster_weight_threshold;
static double expected_wingspan_ratio;
static double low_hom_weight_fraction;
static double low_hom_sample_inflation;
static double starting_cluster_weight;
static size_t max_iter;
static double min_covar;
static double max_covar1;
static double max_covar2;
static double covar_floor_decay;
static double low_covar_threshold;
static double wing_length_delta_penalty;
static double small_cluster_penalty;
static double unbalanced_wingspan_penalty;
static double low_covar_penalty;
static double bic_weight;
// Not static because used by GenotypeCaller
double birdseed::v1::final_weight_min;
static double cluster_variance_regularization_factor;
static double var_mult;
static bool  allow_unlikely_clusters;
static double hom_hom_penalty;
static double mono_het_penalty;
static double two_cluster_low_observation_penalty_factor;
// This is not static because GenotypeCaller uses it.
double birdseed::v1::relative_distance_confidence_weight;
static double std_inflection_point;

static const bool special = true;

#ifdef _MSC_VER
#define log(value) log((double)value)
#endif

#ifdef __sun__
#define log(value) log((double)value)
#endif

// This class translates between the local parameters in this file and PgOptions world.
template<class T>
class OptionSpec
{
  public:
    T *targetVariable;
    PgOpt pgOpt;
    SelfDoc::Opt selfDocOpt;
  OptionSpec(const char *longName, T *targetVariable, T defaultValue, const char *defaultValueAsString, const char *helpDoc, PgOpt::PgOptType_t optType, bool makeSelfDocOpt):
        targetVariable(targetVariable),
        pgOpt()
    {
        pgOpt.m_shortName = "";
        pgOpt.m_longName = longName;
        pgOpt.m_type = optType;
        pgOpt.m_help = helpDoc;
        pgOpt.m_defaultValue = defaultValueAsString;
        *targetVariable = defaultValue;
        if (makeSelfDocOpt) {
            selfDocOpt.name = longName;
            selfDocOpt.value = defaultValueAsString;
            selfDocOpt.defaultVal = defaultValueAsString;
            selfDocOpt.minVal = "NA";
            selfDocOpt.maxVal = "NA";
            selfDocOpt.descript = helpDoc;
            switch (optType) {
            case PgOpt::BOOL_OPT:
                selfDocOpt.type = SelfDoc::Opt::Boolean;
                break;
            case PgOpt::INT_OPT:
                selfDocOpt.type = SelfDoc::Opt::Integer;
                break;
            case PgOpt::DOUBLE_OPT:
                selfDocOpt.type = SelfDoc::Opt::Double;
                break;
            default:
                throw BroadException("Unrecognized option type", __FILE__, __LINE__);
                break;
            }
        }
    }
};

#define DEFINE_OPTION(TYPE, VAR_NAME, DEFAULT_VALUE, HELP_DOC, OPT_TYPE, MAKE_SELF_DOC) \
    OptionSpec<TYPE>(#VAR_NAME, &VAR_NAME, DEFAULT_VALUE, #DEFAULT_VALUE, HELP_DOC, OPT_TYPE, MAKE_SELF_DOC)
#define DEFINE_DOUBLE_OPTION(VAR_NAME, DEFAULT_VALUE, HELP_DOC, MAKE_SELF_DOC) DEFINE_OPTION(double, VAR_NAME, DEFAULT_VALUE, HELP_DOC, PgOpt::DOUBLE_OPT, MAKE_SELF_DOC)

static OptionSpec<double> doubleOptions[] = {
    DEFINE_DOUBLE_OPTION(std_slope, 0.062, "expected slope of cluster standard deviation versus cluster mean intensity.", true),
    DEFINE_DOUBLE_OPTION(epsilon, 0.001, "tolerance at which to stop optimizing cluster locations.", true),
    DEFINE_DOUBLE_OPTION(eps, 0.00000000000000022204, "a very small number.", true),
    DEFINE_DOUBLE_OPTION(var_start, 1.1, "intialize the variances to be var_start times the expected.", true),
    DEFINE_DOUBLE_OPTION(cluster_distance_ratio_cutoff, 0.85, "the ratio of adjacent cluster means in each direction must exceed this value.", true),
    DEFINE_DOUBLE_OPTION(merged_cluster_threshold, .025, "if two cluster means get this close to each other, consider them merged, and stop trying EM.", true),
    DEFINE_DOUBLE_OPTION(small_cluster_weight_threshold, 0.01, "if k==3, any weight<small_cluster_weight_threshold, penalize ll with small_cluster_penalty.", true),
	DEFINE_DOUBLE_OPTION(low_hom_weight_fraction, 0.5, "Hom cluster should not have low weight.", true),
	DEFINE_DOUBLE_OPTION(low_hom_sample_inflation, 100, "Hom cluster should not have low weight.  Sample inflation factor.", true),
	DEFINE_DOUBLE_OPTION(starting_cluster_weight,0.05, "Starting weight for uninitialized clusters.", true),
    DEFINE_DOUBLE_OPTION(small_cluster_penalty, 10.0, "how much to penalize small clusters when k=3.", true),
    DEFINE_DOUBLE_OPTION(expected_wingspan_ratio, 1.15, "penalize ll if ratio of wing lengths is above this number.", true),
    DEFINE_DOUBLE_OPTION(unbalanced_wingspan_penalty, 5.0, "how much to penalize differences from expected_wingspan_ratio.", true),
    DEFINE_DOUBLE_OPTION(min_covar, -0.7, "don't let covar get lower than this.", true),
	DEFINE_DOUBLE_OPTION(max_covar1, 0.9, "covar1 not larger than this.", true),
	DEFINE_DOUBLE_OPTION(max_covar2, 0.95,"covar2 not larger than this.", true),
	DEFINE_DOUBLE_OPTION(covar_floor_decay,8,"Covariance decays over this iteration scale. Default ", true),
    DEFINE_DOUBLE_OPTION(low_covar_threshold, 1.0, "penalize covariances below this number.", true),
    DEFINE_DOUBLE_OPTION(low_covar_penalty, 15.0, "how much to penalize covar below low_covar_threshold.", true),
    DEFINE_DOUBLE_OPTION(wing_length_delta_penalty, 50.0, "how much to penalize differences from the prior.", true),
    DEFINE_DOUBLE_OPTION(bic_weight, 1.0, "how much to penalize higher-order k's.", true),
    DEFINE_DOUBLE_OPTION(final_weight_min, 0.333, "After calculating clusters, ensure all weights are >= this.", true),
    DEFINE_DOUBLE_OPTION(cluster_variance_regularization_factor, 1.0, "How much cluster variances are regularized to look like each other.", true),
    DEFINE_DOUBLE_OPTION(var_mult, 1.2, "Multiply the variance for missing clusters by this value squared.", true),
    DEFINE_DOUBLE_OPTION(hom_hom_penalty, 2.1, "Multiply the average distance between clusters and priors in 2-cluster model by this, when trying to fit clusters to hom priors.", true),

    DEFINE_DOUBLE_OPTION(mono_het_penalty, 10, "Multiply the distance squared between cluster and prior in 1-cluster model by this, when trying to fit single cluster to AB prior.", true),
    DEFINE_DOUBLE_OPTION(two_cluster_low_observation_penalty_factor, 1000000, "When penalizing an alignment of two-cluster model, use this factor to place a floor on penalty when number of prior observations for a prior is low.", true),

    DEFINE_DOUBLE_OPTION(relative_distance_confidence_weight, 0.8, "How much to weight confidence factor determined by comparing probability of best match vs. probability of second best match. Confidence factor determined by measuring distance of sample from cluster center is weighted by 1-this value.", true),
    DEFINE_DOUBLE_OPTION(std_inflection_point, 4.0, "Factor in determination of confidence based on distance of sample from cluster center.", true),

    DEFINE_DOUBLE_OPTION(mean_dist2, 1.25, "if two means get this close, penalize ll by 3k*log(n)/2;.", true),
    DEFINE_DOUBLE_OPTION(lambda3, 15.0, "how much to penalize differences from mean_dist2.", true)
};

static OptionSpec<size_t> size_tOptions[] = {
    DEFINE_OPTION(size_t, max_iter, 50, "if it reaches more than max_iter/k iterations, just stop.", PgOpt::INT_OPT, true)
};

static OptionSpec<bool> boolOptions[] = {
    DEFINE_OPTION(bool, allow_unlikely_clusters, true, "Allow 2-cluster model to match to AA && BB priors, and allow 1-cluster model to match AB prior.", PgOpt::BOOL_OPT, true)
};


static PgOpt *pgOptions[ARRAY_LEN(doubleOptions) + ARRAY_LEN(size_tOptions) + ARRAY_LEN(boolOptions) + 1];

// This class is here merely to get some static initializer code run, so that parameters get
// their default values, and the PgOpt array is initialized.
static class PgOptionCreator
{
  public:
    PgOptionCreator()
    {
        PgOpt **p = pgOptions;
        for (size_t i = 0; i < ARRAY_LEN(doubleOptions); ++i) {
            *p++ = &doubleOptions[i].pgOpt;
        }
        for (size_t i = 0; i < ARRAY_LEN(size_tOptions); ++i) {
            *p++ = &size_tOptions[i].pgOpt;
        }
        for (size_t i = 0; i < ARRAY_LEN(boolOptions); ++i) {
            *p++ = &boolOptions[i].pgOpt;
        }
        *p = NULL;
    }
} theOptionCreator;

PgOpt **birdseed::v1::getAlgorithmicOpts()
{
    return pgOptions;
}

void birdseed::v1::setAlgorithmicParameters(CONSTHACK PgOptions &opts)
{
    for (size_t i = 0; i < ARRAY_LEN(doubleOptions); ++i) {
        *doubleOptions[i].targetVariable = opts.getDouble(doubleOptions[i].pgOpt.m_longName);
    }
    for (size_t i = 0; i < ARRAY_LEN(size_tOptions); ++i) {
        *size_tOptions[i].targetVariable = opts.getInt(size_tOptions[i].pgOpt.m_longName);
    }
    for (size_t i = 0; i < ARRAY_LEN(boolOptions); ++i) {
        *boolOptions[i].targetVariable = opts.getBool(boolOptions[i].pgOpt.m_longName);
    }
}

void birdseed::v1::getSelfDocOptions(vector<SelfDoc::Opt> *opts)
{
    for (size_t i = 0; i < ARRAY_LEN(doubleOptions); ++i) {
        if (!doubleOptions[i].selfDocOpt.name.empty()) {
            opts->push_back(doubleOptions[i].selfDocOpt);
        }
    }
    for (size_t i = 0; i < ARRAY_LEN(size_tOptions); ++i) {
        if (!size_tOptions[i].selfDocOpt.name.empty()) {
            opts->push_back(size_tOptions[i].selfDocOpt);
        }
    }
    for (size_t i = 0; i < ARRAY_LEN(boolOptions); ++i) {
        if (!boolOptions[i].selfDocOpt.name.empty()) {
            opts->push_back(boolOptions[i].selfDocOpt);
        }
    }
}

void birdseed::v1::setSelfDocOptions(SelfDoc *doc, CONSTHACK std::map<std::string,std::string> &param)
{
    for (size_t i = 0; i < ARRAY_LEN(doubleOptions); ++i) {
        if (!doubleOptions[i].selfDocOpt.name.empty()) {
            SelfCreate::fillInValue(*doubleOptions[i].targetVariable, doubleOptions[i].selfDocOpt.name, param, *doc);
        }
    }
    for (size_t i = 0; i < ARRAY_LEN(size_tOptions); ++i) {
        if (!size_tOptions[i].selfDocOpt.name.empty()) {
            int value;
            SelfCreate::fillInValue(value, size_tOptions[i].selfDocOpt.name, param, *doc);
            *size_tOptions[i].targetVariable = value;
        }
    }
    for (size_t i = 0; i < ARRAY_LEN(boolOptions); ++i) {
        if (!boolOptions[i].selfDocOpt.name.empty()) {
            SelfCreate::fillInValue(*boolOptions[i].targetVariable, boolOptions[i].selfDocOpt.name, param, *doc);
        }
    }
}


static const size_t diploidNumGaussians[] = {1, 2, 3, 3};
static const size_t lenDiploidNumGaussians = ARRAY_SIZE(diploidNumGaussians);
static const size_t haploidNumGaussians[] = {1, 2};
static const size_t lenHaploidNumGaussians = ARRAY_SIZE(haploidNumGaussians);

static int verbosity = 0;

void birdseed::v1::setFitSNPGaussiansVerbosity(int verbose)
{
    verbosity = verbose;
}

static double computeStdInterceptPrior(const Priors &priors, size_t whichAllele)
{
    double numerator = 0.0;
    size_t denominator = 0;
    for (size_t i = 0; i < priors.getNumPriors(); ++i) {
        denominator += priors.getPrior(i).m_numObservations;
        numerator += priors.getPrior(i).m_numObservations *
            (sqrt(priors.getPrior(i).m_covarMatrix[whichAllele][whichAllele]) - std_slope * priors.getPrior(i).m_mean[whichAllele]);
    }
    return numerator / denominator;
}

static double determinant(const Matrix2x2 &mat)
{
    return mat[1][1]*mat[0][0] - mat[1][0]*mat[0][1];
}

static Matrix2x2 invert(const Matrix2x2 &inMatrix)
{
    Matrix2x2 outMatrix;
    double det = determinant(inMatrix);
    outMatrix[0][0] = inMatrix[1][1] / det;
    outMatrix[1][1] = inMatrix[0][0] / det;
    outMatrix[0][1] = -inMatrix[0][1] / det;
    outMatrix[1][0] = -inMatrix[1][0] / det;
    return outMatrix;
}

template<class MATRIX, class VECTOR>
static Matrix2x2 covariance(const MATRIX &matrix, const VECTOR &means)
{
    assert(means.size() == NUM_ALLELES);
    assert(matrix.numCols() == NUM_ALLELES);
    Matrix2x2 ret;
    for (size_t i = 0; i < NUM_ALLELES; ++i) {
        ret[i][i] = variance(matrix, means[i], i);
    }
    double acc = 0.0;
    for (size_t i = 0; i < matrix.numRows(); ++i) {
        acc += (matrix[i][0] - means[0]) * (matrix[i][1] - means[1]);
    }
    ret[0][1] = ret[1][0] = acc / (matrix.numRows() - 1);
    return ret;
}

template <class MATRIX>
static MATRIX pow(const MATRIX &matrix, typename MATRIX::value_type exponent)
{
    MATRIX ret;
    if (exponent == 2) {
        for (size_t row = 0; row < matrix.numRows(); ++row) {
            for (size_t col = 0; col < matrix.numCols(); ++col) {
                ret[row][col] = matrix[row][col] * matrix[row][col];
            }
        }
    }
    else {
        for (size_t row = 0; row < matrix.numRows(); ++row) {
            for (size_t col = 0; col < matrix.numCols(); ++col) {
                ret[row][col] = pow(matrix[row][col], exponent);
            }
        }
    }
    return ret;
}

template<class VECTOR>
static VECTOR powVector(const VECTOR &vec, typename VECTOR::value_type exponent)
{
    VECTOR ret(vec.size());
    if (exponent == 2) {
        for (size_t i = 0; i < vec.size(); ++i) {
            ret[i] = vec[i] * vec[i];
        }
    }
    else {
        for (size_t i = 0; i < vec.size(); ++i) {
            ret[i] = pow(vec[i], exponent);
        }
    }
    return ret;
}


template <class MATRIX>
static MATRIX sqrt(const MATRIX &matrix)
{
    MATRIX ret(matrix.numRows());
    for (size_t row = 0; row < matrix.numRows(); ++row) {
        for (size_t col = 0; col < matrix.numCols(); ++col) {
            ret[row][col] = sqrt(matrix[row][col]);
        }
    }
    return ret;
}

template <class VECTOR>
static VECTOR sqrtVector(const VECTOR &vec)
{
    VECTOR ret(vec.size());
    for (size_t i = 0; i < vec.size(); ++i) {
        ret[i] = sqrt(vec[i]);
    }
    return ret;
}

static inline double calculateBestPriorFit(const Matrix2x2 &B, const Priors &priors, size_t otherprior, size_t bestpriormatch)
{
    FixedVector<double, NUM_ALLELES> priorMeanDelta = priors.getPrior(otherprior).m_mean - priors.getPrior(bestpriormatch).m_mean;
    return exp(-((priorMeanDelta * B * priorMeanDelta)/2.0));
}

// Returns the index of the prior that the cluster matched
static size_t monomorphicCluster(Clusters *ret,
                                 const IntensityMatrix &intensities,
                                 const Priors &priors,
                                 const FixedVector<double, NUM_ALLELES> &stdInterceptPrior)
{
    int numSamples = intensities.numRows();

    // number of gaussians
    size_t k = 1;

    FixedVector<double, NUM_ALLELES> means;
    for (size_t i = 0; i < NUM_ALLELES; ++i) {
        means[i] = mean(intensities, i);
    }

    // check out concordance with the priors:
    size_t priorIndexA = Priors::AA_INDEX;
    size_t priorIndexB = Priors::BB_INDEX;
    if (!priors.isDiploid()) {
        priorIndexA = Priors::A_INDEX;
        priorIndexB = Priors::B_INDEX;

    }

     double delta_meanA = priors.getPrior(priorIndexA).m_mean[A_ALLELE_INDEX] - means[A_ALLELE_INDEX];
     double delta_meanB = priors.getPrior(priorIndexA).m_mean[B_ALLELE_INDEX] - means[B_ALLELE_INDEX];
     double aa_distance_sqrd = (delta_meanA * delta_meanA) + (delta_meanB * delta_meanB);
     delta_meanA = priors.getPrior(priorIndexB).m_mean[A_ALLELE_INDEX] - means[A_ALLELE_INDEX];
     delta_meanB = priors.getPrior(priorIndexB).m_mean[B_ALLELE_INDEX] - means[B_ALLELE_INDEX];
     double bb_distance_sqrd = (delta_meanA * delta_meanA) + (delta_meanB * delta_meanB);
     delta_meanA = priors.getPrior(Priors::AB_INDEX).m_mean[A_ALLELE_INDEX] - means[A_ALLELE_INDEX];
     delta_meanB = priors.getPrior(Priors::AB_INDEX).m_mean[B_ALLELE_INDEX] - means[B_ALLELE_INDEX];
     double ab_distance_sqrd = (delta_meanA * delta_meanA) + (delta_meanB * delta_meanB);
    size_t bestpriormatch;
    size_t otherprior;
    if (priors.isDiploid()) {
        bestpriormatch = (aa_distance_sqrd <= bb_distance_sqrd? Priors::AA_INDEX: Priors::BB_INDEX);
        double bestdistance_sqrd = (aa_distance_sqrd <= bb_distance_sqrd? aa_distance_sqrd: bb_distance_sqrd);
        otherprior = Priors::AB_INDEX;
        if (allow_unlikely_clusters) {
            ab_distance_sqrd *= mono_het_penalty;
            if (ab_distance_sqrd < bestdistance_sqrd) {
                bestpriormatch = Priors::AB_INDEX;
            }
        }
    } else if (aa_distance_sqrd <= bb_distance_sqrd) {
        bestpriormatch = Priors::A_INDEX;
        otherprior = Priors::B_INDEX;
    } else {
        bestpriormatch = Priors::B_INDEX;
        otherprior = Priors::A_INDEX;
    }

    // Found best match
	if (verbosity>2)
	{
		cout << "MonomorphicClusterDistances:\t";
		cout << aa_distance_sqrd << "\t";
		cout << ab_distance_sqrd << "\t";
		cout << bb_distance_sqrd << "\t";
		cout << endl;
		cout << "MonomorphicPriorMatch:\t";
		cout << bestpriormatch << "\t";
		cout << endl;
	}

    double vars[NUM_ALLELES];
    for (size_t i = 0; i < NUM_ALLELES; ++i) {
        vars[i] = numSamples/(numSamples) * variance(intensities, means[i], i);
    }
    Matrix2x2 temp = covariance(intensities, means);
    double covar2 = temp[0][1]/sqrt(temp[0][0]*temp[1][1]);
    double covar = min(max_covar2, max(covar2, min_covar));
    Matrix2x2 varmat;
    varmat[0][0] = vars[0];
    varmat[1][0] = varmat[0][1] = covar * sqrt(vars[0]*vars[1]);
    varmat[1][1] = vars[1];
    Matrix2x2 B = invert(varmat);
    double bestpriorfit;
    if (priors.isDiploid() && bestpriormatch == Priors::AB_INDEX) {
        bestpriorfit = max(calculateBestPriorFit(B, priors, Priors::AA_INDEX, Priors::AB_INDEX),
                           calculateBestPriorFit(B, priors, Priors::BB_INDEX, Priors::AB_INDEX));
    } else {
        bestpriorfit = calculateBestPriorFit(B, priors, otherprior, bestpriormatch);
    }
    double bestscale1 = means[0]/priors.getPrior(bestpriormatch).m_mean[0];
    double bestscale2 = means[1]/priors.getPrior(bestpriormatch).m_mean[1];

    double A = (1/sqrt(determinant(varmat)))/(2*pi);
    double sum_log_pxi_zj = 0.0;
    for (int i = 0; i < numSamples; ++i) {
        sum_log_pxi_zj += log(A * exp(-((intensities[i] - means) * B * (intensities[i] - means))/2.0));
    }

    // fill in the missing gaussians
    FixedVector<double, NUM_ALLELES> stdIntercept;
    for (size_t d = 0; d < NUM_ALLELES; ++d) {
        stdIntercept[d] = sqrt(vars[d]) - std_slope*means[d];
    }

    ret->means.resize(priors.getNumPriors());
    for (size_t i = 0; i < priors.getNumPriors(); ++i) {
        ret->means[i][0] = priors.getPrior(i).m_mean[0] * bestscale1;
        ret->means[i][1] = priors.getPrior(i).m_mean[1] * bestscale2;
    }

    double var_mult_sq = var_mult * var_mult;
    ret->vars.resize(priors.getNumPriors());
    for (size_t row = 0; row < priors.getNumPriors(); ++row) {
        for (size_t col = 0; col < NUM_ALLELES; ++col) {
          ret->vars[row][col] = 
            (stdIntercept[col] + std_slope*ret->means[row][col]) * 
            (stdIntercept[col] + std_slope*ret->means[row][col]);

            if (row != bestpriormatch) {
                ret->vars[row][col] *= var_mult_sq;
            }
        }
    }
    ret->weights.resize(priors.getNumPriors());
    ret->weights.setAllElements(starting_cluster_weight);
    if (priors.getNumPriors() == MAX_NUM_CLUSTERS) {
        ret->weights[bestpriormatch] = 1-2*starting_cluster_weight;
    } else {
        ret->weights[bestpriormatch] = 1-starting_cluster_weight;
    }

	if (verbosity>2)
	{
		cout << "MonomorphicClusterLogL:\t";
		cout << sum_log_pxi_zj << "\t";
		cout << max(0.0,low_covar_threshold-covar2) << "\t";
		cout << bestpriorfit << "\t";
		cout << 0.5*(2*k+2+k)*log(static_cast<double>(numSamples)) << "\t";
		cout << endl;
	}
    ret->log_likelihood = sum_log_pxi_zj
        - low_covar_penalty*log(static_cast<double>(numSamples))*max(0.0, low_covar_threshold-covar2)
        - wing_length_delta_penalty*bestpriorfit
        - bic_weight*0.5*(2*k+2+k)*log(static_cast<double>(numSamples));
    ret->covar = covar;

    return bestpriormatch;
}


static void initialize(Clusters *ret,
                       size_t k,
                       size_t count,
                       const IntensityMatrix &intensities,
                       const Priors &priors,
                       const FixedVector<double, NUM_ALLELES> &stdInterceptPrior,
                       const VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> &oldavgvars,
                       size_t bestpriormatch)
{
    assert(k > 1);
    if (k == 2) {
        ret->covar = 0;
        ret->weights.resize(k);
        ret->weights.setAllElements(1.0/k);
        ret->means.resize(2);
        for (size_t col = 0; col < NUM_ALLELES; ++col) {
            ret->means[0][col] = priors.getAPrior().m_mean[col];
            ret->means[1][col] = priors.getBPrior().m_mean[col];
        }
        Matrix2x2 temp;
        temp[0][0] = temp[1][0] = stdInterceptPrior[0];
        temp[0][1] = temp[1][1] = stdInterceptPrior[1];

        ret->vars.resize(2);
        ret->vars = pow(temp + ret->means * std_slope, 2) * var_start;
        if (special && count == 1) {
            for (size_t row = 0; row < k; ++row) {
                for (size_t col = 0; col < NUM_ALLELES; ++col) {
                    ret->vars[row][col] =
                        0.666 * oldavgvars[bestpriormatch][col];
                }
            }
        }
    } else if (k == 3) {
        ret->weights.resize(k);
        ret->weights.setAllElements(1.0/k);
        if (special && count == 2) {
            ret->means.resize(3);
            for (size_t row = 0; row < MAX_NUM_CLUSTERS; ++row) {
                ret->means[row] = priors.getPrior(row).m_mean;
            }
            ret->vars.resize(3);
            for (size_t row = 0; row < MAX_NUM_CLUSTERS; ++row) {
                ret->vars[row][0] = priors.getPrior(row).m_covarMatrix[0][0];
                ret->vars[row][1] = priors.getPrior(row).m_covarMatrix[1][1];
            }
            double acc = 0;
            for (size_t i = 0; i < MAX_NUM_CLUSTERS; ++i) {
                acc += priors.getPrior(i).m_numObservations * priors.getPrior(i).m_covarMatrix[1][0] /
                    sqrt(priors.getPrior(i).m_covarMatrix[0][0] * priors.getPrior(i).m_covarMatrix[1][1]);
            }
            ret->covar = max(min_covar,min(max_covar1, acc/priors.numObservations()));
        } else if (special) {
            ret->means.resize(3);
            ret->means[0][0] = minInColumn(intensities, 0);
            ret->means[0][1] = maxInColumn(intensities, 1);
            ret->means[2][0] = maxInColumn(intensities, 0);
            ret->means[2][1] = minInColumn(intensities, 1);
            ret->means[1] = (ret->means[0] + ret->means[2])/2;

            ret->vars.resize(3);
            for (size_t row = 0; row < ret->means.numRows(); ++row) {
                for (size_t col = 0; col < ret->means.numCols(); ++col) {
                  ret->vars[row][col] = 
                      ((ret->means[row][col] * std_slope + stdInterceptPrior[col]) * 
                       (ret->means[row][col] * std_slope + stdInterceptPrior[col])) * 
                    var_start;
                }
            }
            for (size_t col = 0; col < ret->vars.numCols(); ++col) {
                ret->vars[1][col] /= 20;
            }
            ret->covar = 0.0;
        } else {
            ret->means.resize(3);
            ret->vars.resize(3);
            for (size_t row = 0; row < MAX_NUM_CLUSTERS; ++row) {
                ret->means[row] = priors.getPrior(row).m_mean;
                ret->vars[row][0] = priors.getPrior(row).m_covarMatrix[0][0];
                ret->vars[row][1] = priors.getPrior(row).m_covarMatrix[1][1];
            }
            ret->covar = 0.0;
        }
    } else {
        assert(false);
    }
}

// k=2 case
// Not called for male chrX
static void fillInMissingGaussians(Clusters *clusters, const Priors &priors, size_t numSamples)
{
    // fill in missing Gaussians:
    // check out concordance with the priors:
    // If all else fails, consider clusters to be AA and AB
    size_t bestpriormatch[2] = {0, 1};
    double bestdistance = INFINITY;
    double bestshift1 = 0;
    double bestshift2 = 0;
    double bestscale1 = 1;
    double bestscale2 = 1;

    // Just syntactic sugar
    VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> &means = clusters->means;

    for (size_t firstCluster = Priors::AA_INDEX; firstCluster < Priors::BB_INDEX; ++firstCluster) {
        size_t maxSecondCluster;
        maxSecondCluster = (allow_unlikely_clusters? Priors::BB_INDEX: firstCluster + 1);
        for (size_t secondCluster = firstCluster + 1; secondCluster <= maxSecondCluster; ++secondCluster) {

            // Hom cluster should not have low weight
		if (verbosity>2)
		{
			cout << "2fillCluster:\t";
			cout << firstCluster << "\t";
			cout << clusters->weights[0] << "\t";
			cout << low_hom_weight_fraction*(numSamples/(numSamples+low_hom_sample_inflation)) << "\t";
			cout << secondCluster << "\t";
			cout << clusters->weights[1] << "\t";
			cout << endl;
		}
            if (firstCluster == Priors::AA_INDEX && clusters->weights[0] < low_hom_weight_fraction*numSamples/(numSamples+low_hom_sample_inflation)) {
                continue;
            } else if (secondCluster == Priors::BB_INDEX && clusters->weights[1] < low_hom_weight_fraction*numSamples/(numSamples+low_hom_sample_inflation)) {
                continue;
            }

            Matrix2x2 tempmat;
            tempmat[0][1] = tempmat[1][1] = 1;
            tempmat[0][0] = priors.getPrior(firstCluster).m_mean[0];
            tempmat[1][0] = priors.getPrior(secondCluster).m_mean[0];

            FixedVector<double, 2> temp1 = invert(tempmat) * means.columnSlice(0);

            tempmat[0][1] = tempmat[1][1] = 1;
            tempmat[0][0] = priors.getPrior(firstCluster).m_mean[1];
            tempmat[1][0] = priors.getPrior(secondCluster).m_mean[1];

            FixedVector<double, 2> temp2 = invert(tempmat) * means.columnSlice(1);

            FixedVector<double, 2> priorMeanSlice0;
            priorMeanSlice0[0] = priors.getPrior(firstCluster).m_mean[0];
            priorMeanSlice0[1] = priors.getPrior(secondCluster).m_mean[0];

            FixedVector<double, 2> priorMeanSlice1;
            priorMeanSlice1[0] = priors.getPrior(firstCluster).m_mean[1];
            priorMeanSlice1[1] = priors.getPrior(secondCluster).m_mean[1];

            double distance = meanVector(sqrtVector(powVector(priorMeanSlice0 - means.columnSlice(0), 2) +
                                                    powVector(priorMeanSlice1 - means.columnSlice(1), 2)));

            // Penalize this alignment if the hom cluster has few observations
            if (firstCluster == Priors::AA_INDEX && priors.getPrior(Priors::AA_INDEX).m_numObservations < priors.numObservations()) {
                distance *= sqrt((two_cluster_low_observation_penalty_factor + priors.numObservations())/(two_cluster_low_observation_penalty_factor + priors.getPrior(Priors::AA_INDEX).m_numObservations));
            }
            // There is potentially a double penalty if for the AA-BB alignment, if both have few observations.
            if (secondCluster == Priors::BB_INDEX && priors.getPrior(Priors::BB_INDEX).m_numObservations < priors.numObservations()) {
                distance *= sqrt((two_cluster_low_observation_penalty_factor + priors.numObservations())/(two_cluster_low_observation_penalty_factor + priors.getPrior(Priors::BB_INDEX).m_numObservations));
            }
            if (firstCluster == Priors::AA_INDEX && secondCluster == Priors::BB_INDEX) {
                distance *= hom_hom_penalty;
            }
		if (verbosity>2)
		{
			// decision made on distance
			cout << "2fillClusterDistance:\t";
			cout << distance << "\t";
			cout << temp1[0] << "\t" << temp2[0] << "\t";
			cout << temp1[1] << "\t" << temp2[1] << "\t";
			cout << endl;
		}
            if (distance<bestdistance) {
                bestdistance = distance;
                bestpriormatch[0] = firstCluster;
                bestpriormatch[1] = secondCluster;
                bestscale1 = temp1[0];
                bestscale2 = temp2[0];
                bestshift1 = temp1[1];
                bestshift2 = temp2[1];
            }
        }
    }

    VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> oldMeans(clusters->means);
    VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> oldVars(clusters->vars);
    VarVectorWithReservedLength<double, MAX_NUM_CLUSTERS> oldWeights(clusters->weights);
    // fill in the missing gaussians:
    FixedVector<double, NUM_ALLELES> stdIntercept;
    for (size_t d = 0; d < NUM_ALLELES; ++d) {
        stdIntercept[d] = meanVector(sqrtVector(clusters->vars.columnSlice(d)) -
                              clusters->means.columnSlice(d) * std_slope);
    }
    clusters->means.resize(MAX_NUM_CLUSTERS);
    for (size_t i = 0; i < MAX_NUM_CLUSTERS; ++i) {
        clusters->means[i][0] = bestshift1 + priors.getPrior(i).m_mean[0] * bestscale1;
        clusters->means[i][1] = bestshift2 + priors.getPrior(i).m_mean[1] * bestscale2;
    }
    for (size_t i = 0; i < 2; ++i) {
        clusters->means[bestpriormatch[i]] = oldMeans[i];
    }
    clusters->means = maxMatrixElementwise(clusters->means, 50);
    clusters->vars.resize(MAX_NUM_CLUSTERS);
    double var_mult_sq = var_mult * var_mult;
    for (size_t row = 0; row < MAX_NUM_CLUSTERS; ++row) {
        clusters->vars[row] = powVector(stdIntercept + clusters->means[row]*std_slope, 2) * var_mult_sq;
    }
    for (size_t i = 0; i < 2; ++i) {
        clusters->vars[bestpriormatch[i]] = oldVars[i];
    }
    clusters->weights.resize(MAX_NUM_CLUSTERS);
    clusters->weights.setAllElements(starting_cluster_weight);
    for (size_t i = 0; i < 2; ++i) {
        clusters->weights[bestpriormatch[i]] = min(1 - starting_cluster_weight * (MAX_NUM_CLUSTERS-1),
                                                     max(starting_cluster_weight,
                                                         oldWeights[i] *
                                                         (1 - starting_cluster_weight * (MAX_NUM_CLUSTERS-2))));
    }
}

void birdseed::v1::calculatePXI_ZJ(PXI_ZJMatrix *pxi_zj,
                               const Clusters &clusters,
                               const IntensityMatrix &intensities,
                               size_t k)
{
    assert(pxi_zj->numRows() == intensities.numRows());
    assert(pxi_zj->numCols() == k);
    for (size_t j = 0; j < k; ++j) {
        Matrix2x2 varmat;
        varmat[0][0] = clusters.vars[j][0];
        varmat[0][1] = clusters.covar * sqrt(clusters.vars[j][0] * clusters.vars[j][1]);
        varmat[1][0] = varmat[0][1];
        varmat[1][1] = clusters.vars[j][1];
        double A = clusters.weights[j] * (1/pow(determinant(varmat), .5)) / (2*pi);
        Matrix2x2 B = invert(varmat);
        for (size_t i = 0; i < intensities.numRows(); ++i) {
            FixedVector<double, NUM_ALLELES> meanDelta = intensities[i] - clusters.means[j];
            (*pxi_zj)[i][j] = A * exp(-((meanDelta * B * meanDelta)/2));
        }
    }
}

static void estimation(Clusters *clusters, const IntensityMatrix &intensities, size_t k)
{
    calculatePXI_ZJ(&clusters->pxi_zj, *clusters, intensities, k);
    double theSum = 0.0;
    for (size_t row = 0; row < clusters->pxi_zj.numRows(); ++row) {
        theSum += log(sumRow(clusters->pxi_zj, row));

    }
    clusters->log_likelihood = theSum;
}

static void finishEMLoop(Clusters *clusters, const Priors &priors, size_t k, size_t numSamples)
{
    // Just syntactic sugar
    VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> &means = clusters->means;

    vector<double> meandists(k-1);
    for (size_t i = 0; i < k-1; ++i) {
        meandists[i] = sqrt(sumVector(powVector(means[i+1] - means[i], 2)));
    }
    double avgmeandist = sumVector(meandists)/meandists.size();
    double worstmeandist = maxInVector(meandists)/avgmeandist;
    double tempmeandist = INFINITY;
    for (size_t i = 0; i < k-1; ++i) {
        double temp = (means[i+1][0]/means[i][0] + means[i][1]/means[i+1][1])/2;
        tempmeandist = min(temp, tempmeandist);
    }

	// simplified version of the logic for testing the second cluster
 	tempmeandist = (10.0/3) * (tempmeandist -
                                   (priors.isDiploid()?
                                    mean_dist2:
                                    1.0));

    double tempmeandist2 = INFINITY;
    for (size_t i = 0; i < k-1; ++i) {
        tempmeandist2 = min(tempmeandist2, means[i+1][0]/means[i][0]);
        tempmeandist2 = min(tempmeandist2, means[i][1]/means[i+1][1]);
    }

    if ((k == 2) && priors.isDiploid()) {
        fillInMissingGaussians(clusters, priors, numSamples);
    }

    double bestpriorfit = 0.0;
    for (size_t i = 0; i < priors.getNumPriors() - 1; ++i) {
        double delta_means0 = clusters->means[i][0]-clusters->means[i+1][0];
        double delta_means1 = clusters->means[i][1]-clusters->means[i+1][1];
        double mymeandist = sqrt((delta_means0 * delta_means0) + (delta_means1 * delta_means1));
        delta_means0 = priors.getPrior(i).m_mean[0] - priors.getPrior(i+1).m_mean[0];
        delta_means1 = priors.getPrior(i).m_mean[1] - priors.getPrior(i+1).m_mean[1];
        double priormeandist = sqrt((delta_means0 * delta_means0) + (delta_means1 * delta_means1));
        bestpriorfit += abs(mymeandist-priormeandist)/(mymeandist+priormeandist);
    }

	if (verbosity>2)
	{
		cout << "finishEMlogL:\t";
		cout << clusters->log_likelihood << "\t";
		cout << log(static_cast<double>(numSamples))*max(0.0,low_covar_threshold-clusters->covar) << "\t";
		cout << log(static_cast<double>(numSamples))*max(0.0,worstmeandist-expected_wingspan_ratio) << "\t";
		cout << min(1.0,bestpriorfit) << "\t";
		cout << 0.5*(2*k+2+k)*log(static_cast<double>(numSamples)) << "\t";
		cout << endl;
	}
    // log_likelihood is already sum(log(sum(pxi_zj,2))) ... %ll of data | gaussians
    clusters->log_likelihood = clusters->log_likelihood
        - low_covar_penalty*log(static_cast<double>(numSamples))*(max(0.0,low_covar_threshold-clusters->covar))  // ll of covariance
        - unbalanced_wingspan_penalty*k*log(static_cast<double>(numSamples))*(max(0.0,worstmeandist-expected_wingspan_ratio))  //ll of mean distribution
        - wing_length_delta_penalty*min(1.0, bestpriorfit)  // ll of prior
        - bic_weight*0.5*(2*k+2+k)*log(static_cast<double>(numSamples)); // ll of model (BIC criterion)
    if (tempmeandist < 0) { //make sure the distance between means is enough
        clusters->log_likelihood = clusters->log_likelihood+k*tempmeandist*lambda3*k*log(static_cast<double>(numSamples));
	if (verbosity>3)
	{
		cout << "tempmeandist<0:\t";
		cout << k*tempmeandist*k*log(static_cast<double>(numSamples)) << "\t";
		cout << endl;
	}
    }
    if (tempmeandist2<cluster_distance_ratio_cutoff) { // Make sure the order is correct (by a margin), or throw away:
        clusters->log_likelihood = -INFINITY;
	if (verbosity>3)
	{
		cout << "tempmeandist2<cluster_distance_ratio_cutoff:\t-INFINITY" << endl;
	}
    }
    if ((k==3) && (minInVector(clusters->weights)<small_cluster_weight_threshold)) {
        clusters->log_likelihood = clusters->log_likelihood
            - small_cluster_penalty*k*sqrt((small_cluster_weight_threshold-minInVector(clusters->weights)))*log(static_cast<double>(numSamples));
	if (verbosity>3)
	{
		cout << "Small_cluster_weight_threshold\t";
		cout << k*sqrt((small_cluster_weight_threshold-minInVector(clusters->weights)))*log(static_cast<double>(numSamples)) << "\t";
		cout << endl;
	}
    }
}

// false return means stop E-M loop.
static bool maximization(Clusters *clusters,
                         const IntensityMatrix &intensities,
                         const Priors &priors,
                         const FixedVector<double, NUM_ALLELES> &stdInterceptPrior,
                         size_t numSamples,
                         size_t k,
                         size_t iter)
{
    // Syntactic sugar
    VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> &means = clusters->means;

    VarVectorWithReservedLength<double, MAX_NUM_CLUSTERS> nj(clusters->pxi_zj.numCols());
    for (size_t col = 0; col < clusters->pxi_zj.numCols(); ++col) {
        nj[col] = sumColumn(clusters->pxi_zj, col);
        clusters->weights[col] = nj[col] / numSamples;
    }
    double minNJ = minInVector(nj);
    if ((minNJ < epsilon && iter > 10) || (minNJ < eps && iter > 1)) {
        clusters->log_likelihood = -INFINITY;
        return false;
    }
    for (size_t i = 0; i < nj.size(); ++i) {
        // Floor to be at least smallest non-zero double.
        nj[i] = max(nj[i], numeric_limits<double>::epsilon());
    }

    VarVector<double> covs(k);

    for (size_t j = 0; j < k; ++j) {
        // means(j,:) = sum(repmat(pxi_zj(:,j),1,2).*dataVector) / nj(j);
        for (size_t col = 0; col < NUM_ALLELES; ++col) {
            double meanSum = 0;
            for (size_t row = 0; row < intensities.numRows(); ++row) {
                meanSum += clusters->pxi_zj[row][j] * intensities[row][col];
            }
            means[j][col] = meanSum / nj[j];
        }

        // vars(j,:) = 1/nj(j).*sum(repmat(pxi_zj(:,j),1,2) .* ((dataVector-repmat(means(j,:),n,1)).^2));
        for (size_t col = 0; col < NUM_ALLELES; ++col) {
            double varSum = 0;
            for (size_t row = 0; row < intensities.numRows(); ++row) {
                varSum += clusters->pxi_zj[row][j] * ((intensities[row][col] - means[j][col]) * (intensities[row][col] - means[j][col]));
            }
            clusters->vars[j][col] = varSum/nj[j];
        }

        // covs(j) = 1/nj(j).* (pxi_zj(:,j) .* (dataVector(:,1)-repmat(means(j,1),n,1)))'*(dataVector(:,2)-repmat(means(j,2),n,1));
        double covSum = 0;
        for (size_t row = 0; row < intensities.numRows(); ++row) {
            covSum += clusters->pxi_zj[row][j] *
                (intensities[row][0] - means[j][0]) *
                (intensities[row][1] - means[j][1]);
        }
        covs[j] = covSum / nj[j];
    }
    // covar = min(0.95,max(min_covar,sum(nj'.*covs'./(eps+sqrt(vars(:,1).*vars(:,2))))/sum(nj)));
    clusters->covar = min(max_covar2,max(min_covar,
                                       sumVector(nj.multiplyElementwise(covs).divideElementwise
                                                 (sqrtVector(clusters->vars.columnSlice(0).multiplyElementwise(clusters->vars.columnSlice(1))) + eps))/
                                       sumVector(nj)));

    if ((iter<=covar_floor_decay) && (k==2 || k==3)) {
        clusters->covar = max(max_covar1-(iter-1.0)/covar_floor_decay, clusters->covar);
        clusters->log_likelihood = -INFINITY;
    }

    // make sure the means aren't too close
    // if(max((sqrt(sum((means(2:k,:)-means(1:k-1,:)).^2,2))'./sqrt(sum(means(2:k,:).^2))))<merged_cluster_threshold)
    double maxDistance = -INFINITY;
    for (size_t i = 1; i < k; ++i) {
        double distance = sqrt(((means[i][0] - means[i-1][0]) * (means[i][0] - means[i-1][0])) + ((means[i][1] - means[i-1][1]) * (means[i][1] - means[i-1][1]))) /
            sqrt((means[i][0] * means[i][0]) + (means[i][1] * means[i][1]));
        maxDistance = max(maxDistance, distance);
    }
    if (maxDistance < merged_cluster_threshold) {
        clusters->log_likelihood = -INFINITY;
        return false;
    }

    // regularize the variances
    for (size_t d = 0; d < 2; ++d) {
        // vars(:,d) = nj'./(nj'+np).*vars(:,d) + ...
        //        np./(nj'+np).*((repmat(stdInterceptPrior(d),k,1)+std_slope*means(:,d)).^2); %regularization term 1
        VarVectorWithReservedLength<double, MAX_NUM_CLUSTERS> temp(nj.size());
        temp = nj.divideElementwise(nj).multiplyElementwise(clusters->vars.columnSlice(d)); // regularization term 1
        for (size_t i = 0; i < temp.size(); ++i) {
            clusters->vars[i][d] = temp[i];
        }
    }

    // how much variances are regularized to look like each other
    double cluster_variance_regularization = cluster_variance_regularization_factor * intensities.numRows() * ((double(k-1)/k) * (double(k-1)/k));

    for (size_t d = 0; d < 2; ++d) {
        // expectedvars = (mean(sqrt(vars(:,d))-std_slope*means(:,d)) +std_slope*means(:,d)).^2;
        VarVectorWithReservedLength<double, MAX_NUM_CLUSTERS> expectedvars(clusters->vars.numRows());

        expectedvars = powVector(clusters->means.columnSlice(d)*std_slope +
                                 meanVector(sqrtVector(clusters->vars.columnSlice(d)) - clusters->means.columnSlice(d)*std_slope), 2);

        // vars(:,d) = nj'./(nj'+np+cluster_variance_regularization).*vars(:,d) + ...
        //    cluster_variance_regularization./(nj'+np+cluster_variance_regularization).*expectedvars; %regularization term 2
        VarVectorWithReservedLength<double, MAX_NUM_CLUSTERS> temp(clusters->vars.numRows());
        VarVectorWithReservedLength<double, MAX_NUM_CLUSTERS> temp2(nj.size());
        temp2 = nj + cluster_variance_regularization;
        temp = nj.divideElementwise(temp2).multiplyElementwise(clusters->vars.columnSlice(d)) +
            temp2.divisor(cluster_variance_regularization).multiplyElementwise(expectedvars);
        for (size_t i = 0; i < temp.size(); ++i) {
            clusters->vars[i][d] = temp[i];
        }
    }
    return true;
}


Clusters birdseed::v1::FitSNPGaussianPriors3(const IntensityMatrix &intensities, const Priors &priors)
{
    size_t numSamples = intensities.numRows();

    FixedVector<double, NUM_ALLELES> stdInterceptPrior;
    for (size_t i = 0; i < NUM_ALLELES; ++i) {
        stdInterceptPrior[i] = computeStdInterceptPrior(priors, i);
    }

    const size_t *numGaussians = (priors.isDiploid()? diploidNumGaussians: haploidNumGaussians);
    const size_t lenNumGaussians = (priors.isDiploid()? lenDiploidNumGaussians: lenHaploidNumGaussians);

    // It is required that we start with the 1 gaussian case
    assert(lenNumGaussians > 0);
    assert(numGaussians[0] == 1);
    Clusters bestResult(1,1, 0);
    size_t bestpriormatch = monomorphicCluster(&bestResult, intensities, priors, stdInterceptPrior);

    if (verbosity >= 2) {
        cout << "oneGaussian result: " << bestResult.tostring() << "\n";
    }

    // Use k=1 case to estimate the variances for k=2 case
    VarMatrixWithReservedRows<double, MAX_NUM_CLUSTERS, NUM_ALLELES> oldavgvars(bestResult.vars);


    for (size_t count = 1; count < lenNumGaussians; ++count) {
        size_t k = numGaussians[count];
        assert(k > 1);

        Clusters clusters(intensities.numRows(), k, count);
        initialize(&clusters, k, count, intensities, priors, stdInterceptPrior, oldavgvars, bestpriormatch);

        // EM loop:
        for (size_t iter = 1; ; ++iter) {
            double old_log_likelihood = clusters.log_likelihood;
            estimation(&clusters, intensities, k);

            if ((clusters.log_likelihood - old_log_likelihood < epsilon) || (iter>max_iter)) {
                if (verbosity >= 3) {
                    cout << "Before finishEMLoop " << "count: " << count << "; iter: " << iter << "; " << clusters.tostring() << "\n";
                }
                finishEMLoop(&clusters, priors, k, numSamples);
                if (verbosity >= 3) {
                    cout << "After finishEMLoop " << "count: " << count << "; iter: " << iter << "; " << clusters.tostring() << "\n";
                }
                break;
            }
            // Normalize each row of pxi_zj so that each value is divided by the sum of values in that row
            for (size_t row = 0; row < clusters.pxi_zj.numRows(); ++row) {
                double rowSum = sumRow(clusters.pxi_zj, row);
                for (size_t col = 0; col < clusters.pxi_zj.numCols(); ++col) {
                    if (rowSum > 0) {
                        clusters.pxi_zj[row][col] /= rowSum;
                    } else {
                        clusters.pxi_zj[row][col] = 1.0 / k;
                    }
                }
            }
            if (verbosity >= 3) {
                cout << "Before maximization " << "count: " << count << "; iter: " << iter << "; " << clusters.tostring() << "\n";
            }

            if (!maximization(&clusters, intensities, priors, stdInterceptPrior, numSamples, k, iter)) {
                break;
            }
            if (verbosity >= 3) {
                cout << "After maximization " << "count: " << count << "; iter: " << iter << "; " << clusters.tostring() << "\n";
            }
        }
        if (verbosity >= 2) {
            cout << "After iter loop count: " << count << ": " << clusters.tostring() << "\n";
        }
        if (clusters.log_likelihood > bestResult.log_likelihood) {
            bestResult = clusters;
        }
    }
    if (verbosity >= 2) {
        cout << "final result: " << bestResult.tostring() << "\n";
    }

    // This is no longer necessary so delete it.
    bestResult.pxi_zj.resizeRows(0);

    return bestResult;
}

double birdseed::v1::calculateSNPSpecificCorrectionFactor(const IntensityMatrix &intensities,
                                                      const Priors &priors)
{
    // This could be done in a for loop to avoid memory allocation
    double numerator = meanVector(sqrtVector(powVector(intensities.columnSlice(0), 2) + powVector(intensities.columnSlice(1), 2)));
    double denominator = 0.0;
    for (size_t i = 0; i < priors.getNumPriors(); ++i) {
        double mean0 = priors.getPrior(i).m_mean[0];
        double mean1 = priors.getPrior(i).m_mean[1];
        denominator += sqrt((mean0 * mean0) + (mean1 * mean1));
    }
    denominator /= priors.getNumPriors();
    return numerator / denominator;
}

double birdseed::v1::calculateConfidenceRelativeToCluster(const Matrix2x2 &B,
                                                      const IntensityMatrix::ROW &intensities,
                                                      const IntensityMatrix::ROW &clusterMean)
{
    static double exp_stdInflectionPoint = -1.0;
    if (exp_stdInflectionPoint < 0) {
        exp_stdInflectionPoint = exp(std_inflection_point);
    }
    IntensityMatrix::ROW delta = intensities - clusterMean;
    double stdDist = sqrt((delta * B * delta)/2);
    double confidenceRelativeToCluster = 1/(1+exp(std_inflection_point-stdDist)) -
        1/(1+exp_stdInflectionPoint);
    return confidenceRelativeToCluster;
}

Matrix2x2 birdseed::v1::makeBMatrix(const Matrix2x2::ROW &clusterVars, double covar)
{
    Matrix2x2 varmat;
    varmat[0][0] = clusterVars[0];
    varmat[1][0] = varmat[0][1] = covar * sqrt(clusterVars[0]*clusterVars[1]);
    varmat[1][1] = clusterVars[1];
    return invert(varmat);
}

/******************************************************************/
/**************************[END OF FitSNPGaussiansPriors3.cpp]*************************/
/******************************************************************/

/* Emacs configuration
 * Local Variables:
 * mode: C++
 * tab-width:4
 * End:
 */
