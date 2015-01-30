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
#include "copynumber/CNIntensityAdjustmentMethodHighPassFilter.h"
//
#include "copynumber/DataBlock.h"
#include "copynumber/Indexer.h"
//
#include "util/Verbose.h"
//
#include <new>
#include <utility>
//


using namespace std;

/**
 * @brief Supply a little how/what/why about the algorithms this
 * class performs and what parameters it takes.
 * @return SelfDoc
 */
SelfDoc CNIntensityAdjustmentMethodHighPassFilter::explainSelf() 
{
	CNIntensityAdjustmentMethodHighPassFilter obj;
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
std::vector<SelfDoc::Opt> CNIntensityAdjustmentMethodHighPassFilter::getDefaultDocOptions() 
{
	std::vector<SelfDoc::Opt> opts;

	// SelfDoc::Opt(name, type, value, default, min, max, description)
	SelfDoc::Opt DataBlockRows = {"data-block-rows", SelfDoc::Opt::Integer, "320", "320", "NA", "NA", "Data block rows."};
	opts.push_back(DataBlockRows);

	SelfDoc::Opt DataBlockCols = {"data-block-cols", SelfDoc::Opt::Integer, "2015", "2015", "NA", "NA", "Data block columns."};
	opts.push_back(DataBlockCols);

	SelfDoc::Opt MiniBlockRows = {"mini-block-rows", SelfDoc::Opt::Integer, "8", "8", "NA", "NA", "Mini block rows."};
	opts.push_back(MiniBlockRows);

	SelfDoc::Opt MiniBlockCols = {"mini-block-cols", SelfDoc::Opt::Integer, "8", "8", "NA", "NA", "Mini block columns."};
	opts.push_back(MiniBlockCols);

	SelfDoc::Opt GlobalSmoothWeight = {"global-smooth-weight", SelfDoc::Opt::Double, "256.0", "256.0", "NA", "NA", "Global smooth weight."};
	opts.push_back(GlobalSmoothWeight);

	SelfDoc::Opt LocalSmoothWeight = {"local-smooth-weight", SelfDoc::Opt::Double, "64.0", "64.0", "NA", "NA", "Local smooth weight."};
	opts.push_back(LocalSmoothWeight);

	SelfDoc::Opt Converged = {"converged", SelfDoc::Opt::Double, "0.0001", "0.0001", "NA", "NA", "Converged."};
	opts.push_back(Converged);

	SelfDoc::Opt UseSingleBlock = {"use-single-block", SelfDoc::Opt::Boolean, "false", "false", "NA", "NA", "Determines whether or not the entire chip is analyzed as a single block or not."};
	opts.push_back(UseSingleBlock);

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
SelfCreate* CNIntensityAdjustmentMethodHighPassFilter::newObject(std::map<std::string,std::string>& params) 
{
	SelfDoc doc = explainSelf();
	std::vector<SelfDoc::Opt> opts = getDefaultDocOptions();
	CNIntensityAdjustmentMethodHighPassFilter* pMethod = new CNIntensityAdjustmentMethodHighPassFilter();
	std::string strPrefix = getPrefix();

	pMethod->data_block_rows = setupIntParameter("data-block-rows", strPrefix, params, doc, opts);
	pMethod->data_block_cols = setupIntParameter("data-block-cols", strPrefix, params, doc, opts);
	pMethod->mini_block_rows = setupIntParameter("mini-block-rows", strPrefix, params, doc, opts);
	pMethod->mini_block_cols = setupIntParameter("mini-block-cols", strPrefix, params, doc, opts);
	pMethod->global_smooth_weight = setupDoubleParameter("global-smooth-weight", strPrefix, params, doc, opts);
	pMethod->local_smooth_weight = setupDoubleParameter("local-smooth-weight", strPrefix, params, doc, opts);
	pMethod->converged = setupDoubleParameter("converged", strPrefix, params, doc, opts);
	pMethod->m_bUseSingleBlock = setupBoolParameter("use-single-block", strPrefix, params, doc, opts);
	
	return pMethod;
}

/**
 * @brief Constructor
 */
CNIntensityAdjustmentMethodHighPassFilter::CNIntensityAdjustmentMethodHighPassFilter()
{
	data_block_cols = 0;
	data_block_rows = 0;
	mini_block_cols = 0;
	mini_block_rows = 0;
	global_smooth_weight = 0;
	local_smooth_weight = 0;
	converged = 0;
	m_bUseSingleBlock = false;
}

/**
 * @brief Destructor
 */
CNIntensityAdjustmentMethodHighPassFilter::~CNIntensityAdjustmentMethodHighPassFilter()
{
}

/**
 * @brief Run the calculations.
 */
void CNIntensityAdjustmentMethodHighPassFilter::run()
{
	Verbose::out(1, "CNIntensityAdjustmentMethodHighPassFilter::run(...) start");

	// TODO: Parameterize this.
	// This will currently be set up for both cyto2 full and cyto2 focused
	// arrays.  The algorithm was developed on cyto2 full arrays where the
	// copy number probes are laid out in 4 bands.  If the image is oriented
	// so that the logo is in the top left corner the bands are horizontal
	// spanning the entire image.  In order to save RAM the bands were
	// one at a time.  The function index_to_rc takes a probe ID along with
	// the width of the image to find the row and column orientation used
	// in this module.  This is a right-hand coordinate system.


/*     This is a set of additional comments added when this module was extended to be used on the
 *     snp6 and snp7 chip.  The original use on the 2.7M chip only analyzed the CN probes.  These
 *     probes were located on a set of 4 bands across the chip. The extension to snp6 and snp7
 *     chips treats the whole chip as a single block and takes all SNP and CN probes into analysis. 
 *
 *     The original meaning of the data_block_rows  and data_block_column parameter values, set on 
 *     the command line input, were the number of rows and columns in the blocks of CN probes.  The  
 *     blocks themselves were then defined by hard coding the row value at which the block started.
 *
 *     In order to do a single block analysis these two parameter values are now interpreted as the
 *     number of rows and columns in the entire chip, and the above row value of the starting position
 *     of the block is given as 0.
 *
 *     In order to maintain backwards compatibility a new module option was defined. It is 
 *     "--use-single-block".  It is only used when single block analysis is desired on the 2.7M
 *     chip.  Without this option being defined on the command line input, multiple block analysis
 *     is performed on 2.7M, single block is performed on snp6/snp7 and 310K. 
 */

	vector<int> data_block_starts; 
        // We now determine what type of chip and whether or not we do single or multiple block analysis.
        // Chip is 310K.  Always do single block analysis. 
	if (data_block_rows == data_block_cols) // assume we are running the focused array.
	{
		data_block_starts.push_back(0);  // first row of the image
	}
        if(m_pEngine->getOptBool("cyto2"))
	{
                // Chip is 2.7M and we want multiple block analysis 
                if(! m_bUseSingleBlock)
                {
		        // The 4 cyto2 full bands are described by the start row and the height
		        // of the band.  Each band is called a data block.

		        data_block_starts.push_back(184);  // Top band / data block
		        data_block_starts.push_back(686);  // Second from top
		        data_block_starts.push_back(1196); // Second from bottom
		        data_block_starts.push_back(1694); // Bottom band
                } 
                // Chip is 2.7M and single block analysis requested on command line.
                else
                {
		    data_block_starts.push_back(0);  // Top band / data block
                } 
	}
        // Chip is either snp6 or snp7. Single block analysis always performed.
        else
        {
                data_block_starts.push_back(0);  // Top band / data block
        }


	int image_cols = data_block_cols;      // The width of the image.  Probably this 
	                                       // will always be the same as data_block_cols.

	// To speed up the filtering the background (low frequency component) is
	// computed by small blocks rather than individual probe cells.  For 
	// example an 8x8 block will reduce the computations in the iterative
	// portion of the algorithm by 64.  The blocks need not be square or
	// a power of two or divide into a larger number evenly. Just not too
	// big for smoothness and not too small for speed.
//	int mini_block_rows=8;
//	int mini_block_cols=8;

	// The high pass filter depends on a local and global weight.  The
	// global weight enforces long term trends in the bakcground through
	// an imposed hierarchy in the data. The local weight enforces short
	// term trends by enforcing that the background be similar to the 
	// immediately surrounding region.   Implicit here is a wieght of 1.0
	// for the observed value.  The local and global weights of 8.0 and 8.0
	// are reasonable.  They can be doubled or halved without excessive
	// changes in the background.  These would make good default values.
	// They should never be less than 0.0.  Large values, say greater
	// than 32.0 are may not be harmful but not particularly useful and
	// could cause convergence to slow.  As far as I am concerned, these
	// values could be hard coded in by assigning them to a macro.
//	double local_smooth_weight = 8.0;
//	double global_smooth_weight = 8.0;

	// A reasonable value.  Can make it 0.001 if wanted.  Convergence is
	// fast.  I doubt that the extra iteration or two to get to 0.001
	// would be noticed.
//	double converged = 0.01;

//	For the Cyto2 array. The entire image is treated as a single block
//	and would look like what is below.
/*
	if (m_pEngine->getOptInt("array-size") == 658) // CytogeneticsFocused_Arrayy
	{
		data_block_starts.clear();
		data_block_starts.push_back(0);  // first row of the image
		data_block_rows = 658;       // height of the image
		data_block_cols = 658;       // width of the image

		image_cols = 658;            // width of the image
	}
*/
	// Do the bands / data blocks one at a time
	for (int iblock=0; iblock<data_block_starts.size(); iblock++) {
		// start and stop of the band
		int row_start = data_block_starts[iblock];
		int row_stop = row_start + data_block_rows;

		// For speed we need a fast lookup table to map probe ids to
		// their position in the input probe id vector vProbeIDs.
		// .first=id .second=position
		vector<pair<int,int> > lookup;
		pid_lookup(lookup,row_start,row_stop,image_cols);

		// Make an instance of the smoother.  Resizing or clearing is not
		// implemented so make sure it gets deleted after each use.
		DataBlock * DB = new DataBlock(data_block_rows, data_block_cols,
				mini_block_rows, mini_block_cols);

		// Load up the DataBlock with whatever was observed.  Missing
		// values are okay.  If not entered the entry is missing.  The
		// indexing using the lookup table is not user friendly but
		// it is fast and otherwise would be slow.
		for (int i=0; i<data_block_rows; i++) {
			for (int j=0; j<data_block_cols; j++) {
				int idx = i*image_cols + j;
				if (lookup[idx].first) {
                                        float fIntensity = getProbes()->at(lookup[idx].second)->getIntensity();
                                        float fMedianIntensity = getProbes()->at(lookup[idx].second)->getMedianIntensity();
                                        DB->set_value(i,j, log2(fIntensity/fMedianIntensity));
					}
				}
			}

		// Once the data has been entered load the full hierarchy of
		// the data block.
		DB->load_array_level();

		//  Run the smoother to convergence.  The DataBlock knows what to do.
		DB->mrf_smooth(local_smooth_weight, global_smooth_weight, converged);
		if (iblock == 0) {CNAnalysisMethod::memory("HighPassFilter::run(...)");}

		// Pull the background estimate out of the DataBlock.
		for (int i=0; i<data_block_rows; i++) {
			for (int j=0; j<data_block_cols; j++) {
				int idx = i*image_cols + j;
				if (lookup[idx].first) {
					CNProbe* p = getProbes()->at(lookup[idx].second);
                                        float fIntensity = log2(p->getIntensity());
                                        fIntensity -= DB->background(i,j);
                                        float fTwo= 2.0;
                                        p->setIntensity(pow(fTwo,fIntensity));
				}
			}

		}
                delete DB;
		
	}
        if (m_pEngine->getOptBool("keep-intermediate-data"))
        {
                writeIntensities("highpass-filter");
        }
	Verbose::out(1, "CNIntensityAdjustmentMethodHighPassFilter::run(...) end");
}


void CNIntensityAdjustmentMethodHighPassFilter::pid_lookup(vector<pair<int,int> > & lookup, const int row_start,
		const int row_stop, const int image_width)
{
	int block_sz = (row_stop - row_start)*image_width;
	lookup.resize(block_sz,make_pair<int,int>(0,0));

	for (int i=0; i<getProbes()->size(); i++) {
		pair<int,int> row_col = index_to_rc(getProbes()->at(i)->getProbeID(),image_width);
		if (row_col.first >= row_start && row_col.first < row_stop) {
			int row = row_col.first - row_start;
			int array_index = row*image_width + row_col.second;
			lookup[array_index].first = getProbes()->at(i)->getProbeID();
			lookup[array_index].second = i;
		}
	}
}
