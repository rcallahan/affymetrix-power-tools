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

//
#include "chipstream/QuantLabelZIO.h"
//
#include "util/Convert.h"
#include "util/Guid.h"
#include "util/Fs.h"

using namespace affx;

//////////

///@ brief Defines seven columns for each cluster.
void QuantLabelZ__writeSnpPosteriorTsv_defC_version_2(TsvReport& tsv,const std::string& prefix,int* cidx)
{
  // order is as defined in "snp.label.h"
  // 1d
  tsv.defineColumn(0,(*cidx)++, prefix+"_MeanStrength", affx::FILE5_DTYPE_DOUBLE); // k
  tsv.defineColumn(0,(*cidx)++, prefix+"_Mean", affx::FILE5_DTYPE_DOUBLE); // m
  tsv.defineColumn(0,(*cidx)++, prefix+"_Variance", affx::FILE5_DTYPE_DOUBLE); // ss
  tsv.defineColumn(0,(*cidx)++, prefix+"_VarianceStrength", affx::FILE5_DTYPE_DOUBLE); // v
}

///@ brief Defines seven columns for each cluster.
void QuantLabelZ__writeSnpPosteriorTsv_defC_version_3(TsvReport& tsv,const std::string& prefix,int* cidx)
{
  // order is as defined in "snp.label.h"
  // 1d
  tsv.defineColumn(0,(*cidx)++, prefix+"_Mean", affx::FILE5_DTYPE_DOUBLE); // m
  tsv.defineColumn(0,(*cidx)++, prefix+"_Variance", affx::FILE5_DTYPE_DOUBLE); // ss
  tsv.defineColumn(0,(*cidx)++, prefix+"_MeanStrength", affx::FILE5_DTYPE_DOUBLE); // k
  tsv.defineColumn(0,(*cidx)++, prefix+"_VarianceStrength", affx::FILE5_DTYPE_DOUBLE); // v
  // 2d
  tsv.defineColumn(0,(*cidx)++, prefix+"_YM", affx::FILE5_DTYPE_DOUBLE); // YM
  tsv.defineColumn(0,(*cidx)++, prefix+"_YSS", affx::FILE5_DTYPE_DOUBLE); // YSS
  tsv.defineColumn(0,(*cidx)++, prefix+"_YXSS", affx::FILE5_DTYPE_DOUBLE); // XYSS
}

/// @brief Defines the 12 columns for
void QuantLabelZ__writeSnpPosteriorTsv_defD_version_2(TsvReport& tsv,const std::string& prefix,int* cidx)
{
  // also as defined in snp.label.h
  tsv.defineColumn(0,(*cidx)++, prefix+"_XAH", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_XAB", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_XHB", affx::FILE5_DTYPE_DOUBLE);
}

/// @brief Defines the 12 columns for
void QuantLabelZ__writeSnpPosteriorTsv_defD_version_3(TsvReport& tsv,const std::string& prefix,int* cidx)
{
  // also as defined in snp.label.h
  tsv.defineColumn(0,(*cidx)++, prefix+"_XAH", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_XAB", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_XHB", affx::FILE5_DTYPE_DOUBLE);
  //
  tsv.defineColumn(0,(*cidx)++, prefix+"_YAH", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_YAB", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_YHB", affx::FILE5_DTYPE_DOUBLE);
  //
  tsv.defineColumn(0,(*cidx)++, prefix+"_XYAH", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_XYAB", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_XYHB", affx::FILE5_DTYPE_DOUBLE);
  //
  tsv.defineColumn(0,(*cidx)++, prefix+"_YXAH", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_YXAB", affx::FILE5_DTYPE_DOUBLE);
  tsv.defineColumn(0,(*cidx)++, prefix+"_YXHB", affx::FILE5_DTYPE_DOUBLE);
}

/// @brief Define the columns in a tsv report.
void QuantLabelZ__writeSnpPosteriorTsv(TsvReport& tsv,
                                       const std::string& fileName,
                                       affx::TsvReport::TsvReportFmt_t format_file,
                                       int format_ver)
{
  tsv.setFilename(fileName);
  tsv.setFormat(format_file);

  // record the version...
  tsv.addHeader("SnpPosteriorFormatVer",ToStr(format_ver));
  //
  tsv.addHeader("guid",affxutil::Guid::GenerateNewGuid());

  if (format_ver==1) {
    // the order comma seperated values are written in.
    tsv.addHeader("data-order",QUANTLABELZ__SNPTSV_DATAORDER);
    //
    tsv.defineStringColumn(0,0,"id",TSVREPORT_PROBESET_STRLEN);
    tsv.defineStringColumn(0,1,"BB",200);
    tsv.defineStringColumn(0,2,"AB",200);
    tsv.defineStringColumn(0,3,"AA",200);
    tsv.defineStringColumn(0,4,"CV",200);
  }
  else if (format_ver==2) {
    int cidx=0;
    tsv.defineColumn(0, cidx++, "id", affx::FILE5_DTYPE_STRING, TSVREPORT_PROBESET_STRLEN);
    // the posterior columns for each cluster
    QuantLabelZ__writeSnpPosteriorTsv_defC_version_2(tsv,"Cluster_AA",&cidx);
    QuantLabelZ__writeSnpPosteriorTsv_defC_version_2(tsv,"Cluster_AB",&cidx);
    QuantLabelZ__writeSnpPosteriorTsv_defC_version_2(tsv,"Cluster_BB",&cidx);
    // The cross term clusters
    QuantLabelZ__writeSnpPosteriorTsv_defD_version_2(tsv,"Cross-Term",&cidx);
  }
  else if (format_ver==3) {
    int cidx=0;
    tsv.defineColumn(0, cidx++, "id", affx::FILE5_DTYPE_STRING, TSVREPORT_PROBESET_STRLEN);
    // the posterior columns for each cluster
    QuantLabelZ__writeSnpPosteriorTsv_defC_version_3(tsv,"Cluster_AA",&cidx);
    QuantLabelZ__writeSnpPosteriorTsv_defC_version_3(tsv,"Cluster_AB",&cidx);
    QuantLabelZ__writeSnpPosteriorTsv_defC_version_3(tsv,"Cluster_BB",&cidx);
    // The cross term clusters
    QuantLabelZ__writeSnpPosteriorTsv_defD_version_3(tsv,"Cross-Term",&cidx);
  }
  else	{
    // opps!
    Err::errAbort("Bad format_ver (="+ToStr(format_ver));
  }
  //
  tsv.writeTsv_v1();
}

/// @brief write a cluster data value to the report
void QuantLabelZ__writeSnpPosteriorValue_v2_C(TsvReport& tsv,cluster_data& c,int *cidx)
{
  // the columns must be in the same order as defined in writeSnpPosteriorTsv_defP
  // 1d
  tsv.set_d(0,(*cidx)++, c.k);
  tsv.set_d(0,(*cidx)++, c.m);
  tsv.set_d(0,(*cidx)++, c.ss);
  tsv.set_d(0,(*cidx)++, c.v);
}

/// @brief write a distribution data value to the report.
void QuantLabelZ__writeSnpPosteriorValue_v2_D(TsvReport& tsv,snp_distribution& d,int *cidx)
{
  // the columns must be in the same order as defined in writeSnpPosteriorTsv_defC
  // 1d
  tsv.set_d(0,(*cidx)++, d.xah);
  tsv.set_d(0,(*cidx)++, d.xab);
  tsv.set_d(0,(*cidx)++, d.xhb);
}

/// @brief write a cluster data value to the report
void QuantLabelZ__writeSnpPosteriorValue_v3_C(TsvReport& tsv,cluster_data& c,int *cidx)
{
  // the columns must be in the same order as defined in writeSnpPosteriorTsv_defP
  // 1d
  tsv.set_d(0,(*cidx)++, c.m);
  tsv.set_d(0,(*cidx)++, c.ss);
  tsv.set_d(0,(*cidx)++, c.k);
  tsv.set_d(0,(*cidx)++, c.v);
  // 2d
  tsv.set_d(0,(*cidx)++, c.ym);
  tsv.set_d(0,(*cidx)++, c.yss);
  tsv.set_d(0,(*cidx)++, c.xyss);
}

/// @brief write a distribution data value to the report.
void QuantLabelZ__writeSnpPosteriorValue_v3_D(TsvReport& tsv,snp_distribution& d,int *cidx)
{
  // the columns must be in the same order as defined in writeSnpPosteriorTsv_defC
  // 1d
  tsv.set_d(0,(*cidx)++, d.xah);
  tsv.set_d(0,(*cidx)++, d.xab);
  tsv.set_d(0,(*cidx)++, d.xhb);
  // 2d
  tsv.set_d(0,(*cidx)++, d.yah);
  tsv.set_d(0,(*cidx)++, d.yab);
  tsv.set_d(0,(*cidx)++, d.yhb);
  //
  tsv.set_d(0,(*cidx)++, d.xyah);
  tsv.set_d(0,(*cidx)++, d.xyab);
  tsv.set_d(0,(*cidx)++, d.xyhb);
  //
  tsv.set_d(0,(*cidx)++, d.yxah);
  tsv.set_d(0,(*cidx)++, d.yxab);
  tsv.set_d(0,(*cidx)++, d.yxhb);
}

/// @brief Write a snp_param to the report.
void QuantLabelZ__writeSnpPosteriorValue(TsvReport& tsv,
                                         int format_ver,
                                         const std::string& probeset_id,
                                         snp_param& tsp)
{
  if (format_ver==1) {
    //
    tsv.set_string(0,0,probeset_id);
    // use names as the column order is defined as ("BB","AB","AA") which is odd
    tsv.set_string(0,1,c_to_str(tsp.posterior.bb,",",tsp.clustertype));
    tsv.set_string(0,2,c_to_str(tsp.posterior.ab,",",tsp.clustertype));
    tsv.set_string(0,3,c_to_str(tsp.posterior.aa,",",tsp.clustertype));
    // CV
    std::ostringstream stream;
    stream << tsp.posterior.xah << "," << tsp.posterior.xab << "," << tsp.posterior.xhb;
    if (tsp.clustertype>1) {
      stream << "," << tsp.posterior.yah << "," << tsp.posterior.yab << "," << tsp.posterior.yhb;
      stream << "," << tsp.posterior.xyah << "," << tsp.posterior.xyab << "," << tsp.posterior.xyhb;
      stream << "," << tsp.posterior.yxah << "," << tsp.posterior.yxab << "," << tsp.posterior.yxhb;
    }
    tsv.set_string(0,4,stream.str());
  }
  //
  else if (format_ver==2) {
    int cidx=0;
    tsv.set_string(0, cidx++, probeset_id);
    //
    QuantLabelZ__writeSnpPosteriorValue_v2_C(tsv,tsp.posterior.aa,&cidx);
    QuantLabelZ__writeSnpPosteriorValue_v2_C(tsv,tsp.posterior.ab,&cidx);
    QuantLabelZ__writeSnpPosteriorValue_v2_C(tsv,tsp.posterior.bb,&cidx);
    // cross terms
    QuantLabelZ__writeSnpPosteriorValue_v2_D(tsv,tsp.posterior,&cidx);
  }
  else if (format_ver==3) {
    int cidx=0;
    tsv.set_string(0, cidx++, probeset_id);
    //
    QuantLabelZ__writeSnpPosteriorValue_v3_C(tsv,tsp.posterior.aa,&cidx);
    QuantLabelZ__writeSnpPosteriorValue_v3_C(tsv,tsp.posterior.ab,&cidx);
    QuantLabelZ__writeSnpPosteriorValue_v3_C(tsv,tsp.posterior.bb,&cidx);
    // cross terms
    QuantLabelZ__writeSnpPosteriorValue_v3_D(tsv,tsp.posterior,&cidx);
  }
  else {
    Err::errAbort("writeSnpPosteriorValue: bad ver");
  }
  //
  tsv.writeLevel(0);
}

/**
 * export prior/posterior parameters for a genotype as a string
 *
 * @param d - cluster data
 * @param sep - separation
 */
std::string c_to_str(cluster_data d, const std::string& sep, int dimension)
{
  // some of the values are coming out as "-0".
  // add a test for 0.0 and make it a postive 0.

  // export as string
  std::string st;
  st += ToStr((d.m==0.0)?0.0:d.m);
  st += sep;
  st += ToStr(d.ss);
  st += sep;
  st += ToStr(d.k);
  st += sep;
  st += ToStr(d.v);
  // updated with new posterior
  // adding 2nd dimension
  if (dimension>1)
    {
      st += sep;
      st += ToStr(d.ym);
      st += sep;
      st += ToStr(d.yss);
      st += sep;
      st += ToStr((d.xyss==0.0)?0.0:d.xyss);
    }
  return(st);
}

/**
 * take string and manufacture cluster data for genotype
 *
 * @param d - cluster data
 * @param st - string containing data
 */
void str_to_c(cluster_data &d,const std::string& str) {
  vector<string> words;

  Util::chopString(str, ',', words);

  // APT-550: give a better error message.
  if (!((words.size()==4)||(words.size()==7))) {
    APT_ERR_ABORT("QuantLabelZ__SnpPriorFromStrings: misformatted Prior string."
                  "'"+str+"' does not have 4 or 7 elements.");
  }

  d.m = Convert::toFloat(words[0]);
  d.ss = Convert::toFloat(words[1]);
  d.k = Convert::toFloat(words[2]);
  d.v = Convert::toFloat(words[3]);
  // if this is 2d clustering
  if (words.size()==7){
	  d.ym = Convert::toFloat(words[4]);
	  d.yss = Convert::toFloat(words[5]);
	  d.xyss = Convert::toFloat(words[6]);
	}
}

/**
 * Generate covariances from string
 *
 * @param ah - aa vs ab genotype covariance
 * @param ab - aa vs bb genotype covariance
 * @param hb - ab vs bb genotype covariance
 * @param st - string to be broken up into covariances
 */
void str_to_cv(snp_distribution& Dist,const std::string& str)
{
  vector<string> words;
  Util::chopString(str,',',words);

   // APT-550: give a better error message.
  if (!((words.size()==3)||(words.size()==12))) {
    APT_ERR_ABORT("QuantLabelZ__SnpPsteriorFromStrings: misformatted Posterior string."
                  "'"+str+"' does not have 3 or 12 elements.");
  }
  // match writeSnpPosteriorValue order

  // The 1d...
  Dist.xah = Convert::toFloat(words[0]);
  Dist.xab = Convert::toFloat(words[1]);
  Dist.xhb = Convert::toFloat(words[2]);

  // If there is more data then it is 2d clustering.
  if (words.size()>3) {
    // 2d...
	  Dist.yah = Convert::toFloat(words[3]);
	  Dist.yab = Convert::toFloat(words[4]);
	  Dist.yhb = Convert::toFloat(words[5]);
	  Dist.xyah = Convert::toFloat(words[6]);
	  Dist.xyab = Convert::toFloat(words[7]);
	  Dist.xyhb = Convert::toFloat(words[8]);
	  Dist.yxah = Convert::toFloat(words[9]);
	  Dist.yxab = Convert::toFloat(words[10]);
	  Dist.yxhb = Convert::toFloat(words[11]);
  }
}

/**
 * construct prior from strings
 *
 * @param tsp - snp parameters including prior to be set
 * @param bb - string describing bb genotype
 * @param ab - string describing ab genotype
 * @param aa - string describing aa genotype
 * @param cv - string describing covariances between clusters
 */
void QuantLabelZ__SnpPriorFromStrings(snp_labeled_distribution &tsp,
                                      const std::string& bb,
                                      const std::string& ab,
                                      const std::string& aa,
                                      const std::string& cv)
{
  str_to_c(tsp.Dist.aa,aa);
  str_to_c(tsp.Dist.ab,ab);
  str_to_c(tsp.Dist.bb,bb);
  str_to_cv(tsp.Dist,cv);
}

/**
 * construct prior from strings
 *
 * @param bb - string for bb genotype
 * @param ab - string for ab genotype
 * @param aa - string for aa genotype
 * @param cv - string for covariances between clusters
 * @return snp parameters including prior
 */
snp_labeled_distribution
QuantLabelZ__CreatePriorFromStrings(const std::string& bb,
                                    const std::string& ab,
                                    const std::string& aa,
                                    const std::string& cv)
{
  snp_labeled_distribution tsp;
  // tsp.Initialize();
  tsp.Dist.Clear();
  QuantLabelZ__SnpPriorFromStrings(tsp,bb,ab,aa,cv);
  return(tsp);
}

/**
 * load snp specific prior table mapping snp names to specific priors
 *
 * @param PriorMap - map of names to snp parameters
 * @param fileName - file to read these parameters
 */
void QuantLabelZ__readSnpPriorMap(AffxArray<snp_labeled_distribution>& Priors, const std::string& fileName)
{
  // is this an A5/HDF5 file?
  if (Fs::isHdf5File(fileName)) {
    APT_ERR_ABORT("The file '"+fileName+"' appears to be an A5/HDF5 file.");
  }
  
  affx::TsvFile tsv;
  int cidx_aa=-1;

  // look to see if there is a "AA" column,
  // that means v1 format.
  tsv.open(fileName);
  cidx_aa=tsv.cname2cidx(0,"AA");
  tsv.close();

  //
  if (cidx_aa>=0) {
    QuantLabelZ__readSnpPriorMap_v1(Priors,fileName);
  }
  else {
    QuantLabelZ__readSnpPriorMap_v2(Priors,fileName);
  }
}

void QuantLabelZ__readSnpPriorMap_v1(AffxArray<snp_labeled_distribution>& Priors,
                                     const std::string& fileName)
{
  affx::TsvFile tsv;
  string aa, ab, bb, id,cv;
  Priors.deleteAll();
  Priors.reserve(1000000);
  tsv.open(fileName);

  //  check the order
  std::string data_order;
  if (tsv.getHeader("data-order",data_order)==affx::TSV_OK) {
    /*
    if (data_order!=QUANTLABELZ__SNPTSV_DATAORDER) {
      Err::errAbort("QuantLabelZ__readSnpPriorMap_v1: data-order is '"+data_order+"'."
                    "Expecing '" QUANTLABELZ__SNPTSV_DATAORDER "'.");
    }
    */
  }

  // This is the "v1" format with comma seperated values...
  tsv.bind(0,"id", &id, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"BB", &bb, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"AB", &ab, affx::TSV_BIND_REQUIRED);
  tsv.bind(0,"AA", &aa, affx::TSV_BIND_REQUIRED);

  int fc_cidx=tsv.cname2cidx(0,"CV");
  if (fc_cidx!=affx::TSV_ERR_NOTFOUND) {
    tsv.bind(0,"CV", &cv, affx::TSV_BIND_REQUIRED);
  }
  else {
    cv = "0,0,0";
  }

 int iCount = 0;
  while (tsv.nextLevel(0) == affx::TSV_OK) {
    snp_labeled_distribution* p = new snp_labeled_distribution();
    *p = QuantLabelZ__CreatePriorFromStrings(bb,ab,aa,cv);
    p->probeset_id = id;
    Priors.add(p);
    iCount++;
  }
  tsv.close();
  Priors.resize(iCount);
  Priors.quickSort(0);
  //
  if (iCount == 0) {
	  Verbose::out(1, "WARNING: Cannot find any SNP priors in the input file.");
  }
  else {
    Verbose::out(3, "Found " + ToStr(iCount) + " SNP priors in the input file.");
  }
}

void QuantLabelZ__readSnpPriorMap_tsv_v2_C(affx::TsvFile& tsv,
                                           const std::string& prefix,
                                           cluster_data& c)
{
  //
  c.m=c.ss=c.k=c.v=0.0;
  tsv.get(0,prefix+"_Mean"            , c.m );
  tsv.get(0,prefix+"_Variance"        , c.ss);
  tsv.get(0,prefix+"_MeanStrength"    , c.k );
  tsv.get(0,prefix+"_VarianceStrength", c.v );
  //
  c.ym=c.yss=c.xyss=0.0;
  tsv.get(0,prefix+"_YM"  , c.ym);
  tsv.get(0,prefix+"_YSS" , c.yss);
  tsv.get(0,prefix+"_YXSS", c.xyss);
}

void QuantLabelZ__readSnpPriorMap_tsv_v2_D(affx::TsvFile& tsv,
                                           const std::string& prefix,
                                           snp_distribution& d)
{
  // same as writeSnpPosteriorValue_v2_D
  // 1d
  tsv.get(0,prefix+"_XAH",d.xah);
  tsv.get(0,prefix+"_XAB",d.xab);
  tsv.get(0,prefix+"_XHB",d.xhb);
  // 2d
  tsv.get(0,prefix+"_YAH",d.yah);
  tsv.get(0,prefix+"_YAB",d.yab);
  tsv.get(0,prefix+"_YHB",d.yhb);
  //
  tsv.get(0,prefix+"_XYAH",d.xyah);
  tsv.get(0,prefix+"_XYAB",d.xyab);
  tsv.get(0,prefix+"_XYHB",d.xyhb);
  //
  tsv.get(0,prefix+"_YXAH",d.yxah);
  tsv.get(0,prefix+"_YXAB",d.yxab);
  tsv.get(0,prefix+"_YXHB",d.yxhb);
}

void QuantLabelZ__readSnpPriorMap_v2(AffxArray<snp_labeled_distribution>& Priors,
                                     const std::string& fileName)
{
  Priors.deleteAll();
  Priors.reserve(1000000);
  int iCount = 0;
  int cidx=0;

  affx::TsvFile tsv;
  if (tsv.open(fileName)!=affx::TSV_OK) {
    Err::errAbort("QuantLabelZ__readSnpPriorMap_v2: open failed.");
  }

  while (tsv.nextLevel(0) == affx::TSV_OK) {
    snp_labeled_distribution* p = new snp_labeled_distribution;
    cidx=0;
    tsv.get(0,cidx++,p->probeset_id);
    QuantLabelZ__readSnpPriorMap_tsv_v2_C(tsv,"Cluster_AA",p->Dist.aa);
    QuantLabelZ__readSnpPriorMap_tsv_v2_C(tsv,"Cluster_AB",p->Dist.ab);
    QuantLabelZ__readSnpPriorMap_tsv_v2_C(tsv,"Cluster_BB",p->Dist.bb);
    //
    QuantLabelZ__readSnpPriorMap_tsv_v2_D(tsv,"Cross-Term",p->Dist);
    //
    Priors.add(p);
    iCount++;
  }
  //
  Priors.resize(iCount);
  Priors.quickSort(0);
  //
  if (iCount == 0) {
	  Verbose::out(1, "WARNING: Cannot find any SNP priors in the input file.");
  }
  else {
    Verbose::out(3, "Found " + ToStr(iCount) + " SNP priors in the input file.");
  }
}

//////////


void QuantLabelZ__readSnpPriorMap_tsv5_v1(AffxArray<snp_labeled_distribution>& Priors,
                                          affx::File5_Tsv* tsv5)
{
  string  id, aa, ab, bb, cv;

  // This is the "v1" format with comma seperated values...
  int tsv5_id_cidx=tsv5->getColumnIdx(0,"id");
  if(tsv5_id_cidx < 0)
    Err::errAbort("Unable to find column id in SNP priors file.");

  int tsv5_AA_cidx=tsv5->getColumnIdx(0,"AA");
  if(tsv5_AA_cidx < 0)
    Err::errAbort("Unable to find column AA in SNP priors file.");

  int tsv5_AB_cidx=tsv5->getColumnIdx(0,"AB");
  if(tsv5_AB_cidx < 0)
    Err::errAbort("Unable to find column AB in SNP priors file.");

  int tsv5_BB_cidx=tsv5->getColumnIdx(0,"BB");
  if(tsv5_BB_cidx < 0)
    Err::errAbort("Unable to find column BB in SNP priors file.");

  int tsv5_CV_cidx=tsv5->getColumnIdx(0,"CV");
  if (tsv5_CV_cidx<0)
    cv="0,0,0";

  Priors.deleteAll();
  Priors.reserve(1000000);

  int iCount = 0;
  while (tsv5->nextLevel(0) == affx::FILE5_OK) {
    tsv5->get(0,tsv5_id_cidx,&id);
    tsv5->get(0,tsv5_AA_cidx,&aa);
    tsv5->get(0,tsv5_AB_cidx,&ab);
    tsv5->get(0,tsv5_BB_cidx,&bb);
    if (tsv5_CV_cidx>=0) {
      tsv5->get(0,tsv5_CV_cidx,&cv);
    }
    //
    snp_labeled_distribution* p = new snp_labeled_distribution();
    *p = QuantLabelZ__CreatePriorFromStrings(bb,ab,aa,cv);
    p->probeset_id = id;
    Priors.add(p);
    iCount++;
  }

  Priors.resize(iCount);
  Priors.quickSort(0);
  //
  if (iCount == 0) {
	  Verbose::out(1, "WARNING: Cannot find any SNP priors in the input file.");
  }
  else {
    Verbose::out(3, "Found " + ToStr(iCount) + " SNP priors in the input file.");
  }
}

void QuantLabelZ__readSnpPriorMap_tsv5_v2_C(affx::File5_Tsv* tsv5,
                                            const std::string& prefix,
                                            cluster_data& c)
{

  c.Clear();
  //
  tsv5->get(0,prefix+"_Mean"            , &c.m );
  tsv5->get(0,prefix+"_Variance"        , &c.ss);
  tsv5->get(0,prefix+"_MeanStrength"    , &c.k );
  tsv5->get(0,prefix+"_VarianceStrength", &c.v );
  //
  tsv5->get(0,prefix+"_YM"  , &c.ym);
  tsv5->get(0,prefix+"_YSS" , &c.yss);
  tsv5->get(0,prefix+"_YXSS", &c.xyss);
}

void QuantLabelZ__readSnpPriorMap_tsv5_v2_D(affx::File5_Tsv* tsv5,
                                            const std::string& prefix,
                                            snp_distribution& d)
{
  d.Clear_Cross();
  // 1d
  tsv5->get(0,prefix+"_XAH",&d.xah);
  tsv5->get(0,prefix+"_XAB",&d.xab);
  tsv5->get(0,prefix+"_XHB",&d.xhb);
  // 2d
  tsv5->get(0,prefix+"_YAH",&d.yah);
  tsv5->get(0,prefix+"_YAB",&d.yab);
  tsv5->get(0,prefix+"_YHB",&d.yhb);
  //
  tsv5->get(0,prefix+"_XYAH",&d.xyah);
  tsv5->get(0,prefix+"_XYAB",&d.xyab);
  tsv5->get(0,prefix+"_XYHB",&d.xyhb);
  //
  tsv5->get(0,prefix+"_YXAH",&d.yxah);
  tsv5->get(0,prefix+"_YXAB",&d.yxab);
  tsv5->get(0,prefix+"_YXHB",&d.yxhb);
}

void QuantLabelZ__readSnpPriorMap_tsv5_v2(int iXChromosome, int iYChromosome, CNExperiment& objExperiment, CNProbeSetArray& vProbeSets,
                                          affx::File5_Tsv* tsv5)
{
  // function assumes probe sets are already sorted.
  int cidx=0;
  AffxString str;
  CNProbeSet objSearch;
  int iCopyNumber = 2;
  int iFindString = -1;
  while (tsv5->nextLevel(0) == affx::FILE5_OK) {
    cidx=0;
    tsv5->get(0,cidx++, &str);
	iCopyNumber = 2;
	iFindString = str.find(":");
	if (iFindString != std::string::npos)
	{
		iCopyNumber = ::getInt(str.substring(iFindString + 1));
		str = str.substring(0, iFindString);
	}
	objSearch.setProbeSetName(str);
	int iSearchIndex = vProbeSets.binarySearch(objSearch, 0);
	if (iSearchIndex != -1)
	{
		CNProbeSet* pobjProbeSet = vProbeSets.getAt(iSearchIndex);
		bool bMatch = false;
		bool bSexChromosome = false;
		if ((pobjProbeSet->getChromosome() == iXChromosome) && (!pobjProbeSet->isPseudoAutosomalRegion()))
                {
			bSexChromosome = true;
                        if (objExperiment.getRawIntensityRatioGenderAsInt() == affx::Male)   // Male on X chromosome
                        {
			        if ((iCopyNumber == 1) && pobjProbeSet->getImputedCNState() == 1) 
                                {
                                    bMatch = true;
                                }
			        if ((iCopyNumber == 2) && pobjProbeSet->getImputedCNState() >= 2) 
                                {
                                    bMatch = true;
                                }
			        if ((iCopyNumber == 2) && pobjProbeSet->getImputedCNState() < 1) 
                                {
                                    bMatch = true;
                                }
                        }
                        else  // Female on X chromosome
                        {
		                if (iCopyNumber == 2) 
                                {
                                        bMatch = true;
                                }
                        }
                }
                else  // Y chromosome Autosome or Par
                { 
                        if (pobjProbeSet->getChromosome() == iYChromosome)
                        {
                                bSexChromosome = true;
			        if ((objExperiment.getRawIntensityRatioGenderAsInt() == affx::Male) && (iCopyNumber == 1)) // Male on Y 
                                {
                                        bMatch = true;
                                }
			        if ( !(objExperiment.getRawIntensityRatioGenderAsInt() == affx::Male)  && (iCopyNumber == 0))  // Female on Y 
                                {
                                        bMatch = true;
                                }
		        } 
		        else // Autosome or Par 
                        {
                                if (iCopyNumber == 2) 
                                {
                                        bMatch = true;
                                }
                        }
                } 
		if (bMatch)
		{
			if ((bSexChromosome) || (pobjProbeSet->getSnpDistribution() == NULL))
			{
				snp_distribution* p = NULL;
				if (pobjProbeSet->getSnpDistribution() == NULL) {p = new snp_distribution;}
				else {p = pobjProbeSet->getSnpDistribution();}

				QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5,"Cluster_AA",p->aa);
				QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5,"Cluster_AB",p->ab);
				QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5,"Cluster_BB",p->bb);
				//
				QuantLabelZ__readSnpPriorMap_tsv5_v2_D(tsv5,"Cross-Term",*p);
				//
				pobjProbeSet->setSnpDistribution(p);
			}
		}
    }
  }
}

void QuantLabelZ__readSnpPriorMap_tsv5_v2(AffxArray<snp_labeled_distribution>& Priors,
                                          affx::File5_Tsv* tsv5)
{
  Priors.deleteAll();
  Priors.reserve(1000000);
  int iCount = 0;

  while (tsv5->nextLevel(0) == affx::FILE5_OK) {
    snp_labeled_distribution* p = new snp_labeled_distribution;
    tsv5->get(0,0, &p->probeset_id);
    QuantLabelZ__readSnpPrior_tsv5_v2(tsv5, p->Dist);
    // QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5,"Cluster_AA",p->Dist.aa);
    // QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5,"Cluster_AB",p->Dist.ab);
    // QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5,"Cluster_BB",p->Dist.bb);
    // //
    // QuantLabelZ__readSnpPriorMap_tsv5_v2_D(tsv5,"Cross-Term",p->Dist);
    //
    Priors.add(p);
    iCount++;
  }
  //
  Priors.resize(iCount);
  Priors.quickSort(0);
  //
  if (iCount == 0) {
	  Verbose::out(1, "WARNING: Cannot find any SNP priors in the input file.");
  }
  else {
    Verbose::out(3, "Found " + ToStr(iCount) + " SNP priors in the input file.");
  }
}

void QuantLabelZ__readSnpPriorMap_tsv5(AffxArray<snp_labeled_distribution>& Priors,
                                       affx::File5_Tsv* tsv5)
{
  int cidx;
  cidx=tsv5->getColumnIdx(0,"AA");
  if (cidx>=0) {
    return QuantLabelZ__readSnpPriorMap_tsv5_v1(Priors,tsv5);
  }
  //
  cidx=tsv5->getColumnIdx(0,"Cluster_AA_Mean");
  if (cidx>=0) {
    return QuantLabelZ__readSnpPriorMap_tsv5_v2(Priors,tsv5);
  }
  //
  Err::errAbort("QuantLabelZ__readSnpPriorMap_tsv5");
}

void QuantLabelZ__readSnpPrior_tsv5_v2(File5_Tsv*tsv5, snp_distribution &cluster) {
  cluster.Clear();
  QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5, "Cluster_AA", cluster.aa);
  QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5, "Cluster_AB", cluster.ab);
  QuantLabelZ__readSnpPriorMap_tsv5_v2_C(tsv5, "Cluster_BB", cluster.bb);
  QuantLabelZ__readSnpPriorMap_tsv5_v2_D(tsv5, "Cross-Term", cluster);
}

/// @brief convert a set of priors from one format to another.
void QuantLabelZ__convert_priors_tsv(const std::string& file_in,
                                     const std::string& file_out,
                                     int format_ver,
                                     int verbose)
{
  // where to park the data
  AffxArray<snp_labeled_distribution> Priors;

  // slurp
  if (verbose>=1) {
    printf("### reading priors from '%s'...\n",file_in.c_str());
  }
  QuantLabelZ__readSnpPriorMap(Priors,file_in);

  //
  if (verbose>=1) {
    printf("### writing posteriors to '%s'...\n",file_out.c_str());
  }

  affx::TsvReport tsv;
  QuantLabelZ__writeSnpPosteriorTsv(tsv,file_out,TsvReport::FMT_TSV,format_ver);

  for (int i=0;i<Priors.size();i++) {
    // bleah -- why arent these symmetric?
    snp_param tsp;
    tsp.posterior=Priors[i]->Dist;
    tsp.clustertype=2;
    //
    QuantLabelZ__writeSnpPosteriorValue(tsv,format_ver,Priors[i]->probeset_id,tsp);
  }
  tsv.close();
}

/// @brief convert a set of priors from one format to another.
void QuantLabelZ__convert_priors_tsv5(const std::string& file_in,
                                      const std::string& tsv_name,
                                      const std::string& file_out,
                                      int format_ver,
                                      int verbose)
{
  // where to park the data
  AffxArray<snp_labeled_distribution> Priors;

  // slurp
  if (verbose>=1) {
    printf("### reading priors from A5 '%s:%s'...\n",
           file_in.c_str(),tsv_name.c_str());
  }
  //
  affx::File5_File* file5;
  affx::File5_Tsv* tsv5;
  file5=new affx::File5_File();
  file5->open(file_in,affx::FILE5_OPEN_RO);
  tsv5=file5->openTsv(tsv_name);
  //
  QuantLabelZ__readSnpPriorMap_tsv5(Priors,tsv5);
  //
  tsv5->close();
  delete tsv5;
  file5->close();
  delete file5;

  //
  if (verbose>=1) {
    printf("### writing posteriors to '%s'...\n",file_out.c_str());
  }

  affx::TsvReport tsv;
  QuantLabelZ__writeSnpPosteriorTsv(tsv,file_out,TsvReport::FMT_TSV,format_ver);

  for (int i=0;i<Priors.size();i++) {
    // bleah -- why arent these symmetric?
    snp_param tsp;
    tsp.posterior=Priors[i]->Dist;
    tsp.clustertype=2;
    //
    QuantLabelZ__writeSnpPosteriorValue(tsv,format_ver,Priors[i]->probeset_id,tsp);
  }
  tsv.close();

}

/// @brief convert a set of priors from one format to another.
void QuantLabelZ__convert_priors_tsv_to_tsv5(const std::string& file_in,
                                             const std::string& tsv_name,
                                             const std::string& file_out,
                                             int format_ver,
                                             int verbose)
{
  // where to park the data
  AffxArray<snp_labeled_distribution> Priors;

  // slurp
  if (verbose>=1) {
    printf("### reading priors from A5 '%s:%s'...\n",
           file_in.c_str(),tsv_name.c_str());
  }

  QuantLabelZ__readSnpPriorMap(Priors,file_in);

  if (verbose>=1) {
    printf("### writing posteriors to '%s'...\n",file_out.c_str());
  }

  affx::TsvReport tsv;
  tsv.setGroupname(tsv_name);
  QuantLabelZ__writeSnpPosteriorTsv(tsv,file_out,TsvReport::FMT_A5,format_ver);

  for (int i=0;i<Priors.size();i++) {
    // bleah -- why arent these symmetric?
    snp_param tsp;
    tsp.posterior=Priors[i]->Dist;
    tsp.clustertype=2;
    //
    QuantLabelZ__writeSnpPosteriorValue(tsv,format_ver,Priors[i]->probeset_id,tsp);
  }
  tsv.close();

}
