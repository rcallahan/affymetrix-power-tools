////////////////////////////////////////////////////////////////
//
// Copyright (C) 2011 Affymetrix, Inc.
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
#include "file/TsvFile/ClfFile.h"
#include "file/TsvFile/PgfFile.h"
#include "file/TsvFile/TsvFile.h"
//
#include <cstring>
#include <string.h>
#include <string>
//

using namespace std;
using namespace affx;

void makeMapEpsilonTsv(std::string file_name,
                       const char* mapParameters[],
                       double mapEpsilons[])
{
  affx::TsvFile datasetTsv;
  char string_buf[100];

  int col_cnt = 2;
  int row_cnt = 10;

  //Print Comment
  snprintf(string_buf,sizeof(string_buf),
           "'%s' %d rows, %d cols",
           file_name.c_str(),row_cnt,col_cnt);
  datasetTsv.addHeader("Note",string_buf);

  //Create Columns and their appropriate titles
  snprintf(string_buf,sizeof(string_buf),"MapParameter");
  datasetTsv.defineColumn(0,0,string_buf);
  datasetTsv.setPrecision(0,0,1);

  snprintf(string_buf,sizeof(string_buf),"MapEpsilon");
  datasetTsv.defineColumn(0,1,string_buf);
  datasetTsv.setPrecision(0,1,2);

  //Print out the template
  datasetTsv.writeTsv(file_name);

  //Fill the columns with Data
  int i = 0;
  while (mapParameters[i]!=NULL || mapEpsilons[i]!=0.0)
    {
      datasetTsv.set(0, 0, mapParameters[i]);
      datasetTsv.set(0, 1, mapEpsilons[i]);
      datasetTsv.writeLevel(0);
      i++;
    }

  datasetTsv.close();
}

void makeDatasetTsv(std::string file_name, const char* groupNames[], const char* tsvNames[])
{
  affx::TsvFile datasetTsv;
  char string_buf[100];

  int col_cnt = 2;
  int row_cnt = 10;

  //Print Comment
  snprintf(string_buf,sizeof(string_buf),
           "'%s' %d rows, %d cols",
           file_name.c_str(),row_cnt,col_cnt);
  datasetTsv.addHeader("Note",string_buf);

  //Create Columns and their appropriate titles
  snprintf(string_buf,sizeof(string_buf),"GroupName");
  datasetTsv.defineColumn(0,0,string_buf);
  datasetTsv.setPrecision(0,0,1);

  snprintf(string_buf,sizeof(string_buf),"TsvName");
  datasetTsv.defineColumn(0,1,string_buf);
  datasetTsv.setPrecision(0,1,2);

  //Print out the template
  datasetTsv.writeTsv(file_name);

  //Fill the columns with Data
  int i = 0;
  while(groupNames[i]!= NULL || tsvNames[i]!=NULL)
    {
      datasetTsv.set(0, 0, groupNames[i]);
      datasetTsv.set(0, 1, tsvNames[i]);
      datasetTsv.writeLevel(0);
      i++;
    }

  datasetTsv.close();
}

void makeIgnoreColumnTsv(std::string file_name, const char* ignoreColumns[])
{
  affx::TsvFile ignoreTsv;
  char string_buf[100];

  int col_cnt = 1;
  int row_cnt = 17;

  //Print Comment
  snprintf(string_buf,sizeof(string_buf),
           "'%s' %d rows, %d cols",
           file_name.c_str(),row_cnt,col_cnt);
  ignoreTsv.addHeader("Note",string_buf);

  //Create Columns and their appropriate titles
  snprintf(string_buf,sizeof(string_buf),"IgnoreColumns");
  ignoreTsv.defineColumn(0,0,string_buf);
  ignoreTsv.setPrecision(0,0,1);


  //Print out the template
  ignoreTsv.writeTsv(file_name);

  //Fill the columns with Data
  int i = 0;
  while(ignoreColumns[i] != NULL)
    {
      ignoreTsv.set(0, 0, ignoreColumns[i]);
      ignoreTsv.writeLevel(0);
      i++;
    }

  ignoreTsv.close();
}


void makeIgnoreParamTsv(std::string file_name, const char* ignoreParameters[])
{
  affx::TsvFile ignoreTsv;
  char string_buf[100];

  int col_cnt = 1;
  int row_cnt = 17;

  //Print Comment
  snprintf(string_buf,sizeof(string_buf),
           "'%s' %d rows, %d cols",
           file_name.c_str(),row_cnt,col_cnt);
  ignoreTsv.addHeader("Note",string_buf);

  //Create Columns and their appropriate titles
  snprintf(string_buf,sizeof(string_buf),"IgnoreParameters");
  ignoreTsv.defineColumn(0,0,string_buf);
  ignoreTsv.setPrecision(0,0,1);


  //Print out the template
  ignoreTsv.writeTsv(file_name);

  //Fill the columns with Data
  int i = 0;
  while(ignoreParameters[i] != NULL)
    {
      ignoreTsv.set(0, 0, ignoreParameters[i]);
      ignoreTsv.writeLevel(0);
      i++;
    }

  ignoreTsv.close();
}


int main()
{
  const char* groupNames[] = {"SNPReference", "SNPReference","SNPReference",
                  "Cyto2", "Cyto2", "Cyto2",
                  "Cyto2", "Cyto2", "Cyto2",
                  "Cyto2",
                  NULL
  };

  const char* tsvNames[] = {"EMSNParameters", "EMSNArray", "SNPReference",
                "AntigenomicProbes", "MedianSignals", "ProbeEffects",
                "SNPPosteriors", "SketchCN", "SketchSNP",
                "WaveCorretion",
                NULL
  };

  const char* ignoreParameters[] = {
    "FileCreationTime", "FileIdentifier","program-version",
    "create_date", "create-date","affymetrix-algorithm-param-option-verbose",
    "affymetrix-algorithm-param-option-exec-guid", "affymetrix-algorithm-param-option-program-cvs-id","affymetrix-algorithm-param-option-version-to-report",
    "affymetrix-algorithm-param-option-command-line", "affymetrix-algorithm-param-option-mem-usage", "affymetrix-algorithm-param-option-run-probeset-genotype",
    "affymetrix-algorithm-param-option-cels", "affymetrix-algorithm-param-option-out-dir", "affymetrix-algorithm-param-state-time-start",
    "affymetrix-algorithm-param-state-free-mem-at-start", "affymetrix-algorithm-param-option-temp-dir",
    NULL
  };

  const char* mapParameters[] = {"Segments.CN.MeanMarkerDistance",
                 NULL
  };
  double mapEpsilons[] = {376.0f,
                   NULL
  };

  makeDatasetTsv("datasetFile.tsv", groupNames, tsvNames);
  makeIgnoreParamTsv("ignoreParamFile.tsv", ignoreParameters);
  makeIgnoreColumnTsv("ignoreColumnFile.tsv", ignoreParameters);
  makeMapEpsilonTsv("qt_mapEpsilonFile.tsv", mapParameters, mapEpsilons);
}

