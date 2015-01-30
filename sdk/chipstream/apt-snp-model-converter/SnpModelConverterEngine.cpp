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

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "chipstream/apt-snp-model-converter/SnpModelConverterEngine.h"
//
#include "chipstream/SnpModelDb.h"
#include "chipstream/SnpModelConverter.h"
//
#include "util/PgOptions.h"
//
#include <map>
#include <string>
#include <sstream>


SnpModelConverterEngine::Reg SnpModelConverterEngine::reg;

SnpModelConverterEngine * SnpModelConverterEngine::FromBase(BaseEngine *engine)
{
  if(engine != NULL && engine->getEngineName() == SnpModelConverterEngine::EngineName())
    return (SnpModelConverterEngine *)engine;
  return NULL;
}

SnpModelConverterEngine::SnpModelConverterEngine()
{
  defineOptions(this);
}

SnpModelConverterEngine::~SnpModelConverterEngine()
{
}

void SnpModelConverterEngine::defineOptions(SnpModelConverterEngine *engine)
{

  // Program Control
  engine->defineOptionSection("SnpModelConverter Program Control");
  defineOption("", "convert-to-sqlite", PgOpt::BOOL_OPT,
               "Convert a TSV or HDF5 file to a SQLite database.",
               "false");
  defineOption("d", "dump-headers", PgOpt::STRING_OPT,
               "Dump both the physical and derived headers with comments. "
               "The possible derived headers require one of the types [none,GTC4.1].",
               "");
  defineOption("", "header-delimiter", PgOpt::STRING_OPT,
               "Delimiter character used when more than one header exists in a file,  "
               "or a header is ambiguous. If not supplied then only the first value "
               "is returned when multiple values are detected.",
               "");
  defineOption("", "filter-headers-file", PgOpt::STRING_OPT,
               "A file containing sample TSV headers to filter on when --dump-headers is used.  "
               "NOTE: no validation is done on the headers in as possible model headers.",
               "");
  defineOption("", "validation-exit", PgOpt::BOOL_OPT,
               "The program exits with code 0 if validation is the only error. "
               "Use this option to have the program exit -1 when one or more of the "
               "the input model files does not validate.",
               "false");


  // Input Section
  engine->defineOptionSection("Input Options");
  defineOption("", "input-file", PgOpt::STRING_OPT,
               "Required input models file of type TSV, A5, SQLite or txt.",
               "");
  defineOption("", "model-files", PgOpt::STRING_OPT,
               "Input file listing a column of models files of type TSV, A5, SQLite or txt."
               "File must have a column header of \"model-files\".",
               "");

  // Output Section
  engine->defineOptionSection("Output Options");
  defineOption("", "out-file", PgOpt::STRING_OPT,
               "The name of an optional output file.",
               "");

}


/**
 * Make sure that our options are sane. Call Err::errAbort if not.
 */
void SnpModelConverterEngine::checkOptionsImp()
{

  //DUMP HEADERS
  std::string derivedHeaderType = getOpt("dump-headers");
  if(! derivedHeaderType.empty()) {
    bool optFound = false;
    for(int i = 0; i < SnpModelConverter::LAST_DERIVED_TYPE; i++) {
      if(derivedHeaderType == SNP_MODEL_DERIVED_HEADER_TYPE_STRINGS[i]) {
        optFound = true;
        break;
      }
    }

    if(! optFound) {
      std::stringstream msg;
      msg << " " << derivedHeaderType << " : ";
      msg << "--dump-headers is an invalid derived type. " << std::endl;
      msg << "The derived type must be one of [";
      for(int i = 0; i < SnpModelConverter::LAST_DERIVED_TYPE; i++) {
        msg << SNP_MODEL_DERIVED_HEADER_TYPE_STRINGS[i] << " ";
      }
      msg << "].";

      APT_ERR_ABORT(msg.str());
    }
  }


  // DELIMITER
  if(!getOpt("header-delimiter").empty()) {
    if(getOpt("header-delimiter").length() > 1) {
      APT_ERR_ABORT(" --header delimiter must be a single character.");
    }
  }

  // INPUT FILE[S]
  std::string inputFile = getOpt("input-file");
  std::string modelFiles = getOpt("model-files");
  
  if ( inputFile.empty() && !modelFiles.empty() ) {
    // MODEL FILES LISTING
    APT_ERR_ASSERT( Fs::fileExists(modelFiles),
                    modelFiles + std::string(" file not found."));
    APT_ERR_ASSERT( affx::TsvFile::extractColToVec(modelFiles, "model-files",
                                                 &m_modelFiles ) == affx::TSV_OK,
                    modelFiles + std::string(" invalid file column header 'model-files' is required."));
  }  
  else {
    if (inputFile.empty()) {
      APT_ERR_ABORT(" --input-file or --model-files is required.");
    }
    else if (!Fs::fileExists(inputFile)) {
      APT_ERR_ABORT(inputFile + " file not found.");
    }
    m_modelFiles.push_back( inputFile );
  }

  // FILTER HEADERS FILE
  std::string filterHeadersFile = getOpt("filter-headers-file");
  if ( !filterHeadersFile.empty() ) {
    if (!Fs::fileExists(filterHeadersFile)) {
      APT_ERR_ABORT(filterHeadersFile + " file not found.");
    }    
  }
  
}

void SnpModelConverterEngine::runImpDumpHeaders(const std::map< std::string , std::string *> & headers)
{

  std::map< std::string , std::string *>::const_iterator it;
  std::map< std::string , std::string > filterHeaders;
  std::string filterHeadersFile = getOpt("filter-headers-file");

  try {
    if ( !filterHeadersFile.empty() ) {
      affx::TsvFile tsv;
      tsv.open(filterHeadersFile);
      tsv.headersBegin();
      std::string key, value;
      while ( tsv.headersNext( key, value ) == affx::TSV_OK ) {
        filterHeaders[key] = value;
      }
      tsv.close();
    }

    for(it = headers.begin(); it != headers.end(); it++) {
      if ( filterHeaders.size() &&
           (filterHeaders.count(it->first) == 0) ) {
        continue;
      }
      cout << it->first << "=";
      if(it->second == NULL) {
        cout << "NULL";
      } else {
        cout << *(it->second) ;
      }
      cout << endl;
    }
  }
  catch (...) {
  }
}

void SnpModelConverterEngine::runImpConvertToDb( const std::string & modelFilePath ) {
        
  if ( SnpModelDb::isSqlModelFile( modelFilePath ) ) {
    Verbose::warn(1, " --convert-to-sqlite option ignored, file is already SQLite");
    
  }
  else {
    std::string sqlFile = modelFilePath;
    if ( !getOpt("out-file").empty() ) {
      sqlFile = getOpt("out-file");
    }
    else {
      std::string::size_type pos;
      std::string ext;
      pos  = modelFilePath.rfind(".");
      if ( pos != std::string::npos ) {
        ext = modelFilePath.substr(pos);
      }
      if ( (ext == ".a5")  || ( ext == ".txt" ) || (ext == ".tsv" ) ) {
        sqlFile = Fs::noextname1(modelFilePath) + ".db";
      }
      else {
        sqlFile = modelFilePath + ".db";
      }
    }
    if ( Fs::touch(sqlFile, false) != APT_OK ) {
      Verbose::warn(1, " --convert-to-sqlite permission denied for file " + sqlFile );
    }
    else {
      SnpModelConverter dummy;
      dummy.convertToDbModel(modelFilePath, sqlFile);
      Verbose::out(1, sqlFile + " sqlite file created.");
    }
  }
}

/**
   This is the "main()" equivalent of the engine.
*/
void SnpModelConverterEngine::runImp()
{

  
  bool exitOk = true;
  std::string src;
  std::string test = getOpt("header-delimiter");
  std::string headerType = getOpt("dump-headers");
  
  for ( int i = 0; i < m_modelFiles.size(); i++ ) {
    src = Fs::Unc(m_modelFiles[i]);

    if ( ! Fs::fileExists( src ) ) {
      Verbose::warn(1, src + " file not found." );
      exitOk = false;
      continue;
    }
    std::map< std::string, std::string *> headers;

    bool isModel =
      SnpModelConverter::GetModelHeaders(src, headerType , headers,
                                         test.empty() ? NULL : &(test[0]));
    // DUMP HEADERS
    if(! headerType.empty()) {
      runImpDumpHeaders(headers);
    }

    if(isModel) {
      Verbose::out(1, src + ToStr(" valid model file detected."));
      if ( getOptBool("convert-to-sqlite") ) {
        runImpConvertToDb(src);
      }
    }
    else if ( !getOpt("model-files").empty() ) {
      Verbose::warn(1,src + " is not a recognized model file.");
    }

    exitOk = exitOk && isModel;
  }

  if ( ! exitOk ) {
    if ( ! getOpt("input-file").empty() ) {
      Verbose::warn(1, src + " is not a recognized model file.");
    }
    else {
      Verbose::warn(1, " Invalid model files detected.");
    }
    if ( getOptBool("validation-exit") ) {
      APT_ERR_ABORT(" validation error exit requested.");
    }
  }
  return;
}
