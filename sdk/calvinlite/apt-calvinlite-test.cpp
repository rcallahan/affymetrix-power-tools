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
// 
// CalvinLite/apt-calvinlite-test.cpp ---
// 
// $Id: apt-calvinlite-test.cpp,v 1.2 2009-10-29 22:28:59 harley Exp $
// 

//
#include "calvinlite/CalvinLite.h"
//
#include "calvinlite/CL_File.h"
#include "calvinlite/CL_Gdh.h"
#include "calvinlite/CL_string.h"
#include "calvinlite/CL_util.h"
//
#include <assert.h>
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

//////////

void assert_str_state(CL_string& str,
                      bool n_valid,const std::string& n_str,
                      bool w_valid,const std::string& w_str)
{
  str.dump();
  //
  assert(str.m_nstr_valid==n_valid);
  if (str.m_nstr_valid) {
    assert(str.m_nstr==n_str);
  }
  assert(str.m_wstr_valid==w_valid);
}
                  
void test_cl_string()
{
  CL_string str;

  str.setNstr("abc");
  assert_str_state(str,true,"abc",false,"");
  str.ensureWstrValid();
  assert_str_state(str,true,"abc",true,"");

  str.setNstr("01234");
  assert_str_state(str,true,"01234",false,"");
  str.ensureWstrValid();
  assert_str_state(str,true,"01234",true,"");

}

//////////

void test_err_handling() 
{
  bool caught_an_error=false;
  CL_File* clfile=new CL_File();

  // by default errors should be thrown.
  try {
    printf("Expected error: ");
    clfile->setErr(CL_ERR,"foo");
    assert(0); // should not be reached.
  }
  catch (CalvinLiteException e) {
    caught_an_error=true;
  }
  catch (...) {
    // caught something but it wasnt what we wanted... throw it on.
    assert(0);
  }
  if (caught_an_error==false) {
    throw std::runtime_error("whoops!");
  }

  // But if we call "throwOnErr(false)" we shouldnt have a throw.
  caught_an_error=false;
  try {
    clfile->throwOnErr(false);
    assert(clfile->setErr(CL_ERR,"foo")==CL_ERR);
    assert(clfile->errNum()==CL_ERR);
  }
  catch (CalvinLiteException e) {
    caught_an_error=true;
  }
  catch (...) {
    // caught something but it wasnt what we wanted... throw it on.
    assert(0);
  }
  if (caught_an_error==true) {
    throw std::runtime_error("whoops!");
  }

  delete clfile;
}
//////////

void test_small_file_1(const std::string& filename)
{
  CL_File* clfile;

  // the smallest calvin file.
  clfile=new CL_File();
  clfile->write_File(filename);
  delete clfile;

  // open the file
  clfile=new CL_File();
  clfile->read_File(filename);
  clfile->dumpSegs();
  delete clfile;
}

void test_small_file_2(const std::string& filename)
{
  char buf[100];
  CL_File* clfile;

  clfile=new CL_File();

  CL_Gdh* gdh_0=clfile->getGdh();
  gdh_0->defaultValues();
  CL_Gdh* gdh_1=gdh_0->newLeafGdh();
  gdh_1->defaultValues();
  
  CL_DataGroup* dg=clfile->newDataGroup("ABC");

  const int dset_cnt=3;
  const int drow_cnt=0x10;
  CL_DataSet* ds;

  for (int dsi=0;dsi<dset_cnt;dsi++) {
    sprintf(buf,"set-%02d",dsi);
    ds=dg->newDataSet(buf);
    int col_11_cidx=ds->newColumn("col-11",CL_TC_BYTE);
    //int col_12_cidx=ds->newColumn("col-12",CL_TC_BYTE);
    //int col_13_cidx=ds->newColumn("col-13",CL_TC_BYTE);
    //int col_14_cidx=ds->newColumn("col-14",CL_TC_INT);
    //const int r_cnt=1;
    ds->setRowCount(drow_cnt);
    
    for (int r_i=0;r_i<drow_cnt;r_i++) {
      ds->set(r_i,col_11_cidx,r_i);
      //ds->set(r_i,col_11_cidx,0xFF);
      // ds->set(r_i,col_12_cidx,0x22);
      // ds->set(r_i,1,0x11);
      // ds->set(r_i,2,0x22);
      // ds->set(r_i,col_14_cidx,0x8899aabb);
    }

    ds->dump("ds");
  }

  clfile->write_File(filename+"-1");
  
  int col_21_cidx=ds->newColumn("col-21",CL_TC_BYTE);
  //int col_22_cidx=ds->newColumn("col-22",CL_TC_BYTE);
  for (int r_i=0;r_i<drow_cnt;r_i++) {
    ds->set(r_i,col_21_cidx,0xFF);
    //ds->set(r_i,col_22_cidx,0x44);
  }

  //ds->newColumn("col-21",CL_TC_FLOAT);
  //for (int r_i=0;r_i<r_cnt;r_i++) {
  //  ds->setFloat(r_i,2,100.0*r_i);
  //}

  //ds->dump("");
  clfile->write_File(filename+"-2");

  delete clfile;
}

void test_small_file_3(const std::string& filename)
{
  char buf[100];
  CL_File* clfile;
  
  clfile=new CL_File();

  CL_Gdh* gdh_0=clfile->getGdh();
  gdh_0->defaultValues();
  
  CL_DataGroup* dg=clfile->newDataGroup("ABC");

  CL_DataSet*   ds=dg->newDataSet("DEF");
  int col_00_cidx=ds->newColumn("slen-0",CL_TC_TEXT_ASCII8,0);
  int col_01_cidx=ds->newColumn("slen-1",CL_TC_TEXT_ASCII8,1);
  int col_02_cidx=ds->newColumn("slen-2",CL_TC_TEXT_ASCII8,2);
  int col_03_cidx=ds->newColumn("slen-10",CL_TC_TEXT_ASCII8,10);
  //int col_04_cidx=ds->newColumn("str-len-100",CL_TC_STRING,100);
  ds->setRowCount(10);

#define GENDATA(_ridx,_cidx) { sprintf(buf,"r%03d-c%03d",_ridx,_cidx); ds->set(_ridx,_cidx,buf); }

  for (int r_i=0;r_i<ds->rowCount();r_i++) {
    GENDATA(r_i,col_00_cidx);
    GENDATA(r_i,col_01_cidx);
    GENDATA(r_i,col_02_cidx);
    GENDATA(r_i,col_03_cidx);
    //GENDATA(r_i,col_04_cidx);
  }

  clfile->write_File(filename);
  delete clfile;
}

// This is on the wiki page.
// http://infowiki.ev.affymetrix.com/index.php/CalvinLite_Examples
void test_small_file_4(const std::string& filename)
{
  // a new yet-unamed file
  CL_File* clfile=new CL_File();

  // a new group and dataset
  CL_DataGroup* dg=clfile->newDataGroup("MyDataGroup");
  CL_DataSet* ds=dg->newDataSet("MyDataSet");

  // three columns in our dataset.
  int col_1_cidx=ds->newColumn("chr",CL_TC_BYTE);
  int col_2_cidx=ds->newColumn("base-pos",CL_TC_INT);
  int col_3_cidx=ds->newColumn("probe-id",CL_TC_TEXT_ASCII8,16);

  const int drow_cnt=0x10;
  ds->setRowCount(drow_cnt);

  // put some data in the dataset
  // can index by column idx or name
  for (int r_i=0;r_i<drow_cnt;r_i++) {
    ds->set(r_i,col_1_cidx,r_i);
    ds->set(r_i,col_2_cidx,r_i*10);
    ds->set(r_i,col_3_cidx,"abc"+ToStr(r_i));
  }

  // all done... 
  clfile->write_File(filename);
  delete clfile;
}

void test_rw_file(const std::string& filename,int opt_write)
{
  CL_File* cl_file;
  cl_file=new CL_File();

  assert(cl_file->open(filename)==CL_OK);

  CL_DataSet* ds=cl_file->getDataSet(0,0);
  if (ds!=NULL) {
    float intensity;
    for (int r=0;((r<10)&&(r<ds->rowCount()));r++) {
      ds->get(r,0,&intensity);
      printf("%3d : %f\n",r,intensity);
    }
  }

  if (opt_write==1) {
    cl_file->write(filename+".out");
  }

  delete cl_file;
}

void test_uuid()
{
  for (int i=0;i<10;i++) {
    std::string uuid=CL_gen_uuid();
    printf("%2d : %s\n",i,uuid.c_str());
  }
}

//////////

#define TEST_WRITEARRAY_1_LEN 256

void test_writearray_1()
{
  char  arrayBuf_char[TEST_WRITEARRAY_1_LEN];
  int   tmp_int;
  int   arrayBuf_int[TEST_WRITEARRAY_1_LEN];
  int   arrayBuf_zeros[TEST_WRITEARRAY_1_LEN];
  float arrayBuf_float[TEST_WRITEARRAY_1_LEN];
  
  // set some values
  for (int i=0;i<TEST_WRITEARRAY_1_LEN;i++) {
    arrayBuf_char[i]=i;
    arrayBuf_int[i]=10*i;
    arrayBuf_float[i]=100.0*i;
    arrayBuf_zeros[i]=0;
  }

  std::string filename="test-writearray-1";

  CL_File* clfile=new CL_File();
	CL_DataGroup* dg;
  CL_DataSet* ds;

  dg=clfile->newDataGroup("DG-1");
  ds=dg->newDataSet("DS-1-1");
  ds->newColumn("col-1",CL_TC_UBYTE);
  ds->newColumn("col-2",CL_TC_INT);
  ds->newColumn("col-3",CL_TC_FLOAT);
  ds->setRowCount(TEST_WRITEARRAY_1_LEN);

  //
  ds->setFromArray(0,0,arrayBuf_char,TEST_WRITEARRAY_1_LEN);
  ds->get(TEST_WRITEARRAY_1_LEN,0,&tmp_int);
  //
  ds->setFromArray(0,1,arrayBuf_int,TEST_WRITEARRAY_1_LEN);
  ds->get(TEST_WRITEARRAY_1_LEN,0,&tmp_int);
  //
  ds->setFromArray(0,2,arrayBuf_float,TEST_WRITEARRAY_1_LEN);
  ds->get(TEST_WRITEARRAY_1_LEN,0,&tmp_int);
  
  ds=dg->newDataSet("DS-1-2");
  ds->newColumn("col-1",CL_TC_UBYTE);
  ds->newColumn("col-2",CL_TC_INT);
  ds->newColumn("col-3",CL_TC_FLOAT);
  ds->setRowCount(TEST_WRITEARRAY_1_LEN);

  // rotate the cols to test conversions.
  ds->setFromArray(0,1,arrayBuf_char,TEST_WRITEARRAY_1_LEN);
  ds->get(TEST_WRITEARRAY_1_LEN,0,&tmp_int);
  //
  ds->setFromArray(0,2,arrayBuf_int,TEST_WRITEARRAY_1_LEN);
  ds->get(TEST_WRITEARRAY_1_LEN,0,&tmp_int);
  //
  ds->setFromArray(0,0,arrayBuf_float,TEST_WRITEARRAY_1_LEN);
  ds->get(TEST_WRITEARRAY_1_LEN,0,&tmp_int);

  // These shouldnt overwrite the data; they are out of bounds.
  assert(ds->setFromArray(-1,1,arrayBuf_zeros,TEST_WRITEARRAY_1_LEN)==CL_ERR);
  assert(ds->setFromArray(100,1,arrayBuf_zeros,TEST_WRITEARRAY_1_LEN)==CL_ERR);

  // exercise the resizing...
  dg=clfile->newDataGroup("DG-2");
  ds=dg->newDataSet("DS-2-1");

  int cidx=0;
  for (int i=0;i<50;i++) {
    ds->newColumn("col-"+ToStr(cidx++),CL_TC_UBYTE);
    ds->setFromArray(0,cidx-1,arrayBuf_char,TEST_WRITEARRAY_1_LEN);
    ds->newColumn("col-"+ToStr(cidx++),CL_TC_UINT);
    ds->setFromArray(0,cidx-1,arrayBuf_int,TEST_WRITEARRAY_1_LEN);
    ds->newColumn("col-"+ToStr(cidx++),CL_TC_TEXT_ASCII8,i*10);
    //ds->setRowCount(0);
    ds->setRowCount(TEST_WRITEARRAY_1_LEN);
  }

  ds=dg->newDataSet("DG-3");

  ds=dg->newDataSet("DS-3-1");
  ds->newColumn("col-1",CL_TC_UBYTE);
  ds->newColumn("col-2",CL_TC_INT);
  ds->newColumn("col-3",CL_TC_FLOAT);
  ds->setRowCount(0);

  ds=dg->newDataSet("DG-4");
  ds=dg->newDataSet("DG-5");
  ds=dg->newDataSet("DG-6");

  clfile->write_File(filename+".calvin");
  clfile->close();
  delete clfile;
}

//////////

void test_writecychpfile_1()
{
  std::string filename="foo.cychp";
  CL_File* clfile;

  clfile=new CL_File();

  CL_Gdh* gdh_0 = clfile->getGdh();
  gdh_0->defaultValues();
  // in order
  gdh_0->setDataType("affymetrix-multi-data-type-analysis");
  gdh_0->setUuid("3925ed53-077a-48ca-2189-8e84f4541dbb");
  gdh_0->setDateTime("2009-07-15T17:58:34Z");
  gdh_0->setLocale("en-US");
  //
  gdh_0->addParam("affymetrix-algorithm-name", "CYTO2");

  //
  clfile->write_File(filename);
  
  if (clfile->errNum()!=CL_OK) {
    clfile->printErrMsg();
    printf("Could not write calvin file: '%s'\n",filename.c_str());
    assert(0);
  }
  
  delete clfile;
  return;
}

//////////

void test_writecychpfile_2()
{
	//int outputfileLen = (int)mxGetN(OUTPUTFILE) + 1;
	//char *outputfile = new char [outputfileLen];
	//if (mxGetString(OUTPUTFILE, outputfile, outputfileLen) != 0) {
	//	mexErrMsgTxt("Failed to copy output file name. Error in WriteCyChpFile");
	//	return;
	//}
	//std::string filename(outputfile);
	//delete outputfile;
	std::string filename="test_writecychpfile_2.cychp";
  char buf[100];

	CL_File* clfile;

	clfile=new CL_File();
	CL_Gdh* gdh_0 = clfile->getGdh();
	//gdh_0->defaultValues();
	gdh_0->setDataType("affymetrix-multi-data-type-analysis");

  //gdh_0->dump();

	std::string stringval = "CYTO2";
	gdh_0->addParam("affymetrix-algorithm-name", stringval);

#define getString(_x) _x

	std::string stringval0 = getString("ALGORITHMVERSION");
	gdh_0->addParam("affymetrix-algorithm-version", stringval0);
	std::string stringval1 = getString("CHIPTYPE");
	gdh_0->addParam("affymetrix-array-type", stringval1);
	gdh_0->addParam("affymetrix-algorithm-param-ArraySet", stringval1);
	std::string stringval2 = getString("PROGRAMNAME");
	gdh_0->addParam("program-name", stringval2);
	std::string stringval3 = getString("PROGRAMVERSION");
	gdh_0->addParam("program-version", stringval3);
	std::string stringval4 = getString("PROGRAMCOMPANY");
	gdh_0->addParam("program-company", stringval4);
	std::string stringval5 = getString("DATECREATED");
	gdh_0->addParam("create-date", stringval5);

	/* CL_Gdh* gdh_1= */ gdh_0->newLeafGdh();
	//gdh_1->defaultValues();

	CL_DataGroup* dg_1;
	CL_DataSet* dg_1_ds_1;
	CL_DataSet* dg_1_ds_2;
  
	int drow_cnt = 10;

	// Chromosomes group
	//   Summary: Chromosome Ubyte, Display ASCII-4, StartIndex UInt, MarkerCount UInt,
	//            MinSignal Float, MaxSignal Float, MedianCNState Float, HomFrequency Float, HetFrequency Float,
	//            Mosaicism Float, LOH Float
	dg_1=clfile->newDataGroup("Chromosomes");
	dg_1_ds_1=dg_1->newDataSet("Summary");
	int col1 = dg_1_ds_1->newColumn("Chromosome", CL_TC_UBYTE);
	//int col2 = dg_1_ds_1->newColumn("Display", CL_TC_TEXT_ASCII8, 4);
	/* int col3  = */ dg_1_ds_1->newColumn("StartIndex", CL_TC_UINT);
	/* int col4  = */ dg_1_ds_1->newColumn("MarkerCount", CL_TC_UINT);
	/* int col5  = */ dg_1_ds_1->newColumn("MinSignal", CL_TC_FLOAT);
	/* int col6  = */ dg_1_ds_1->newColumn("MaxSignal", CL_TC_FLOAT);
	/* int col7  = */ dg_1_ds_1->newColumn("MedianCNState", CL_TC_FLOAT);
	/* int col8  = */ dg_1_ds_1->newColumn("HomFrequency", CL_TC_FLOAT);
	/* int col9  = */ dg_1_ds_1->newColumn("HetFrequency", CL_TC_FLOAT);
	/* int col10 = */ dg_1_ds_1->newColumn("Mosaicism", CL_TC_FLOAT);
	/* int col11 = */ dg_1_ds_1->newColumn("LOH", CL_TC_FLOAT);
  //
  int col12 = dg_1_ds_1->newColumn("Display", CL_TC_TEXT_ASCII8, 40);
  //
	dg_1_ds_1->setRowCount(drow_cnt);

  //
  for (int i=0;i<drow_cnt;i++) {
    dg_1_ds_1->set(i,col1,i);
    //
    sprintf(buf,"string-%03d-gnirts",i);
    dg_1_ds_1->set(i,col12,buf);
  }
  
  if (1) {
    dg_1_ds_2=dg_1->newDataSet("AllelePeaks");
    int col_2_1 = dg_1_ds_2->newColumn("Chromosome", CL_TC_UBYTE);
    int col_2_2 = dg_1_ds_2->newColumn("Display", CL_TC_TEXT_ASCII8, 4);
    /* int col_2_3 = */ dg_1_ds_2->newColumn("StartIndex", CL_TC_UINT);
    
    dg_1_ds_2->setRowCount(drow_cnt);
    
    for (int i=0;i<drow_cnt;i++) {
      dg_1_ds_2->set(i,col_2_1,i);
      //
      sprintf(buf,"%03d",i);
      dg_1_ds_2->set(i,col_2_2,buf);
    }
  }

	// ProbeSets group
	//   CopyNumber table
	//   AllelePeaks table
	// AlgorithmData
	//   MarkerABSignal table
	// Segments
	//   CN table
	//   LOH table
	//   CNNeutralLOH table
	//   NormalDiploid table
	//   Mosaicism table

	clfile->write_File(filename);

	if (clfile->errNum()!=CL_OK) {
    clfile->printErrMsg();
	}
	delete clfile;
}

void test_writecychpfile_3()
{
  CL_File* clfile;

  clfile=new CL_File();

  CL_Gdh* gdh_0 = clfile->getGdh();
  gdh_0->setDataType("test_writecychpfile_3");
  //
  gdh_0->addParam("int-00", 0x00);
  gdh_0->addParam("int-0F", 0x0F);
  gdh_0->addParam("int-F0", 0xF0);
  gdh_0->addParam("int-FF", 0xFF);
  gdh_0->addParam("int-FFFF", 0xFFFF);
  gdh_0->addParam("int-FFFFFFFF",(int)0xFFFFFFFF);
  //
  gdh_0->addParam("string-null","");
  gdh_0->addParam("string-A","A");
  gdh_0->addParam("string-AB","AB");
  gdh_0->addParam("string-ABC","ABC");
  gdh_0->addParam("string-ABCDEFGHJ","ABCDEFGHJ");
  //
  gdh_0->addParam("double-0", 0.0);
  gdh_0->addParam("double-100", 100.0);
  //

  // now that the headers are there, lets add some data.
  CL_DataGroup* dg_1=clfile->newDataGroup("AlgorithmData");
  CL_DataSet*   dg_1_ds_1=dg_1->newDataSet("MarkerABSignal");

  /* int col_1_1 = */ dg_1_ds_1->newColumn("Index", CL_TC_UINT); //1
  /* int col_1_2 = */ dg_1_ds_1->newColumn("ASignal", CL_TC_FLOAT); //1
  /* int col_1_3 = */ dg_1_ds_1->newColumn("BSignal", CL_TC_FLOAT); //1
  /* int col_1_4 = */ dg_1_ds_1->newColumn("SCAR", CL_TC_FLOAT); //1
  int col_1_5 = dg_1_ds_1->newColumn("Strings", CL_TC_TEXT_ASCII8,64); //1

  /* CL_DataSet*   dg_1_ds_2= */ dg_1->newDataSet("MarkerABSignal");

  //
  /* CL_DataGroup* dg_2= */ clfile->newDataGroup("Null-2-DS");

  //
  CL_DataGroup* dg_3=clfile->newDataGroup("Null-3-DG");
  /* CL_DataSet*   dg_3_ds_1= */ dg_3->newDataSet("Null-3-DS");

  //
  for (int i=0;i<5;i++) {
    int rowcount=i*10*1000;
    //dg_1_ds_1->setRowCount(0);
    dg_1_ds_1->setRowCount(rowcount);
    for (int r=0;r<rowcount;r++) {
      dg_1_ds_1->set(r,0,r);
      dg_1_ds_1->set(r,col_1_5,ToStr(r));
    }
    //
    clfile->write_File("test_writecychpfile_3-"+ToStr(rowcount)+".cychp");
    clfile->dumpSegs();
  }

  delete clfile;
}

//////////

int main(int argc,const char* argv[])
{
  const char** args=argv;
  int opt_write=0;

  test_err_handling();

  // return 0;
  // CL_File_debug_flags=1;

  test_writearray_1();
  //return 0;

  // tests written while getting cychip files to work.
  // test_writecychpfile_1();
  // test_writecychpfile_2();
  // test_writecychpfile_3();

  //
  //return 0;

  // skip argv[0], the program name
  args++;
  // done?
  if (*args==NULL) {
    return 0;
  }

  //
  if (strcmp(*args,"--dump")==0) {
    args++;
    return CL_File::dumpFile(*args);
  }

  //
  if (strcmp(*args,"--rw")==0) {
    opt_write=1;
    args++;
    //
    while (*args!=NULL) {
      test_rw_file(*args,opt_write);
      args++;
    }
    return 0;
  }

  //
  if ((*args!=NULL)&&(strcmp(*args,"--test-small-1")==0)) {
    args++;
    test_small_file_1(*args);
  }
  if ((*args!=NULL)&&(strcmp(*args,"--test-small-2")==0)) {
    args++;
    test_small_file_2(*args);
  }
  if ((*args!=NULL)&&(strcmp(*args,"--test-small-3"))==0) {
    args++;
    test_small_file_3(*args);
  }

  if ((*args!=NULL)&&(strcmp(*args,"--uuid")==0)) {
    args++;
    test_uuid();
  }

  //
  return 0;
}
