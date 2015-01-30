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
// affy/sdk/bboard/dao/DaoUtil.cpp ---
//
// $Id: DaoUtil.cpp,v 1.1 2009-11-05 20:41:39 harley Exp $
//

#ifdef _MSC_VER
#define _CRT_SECURE_NO_WARNINGS
#endif

//
#include "bboard/dao/DaoUtil.h"
//
//#include "bboard/dao/DaoTable.h"
#include "bboard/dao/DaoTypes.h"
#include "bboard/BboardTypes.h"
#include "portability/affy-system-api.h"
//
#include <stdio.h>

//////////

#define DAOUTIL_DEFINE_OBJ(_type)                                       \
  AptErr_t DaoUtil::readFromFile(const std::string& path,_type* obj) {  \
    return readFromFile(FsPath(path),obj);                             \
  };                                                                    \
  AptErr_t DaoUtil::readFromFile(const FsPath& dao_p,_type* obj) {     \
    Dao_File*   dao_f=new Dao_File();                                   \
    dao_f->open(dao_p,0);                                               \
    Dao_Group*  dao_g=dao_f->openGroup(dao_p.getInternalGroupName(),0); \
    Dao_Table*  dao_t=dao_g->openTable(dao_p.getInternalTableName(),0); \
                                                                        \
    AptErr_t rv=readFromTable(dao_t,obj);                               \
                                                                        \
    delete dao_t;                                                       \
    dao_g->close();                                                     \
    delete dao_g;                                                       \
    dao_f->close();                                                     \
    delete dao_f;                                                       \
    return rv;                                                          \
  }                                                                     \
  AptErr_t DaoUtil::writeToFile(const std::string& path,_type* obj) {   \
    return writeToFile(FsPath(path),obj);                              \
  }                                                                     \
  AptErr_t DaoUtil::writeToFile(const FsPath& dao_p,_type* obj) {      \
    Dao_File*   dao_f=new Dao_File();                                   \
    dao_f->create(dao_p,0);                                             \
    Dao_Group*  dao_g=dao_f->createGroup(dao_p.getInternalGroupName(),0); \
    Dao_Table*  dao_t=dao_g->createTable(dao_p.getInternalTableName(),0); \
                                                                        \
    AptErr_t rv=writeToTable(dao_t,obj);                                \
                                                                        \
    delete dao_t;                                                       \
    dao_g->close();                                                     \
    delete dao_g;                                                       \
    dao_f->close();                                                     \
    delete dao_f;                                                       \
    return rv;                                                          \
  }

DAOUTIL_DEFINE_OBJ(Bboard);
DAOUTIL_DEFINE_OBJ(PgOptions);

//////////

std::string gen_sub_table_name(const std::string& prefix,int& cnt)
{
  char buf[1000];
  snprintf(buf,sizeof(buf),"%02d",cnt++);
  std::string sub_table_name=prefix;
  if (sub_table_name!="") {
    sub_table_name+="-";
  }
  sub_table_name+=buf;
  return sub_table_name;
}

//////////

AptErr_t DaoUtil::readFromTable(Dao_Table* dao_t,Bboard *bb)
{
  std::string key_str;
  std::string type_str;
  std::string val_str;

  dao_t->rewind();

  while (dao_t->nextRow()==APT_OK) {
    dao_t->get(0,"key",&key_str);
    dao_t->get(0,"type",&type_str);
    dao_t->get(0,"value",&val_str);
    //
    BboardType_t bbt=BboardTypeFromString(type_str);
    //
    if (0) {
      printf("readFromTsvFile(): '%s','%s'==%d,'%s'\n",
             key_str.c_str(),type_str.c_str(),bbt,val_str.c_str());
    }
    //
    switch (bbt) {
    case BBT_CHAR:
      {
        char val_char;
        sscanf(val_str.c_str(),"%c",&val_char);
        bb->set(key_str,val_char);
      }
      break;
    case BBT_DOUBLE:
      {
        float val_float;
        sscanf(val_str.c_str(),"%f",&val_float);
        double val_double=val_float;
        bb->set(key_str,val_double);
      }
      break;
    case BBT_FLOAT:
      {
        float val_float;
        sscanf(val_str.c_str(),"%f",&val_float);
        bb->set(key_str,val_float);
      }
      break;
    case BBT_INT:
      {
        int val_int;
        sscanf(val_str.c_str(),"%d",&val_int);
        bb->set(key_str,val_int);
      }
      break;
    case BBT_STRING:
      bb->set(key_str,val_str);
      break;
      //
    case BBT_VEC_INT:
      {
        std::vector<int>* vecptr=new std::vector<int>();
        bb->setPtrAndType(key_str,vecptr,BBT_VEC_INT);
        // bleah -- should "#import sdk/util/StringUtil.h"
        int tmp_int;
        int scan_cnt;
        const char* ptr=val_str.c_str();
        while (*ptr!=0) {
          scan_cnt=0;
          if (sscanf(ptr,"%d%n",&tmp_int,&scan_cnt)!=0) {
            //printf("'%d' %d\n",tmp_int,scan_cnt);
            vecptr->push_back(tmp_int);
          }
          ptr+=scan_cnt;
          if (*ptr==0) {
            break;
          }
          ptr++; // ",";
        }
      }
      break;
      //
    case BBT_BBOARD:
      {
        Bboard* sub_bb=new Bboard(val_str);
        Dao_Table* sub_dao_t=dao_t->getGroup()->openTable(val_str,0);
        DaoUtil::readFromTable(sub_dao_t,sub_bb);
        bb->set(key_str,sub_bb);
        sub_dao_t->close();
        delete sub_dao_t;
      }
      break;
      //
    case BBT_DAO_FILE:
      {
        printf("readFromDaoGroup():BBT_DAO_FILE: '%s'\n",val_str.c_str());
        Dao_File* sub_dao_f=new Dao_File();
        sub_dao_f->m_dao_path.setPath(val_str);
        bb->set(key_str,sub_dao_f);
      }
      break;
    case BBT_PGOPTIONS:
      {
        PgOptions* sub_pgopts;
        bb->get(key_str,&sub_pgopts);
        if (sub_pgopts==NULL) {
          printf("the bboard needs an options obj called '%s' to read options.\n",key_str.c_str());
        }
        else {
          Dao_Table* sub_dao_t=dao_t->getGroup()->openTable(val_str,0);
          readFromTable(sub_dao_t,sub_pgopts);
          delete sub_dao_t;
        }
        // Well... Here we have a bit of a wrinkle...
        // We only want to read the option/value pairs, we dont want to 
        // construct the entire PgOption object.  The program should have it prebuilt
        // and passed in.
        // for now we ignore it.
      }
      break;
    default:
      printf("unhandled: '%s','%s'==%d,'%s'\n",
             key_str.c_str(),type_str.c_str(),bbt,val_str.c_str());
      break;
    };
  }

  dao_t->close();
  //
  return APT_OK;
}

AptErr_t DaoUtil::writeToTable(Dao_Table* dao_t,Bboard* bb)
{
  // debug
  bb->dump();
  //
  dao_t->clearSchema();
  //
  dao_t->addHeader("dao-table-type","blackboard");
  dao_t->addHeader("dao-version",1);
  //
  dao_t->defineColumn(0,0,"key",DAO_STRING);
  dao_t->defineColumn(0,1,"type",DAO_STRING);
  dao_t->defineColumn(0,2,"value",DAO_STRING);

  dao_t->endHeaders();

  //
  int sub_table_cnt=0;
  for (Bboard::NameRefMap_t::iterator i=bb->m_name_ref_map.begin();
       i!=bb->m_name_ref_map.end();
       ++i) {
    //
    BboardBox* bbbox=i->second.boxPtr();
    //
    dao_t->set(0,0,i->first);
    dao_t->set(0,1,BboardTypeToString(bbbox->getBbType()));

    switch (bbbox->getBbType()) {
    case BBT_CHAR:
    case BBT_DOUBLE:
    case BBT_FLOAT:
    case BBT_INT:
    case BBT_STRING:
      dao_t->set(0,2,bbbox->asString());
      break;
      //
    case BBT_VEC_INT:
      // There should be some smarts here.
      // If the vector is small put it in as a comma list, "vec-int"
      // otherwise write it as a "vec-int-file" and put the filename here.
      // for now just put it all here.
      {
        std::string outstr;
        char buf[100];
        void* voidptr;
        std::vector<int>* vecptr;
        //
        bbbox->getPtrIfType(voidptr,BBT_VEC_INT);
        vecptr=(std::vector<int>*)voidptr;
        for (int i=0;i<vecptr->size();i++) {
          if (i!=0) {
            outstr+=",";
          }
          sprintf(buf,"%d",(*vecptr)[i]);
          outstr+=buf;
        }
        //
        dao_t->set(0,2,outstr);
      }
      break;

      // For complex types, we use an external file.
      // For example if we have a sub blackboard, we dump that to another file.
    case BBT_BBOARD:
      {
        Bboard* sub_bb;
        bbbox->getData(&sub_bb);
        // loop detection
        assert(sub_bb!=bb);
        assert(sub_bb!=NULL);
        sub_bb->dump();
        

        // generate the sub table_name
        std::string sub_table_name=gen_sub_table_name(dao_t->getName(),sub_table_cnt);
        // printf("BBT_BBOARD: '%s'\n",sub_table_name.c_str());
        dao_t->set(0,2,sub_table_name);
        //
        Dao_Table* sub_dao_t=dao_t->getGroup()->createTable(sub_table_name,0);
        writeToTable(sub_dao_t,sub_bb);
        delete sub_dao_t;
      }
      break;
      //
    case BBT_PGOPTIONS:
      {
        //
        PgOptions* sub_pgopts;
        bbbox->getData(&sub_pgopts);
        //
        std::string sub_table_name=gen_sub_table_name(dao_t->getName(),sub_table_cnt);
        dao_t->set(0,2,sub_table_name);
        //
        Dao_Table* sub_dao_t=dao_t->getGroup()->createTable(sub_table_name,0);
        writeToTable(sub_dao_t,sub_pgopts);
        delete sub_dao_t;
      }
      break;
      //
    case BBT_DAO_FILE:
      {
        void* tmp_voidptr;
        Dao_File* tmp_dao_file;
        bbbox->getPtrIfType(tmp_voidptr,BBT_DAO_FILE);
        tmp_dao_file=(Dao_File*)tmp_voidptr;
        dao_t->set(0,2,tmp_dao_file->m_dao_path.asUnixPath());
      }
      break;
      //
    default:
      dao_t->set(0,2,"*opps! Bad type*");
      break;
    }
    //
    dao_t->writeRow();
    // for debugging.
    dao_t->flush();
  }

  //
  dao_t->close();
  return APT_OK;
}

//////////

AptErr_t DaoUtil::readFromTable(Dao_Table* dao_t,PgOptions* pgopts)
{
  std::string opt_name;
  std::string opt_type;
  std::string opt_value;

  assert(pgopts!=NULL);

  dao_t->rewind();

  while (dao_t->nextRow()==APT_OK) {
    dao_t->get(0,0,&opt_name);
    dao_t->get(0,1,&opt_type);
    dao_t->get(0,2,&opt_value);
    //
    PgOpt* opt=pgopts->findOpt(opt_name);
    if (opt==NULL) {
      printf("readFromTable(): Option '%s' not found!\n",opt_name.c_str());
      //assert(0);
      continue;
    }
    //
    printf("setting option '%s' to '%s'\n",opt_name.c_str(),opt_value.c_str());
    opt->setValue(opt_value);
  }
  dao_t->close();
  //
  return APT_OK;
}

AptErr_t DaoUtil::writeToTable(Dao_Table* dao_t,PgOptions* pgopts)
{
  dao_t->addHeader("dao-table-type","pgoptions");
  dao_t->addHeader("dao-version",1);

  dao_t->defineColumn(0,0,"option",DAO_STRING);
  dao_t->defineColumn(0,1,"type",DAO_STRING);
  dao_t->defineColumn(0,2,"value",DAO_STRING);

  dao_t->endHeaders();

  // no iterator?
  for (int i=0;i<pgopts->m_option_vec.size();i++) {
    PgOpt* opt=pgopts->m_option_vec[i];
    //
    if (opt->m_type==PgOpt::INVALID_OPT) {
      continue;
    }
    //
    dao_t->set(0,0,opt->m_longName);
    //
    switch (opt->m_type) {
    case PgOpt::INVALID_OPT:
      assert(0);
      break;
    case PgOpt::BOOL_OPT:
      dao_t->set(0,1,"bool");
      break;
    case PgOpt::DOUBLE_OPT:
      dao_t->set(0,1,"double");
      break;
    case PgOpt::INT_OPT:
      dao_t->set(0,1,"int");
      break;
    case PgOpt::STRING_OPT:
      dao_t->set(0,1,"string");
      break;
    default:
      assert(0);
    }
    //
    dao_t->set(0,2,opt->getValue());
    //
    dao_t->writeRow();
  }

  //
  dao_t->close();
  return APT_OK;
}
