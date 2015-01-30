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
// affy/sdk/chipstream/KitAODb.cpp ---
//

#include "chipstream/KitAODb.h"
//
#include <set>
//
#include "file/TsvFile/TsvFile.h"
#include "util/Err.h"
#include "util/Verbose.h"

KitAODbEntry::KitAODbEntry() {
  m_linenum=0;
  m_pc1=0.0;
  m_means=0.0;
}
KitAODbEntry::KitAODbEntry(int linenum,float pc1,float means) {
  m_linenum=0;
  m_pc1=pc1;
  m_means=means;
}

//////////

KitAODb::KitAODb() {
  clear();
}

KitAODb::~KitAODb() {
  clear();
}

void KitAODb::clear() {
  m_classiferAValue=-9999.9;
  m_classiferOValue=-9999.9;
  //
  m_filename.clear();

  // delete the entries in the set.
  for (entry_iter_t i=m_entries.begin();
       i!=m_entries.end();
       i++) {
    delete *i;
  }
  m_entries.clear();
  //
  m_name2entry.clear();
  //
  m_added_count=0;
  m_file_probeset_count = 0;
}

//////////

double KitAODb::getClassiferAValue() const {
  return m_classiferAValue;
}
double KitAODb::getClassiferOValue() const {
  return m_classiferOValue;
}

KitAODbEntry* KitAODb::findEntry(const std::string& e_name) {
  KitAODb::name2entry_iter_t i=m_name2entry.find(e_name);
  if (i==m_name2entry.end()) {
    return NULL;
  }
  return i->second;
}

void KitAODb::addEntry(const std::string& e_name,KitAODbEntry* e_ptr) {
  // bump the number of additions made.
  m_added_count++;
  //
  KitAODbEntry* found_ptr=findEntry(e_name);
  //printf("%d: %s => %p\n",m_entry_count,e_name.c_str(),e_ptr);
  APT_ERR_ASSERT(((found_ptr==NULL)||(found_ptr==e_ptr)),"duplicate name '"+e_name+"' found.");
  // done if it is already in the map...
  if (found_ptr!=NULL) {
    return;
  }
  // ...not in the maps, add it.
  name2entry_val_t e_pair(e_name,e_ptr);
  m_name2entry.insert(e_pair);
  //
  m_entries.insert(e_ptr);
}

//////////

int KitAODb::addedCount() const {
  return m_added_count;
}

int KitAODb::nameCount() const {
  return m_name2entry.size();
}

int KitAODb::entryCount() const {
  return m_entries.size();
}

//////////

void KitAODb::readTsv(const std::string& filename) {
  m_filename=filename;

  Verbose::out(1,"Reading KitAODb: '"+m_filename+"'...");

  affx::TsvFile tsv;
  APT_ERR_ASSERT(tsv.open(m_filename)==affx::TSV_OK,"open of '"+m_filename+"' failed.");

  // check the version
  int kit_ao_format_version=-1;
  if (tsv.getHeader("db-format-version",kit_ao_format_version)==affx::TSV_OK) {
    APT_ERR_ASSERT(kit_ao_format_version==1,"we only read v1 format files.");
  }

  //
  m_classiferAValue=0.0;
  APT_ERR_ASSERT((tsv.getHeader("classifier-cutoff-a-value",m_classiferAValue)==affx::TSV_OK),
                 "The header key 'classifer-cutoff-a-value' is required.");
  m_classiferOValue=0.0;
  APT_ERR_ASSERT((tsv.getHeader("classifier-cutoff-o-value",m_classiferOValue)==affx::TSV_OK),
                 "The header key 'classifer-cutoff-o-value' is required.");
  // 
  if (tsv.getHeader("probeset-count",m_file_probeset_count)!=affx::TSV_OK) {
    m_file_probeset_count=0;
  }

  // we require the these two columns (PC1 & means)
  int pc1_cidx=tsv.cname2cidx(0,"pc1","PC1");
  APT_ERR_ASSERT((pc1_cidx>=0),"the 'pc1' column was not found.");

  int means_cidx=tsv.cname2cidx(0,"means");
  APT_ERR_ASSERT((means_cidx>=0),"the 'means' column was not found.");

  // for compatiablity with other versions of this file.
  int probeset_id_cidx=tsv.cname2cidx(0,"probeset_id","GS_probeset_id","CEU_probeset_id");
  // and we require a column of probeset names.
  APT_ERR_ASSERT(probeset_id_cidx!=-1,"File did not have a probeset id column.");

  //
  std::string entry_name;
  while (tsv.nextLevel(0)==affx::TSV_OK) {
    KitAODbEntry* entry_ptr=new KitAODbEntry();
    APT_ERR_ASSERT(entry_ptr!=NULL,"new KitAODbEntry() failed.");
    entry_ptr->m_linenum=tsv.lineNumber();
    tsv.get(0,pc1_cidx,entry_ptr->m_pc1);
    tsv.get(0,means_cidx,entry_ptr->m_means);
    tsv.get(0,probeset_id_cidx,entry_name);
    //
    this->addEntry(entry_name,entry_ptr);
  }

  tsv.close();

  // check that we read the correct number.
  if (m_file_probeset_count>0) {
    APT_ERR_ASSERT(m_file_probeset_count==addedCount(),"Didnt read the correct number of items from the input file.");
  }
}

void KitAODb::writeTsv(const std::string& filename) {
  affx::TsvFile tsv;

  tsv.addHeaderComment("## This is a KitAO DB file.");
  tsv.addHeaderComment("");
  tsv.addHeaderComment("The version of the file format.");
  tsv.addHeader("db-format-version",1);

  tsv.addHeaderComment("");
  tsv.addHeaderComment("The file which was read in.");
  tsv.addHeader("original-filename",m_filename);

  tsv.addHeaderComment("");
  tsv.addHeaderComment("The cutoff value.");
  tsv.addHeader("classifier-cutoff-a-value",m_classiferAValue);
  tsv.addHeader("classifier-cutoff-o-value",m_classiferOValue);

  tsv.addHeaderComment("");
  tsv.addHeaderComment("The number of probesets.");
  tsv.addHeader("probeset-count",addedCount());

  tsv.addHeaderComment("");
  tsv.defineColumn(0,0,"probeset_id");
  tsv.defineColumn(0,1,"PC1");
  tsv.defineColumn(0,2,"means");
  
  tsv.writeTsv_v1(filename);

  for (name2entry_iter_t i=m_name2entry.begin();
       i!=m_name2entry.end();
       i++) {
    //
    tsv.set(0,0,i->first);
    tsv.set(0,1,i->second->m_pc1);
    tsv.set(0,2,i->second->m_means);
    tsv.writeLevel(0);
  }
  tsv.close();
}
