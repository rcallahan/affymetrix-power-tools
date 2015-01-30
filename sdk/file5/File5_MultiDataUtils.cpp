////////////////////////////////////////////////////////////////
//
// Copyright (C) 1989, 1991 Free Software Foundation, Inc.
//
// This library is free software; you can redistribute it and/or modify
// it under the terms of the GNU Lesser General Public License
// (version 2.1) as published by the Free Software Foundation.
//
// This library is distributed in the hope that it will be useful, but
// WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
// or FITNESS FOR A PARTICULAR PURPOSE. See the GNU Lesser General Public License
// for more details.
//
// You should have received a copy of the GNU Lesser General Public License
// along with this library; if not, write to the Free Software Foundation, Inc.,
// 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
//
////////////////////////////////////////////////////////////////

#include "file5/File5_MultiDataUtils.h"

using namespace affx;

#define TOC_GROUP_NAME "TOC"

File5_MultiDataUtils::File5_MultiDataUtils()
{
	file5 = NULL;
	tsv5 = NULL;
	colIdx = -1;
}

File5_MultiDataUtils::~File5_MultiDataUtils()
{
	Close();
}

void File5_MultiDataUtils::Close()
{
	CloseTsv();
	
	if (file5 != NULL && file5->is_open() == 1)
		file5->close();
	delete file5;
	file5 = NULL;
}

void File5_MultiDataUtils::CloseTsv()
{
	colIdx = -1;
	if (tsv5 != NULL && tsv5->is_open() == 1)
		tsv5->close();
	delete tsv5;
	tsv5 = NULL;
  // it is safe to call close() on a closed tsv.
  /*
  if (tsv5!=NULL) {
		tsv5->close();
    delete tsv5;
    tsv5 = NULL;
  }
  */
}

void File5_MultiDataUtils::Open(const std::string &fileName)
{
	try
	{
		Close();
		file5 = new File5_File;
		file5->open(fileName, FILE5_OPEN|FILE5_CREATE);
	}
	catch (std::exception e)
	{
		Close();
		throw e;
	}
}

/////////////////////////////////////////WRITE FUNCTIONS///////////////////////////////////////////////////////

void File5_MultiDataUtils::CreateTsv(const std::string &gname, const std::string &dname)
{
	CloseTsv();
	tsv5 = file5->openTsv(gname, FILE5_OPEN|FILE5_CREATE);
	File5_Group *toc5 = file5->openGroup(TOC_GROUP_NAME, FILE5_OPEN|FILE5_CREATE);
	toc5->headersBegin();
	std::string key;
	std::string val;
	bool found = false;
	while (toc5->headersNext(&key, &val) == FILE5_OK)
	{
		if (gname == val)
		{
			found = true;
			break;
		}
	}
	if (found == false)
		toc5->addHeader(dname, gname);
	toc5->close();
	delete toc5;
}

void File5_MultiDataUtils::DefineColumn(const std::string& cname, const File5_dtype_t ctype, int csize)
{
	colIdx = tsv5->getColumnCount(0);
	if (colIdx == -1)
		colIdx = 0;
	tsv5->defineColumn(0, colIdx, cname, ctype, csize);
}

void File5_MultiDataUtils::AddEntry(int ival)
{
	tsv5->set_i(0, colIdx, ival);
	tsv5->writeLevel(0, colIdx);
}

void File5_MultiDataUtils::AddEntry(float fval)
{
	tsv5->set_f(0, colIdx, fval);
	tsv5->writeLevel(0, colIdx);
}

void File5_MultiDataUtils::AddEntry(const std::string &sval)
{
	tsv5->set_string(0, colIdx, sval);
	tsv5->writeLevel(0, colIdx);
}

void File5_MultiDataUtils::AddParameter(const std::string &key, const std::string &val)
{
	if (tsv5 != NULL)
		tsv5->addHeader(key, val);
	else if (file5 != NULL)
		file5->addHeader(key, val);
}

/////////////////////////////////////////READ FUNCTIONS///////////////////////////////////////////////////////

int File5_MultiDataUtils::TsvLineCount()
{
  // @todo assert(tsv5!=NULL);
	return tsv5->getLineCount();
}

template<class T>
static bool NextLine(File5_Tsv *tsv5, std::string& name, std::vector<T> &data)
{
	int n = tsv5->getColumnCount(0);
	if (tsv5->nextLine() != FILE5_OK)
		return false;
	tsv5->get(0, 0, &name);
	data.resize(n-1);
	for (int i=1; i<n; i++)
		tsv5->get(0, i, &data[i-1]);
	return true;
}

bool File5_MultiDataUtils::NextLine(std::string &name, std::vector<int> &data)
{
	return ::NextLine(tsv5, name, data);
}

bool File5_MultiDataUtils::NextLine(std::string &name, std::vector<float> &data)
{
	return ::NextLine(tsv5, name, data);
}

bool File5_MultiDataUtils::NextLine(std::string &name, std::vector<std::string> &data)
{
	return ::NextLine(tsv5, name, data);
}

bool File5_MultiDataUtils::NextLineAsString(std::string &name, std::vector<std::string> &data)
{
	int n = tsv5->getColumnCount(0);
	if (tsv5->nextLine() != FILE5_OK)
		return false;
	tsv5->get(0, 0, &name);
	data.resize(n-1);
	for (int i=1; i<n; i++)
		tsv5->getAsStr(0, i, &data[i-1]);
	return true;
}

void File5_MultiDataUtils::GetTOC(std::vector<TableOfContentsEntry> &toc)
{
	toc.clear();
	if (file5->name_exists(TOC_GROUP_NAME) == false)
		return;
	File5_Group *toc5 = file5->openGroup(TOC_GROUP_NAME, FILE5_OPEN);
	toc5->headersBegin();
	TableOfContentsEntry entry;
	while (toc5->headersNext(&entry.displayName, &entry.groupName) == FILE5_OK)
	{
		if (entry.displayName.empty() == false)
			toc.push_back(entry);
	}
	toc5->close();
	delete toc5;
}

void File5_MultiDataUtils::OpenTsv(const std::string &gname)
{
	CloseTsv();
	tsv5 = file5->openTsv(gname, FILE5_OPEN|FILE5_CREATE);
}

File5_return_t File5_MultiDataUtils::DeleteTsv(const std::string &name)
{
	File5_Group *toc5 = file5->openGroup(TOC_GROUP_NAME, FILE5_OPEN|FILE5_CREATE);
	toc5->headersBegin();
	std::string key;
	std::string val;
	bool found = false;
	int idx = 0;
	while (toc5->headersNext(&key, &val) == FILE5_OK)
	{
		if (name == val)
		{
			found = true;
			break;
		}
		++idx;
	}
	if (found == false)
		return FILE5_ERR_NOTFOUND;
	toc5->usermeta_writeidx(idx, "", "");
	toc5->close();
	delete toc5;
	return file5->deleteTsv(name);
}

bool File5_MultiDataUtils::headersBegin()
{
	bool status = true;
	if (tsv5 != NULL)
		tsv5->headersBegin();
	else if (file5 != NULL)
		file5->headersBegin();
	else
		status = false;
	return status;
}

bool File5_MultiDataUtils::headersNext(std::string* key, std::string* val)
{
	bool status = false;
	if (tsv5 != NULL)
		status = (tsv5->headersNext(key, val) == FILE5_OK);
	else if (file5 != NULL)
		status = (file5->headersNext(key, val) == FILE5_OK);
	return status;
}

template<class T>
static void ExtractColumn(File5_Tsv *tsv5, const std::string& cname, std::list<T> &col)
{
	col.clear();
	int idx = tsv5->getColumnIdx(0, cname);
	T val;
	tsv5->rewind();
	while (tsv5->nextLine() == FILE5_OK)
	{
		tsv5->get(0, idx, &val);
		col.push_back(val);
	}
}

void File5_MultiDataUtils::ExtractColumn(const std::string& cname, std::list<int> &col)
{
	::ExtractColumn(tsv5, cname, col);
}

void File5_MultiDataUtils::ExtractColumn(const std::string& cname, std::list<float> &col)
{
	::ExtractColumn(tsv5, cname, col);
}

void File5_MultiDataUtils::ExtractColumn(const std::string& cname, std::list<std::string> &col)
{
	::ExtractColumn(tsv5, cname, col);
}

void File5_MultiDataUtils::ExtractColumnAsString(const std::string& cname, std::list<std::string> &col)
{
	col.clear();
	int idx = tsv5->getColumnIdx(0, cname);
	std::string val;
	tsv5->rewind();
	while (tsv5->nextLine() == FILE5_OK)
	{
		tsv5->getAsStr(0, idx, &val);
		col.push_back(val);
	}
  /*
  col.clear();
  affx::File5_TsvColumn* colptr=tsv5->getColumnPtr(0,cname);
  if (colptr==NULL) {
    return;
  }
  colptr->getColumnAsStringVector(&col);
  */
}

void File5_MultiDataUtils::ExtractColumnNames(std::list<ColumnNameType>& colInfo)
{
	colInfo.clear();
	int ncols = tsv5->getColumnCount(0);
	ColumnNameType col;	
	for (int i=0; i<ncols; i++)
	{
		tsv5->getColumnName(0, i, &col.cname);
		col.ctype = tsv5->getColumnDtype(0, i);
		colInfo.push_back(col);
	}
}
