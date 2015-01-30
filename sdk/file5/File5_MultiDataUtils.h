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

#ifndef _FILE5_MULTIDATAUTILS_H_
#define _FILE5_MULTIDATAUTILS_H_

//
#include "file5/File5_types.h"
#include "file5/File5_File.h"
#include "file5/File5_Tsv.h"
//
#include <string>
#include <list>

namespace affx
{

/*! The name and column type. */
typedef struct _ColumnNameType
{
	std::string cname;		/// The column name
	File5_dtype_t ctype;	/// The column type
} ColumnNameType;

/*! Stores a display and associated group name. */
typedef struct _TableOfContentsEntry
{
	std::string displayName;	/// A name to use for user display.
	std::string groupName;		/// The group name
} TableOfContentsEntry;

/*! A utility class for handling a file with multiple data tables. */
// Normally we use "m_" as a prefix for member names.
class File5_MultiDataUtils
{
private:
	/*! The pointer to the file. */
	File5_File *file5;

	/*! A pointer to the data group. */
	File5_Tsv *tsv5;

	/*! The column index. */
	int colIdx;

public:
	/*! The constructor. */
	File5_MultiDataUtils();
	
	/*! The destructor. */
	virtual ~File5_MultiDataUtils();

	/*! Open an existing, or create a new, file.
	 * @param fileName The name of the file.
	 */
	void Open(const std::string &fileName);

	/*! Close the file. */
	void Close();

/////////////////////////////////////////WRITE FUNCTIONS///////////////////////////////////////////////////////

	/*! Create a new tsv group.
	 * @param gname The name of the group.
	 * @param dname The display name for the group.
	 */
	void CreateTsv(const std::string &gname, const std::string &dname);

	/*! Define a new column
	 * @param cname The name of the column.
	 * @param ctype The data type for the column.
	 * @param csize The size of the column.
	 */
	void DefineColumn(const std::string& cname, const File5_dtype_t ctype, int csize=-1);

	/*! Add an entry to the column.
	 * @param ival The value to add.
	 */
	void AddEntry(int ival);

	/*! Add an entry to the column.
	 * @param fval The value to add.
	 */
	void AddEntry(float fval);

	/*! Add an entry to the column.
	 * @param sval The value to add.
	 */
	void AddEntry(const std::string &sval);

	/*! Add a parameter to the open group, or file header before DefineColumn is called.
	 * @param key The parameter key name.
	 * @param val The parameter value.
	 */
	void AddParameter(const std::string &key, const std::string &val);

/////////////////////////////////////////READ FUNCTIONS///////////////////////////////////////////////////////

	/*! Open an existing tsv group.
	 * @param gname The name of the group.
	 */
	void OpenTsv(const std::string &gname);
	
	/*! Close the tsv group. */
	void CloseTsv();

	/*! Delete a tsv group. */
	File5_return_t DeleteTsv(const std::string &name);

	/*! Get the table of contents.
	 * @param toc The table of contents.
	 */
	void GetTOC(std::vector<TableOfContentsEntry> &toc);

	/*! Initialize retrieval of the header.
	 * @param groupType The type of group.
	 * @return The status
	 */
	bool headersBegin();

	/*! Get the next header.
	 * @param key The output header key.
	 * @param val The output header value.
	 * @return The status
	 */
	bool headersNext(std::string* key, std::string* val);

	/*! Extract a column
	 * @param cname The name of the column.
	 * @param col The column data.
	 */
	void ExtractColumn(const std::string& cname, std::list<int> &col);

	/*! Extract a column
	 * @param cname The name of the column.
	 * @param col The column data.
	 */
	void ExtractColumn(const std::string& cname, std::list<float> &col);

	/*! Extract a column
	 * @param cname The name of the column.
	 * @param col The column data.
	 */
	void ExtractColumn(const std::string& cname, std::list<std::string> &col);

	/*! Extract a column as string values
	 * @param cname The name of the column.
	 * @param col The column data.
	 */
	void ExtractColumnAsString(const std::string& cname, std::list<std::string> &col);

	/*! Extract the column names.
	 * @param colInfo The names and data types of the columns.
	 */
	void ExtractColumnNames(std::list<ColumnNameType>& colInfo);

	/*! The line count of the open TSV group. */
	int TsvLineCount();

	/*! Get the next lines data. Assumption is the first column is a name (string) value.
	 * @param name The row name (first column)
	 * @param data The remaining columns data.
	 * @return True if success.
	 */
	bool NextLine(std::string &name, std::vector<int> &data);

	/*! Get the next lines data. Assumption is the first column is a name (string) value.
	 * @param name The row name (first column)
	 * @param data The remaining columns data.
	 * @return True if success.
	 */
	bool NextLine(std::string &name, std::vector<float> &data);

	/*! Get the next lines data. Assumption is the first column is a name (string) value.
	 * @param name The row name (first column)
	 * @param data The remaining columns data.
	 * @return True if success.
	 */
	bool NextLine(std::string &name, std::vector<std::string> &data);

	/*! Get the next lines data as string values. Assumption is the first column is a name (string) value.
	 * @param name The row name (first column)
	 * @param data The remaining columns data.
	 * @return True if success.
	 */
	bool NextLineAsString(std::string &name, std::vector<std::string> &data);
};

}

#endif // _FILE5_MULTIDATAUTILS_H_
