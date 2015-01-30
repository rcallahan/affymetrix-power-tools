////////////////////////////////////////////////////////////////
//
// Copyright (C) 2004 Affymetrix, Inc.
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

/*
>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

The PLIER (Probe Logarithmic Error Intensity Estimate) method produces
an improved signal by accounting for experimentally observed patterns 
in probe behavior and handling error at the appropriately at low and 
high signal values.

>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
*/


/*
 * affy_ptr Template Class
 * used as a smart pointer
 */
#ifndef __AFFY_PTR_H_
#define __AFFY_PTR_H_

template <class T> class affy_ptr {
public:
	affy_ptr() : p(0) {}
	affy_ptr(T* p_) : p(p_) { if(p) p->addref(); }
	~affy_ptr() { if(p) p->release(); }
	operator T*() { return p; }
	T** operator&() { return &p; }
	T& operator*() { return *p; }
	T* operator->() { return p; } 
	bool operator==(T* p_) { return p == p_; }
	bool operator==(affy_ptr<T>& p_) {	return operator==(p_.p); }
	bool operator!=(T* p_) { return !(operator==(p_)); }
	bool operator!=(affy_ptr<T>& p_) {	return operator!=(p_.p); }
	affy_ptr& operator=(affy_ptr<T> &p_) { return operator=((T*) p_);}
	affy_ptr& operator=(T* p_)
	{
		if(p == p_)	return *this;
		if(p) p->release();
		p = p_;
		if(p)	p->addref();
		return *this;
	}
	T* get() const { return p; }
private:
	T* p;
};

#endif // __AFFY_PTR_H_
