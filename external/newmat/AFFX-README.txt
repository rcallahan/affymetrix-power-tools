Newmat is a publically available matrix library provided by
Robert Davies.  We provide a copy of the newmat code
(version 10d) here for use by various components of the
Affymetrix C++ SDK.

Urls:
 http://www.robertnz.net/nm_intro.htm
 http://www.robertnz.net/index.html
 
Newmat is distributed without a license, but there 
is a "Conditions of Use" clause in the documentation.
Specifically:

   1.1 Conditions of use

   There are no restrictions on the use of newmat 
   except that I take no liability for any problems 
   that may arise from this use.

   I welcome its distribution as part of low cost 
   CD-ROM collections.

   You can use it in your commercial projects. 
   However, if you distribute the source, please make 
   it clear which parts are mine and that they are 
   available essentially for free over the Internet.

   Please understand that there may still be bugs and 
   errors. Use at your own risk. I take no 
   responsibility for any errors or omissions in this 
   package or for any misfortune that may befall you 
   or others as a result of its use.

   Please report bugs to me at robert (at) 
   statsresearch.co.nz

   When reporting a bug please tell me which C++ 
   compiler you are using, and what version. Also give 
   me details of your computer. And tell me which version 
   of newmat (e.g. newmat03 or newmat04) you are using. 
   Note any changes you have made to my code. If at 
   all possible give me a piece of code illustrating 
   the bug. See the problem report form.

   Please do report bugs to me.

Also note that the readme states:

   This library is freeware and may be freely used and distributed.

The Newmat code provided here may contain modifications
relative to the Newmat code provided by Robert Davies. In
generally these modifications will be kept to a minimum. 
We will try and keep track of the modifications here:

 - Minor modifications, documented on the newmat
   bug list and notes web page, for compatibility
   with Microsoft Visual C++, including the removal
   of cout use, which causes problems when
   linking into a windows form/gui application.

----------

To update the affy cvs version of newmat:

* Fetch the latest newmat source and unpack it.

wget http://www.robertnz.net/ftp/newmat10.tar.gz
mkdir newmat
cd newmat10
tar -zxvf ../newmat10.tar.gz

* import it into cvs

cvs import -m "* Import of newmat-10d." -I '!' robertnz/newmat RobertDavies NEWMAT_$NEWVERSION

* merge the updates into the current head
  
cvs checkout -jNEWMAT_$LASTVERSION -jNEWMAT_$NEWVERSION robertnz/newmat

* resolve the merges

emacs

* commit the results

cvs commit -m "* Merge with NEWMAT_$LASTVERSION"
