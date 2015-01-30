AFFYMETRIX POWER TOOLS - Version 1.16.1, Binary/Source distributions.
Please review the license before using these tools.


CONTENTS 
========

LINUX/OS X (BINARY DISTRIBUTION):
---------------------------------
- bin - Affymetrix Power Tools (APT) executable programs.

- doc/dox-html - Documentation, help, vignettes and release notes.


WINDOWS:
--------
- C:\Program Files\Affymetrix Power Tools\apt-1.16.0\bin - Affymetrix Power Tools (APT) executable programs for windows

- C:\Program Files\Affymetrix Power Tools\apt-1.16.0\doc - Documentation, help, vignettes and release notes

SOURCE DISTRIBUTION:
--------------------
- sdk - Affymetrix source code

- externals - External source code

- sdk/output/<PLATFORM>/bin - Affymetrix Power Tools (APT) executable programs (after make)

- sdk/output/<PLATFORM>/lib - Common libraries used by APT (after make)

- sdk/apt/doc - Documentation, help, vignettes and release notes.

(where <PLATFORM> denotes the architecture and OS of your system)




INSTALLATION FROM BINARY DISTRIBUTIONS
======================================

LINUX/OS X: 
-----------
If installing from a binary distribution, copy the files to a dir on your $PATH or add the bin directory to your path:

   export PATH=apt-1.16.0/bin:$PATH

WINDOWS 64-bit: 
---------------
unzip APT-1.16.0-i636-intel-win64.zip, click APT-1.16.0-i636-intel-win64.msi file and follow the instructions to install windows 64-bit.


NOTE FOR WINDOWS USERS: These executables should be run from an account with Administrator privileges on your system.


INSTALLATION FROM SOURCE
========================

LINUX/OS X:
-------------
Download the source distribution.

  unzip apt-1.16.0-src.zip
  cd apt-1.16.0/sdk
  ./configure
  make

NOTE: By default the configuration requires CPPUnit to be installed. This can be found at http://sourceforge.net/projects/cppunit/.
To compile without CPPUnit do:

  unzip apt-1.16.0-src.zip
  cd apt-1.16.0/sdk
  ./configure --without-cppunit
  make

If you are planning to run integration or regression tests you will need to install CPPUnit, and
to specify a path to regression/integration test data:

  unzip apt-1.16.0-src.zip
  cd apt-1.16.0/sdk
  ./configure --regression-data=<REGRESSION-DATA-PATH>
  make

Once you have built the binary executables, you may wish to put them in your PATH. Either

  PATH=[PATH-TO-SRC-DIR]/sdk/output/[PLATFORM]/bin:$PATH

or copy the relevant executables to a directory already in your PATH.


WINDOWS:
--------
Note: The following procedure has been tested using Visual Studio 2010. Other versions of Visual Studio have not been tested and 
are not currently recommended..

- unzip apt-1.16.0-src.zip
- Open the solution file <DIR>/apt-1.16.0-src/sdk/apt/apt.sln in Visual Studio 2010
- Select "Release" in the Solution Configurations dropdown and your architecture (Win32 or x64) in the Solution Platforms dropdown.
- Select "Build Solution (F7)"

REGRESSION AND INTEGRATION TESTING
==================================

Binary distributions have been tested using Affymetrix regression test set, and further 
regression tests should not normally be necessary.

When compiling from source for critical applications, we recommend running integration tests on
the complete APT distribution, and, in addition, running full regretion tests on those
applications you intend to use. Depending on your hardware and architecture, running full
regression tests on the complete APT distribution may take over 24 hours is not normally
recommended.

For details of how to download Affymetrix regression data and run regression and/or integration
tests see http://www.affymetrix.com/partners_programs/programs/developer/tools/powertools.affx.

