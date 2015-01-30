////////////////////////////////////////////////////////////////
//
// Copyright (C) 2013 Affymetrix, Inc.
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
#ifdef WIN32
#include <winsock2.h>
#endif

#include "util/SocketEngine.h"
#include "util/Verbose.h"
int main(int argc, const char* argv[])
{
  try
  {
    if ( argc > 3 ) {
      std::vector<std::string> strArgs;
      for ( int i = 2 ; i < argc; i++ ) {
        strArgs.push_back(argv[i]);
      }
      SocketEngine engine(argv[1], strArgs);
      engine.run();
    }
    else {
      SocketEngine engine(argv[1],argv[2]);
      engine.run();
    }
  }
  catch(...)
  {     
    Verbose::out(1,"Unexpected Error: Uncaught exception in SocketEngine."); 
  }

  return 0;
}
