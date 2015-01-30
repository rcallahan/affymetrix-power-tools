////////////////////////////////////////////////////////////////
//
// Copyright (C) 2005 Affymetrix, Inc.
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

#include "util/Err.h"
#include "util/Except.h"
#include "util/Util.h"
//

class MyErrHandler : public ErrHandler {

public:
  // Our broken error handler.
  void handleError(const std::string & msg) {
    Verbose::out(0, "Need to fix this error handler to abort. Message is: '" + msg + "'");
  }
  
};

int main(int argc, const char *argv[]) {
  Verbose::out(1, "Welcome to error handling demostration.");
  // use a throwing error handler.
  Err::setThrowStatus(true); 
  // push on a custom error handler
  MyErrHandler *handler = new MyErrHandler();
  Err::pushHandler(handler);
  Err::errAbort("Simulating first dummy error with custom error handler.");

  // remove custom error handler leaving throwing error handler.
  ErrHandler *errHand = Err::popHandler();
  delete errHand;
  try {
    Err::errAbort("Simulating first dummy error with default throwing error handler.");
  }
  catch(Except &e) {
    Verbose::out(1, "Caught exception: '" + ToStr(e.what()) + "'");
  }
  // add a handler that just exit(1) at end of message.
  Err::setThrowStatus(false);
  Err::errAbort("Simulating first dummy error with default exit(1) error handler.");
  return 0;
}
