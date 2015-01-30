////////////////////////////////////////////////////////////////
//
// Copyright (C) 2010 Affymetrix, Inc.
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
#include <cassert>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iostream>
#include <string>
#ifdef WIN32
#include <winsock2.h>
#include <ws2tcpip.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <sys/socket.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <netdb.h>
#include <signal.h>
#endif /* WIN32 */
#include "portability/affy-system-api.h"
#include "util/SocketServer.h"
#include "util/MsgSocketHandler.h"
#include "util/SocketTextHandler.h"
#include "util/SocketEngine.h"
#include "util/Verbose.h"

using namespace std;

void help() {
  cout << "engine-runner - Basic program that launches specified program in a seperate process." << endl;
  cout << endl;
  cout << "usage:" << endl;
  cout << "   engine-runner port [engine args]" << endl;
  exit(1);
}

/* Run an engine. */
int main(int argc, char *argv[]) {
  if (argc < 2) {
    help();
  }
  vector<string> progArgs;
  for (int i = 3; i < argc; i++) {
      progArgs.push_back(argv[i]);
  }
  string s = argv[1];
  unsigned int port = Convert::toInt(s);
  SocketEngine engine(argv[2], progArgs, "127.0.0.1", port);
  engine.run();
  return 0;
}

const char *testCommand[] = {
  "apt-probeset-genotype",
  "--table-output",
  "--chrX-snps", "../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.chrx",
  "--spf-file", "../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf",
  "-o", "doSpf",
  "--cel-files", "../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt",
  "--use-socket", "127.0.0.1:9035",
  "--console-off",
  "true",
  NULL};

const char *testCom=
  "apt-probeset-genotype"
  " --table-output"
  " --chrX-snps" " ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.chrx"
  " --spf-file" " ../../../regression-data/data/idata/lib/Mapping250K_Sty/Mapping250K_Sty.spf"
  " -o" "doSpf"
  " --cel-files" " ../../../regression-data/data/idata/p-geno/Mapping250K_Sty/cel-files.txt"
  " --use-socket" " 127.0.0.1:9035"
  " &> /dev/null";

