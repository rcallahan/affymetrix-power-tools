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

// HACK: resolve winsock2.h compilation errors by including locally and not from SocketBase.h
#ifdef WIN32
// Windows includes are different for sockets than unix.
#include <winsock2.h>  //QZ: windows.h will load winsock.h automatically. It seems we do not use any winsock2.h specific functions. Ellen has tested out the changes.
#endif

#include "util/SocketServer.h"
#include "util/MsgSocketHandler.h"

using namespace std;

string com = "./apt-probeset-genotype -a artifact-reduction.ResType=2.Clip=0.4.Close=2.Open=2.Fringe=4.CC=2,quant-norm.target=1000.sketch=5000,pm-only,brlmm-p.CM=1.bins=100.mix=1.bic=2.lambda=1.0.HARD=3.SB=0.75.transform=MVA.copyqc=0.00000.wobble=0.05.MS=0.15.copytype=-1.clustertype=2.ocean=0.00001.CSepPen=0.1.CSepThr=4 --qmethod-spec med-polish.expon=true --read-models-brlmmp../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.models --spf-file../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.spf --special-snps../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.specialSNPs --chrX-probes../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.chrXprobes --chrY-probes../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.chrYprobes --target-sketch../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.sketch --use-feat-eff../../regression-data/data/idata/lib/Axiom_GW_Hu_SNP/Axiom_GW_Hu_SNP.r2.AxiomGT1.feature-effects --set-gender-method cn-probe-chrXY-ratio --em-gender false --female-thresh 0.54 --male-thresh 1.0 --cc-chp-output --chip-type Axiom_GW_Hu_SNP --chip-type Axiom_GW_Hu_SNP.r2 --set-analysis-name AxiomGT1 --out-dir test-generated/doMultiChannelAxiom --force --use-socket localhost:9034 ../../regression-data/data/idata/cel/Axiom_GW_Hu_SNP/*.CEL";

void help() {
  cout << "socket-listener-demo - Basic program to listen for messages from engine." << endl;
  cout << endl;
  cout << "usage:" << endl;
  cout << "   socket-listener-demo port" << endl;
  exit(1);
}

int main(int argc, char *argv[]) {
 
 SocketServer server;
 string host = "127.0.0.1";
 if (argc < 2) {
   help();
 }
 string port = argv[1];
 cout << "Trying to open socket on host: " << host <<  " on port: " << port << endl;
 server.socketOpen(host, port);
// sleep(2);
 std::vector<std::string> messages;
 cout << "Entering loop" << endl;
 while(1) {
   if (!server.isConnected()) {
	   cout << ".";
	   cout.flush();
	 }
   int numFound = server.acceptNewConnection();
   if (numFound > 0) {
     cout << "Found " << numFound << " new connections." << endl;
   }
   if (server.isConnected()) {
	   server.checkForMsgs(messages);
	   for(int i = 0; i < messages.size(); i++) {
          cout << messages[i];
          cout.flush();
       }
   }
#ifdef WIN32
   Sleep(100);
#else
   sleep(1);
#endif
 }
 return 0;
}
