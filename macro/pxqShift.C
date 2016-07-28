#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TGraph.h"
#include "TSystem.h"
#include "TStyle.h"
#include "TROOT.h"
#include "TGraphErrors.h"

#include <iomanip>
#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <algorithm> 
#include <map>

#include "fitedge3.C"

using namespace std;

struct shiftset{
  int shift[4];
};


int getShiftLines(map<int,shiftset> &shiftLines){

  shiftLines.clear();

  ifstream shiftFile( "shiftParameters.dat" );

  if(shiftFile.bad() || !shiftFile.is_open() ){
    cout << "Shift file could not be found!" << endl;
    return -1;
  }
  else{

    int runnr;
    shiftset shifts;
    for(int mod = 0; mod < 4; mod++){
      shifts.shift[mod] = 0;
    }
    vector<int> runnumbers;

    while( !shiftFile.eof() ) {
      
      string line;
      getline(shiftFile, line);

      if(line.empty()) continue;
      if(line.at(0) == '#') continue;
      
      stringstream thisline(line);
      thisline >> runnr;
      for(int mod = 0; mod < 4; mod++){
	thisline >> shifts.shift[mod];
      }

      vector<int>::iterator it2;
      it2 = find (runnumbers.begin(), runnumbers.end(), runnr);
      
      if (it2 != runnumbers.end()){
	cout << "Run " << runnr << " occurs at least twice. Please correct that! Using the first occurence." << endl;
	continue;
      }
      runnumbers.push_back(runnr);

      std::map<int,shiftset>::iterator it3;
      it3 = shiftLines.find(runnr);

      if(it3 == shiftLines.end()){
	cout << "Did not find them. Of course" << endl;
	shiftLines.insert(pair<int,shiftset>(runnr,shifts));
      }

    }

    shiftFile.close();

    cout << "Successfully read shifts" << endl;
    return 1;

  }

}


void pxqShift() {
  std::cout << "Run ncols(histogram directory, startrun, stoprun)" << std::endl;
}

void pxqShift(const char* inputdir, int startrun, int stoprun) {

  TCanvas *c2 = new TCanvas("c2","pxq",600,600);
  
  TGraph *gShift[4];
  for(int b = 0; b < 4; b++){
    gShift[b] = new TGraph();
  }

  gStyle->SetOptStat(0);

  int nruns = 0;
  int nclusters = 0;

  int minShift = 255;
  int maxShift = -255;

  int shift = 0;

  // Open shift file and store every line in a vector

  map<int,shiftset> shiftLines;

  if( getShiftLines(shiftLines) < -1){
    exit(1);
  }
  
  for(int run = startrun; run <= stoprun; run++) {


    shiftset shifts;
    for(int mod = 0; mod < 4; mod++){
      shifts.shift[mod] = 0;
    }

    std::map<int,shiftset>::iterator it;
    it = shiftLines.find(run);

    if(it == shiftLines.end()){

      for(int mod = 0; mod < 4; mod++){
	shifts.shift[mod] = 0;
      }

      cout << "Did not find shifts. Create them as 0." << endl;

    }else{

      shifts = it->second;
      cout << "Found shifts for run";
      cout << endl;
    }

    
    // Open histogram

    TString fileName;
    TString fileFull;
    fileFull += inputdir;
    if( !fileFull.EndsWith("/") ) fileFull += "/";

    fileName.Form("quad-%d.root",run);
    fileFull.Append(fileName);

    TFile *source;
    if (!gSystem->AccessPathName( fileFull )) source = TFile::Open(fileFull);
    else continue;

    TH1 *h;

    for(int mod=0; mod < 4; mod++){

      TString plotname = "pxq";
      switch(mod){
      case 0:
	plotname += "A";
	break;
      case 1:
	plotname += "B";
	break;
      case 2:
	plotname += "C";
	break;
      case 3:
	plotname += "D";
	break;
      }
      
      gDirectory->GetObject(plotname,h);
      if(!h){
	cout << plotname << " not found." << endl;
	continue;
      }
      // 5 counts per ke
      double thresh = fitedge3(h);
      shift = (1.8-thresh)*5;

      if(shift > 9) shift = 9;
      if(shift < -15) shift = 9; // Fit is probably garbage because of too low threshold

      // Collect statistics:
      nclusters += h->GetEntries();

      cout << "run" << run << "\t mod " << mod << "\t thresh " << thresh << "\t shift " << shift << endl;

      // Update shift value
      shifts.shift[mod] += shift;
      
      gShift[mod]->SetPoint(gShift[mod]->GetN(),run,shifts.shift[mod]);

      if(shifts.shift[mod] > maxShift) maxShift = shifts.shift[mod];
      if(shifts.shift[mod] < minShift) minShift = shifts.shift[mod];

    }

    // Replace shift values in map

    shiftLines[run] = shifts;
    
    nruns++;
    delete source;

  }

  cout << nruns << " runs analyzed with " << nclusters << " clusters in total." << endl;

  // plotting
  gShift[0]->SetMarkerColor(kBlack);
  gShift[1]->SetMarkerColor(kBlue);
  gShift[2]->SetMarkerColor(kGreen);
  gShift[3]->SetMarkerColor(kRed);

  for(int b = 0; b < 4; b++){
    gShift[b]->SetTitle("Pixel charge offset;runnumber;shift charge");
    gShift[b]->SetMarkerStyle(20);
    if(b==0){
      gShift[b]->Draw("AP");
      gShift[b]->GetYaxis()->SetRangeUser(minShift-0.02, maxShift+0.02);
    }
    else gShift[b]->Draw("sameP");
  }


  // write out the new shifts

  ofstream shiftFile("shiftParameters.dat");

  stringstream shiftStream;

  map<int,shiftset>::iterator it;
  int run;
  for ( it = shiftLines.begin(); it != shiftLines.end(); ++it ){

    run = it->first;
    shiftStream << run << "\t";

    shiftset shifty = it->second;
    for(int mod = 0; mod < 4; mod++){
      shiftStream << shifty.shift[mod] << "\t";
    }

    shiftStream << endl;
    shiftFile << shiftStream.str();
    shiftStream.str("");
  }

  shiftFile.close();

  cout << "\nRun quad for these runs in order to apply the shift and get new pxq distributions. You can do a next iteration afterwards." << endl << endl;

}
