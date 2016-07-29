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


using namespace std;

int getRunParameters(int runnr, double &momentum, double &Bfield, double &alpha){

  ifstream parameterFile( "BScanParameters.dat" );

  cout << endl;
  if( parameterFile.bad() || !parameterFile.is_open() ) {
    cout << "File could not be found." << endl;
    return -1;
  }
  else {

    cout << "Reading parameters." << endl;

    int currentRunnr;
    double currentMomentum;
    double currentAlpha;
    double currentBcurrent;
    
    while( ! parameterFile.eof() ) {

      string line;
      getline( parameterFile, line );
      
      if( line.empty() ) continue;
      if( line.at(0) == '#' ) continue;

      stringstream thisline(line);
      
      thisline >> currentRunnr >> currentMomentum >> currentBcurrent >> currentAlpha;

      if(!(currentRunnr && currentMomentum)){
	continue; // No correct data in runlist
      }
      
      if(currentRunnr == runnr){
	momentum = currentMomentum;
	Bfield = currentBcurrent*0.00095; // convert current to field.
	alpha = currentAlpha;
	parameterFile.close();
	return 1;
      }
      
    }
    cout << "Run " << runnr << " not found in list. Please add it." << endl;
    parameterFile.close();
    return -2;
  }
}

void ncolsB() {
  std::cout << "Run ncols(histogram dir, startrun, stoprun)" << std::endl;
}

void ncolsB(const char* inputdir, int startrun, int stoprun) {

  //  TCanvas *c1 = new TCanvas("c1","ncols",600,600);
  //  TProfile *ncols = new TProfile("ncols"," ",130,0,85,0,60,"");
  TProfile2D *ncolsBscan = new TProfile2D("BAlphaScan","Scan of B and alpha", 30, -9.5, 20.5, 15, -0.05, 1.45);

  TCanvas *c2 = new TCanvas("c2","ncolstan",600,600);
  //  TProfile *ncolstan = new TProfile("ncolstan"," ",130,0,12,0,60,"");
  
  TGraphErrors *gNcol[4];
  for(int b = 0; b < 4; b++){
    gNcol[b] = new TGraphErrors();
  }

  gStyle->SetOptStat(0);

  int nruns, nclusters;
  double Bfield;
  double alpha;
  double momentum;

  double minNcol = 10.;
  double maxNcol = 0.;
  double minAlpha = 90.;
  double maxAlpha = -90.;

  for(int run = startrun; run <= stoprun; run++) {

    if(getRunParameters(run, momentum, Bfield, alpha) < 0) continue;

    TString fileName;
    TString fileDir;
    fileDir += inputdir;
    if( !fileDir.EndsWith("/") ) fileDir += "/";

    fileName.Form("quad-%d.root",run);
    fileDir.Append(fileName);


    TFile *source;
    if (!gSystem->AccessPathName( fileDir )) source = TFile::Open(fileDir);
    else continue;

    //cout << "File: " << fileName << endl;

    TH1 *h;
    gDirectory->GetObject("ncolqf4C",h);
    if(!h) continue;
    if(h->GetEntries() < 100){
      cout << "Only " << h->GetEntries() << " entries. Skip this." << endl;
      continue;
    }

    Double_t ncol = h->GetMean();
    Double_t ncolRMS = h->GetRMS();

    // Collect statistics:
    nclusters += h->GetEntries();

    ncolsBscan->Fill(alpha, Bfield, ncol, 1);

    /*
    int bplot;
    if(Bfield > 1.2) bplot=3;
    else if(Bfield <= 1.2 && Bfield > 0.8 ) bplot=2;
    else if(Bfield <= 0.8 && Bfield > 0.3 ) bplot=1;
    else if(Bfield <= 0.3) bplot=0;
    */

    double alphacorr = 2.0;
    int bplot;
    if(Bfield > 1.2){
      bplot=3;
      alpha += alphacorr * 1.;
    }else if(Bfield <= 1.2 && Bfield > 0.8 ){
      bplot=2;
      alpha += alphacorr * 0.75;
    }else if(Bfield <= 0.8 && Bfield > 0.3 ){
      bplot=1;
      alpha += alphacorr * 0.5;
    }else if(Bfield <= 0.3){
      bplot=0;
      alpha += alphacorr * 0.;
    }

    gNcol[bplot]->SetPoint(gNcol[bplot]->GetN(),alpha,ncol);
    gNcol[bplot]->SetPointError(gNcol[bplot]->GetN(),0.1,ncolRMS);

    cout << "run" << run << "\t ncol " << ncol << "\t alpha " << alpha << endl;

    if(ncol > maxNcol) maxNcol = ncol;
    if(ncol < minNcol) minNcol = ncol;
    if(alpha > maxAlpha) maxAlpha = alpha;
    if(alpha < minAlpha) minAlpha = alpha;

    nruns++;
    delete source;
  }

  cout << nruns << " runs analyzed with " << nclusters << " clusters in total." << endl;

  /*
  c1->cd();
  ncolsBscan->SetTitle("2D B alpha Scan;turn angle [#deg];B [T];columns per cluster");
  ncolsBscan->SetMarkerStyle(20);
  ncolsBscan->SetMarkerColor(1);
  ncolsBscan->Draw("colz");
  ncolsBscan->GetZaxis()->SetRangeUser(minNcol-0.02, maxNcol+0.02);
  */
  c2->cd();

  gNcol[0]->SetMarkerColor(kBlack);
  gNcol[1]->SetMarkerColor(kBlue);
  gNcol[2]->SetMarkerColor(kGreen);
  gNcol[3]->SetMarkerColor(kRed);

  for(int b = 0; b < 4; b++){
    gNcol[b]->SetTitle("B alpha scan;turn angle [#deg];columns per cluster");
    gNcol[b]->SetMarkerStyle(20);
    if(b==0){
      gNcol[b]->Draw("AP");
      gNcol[b]->GetYaxis()->SetRangeUser(minNcol-0.02, maxNcol+0.02);
      gNcol[b]->GetXaxis()->SetLimits(minAlpha-1, maxAlpha+1);
    }
    else gNcol[b]->Draw("sameP");
  }

  TLegend *leg = new TLegend(0.4,0.70,0.6,0.88);
  leg->AddEntry(gNcol[0],"   0 A","p");
  leg->AddEntry(gNcol[1]," 700 A","p");
  leg->AddEntry(gNcol[2],"1050 A","p");
  leg->AddEntry(gNcol[3],"1400 A","p");
  leg->Draw();

}
