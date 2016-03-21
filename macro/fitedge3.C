
// Daniel Pitzl, Apr 2014
// fit left edge of pixel charge distribution by tanh
// .x fitedge3.C("cmspxq2")

#include "TH1.h"
#include "TMath.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"

Double_t fitTanh3( Double_t *x, Double_t *par )
{
  return 0.5*par[2] * ( 1 + TMath::TanH( ( x[0] - par[0] ) / par[1] ) );
}

//----------------------------------------------------------------------
void fitedge3( string hs )
{
  TH1 *h = (TH1*)gDirectory->Get(hs.c_str());

  if( h == NULL ) {
    cout << hs << " does not exist\n";
    return;
  }

  h->SetMarkerStyle(21);
  h->SetMarkerSize(0.8);
  h->SetStats(1);
  gStyle->SetOptFit(101);
  gStyle->SetOptStat(11);

  gROOT->ForceStyle();

  // set bin error:

  int nb = h->GetNbinsX();
  //for( int ii = 1; ii <= nb; ++ii ) h->SetBinError( ii, 1.1 );

  // find a few bins:

  int imax = h->GetMaximumBin();
  double amax = h->GetBinContent(imax);

  double n50 = 0.5*amax;

  double x50 = 0;

  for( int ii = 1; ii < nb; ++ii ) {

    int jj = ii + 1;

    if( h->GetBinContent(ii) <= n50 &&
	h->GetBinContent(jj) >= n50 ) {

      double n1 = h->GetBinContent(ii);
      double n2 = h->GetBinContent(jj);
      double x1 = h->GetBinCenter(ii);
      double x2 = h->GetBinCenter(jj);
      double dx = x2 - x1;
      double dn = n2 - n1;
      x50 = x1 + (n50-n1)/dn * dx;
      break;
    }

  } // ii

  cout << "x50 " << x50 << endl;

  int ib0 = 1; // bin couting starts at 1 (0 is underflows)
  for( ; ib0 < nb; ++ib0 )
    if( h->GetBinContent(ib0) > 0 ) break; // first bin above 0

  double x0 = h->GetBinLowEdge(ib0);
  double x9 = h->GetBinLowEdge(imax) + h->GetBinWidth(imax);

  // create a TF1 with the range from x0 to x9 and 3 parameters

  TF1 *tanhFcn3 = new TF1( "tanhFcn3", fitTanh3, x0, x9, 3 );

  tanhFcn3->SetParName( 0, "edge" );
  tanhFcn3->SetParName( 1, "width" );
  tanhFcn3->SetParName( 2, "scale" );

  // set start values:

  tanhFcn3->SetParameter( 0, x50 ); // edge [ke]
  tanhFcn3->SetParameter( 1, 0.8 ); // width [ke]
  tanhFcn3->SetParameter( 2, amax ); // height

  tanhFcn3->SetNpx(500);
  tanhFcn3->SetLineWidth(4);
  tanhFcn3->SetLineColor(kMagenta);

  h->Fit( "tanhFcn3", "R", "p" ); // R = range from tanhFcn3

  h->Draw( "histpsame" );  // data again on top

  cout << "Ndata = " << tanhFcn3->GetNumberFitPoints() << endl;
  cout << "Npar  = " << tanhFcn3->GetNumberFreeParameters() << endl;
  cout << "NDoF  = " << tanhFcn3->GetNDF() << endl;
  cout << "chisq = " << tanhFcn3->GetChisquare() << endl;
  cout << "prob  = " << tanhFcn3->GetProb() << endl;

}
