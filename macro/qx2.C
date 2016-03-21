
// Daniel Pitzl, Jan 2016
// transform 2D Profile exp(-q) to q

// .x qx2.C("cmsqxvsxm")

#include <iostream> // cout

#include "TProfile2D.h"
#include "TH2D.h"
#include "TMath.h"

//------------------------------------------------------------------------------
void qx2( string hs, double wid = 3.5 )
{
  TProfile2D * p = (TProfile2D*)gDirectory->Get(hs.c_str());

  if( p == NULL ) {
    cout << hs << " does not exist\n";
    return;
  }

  cout << p->GetTitle() << endl;

  TH2D * hq = new
    TH2D( "qx2", p->GetTitle(),
	  p->GetNbinsX(),
	  p->GetXaxis()->GetXmin(),
	  p->GetXaxis()->GetXmax(),
	  p->GetNbinsY(),
	  p->GetYaxis()->GetXmin(),
	  p->GetYaxis()->GetXmax()
	  );
  hq->GetXaxis()->SetTitle( p->GetXaxis()->GetTitle( ) );
  hq->GetYaxis()->SetTitle( p->GetYaxis()->GetTitle( ) );
  hq->GetZaxis()->SetTitle( p->GetZaxis()->GetTitle( ) );

  for( int ii = 1; ii <= p->GetNbinsX(); ++ii )

    for( int jj = 1; jj <= p->GetNbinsY(); ++jj ) {

      double qx = p->GetBinContent(ii,jj);
      double q = -wid*log(qx);
      hq->SetBinContent( ii, jj, q );

    } // jj

  hq->Draw("colz");

}
