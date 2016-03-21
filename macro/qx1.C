
// Daniel Pitzl, Jan 2016
// transform 1D Profile exp(-q) to q

// .x qx1.C("cmsqxvsxm")

#include <iostream> // cout

#include "TProfile.h"
#include "TH1D.h"
#include "TMath.h"

//------------------------------------------------------------------------------
void qx1( string hs, double wid = 3.5 )
{
  TProfile * p = (TProfile*)gDirectory->Get(hs.c_str());

  if( p == NULL ) {
    cout << hs << " does not exist\n";
    return;
  }

  cout << p->GetTitle() << endl;

  TH1D * hq = new
    TH1D( "qx1", p->GetTitle(),
	  p->GetNbinsX(),
	  p->GetXaxis()->GetXmin(),
	  p->GetXaxis()->GetXmax()
	  );
  hq->GetXaxis()->SetTitle( p->GetXaxis()->GetTitle( ) );
  hq->GetYaxis()->SetTitle( p->GetYaxis()->GetTitle( ) );

  for( int ii = 1; ii <= p->GetNbinsX(); ++ii ) {

    double qx = p->GetBinContent(ii);
    double q = -wid*log(qx);
    hq->SetBinContent( ii, q );
    //cout << ii << "  " <<  p->GetBinCenter(ii) << "  " <<  qx << "  " <<  q << endl;

  } // ii

  hq->Draw();

}
