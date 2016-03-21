
// Daniel Pitzl, Sep 2011
// fit Landau function x Gaussian to energy loss in silicon
// .x fitlang.C("h032")

//------------------------------------------------------------------------------
Double_t fitLandauGauss( Double_t *x, Double_t *par )
{
  static int nn = 0;
  nn++;
  static double xbin = 1;
  static double b1 = 0;
  if( nn == 1 ) {
    b1 = x[0];
    cout << "b1 = " << b1 << endl;
  }
  if( nn == 2 ) {
    cout << "b2 = " << x[0] << endl;
    xbin = x[0] - b1;// bin width needed for normalization
    cout << "xbin = " << xbin << endl;
  }

  // Landau:

  Double_t invsq2pi = 0.3989422804014;   // (2 pi)^(-1/2)
  Double_t mpshift  = -0.22278298;       // Landau maximum location

  // MP shift correction:

  double mpc = par[0] - mpshift * par[1]; //most probable value (peak pos)

  //Fit parameters:
  //par[0] = Most Probable (MP, location) parameter of Landau density
  //par[1] = Width (scale) parameter of Landau density
  //par[2] = Total area (integral -inf to inf, normalization constant)
  //par[3] = Gaussian smearing

  // Control constants
  Double_t np = 100.0;      // number of convolution steps
  Double_t sc =   5.0;      // convolution extends to +-sc Gaussian sigmas

  // Range of convolution integral
  double xlow = x[0] - sc * par[3];
  double xupp = x[0] + sc * par[3];

  double step = (xupp-xlow) / np;

  // Convolution integral of Landau and Gaussian by sum

  double sum = 0;
  double xx;
  double fland;

  for( int i = 1; i <= np/2; i++ ) {

    xx = xlow + ( i - 0.5 ) * step;
    fland = TMath::Landau( xx, mpc, par[1] ) / par[1];
    sum += fland * TMath::Gaus( x[0], xx, par[3] );

    xx = xupp - ( i - 0.5 ) * step;
    fland = TMath::Landau( xx, mpc, par[1] ) / par[1];
    sum += fland * TMath::Gaus( x[0], xx, par[3] );
  }

  return( par[2] * invsq2pi * xbin * step * sum / par[3] );
}

//----------------------------------------------------------------------
void fitlang( string hs )
{
  TH1 *h = (TH1*)gDirectory->Get(hs.c_str());

  if( h == NULL ) {
    cout << hs << " does not exist\n";
    return;
  }

  h->SetMarkerStyle(21);
  h->SetMarkerSize(0.8);
  h->SetStats(1);
  //gStyle->SetOptFit(111);//111 = chisq err par
  gStyle->SetOptFit(101);//101 = chisq and par

  gStyle->SetOptStat(11);

  gROOT->ForceStyle();

  double aa = h->GetEntries();//normalization

  // find peak:

  int ipk = h->GetMaximumBin();

  int nn = h->GetNbinsX();
  double ymax = 0;

  for( int ii = 1; ii <= nn; ++ii ){
    if( h->GetBinCenter(ii) < 2.2 ) continue; // avoid peak at 1
    if( h->GetBinContent(ii) > ymax ) {
      ymax = h->GetBinContent(ii);
      ipk = ii;
    }
  }

  double xpk = h->GetBinCenter(ipk);
  double sm = xpk / 9; // sigma
  if( sm < 1.5 ) sm = 1.5;
  double ns = sm; // noise

  cout << "peak at " << xpk << endl;
  cout << "sigma   " << sm << endl;

  // fit range:

  double x0, x9;
  x0 = xpk - 2.5*sm;
  x9 = xpk + 3.0*sm;
  //x9 = xpk + 8.0*sm; // 203i
  //x9 = xpk + 9.0*sm;
  //x9 = 99;
  //x9 = 150; // poor fit in tail

  cout << "fit from " << x0 << " to " << x9 << endl;

  // create a TF1 with the range from x0 to x9 and 4 parameters

  TF1 *fitFcn = new TF1( "fitFcn", fitLandauGauss, x0, x9, 4 );

  fitFcn->SetParName( 0, "peak" );
  fitFcn->SetParName( 1, "sigma" );
  fitFcn->SetParName( 2, "area" );
  fitFcn->SetParName( 3, "smear" );

  fitFcn->SetNpx(500);
  fitFcn->SetLineWidth(4);
  //fitFcn->SetLineColor(kMagenta);
  fitFcn->SetLineColor(kGreen);

  // set start values:

  fitFcn->SetParameter( 0, xpk ); // peak position, defined above
  fitFcn->SetParameter( 1, sm ); // width
  fitFcn->SetParameter( 2, aa ); // area
  fitFcn->SetParameter( 3, ns ); // noise

  h->Fit( "fitFcn", "R", "ep" ); // R = range from fitFcn

  h->Draw( "histepsame" ); // data again on top

  cout << "Ndata = " << fitFcn->GetNumberFitPoints() << endl;
  cout << "Npar  = " << fitFcn->GetNumberFreeParameters() << endl;
  cout << "NDoF  = " << fitFcn->GetNDF() << endl;
  cout << "chisq = " << fitFcn->GetChisquare() << endl;
  cout << "prob  = " << fitFcn->GetProb() << endl;

}
