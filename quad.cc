
// Daniel Pitzl, DESY, 2016
// read, analyze, and plot eudaq data from 4 CMS pixel modules
// eudecoder

// quad -l 59999 185

// D C B A <-- beam

#include <sstream> // stringstream
#include <fstream> // filestream
#include <iomanip>

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1D.h>
#include <TH2D.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include "GblTrajectory.h"
#include <TMatrixD.h>
#include <TMath.h>
#include "MilleBinary.h"

using namespace std;
using namespace gbl;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int adc;
  double cal;
  bool big;
};

struct cluster {
  vector <pixel> vpix;
  int size;
  int sumA; // DP
  double charge;
  double col,row;
  bool big;
  int ncol, nrow;
};

// globals:

pixel pb[66560]; // global declaration: vector of pixels with hit
int fNHit; // global

//------------------------------------------------------------------------------
// inverse decorrelated Weibull PH -> large Vcal DAC
double PHtoVcal( double ph, double a0, double a1, double a2, double a3, double a4, double a5, int mod )
{
  // modph2ps decorrelated: ph = a4 - a3*exp(-t^a2), t = a0 + q/a1

  //return ph; // test !!

  double Ared = ph - a4; // a4 is asymptotic maximum

  if( Ared >= 0 ) {
    Ared = -0.1; // avoid overflow
  }

  // large Vcal = ( (-ln(-(A-a4)/a3))^1/a2 - a0 )*a1

  if( a3 < 1E-9 ) {
    //cout << "PHtoVcal zero a3  " << a3 << endl;
    return ph;
  }
  else if( -Ared > a3 ) {
    //cout << "PHtoVcal small a3  " << a3 << "  " << -Ared << " mod " << mod << endl;
    return ph;
  }

  double vc =
    a1 * ( pow( -log( -Ared / a3 ), 1/a2 ) - a0 );

  if( vc > 999 )
    cout << "overflow " << vc << ", Ared " << Ared << ", a3 " << a3 << endl;

  if( vc != vc ) {

    cout << "PHtoVcal NaN at "
	 << "  " << a0
	 << "  " << a1
	 << "  " << a2
	 << "  " << a3
	 << "  " << a4
	 << "  " << a5
	 << endl;

    return ph;

  }

  return vc * a5; // small Vcal
  //return vc; // large Vcal
}

//------------------------------------------------------------------------------
vector<cluster> getClus()
{
  // returns clusters with local coordinates
  // decodePixels should have been called before to fill pixel buffer pb 
  // simple clusterization
  // cluster search radius fCluCut ( allows fCluCut-1 empty pixels)

  const int fCluCut = 1; // clustering: 1 = no gap (15.7.2012)
  //const int fCluCut = 2;

  vector<cluster> v;
  if( fNHit == 0 ) return v;

  int* gone = new int[fNHit];

  for( int i = 0; i < fNHit; ++i )
    gone[i] = 0;

  int seed = 0;

  while( seed < fNHit ) {

    // start a new cluster

    cluster c;
    c.vpix.push_back( pb[seed] );
    gone[seed] = 1;

    // let it grow as much as possible:

    int growing;
    do{
      growing = 0;
      for( int i = 0; i < fNHit; ++i ) {
        if( !gone[i] ){ // unused pixel
          for( unsigned int p = 0; p < c.vpix.size(); ++p ) { // vpix in cluster so far
            int dr = c.vpix.at(p).row - pb[i].row;
            int dc = c.vpix.at(p).col - pb[i].col;
            if( (   dr>=-fCluCut) && (dr<=fCluCut) 
		&& (dc>=-fCluCut) && (dc<=fCluCut) ) {
              c.vpix.push_back(pb[i]);
	      gone[i] = 1;
              growing = 1;
              break; // important!
            }
          } // loop over vpix
        } // not gone
      } // loop over all pix
    }
    while( growing );

    // added all I could. determine position and append it to the list of clusters:

    c.sumA = 0;
    c.charge = 0;
    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    c.big = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      c.sumA += p->adc; // Aout
      double Qpix = p->cal; // calibrated [Vcal]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
      if( p->big ) c.big = 1;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( ! c.charge == 0 ) {
      c.col /= sumQ;
      c.row /= sumQ;
    }
    else {
      c.col = (*c.vpix.begin()).col;
      c.row = (*c.vpix.begin()).row;
      cout << "GetClus: cluster with zero charge" << endl;
    }

    c.ncol = maxx-minx+1;
    c.nrow = maxy-miny+1;

    v.push_back(c); // add cluster to vector

    // look for a new seed = used pixel:

    while( (++seed < fNHit) && gone[seed] );

  } // while over seeds

  // nothing left, return clusters

  delete gone;
  return v;
}

//------------------------------------------------------------------------------
TMatrixD Jac5( double ds ) // for GBL
{
  /*
    straight line, no B-field
    track = 
    q/p, x', y', x, y
    0,   1,  2,  3, 4
  */
  TMatrixD jac(5, 5);
  jac.UnitMatrix();
  jac[3][1] = ds; // x = xp * ds
  jac[4][2] = ds; // y = yp * ds
  return jac;
}

//------------------------------------------------------------------------------
bool isFiducial( double x, double y)
{
  bool ffiducial = true;
  if(y < -(8.1-0.3-0.06) || y > (8.1-0.3-0.06)
     || x < -(32.4-0.45-0.06) || x > (32.4-0.45-0.06)) ffiducial = false;
  return ffiducial;
}

//------------------------------------------------------------------------------
int main( int argc, char* argv[] )
{
  cout << "main " << argv[0] << " called with " << argc << " arguments" << endl;

  if( argc == 1 ) {
    cout << "give run number" << endl;
    return 1;
  }

  // run number = last arg

  string runnum( argv[argc-1] );
  int run = atoi( argv[argc-1] );

  cout << "run " << run << endl;
  FileReader * reader;
  if( run < 100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X" );
  else if( run < 1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X" );
  else if( run < 10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X" );
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X" );
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X" );

  // further arguments:

  int lev = 999222111; // last event

  double p = 5.2; // [GeV]

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-p" ) )
      p = atof( argv[++i] ); // momentum

  } // argc

  // alignments:

  double dz = 22.5; // [mm] projected z spacing
  if( run > 150 )
    dz = 40.0;

  int aligniteration = 0;
  double alignx[4];
  double aligny[4];
  double fx[4];
  double fy[4];
  double tx[4];
  double ty[4];

  for( int ipl = 0; ipl < 4; ++ipl ) {
    alignx[ipl] = 0;
    aligny[ipl] = 0;
    fx[ipl] = 0;
    fy[ipl] = 0;
    tx[ipl] = 0;
    ty[ipl] = 0;
  }

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;
  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "no " << alignFileName.str() << ", will bootstrap" << endl;
    cout << endl;
  }
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string Hash( "#" );
    string Iteration( "iteration" );
    string Plane( "plane" );
    string Alignx( "alignx" );
    string Aligny( "aligny" );
    string Rx( "fx" );
    string Ry( "fy" );
    string Tx( "tx" );
    string Ty( "ty" );

    int ipl = 0;

    while( ! ialignFile.eof() ) {

      string line;
      getline( ialignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == Hash ) // comments start with #
	continue;

      if( tag == Iteration ) 
	tokenizer >> aligniteration;

      if( tag == Plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 4 ) {
	//cout << "wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == Alignx )
	alignx[ipl] = val;
      else if( tag == Aligny )
	aligny[ipl] = val;
      else if( tag == Rx )
	fx[ipl] = val;
      else if( tag == Ry )
	fy[ipl] = val;
      else if( tag == Tx )
	tx[ipl] = val;
      else if( tag == Ty )
	ty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  // for GBL:

  double resx = 9.9E-3; // [mm] col hit resolution
  double resy = 7.7E-3; // [mm] row hit resolution

  // X0 Si = 21.82/2.33 = 9.365 cm
  // X0 Al = 24.01/2.70 = 8.89 cm
  // X0 Cu = 12.86/8.96 = 1.435 cm
  // X0 air = 36.62/1.204E-3 = 304 m

  double X0Si = ( 0.3 + 0.175 + 3.0 ) / 94; // Sensor + ROC + Alu

  // measurement = residual
  TVectorD meas(2);
  meas.Zero(); // ideal

  TVectorD measPrec(2); // precision = 1/resolution^2
  measPrec[0] = 1.0 / resx / resx;
  measPrec[1] = 1.0 / resy / resy;

  TMatrixD proL2m( 2, 2 ); // measurement projection matrix
  proL2m.UnitMatrix();

  // scatter:
  TVectorD scat(2);
  scat.Zero(); // mean is zero

  double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*log(X0Si) );

  // pair matching cuts:

  double bicutx = 5E-3*dz; // [mm]
  double bicuty = 5E-3*dz; // [mm]

  // triplet linking cuts:

  double tricutx = 0.06; // [mm]
  double tricuty = 0.06; // [mm]

  if( tricutx < 3*dz*tetSi ) {
    tricutx = 3*dz*tetSi;
    tricuty = 3*dz*tetSi;
  }

  cout << "\nlinking cuts " << tricutx << ", " << tricuty << " mm\n" << endl;

  TVectorD wscatSi(2);
  wscatSi[0] = 1.0 / ( tetSi * tetSi ); // weight
  wscatSi[1] = 1.0 / ( tetSi * tetSi );

  // global labels for Pede:

  vector<int> labelsA( 6 );
  labelsA[0] =  1; // dx
  labelsA[1] =  2; // dy
  labelsA[2] =  3; // fx
  labelsA[3] =  4; // fy
  labelsA[4] =  5; // tx
  labelsA[5] =  6; // ty

  vector<int> labelsB( 6 );
  labelsB[0] =  7; // dx
  labelsB[1] =  8; // dy
  labelsB[2] =  9; // fx
  labelsB[3] = 10; // fy
  labelsB[4] = 11; // tx
  labelsB[5] = 12; // ty

  vector<int> labelsC( 6 );
  labelsC[0] = 13; // dx
  labelsC[1] = 14; // dy
  labelsC[2] = 15; // fx
  labelsC[3] = 16; // fy
  labelsC[4] = 17; // tx
  labelsC[5] = 18; // ty

  vector<int> labelsD( 6 );
  labelsD[0] = 19; // dx
  labelsD[1] = 20; // dy
  labelsD[2] = 21; // fx
  labelsD[3] = 22; // fy
  labelsD[4] = 23; // tx
  labelsD[5] = 24; // ty

  MilleBinary * mille = new MilleBinary( "mille.bin" );

  // Landau peak cuts: Mon 27.7.2015

  double qL = 20; // [ke]
  double qR = 36; // [ke]

  string gainFileName[4];
  double ke[4];

  const int A = 0;
  const int B = 1;
  const int C = 2;
  const int D = 3;

  if( run >= 435 ) { // 2016 May

    gainFileName[A] = "D4122-tb24-gaincal.dat";
    ke[0] = 0.0434; // small Vcal -> ke at 26.4 ke with tilt and turn

    gainFileName[B] = "D4139-tb24-gaincal.dat";
    ke[1] = 0.0437; // small Vcal -> ke

    gainFileName[C] = "D4294-tb24-gaincal.dat";
    ke[2] = 0.0428; // small Vcal -> ke

    gainFileName[D] = "D4329-tb24-gaincal.dat";
    ke[3] = 0.0399; // small Vcal -> ke

  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // gain parameters for mod roc col row:

  double ****p0 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p0[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p0[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p0[i][j][k] = new double[80];

  double ****p1 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p1[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p1[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p1[i][j][k] = new double[80];

  double ****p2 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p2[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p2[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p2[i][j][k] = new double[80];

  double ****p3 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p3[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p3[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p3[i][j][k] = new double[80];

  double ****p4 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p4[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p4[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p4[i][j][k] = new double[80];

  double ****p5 = new double ***[4];
  for( int i = 0; i < 4; ++i )
    p5[i] = new double **[16];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      p5[i][j] = new double *[52];
  for( int i = 0; i < 4; ++i )
    for( int j = 0; j < 16; ++j )
      for( int k = 0; k < 52; ++k )
        p5[i][j][k] = new double[80];

  bool haveGain[4] = {0};

  for( int mod = 0; mod < 4; ++mod ) {

    if( gainFileName[mod].length(  ) > 0 ) {

      ifstream gainFile( gainFileName[mod].c_str() );

      if( gainFile ) {

	haveGain[mod] = 1;
	cout << "gain " << mod << ": " << gainFileName[mod] << endl;

	int roc;
	int col;
	int row;
	double a0, a1, a2, a3, a4, a5;

	while( gainFile >> roc ) {
	  gainFile >> col;
	  gainFile >> row;
	  gainFile >> a0;
	  gainFile >> a1;
	  gainFile >> a2;
	  gainFile >> a3;
	  gainFile >> a4;
	  gainFile >> a5;
	  p0[mod][roc][col][row] = a0;
	  p1[mod][roc][col][row] = a1;
	  p2[mod][roc][col][row] = a2;
	  p3[mod][roc][col][row] = a3;
	  p4[mod][roc][col][row] = a4;
	  p5[mod][roc][col][row] = a5;
	}

      } // gainFile open

    } // gainFileName

  } // mod

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream fname; // output string stream

  fname << "quad-" << run << ".root";

  TFile* histoFile = new TFile( fname.str(  ).c_str(  ), "RECREATE" );

  // book histos:
 
  TH1D hcol[4];
  TH1D hrow[4];
  TH1D hpxq[4];
  TH2D * hmap[4];
  TH1D hnpx[4];
  TH1D hsiz[4];
  TH1D hclq[4];
  TH1D hncol[4];
  TH1D hnrow[4];

  TH1D hncl[4];
  for( int mod = 0; mod < 4; ++mod ) {
    char modtos;
    switch( mod ) {
    case 0: modtos = 'A'; break;
    case 1: modtos = 'B'; break;
    case 2: modtos = 'C'; break;
    case 3: modtos = 'D'; break;
    default:modtos = 'X';
    }

    hncl[mod] = TH1D( Form( "ncl%c", modtos ),
		      Form( "plane %c cluster per event;cluster;plane %c events", modtos, modtos ),
		      51, -0.5, 50.5 );
    hcol[mod] = TH1D( Form("col%c", modtos),
		      Form("%c col;col;%c pixels", modtos, modtos), 
		      416, -0.5, 415.5 );
    hrow[mod] = TH1D( Form("row%c",modtos),
		      Form("%c row;row;%c pixels",modtos,modtos),
		      160, -0.5, 159.5 );
    hpxq[mod] = TH1D( Form("pxq%c",modtos),
		      Form("%c pixel charge;pixel q [ke];%c pixels",modtos,modtos),
		      100, 0, 25 );
    hmap[mod] = new  TH2D( Form("pxmap%c",modtos),
			   Form("%c pixel map;column;row;%c pixels",modtos,modtos),
			   416, -0.5, 415.5, 160, -0.5, 159.5 );
    hnpx[mod] = TH1D( Form("npx%c",modtos),
		      Form("%c pixel per event;pixels;%c events",modtos,modtos),
		      51, -0.5, 50.5 );
    hsiz[mod] = TH1D( Form("clsz%c",modtos),
		      Form("%c cluster size;pixels/cluster;%c clusters",modtos,modtos),
		      51, -0.5, 50.5 );
    hclq[mod] = TH1D( Form("clq%c",modtos),
		      Form("%c cluster charge;cluster charge [ke];%c clusters",modtos,modtos),
		      100, 0, 100 );
    hncol[mod]= TH1D( Form("ncol%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hnrow[mod]= TH1D( Form("nrow%c",modtos),
		      Form("%c cluster size;rows/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );

  } // module planes

  TH2D hxxAB( "xxAB", "A vs B;col B;col A;clusters",
	      432, -32.4, 32.4, 432, -32.4, 32.4 );
  TH2D hyyAB( "yyAB", "A vs B;row B;row A;clusters",
	      162, -8.1, 8.1, 162, -8.1, 8.1 );
  TH1D hdxAB( "dxAB", "Ax-Bx;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyAB( "dyAB", "Ay-By;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcAB( "dxcAB", "Ax-Bx;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycAB( "dycAB", "Ay-By;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dxvsxAB( "dxvsxAB", "A-B dx vs x;x [mm];A-B <dx>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyAB( "dxvsyAB", "A-B dx vs y;y [mm];A-B <dx>",
		    81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxAB( "dyvsxAB", "A-B dy vs x;x [mm];A-B <dy>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyAB( "dyvsyAB", "A-B dy vs y;y [mm];A-B <dy>",
		    81, -8.1, 8.1, -1, 1 );

  TProfile mABvst( "mABvst", "AB matches vs time;trigger;AB matches",
		   200, 0E6, 2E6, -1, 2 );

  TH2D hxxCB( "xxCB", "C vs B;col B;col C;clusters",
	      432, -32.4, 32.4, 432, -32.4, 32.4 );
  TH2D hyyCB( "yyCB", "C vs B;row B;row C;clusters",
	      162, -8.1, 8.1, 162, -8.1, 8.1 );
  TH1D hdxCB( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyCB( "dyCB", "Cy-By;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcCB( "dxcCB", "Cx-Bx;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycCB( "dycCB", "Cy-By;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dxvsxCB( "dxvsxCB", "C-B dx vs x;x [mm];C-B <dx>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyCB( "dxvsyCB", "C-B dx vs y;y [mm];C-B <dx>",
		    81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxCB( "dyvsxCB", "C-B dy vs x;x [mm];C-B <dy>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyCB( "dyvsyCB", "C-B dy vs y;y [mm];C-B <dy>",
		    81, -8.1, 8.1, -1, 1 );

  TProfile mCBvst( "mCBvst", "CB matches vs time;trigger;CB matches",
		   200, 0E6, 2E6, -1, 2 );

  TH2D hxxDB( "xxDB", "D vs B;col B;col D;clusters",
	      432, -32.4, 32.4, 432, -32.4, 32.4 );
  TH2D hyyDB( "yyDB", "D vs B;row B;row D;clusters",
	      162, -8.1, 8.1, 162, -8.1, 8.1 );
  TH1D hdxDB( "dxDB", "Dx-Bx;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyDB( "dyDB", "Dy-By;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcDB( "dxcDB", "Dx-Bx;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycDB( "dycDB", "Dy-By;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dxvsxDB( "dxvsxDB", "D-B dx vs x;x [mm];D-B <dx>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyDB( "dxvsyDB", "D-B dx vs y;y [mm];D-B <dx>",
		    81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxDB( "dyvsxDB", "D-B dy vs x;x [mm];D-B <dy>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyDB( "dyvsyDB", "D-B dy vs y;y [mm];D-B <dy>",
		    81, -8.1, 8.1, -1, 1 );

  TProfile mDBvst( "mDBvst", "DB matches vs time;trigger;DB matches",
		   200, 0E6, 2E6, -1, 2 );

  // triplets ACB:

  TH1D hdxACB( "dxACB", "ACB dx;x-x [mm];ACBs", 200, -1, 1 );
  TH1D hdyACB( "dyACB", "ACB dy;y-y [mm];ACBs", 200, -1, 1 );
  TH1D hdxcACB( "dxcACB", "ACB dx;x-x [um];ACBs", 200, -200, 200 );
  TH1D hdycACB( "dycACB", "ACB dy;y-y [um];ACBs", 200, -200, 200 );
  TH1D hdxciACB( "dxciACB", "ACB dx;x-x [um];isolated ACBs",
		 200, -200, 200 );
  TH1D hdyciACB( "dyciACB", "ACB dy;y-y [um];isolated ACBs",
		 200, -200, 200 );
  TH1D hdycfACB( "dycfACB", "ACB dy;y-y [um];inner ACBs",
		 200, -200, 200 );
  TH1D hdycfqACB( "dycfqACB", "ACB dy;y-y [um];Landau peak inner ACBs",
		  200, -200, 200 );

  TProfile dxvsxACB( "dxvsxACB",
		     "ACB dx vs x;x [mm];ACB <dx>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyACB( "dxvsyACB",
		     "ACB dx vs y;y [mm];ACB <dx>",
		     81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxACB( "dyvsxACB",
		     "ACB dy vs x;x [mm];ACB <dy>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyACB( "dyvsyACB",
		     "ACB dy vs y;y [mm];ACB <dy>",
		     81, -8.1, 8.1, -1, 1 );
  TProfile madyvsyACB( "madyvsyACB",
		       "ACB mady vs y;y [mm];ACB MAD(y) [um]",
		       81, -8.1, 8.1, 0, 100 );
  TH2D * hmapACB;
  hmapACB = new TH2D( "mapACB",
		      "ACBplet map;ACBplet col;ACBplet row;ACBplets",
		      8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1 );

  TH1D htxACB( "txACB", "tri angle ACB x;tri angle ACB x [mrad];ACB triplets", 100, -10, 10 );
  TH1D htyACB( "tyACB", "tri angle ACB y;tri angle ACB y [mrad];ACB triplets", 100, -10, 10 );

  // linked D:

  TH1D hsizD3( "clszD3", "D cluster size;pixels/cluster;D3 clusters",
	       51, -0.5, 50.5 );
  TH1D hclqD3( "clqD3", "D cluster charge;cluster charge [ke];D3 clusters",
	       100, 0, 100 );
  TH1D hncolD3( "ncolD3", "D cluster size;columns/cluster;D3 clusters",
		21, -0.5, 20.5 );
  TH1D hnrowD3( "nrowD3", "D cluster size;rows/cluster;D3 clusters",
		21, -0.5, 20.5 );

  // eff D vs BCA:

  TProfile effDvst1( "effDvst1", "effD vs time;trigger;eff D",
		     500, 0, 1E6, -1, 2 );
  TProfile effDvst5( "effDvst5", "effD vs time;trigger;eff D",
		     500, 0, 5E6, -1, 2 );
  TProfile effDvst40( "effDvst40", "effD vs time;trigger;eff D",
		      1000, 0, 40E6, -1, 2 );
  TProfile effDivst1( "effDivst1", "effD vs time;trigger;iso eff D",
		      500, 0, 1E6, -1, 2 );
  TProfile effDivst2( "effDvst21", "effD vs time;trigger;iso eff D",
		      400, 0, 2E6, -1, 2 );
  TProfile effDivst8( "effDivst8", "effD vs time;trigger;iso eff D",
		      500, 0, 8E6, -1, 2 );
  TProfile effDivst10( "effDivst10", "effD vs time;trigger;iso eff D",
		       500, 0, 10E6, -1, 2 );

  TProfile effDvsx0( "effDvsx0", "effD vs lower x;lower ACBplet x [mm];eff D",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effDvsx1( "effDvsx1", "effD vs upper x;upper ACBplet x [mm];eff D",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effDvsy( "effDvsy", "effD vs y;ACBplet y [mm];eff D",
		    81, -8.1, 8.1, -1, 2 );

  TProfile effDvsw( "effDvsw", "effD vs window;link window [mm];eff D",
		     99, 0.025, 4.975, -1, 2 );
  TProfile2D * effDmap1;
  effDmap1 = new TProfile2D( "effDmap1",
			     "D efficiency map;col;row;eff D",
			     8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1, -1, 2 );
  TProfile2D * effDmap4;
  effDmap4 = new TProfile2D( "effDmap4",
			     "D efficiency map;col;row;eff D",
			     4*54, -4*8.1, 4*8.1, 1*81, -8.1, 8.1, -1, 2 );

  // triplets BDC:

  TH1D hdxBDC( "dxBDC", "BDC dx;x-x [mm];BDCs", 200, -1, 1 );
  TH1D hdyBDC( "dyBDC", "BDC dy;y-y [mm];BDCs", 200, -1, 1 );
  TH1D hdxcBDC( "dxcBDC", "BDC dx;x-x [um];BDCs", 200, -200, 200 );
  TH1D hdycBDC( "dycBDC", "BDC dy;y-y [um];BDCs", 200, -200, 200 );
  TH1D hdxciBDC( "dxciBDC", "BDC dx;x-x [um];isolated BDCs",
		 200, -200, 200 );
  TH1D hdyciBDC( "dyciBDC", "BDC dy;y-y [um];isolated BDCs",
		 200, -200, 200 );
  TH1D hdycfBDC( "dycfBDC", "BDC dy;y-y [um];inner BDCs",
		 200, -200, 200 );
  TH1D hdycfqBDC( "dycfqBDC", "BDC dy;y-y [um];Landau peak inner BDCs",
		  200, -200, 200 );

  TProfile dxvsxBDC( "dxvsxBDC",
		     "BDC dx vs x;x [mm];BDC <dx>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyBDC( "dxvsyBDC",
		     "BDC dx vs y;y [mm];BDC <dx>",
		     81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxBDC( "dyvsxBDC",
		     "BDC dy vs x;x [mm];BDC <dy>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyBDC( "dyvsyBDC",
		     "BDC dy vs y;y [mm];BDC <dy>",
		     81, -8.1, 8.1, -1, 1 );
  TProfile madyvsyBDC( "madyvsyBDC",
		       "BDC mady vs y;y [mm];BDC MAD(y) [um]",
		       81, -8.1, 8.1, 0, 100 );
  TProfile madyvsxBDC( "madyvsxBDC",
		       "BDC mady vs x;x [mm];BDC MAD(y) [um]",
		       216, -32.4, 32.4, 0, 100 );
  TH2D * hmapBDC;
  hmapBDC = new TH2D( "mapBDC",
		      "BDCplet map;BDCplet col;BDCplet row;BDCplets",
		      8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1 );

  // eff A vs BCD:

  TProfile effAvst1( "effAvst1", "effA vs time;trigger;eff A",
		     500, 0, 1E6, -1, 2 );
  TProfile effAvst5( "effAvst5", "effA vs time;trigger;eff A",
		     500, 0, 5E6, -1, 2 );
  TProfile effAvst40( "effAvst40", "effA vs time;trigger;eff A",
		      1000, 0, 40E6, -1, 2 );
  TProfile effAivst1( "effAivst1", "effA vs time;trigger;iso eff A",
		      500, 0, 1E6, -1, 2 );
  TProfile effAivst2( "effAvst21", "effA vs time;trigger;iso eff A",
		      400, 0, 2E6, -1, 2 );
  TProfile effAivst8( "effAivst8", "effA vs time;trigger;iso eff A",
		      500, 0, 8E6, -1, 2 );
  TProfile effAivst10( "effAivst10", "effA vs time;trigger;iso eff A",
		       500, 0, 10E6, -1, 2 );

  TProfile effAvsx0( "effAvsx0", "effA vs lower x;lower BDCplet x [mm];eff A",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effAvsx1( "effAvsx1", "effA vs upper x;upper BDCplet x [mm];eff A",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effAvsy( "effAvsy", "effA vs y;BDCplet y [mm];eff A",
		    81, -8.1, 8.1, -1, 2 );

  TProfile effAvsw( "effAvsw", "effA vs window;link window [mm];eff A",
		     99, 0.025, 4.975, -1, 2 );
  TProfile2D * effAmap1;
  effAmap1 = new TProfile2D( "effAmap1",
			     "A efficiency map;col;row;eff A",
			     8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1, -1, 2 );
  TProfile2D * effAmap4;
  effAmap4 = new TProfile2D( "effAmap4",
			     "A efficiency map;col;row;eff A",
			     4*54, -4*8.1, 4*8.1, 1*81, -8.1, 8.1, -1, 2 );

  // linked A:

  TH1D hsizA3( "clszA3", "A cluster size;pixels/cluster;A3 clusters",
	       51, -0.5, 50.5 );
  TH1D hclqA3( "clqA3", "A cluster charge;cluster charge [ke];A3 clusters",
	       100, 0, 100 );
  TH1D hncolA3( "ncolA3", "A cluster size;columns/cluster;A3 clusters",
		21, -0.5, 20.5 );
  TH1D hnrowA3( "nrowA3", "A cluster size;rows/cluster;A3 clusters",
		21, -0.5, 20.5 );

  // A-D tracks:

  TH2D hxxDA( "xxDA", "D vs A;col A;col D;clusters",
	      432, -32.4, 32.4, 432, -32.4, 32.4 );
  TH2D hyyDA( "yyDA", "D vs A;row A;row D;clusters",
	      162, -8.1, 8.1, 162, -8.1, 8.1 );
  TH1D hdxDA( "dxDA", "Dx-Ax;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyDA( "dyDA", "Dy-Ay;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxcDA( "dxcDA", "Dx-Ax;x-x [mm];cluster pairs", 200, -1, 1 );
  TH1D hdycDA( "dycDA", "Dy-Ay;y-y [mm];cluster pairs", 200, -1, 1 );
  TProfile dxvsxDA( "dxvsxDA", "D-A dx vs x;x [mm];D-A <dx>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyDA( "dxvsyDA", "D-A dx vs y;y [mm];D-A <dx>",
		    81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxDA( "dyvsxDA", "D-A dy vs x;x [mm];D-A <dy>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyDA( "dyvsyDA", "D-A dy vs y;y [mm];D-A <dy>",
		    81, -8.1, 8.1, -1, 1 );

  // triplets ADC:

  TH1D hdxADC( "dxADC", "ADC dx;x-x [mm];ADCplets", 200, -1, 1 );
  TH1D hdyADC( "dyADC", "ADC dy;y-y [mm];ADCplets", 200, -1, 1 );
  TH1D hdxcADC( "dxcADC", "ADC dx;x-x [um];ADCplets", 200, -200, 200 );
  TH1D hdycADC( "dycADC", "ADC dy;y-y [um];ADCplets", 200, -200, 200 );
  TH1D hdxciADC( "dxciADC", "ADC dx;x-x [um];isolated ADCplets",
		 200, -200, 200 );
  TH1D hdyciADC( "dyciADC", "ADC dy;y-y [um];isolated ADCplets",
		 200, -200, 200 );

  TProfile dxvsxADC( "dxvsxADC", "ADCplet dx vs x;x [mm];ADCplet <dx>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyADC( "dxvsyADC", "ADCplet dx vs y;y [mm];ADCplet <dx>",
		     81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxADC( "dyvsxADC", "ADCplet dy vs x;x [mm];ADCplet <dy>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyADC( "dyvsyADC", "ADCplet dy vs y;y [mm];ADCplet <dy>",
		     81, -8.1, 8.1, -1, 1 );
  TH1D hxADC( "xADC", "ADCplets;col;ADCplets", 216, -32.4, 32.4 );
  TH1D hyADC( "yADC", "ADCplets;row;ADCplets",  81, -8.1, 8.1 );
  TH2D * hmapADC;
  hmapADC = new TH2D( "mapADC",
		      "ADCplet map;ADCplet col;ADCplet row;ADCplets",
		      8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1 );
  TH1D htxADC( "txADC", "ADCplet angle x;ADCplet angle x;ADCplets",
	       100, -1, 1 );
  TH1D htyADC( "tyADC", "ADCplet angle y;ADCplet angle y;ADCplets",
	       100, -1, 1 );

  // linked B:

  TH1D hsizB4( "clszB4", "B cluster size;pixels/cluster;B4 clusters",
	       51, -0.5, 50.5 );
  TH1D hclqB4( "clqB4", "B cluster charge;cluster charge [ke];B4 clusters",
	       100, 0, 100 );
  TH1D hncolB4( "ncolB4", "B cluster size;columns/cluster;B4 clusters",
		21, -0.5, 20.5 );
  TH1D hnrowB4( "nrowB4", "B cluster size;rows/cluster;B4 clusters",
		21, -0.5, 20.5 );

  // B vs ADC:

  TProfile effBvsx0( "effBvsx0", "effB vs lower x;lower ADCplet x [mm];eff B",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effBvsx1( "effBvsx1", "effB vs upper x;upper ADCplet x [mm];eff B",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effBvsy( "effBvsy", "effB vs y;ADCplet y [mm];eff B",
		    81, -8.1, 8.1, -1, 2 );
  TProfile effBvst1( "effBvst1", "effB vs time;trigger;eff B",
		     500, 0, 1E6, -1, 2 );
  TProfile effBvst5( "effBvst5", "effB vs time;trigger;eff B",
		     500, 0, 5E6, -1, 2 );
  TProfile effBvst40( "effBvst40", "effB vs time;trigger;eff B",
		      1000, 0, 40E6, -1, 2 );
  TProfile effBivst1( "effBivst1", "effB vs time;trigger;iso eff B",
		      500, 0, 1E6, -1, 2 );
  TProfile effBivst8( "effBivst8", "effB vs time;trigger;iso eff B",
		      500, 0, 8E6, -1, 2 );
  TProfile effBivst2( "effBvst21", "effB vs time;trigger;iso eff B",
		      400, 0, 2E6, -1, 2 );
  TProfile effBivst10( "effBivst10", "effB vs time;trigger;iso eff B",
		       500, 0, 10E6, -1, 2 );

  TProfile effBvsw( "effBvsw", "effB vs window;link window [mm];eff B",
		     40, 0.025, 2.025, -1, 2 );
  TProfile2D * effBmap1;
  effBmap1 = new TProfile2D( "effBmap1",
			     "B efficiency map;col;row;eff B",
			     8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1, -1, 2 );
  TProfile2D * effBmap4;
  effBmap4 = new TProfile2D( "effBmap4",
			     "B efficiency map;col;row;eff B",
			     4*54, -4*8.1, 4*8.1, 1*81, -8.1, 8.1, -1, 2 );
  // triplets ADB:

  TH1D hdxADB( "dxADB", "ADB dx;x-x [mm];ADBplets", 200, -1, 1 );
  TH1D hdyADB( "dyADB", "ADB dy;y-y [mm];ADBplets", 200, -1, 1 );
  TH1D hdxcADB( "dxcADB", "ADB dx;x-x [um];ADBplets", 200, -200, 200 );
  TH1D hdycADB( "dycADB", "ADB dy;y-y [um];ADBplets", 200, -200, 200 );
  TH1D hdxciADB( "dxciADB", "ADB dx;x-x [um];isolated ADBplets",
		 200, -200, 200 );
  TH1D hdyciADB( "dyciADB", "ADB dy;y-y [um];isolated ADBplets",
		 200, -200, 200 );

  TProfile dxvsxADB( "dxvsxADB", "ADBplet dx vs x;x [mm];ADBplet <dx>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyADB( "dxvsyADB", "ADBplet dx vs y;y [mm];ADBplet <dx>",
		     81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxADB( "dyvsxADB", "ADBplet dy vs x;x [mm];ADBplet <dy>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyADB( "dyvsyADB", "ADBplet dy vs y;y [mm];ADBplet <dy>",
		     81, -8.1, 8.1, -1, 1 );
  TH1D hxADB( "xADB", "ADBplets;col;ADBplets", 216, -32.4, 32.4 );
  TH1D hyADB( "yADB", "ADBplets;row;ADBplets",  81, -8.1, 8.1 );
  TH2D * hmapADB;
  hmapADB = new TH2D( "mapADB", "ADBplet map;ADBplet col;ADBplet row;ADBplets",
		      8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1 );
  TH1D htxADB( "txADB", "ADBplet angle x;ADBplet angle x;ADBplets",
	       100, -1, 1 );
  TH1D htyADB( "tyADB", "ADBplet angle y;ADBplet angle y;ADBplets",
	       100, -1, 1 );

  // linked C:

  TH1D hsizC4( "clszC4", "C cluster size;pixels/cluster;C4 clusters",
	       51, -0.5, 50.5 );
  TH1D hclqC4( "clqC4", "C cluster charge;cluster charge [ke];C4 clusters",
	       100, 0, 100 );
  TH1D hncolC4( "ncolC4", "C cluster size;columns/cluster;C4 clusters",
		21, -0.5, 20.5 );
  TH1D hnrowC4( "nrowC4", "C cluster size;rows/cluster;C4 clusters",
		21, -0.5, 20.5 );
  TH1D hminxC4( "minxC4", "C first pixel;first pixel mod 2;C4 clusters",
		2, -0.5, 1.5 );
  TH1D hmaxxC4( "maxxC4", "C last pixel;last pixel mod 2;C4 clusters",
		2, -0.5, 1.5 );

  // C vs ADB:

  TProfile effCvsx0( "effCvsx0", "effC vs lower x;lower ADBplet x [mm];eff C",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effCvsx1( "effCvsx1", "effC vs upper x;upper ADBplet x [mm];eff C",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effCvsy( "effCvsy", "effC vs y;ADBplet y [mm];eff C",
		    81, -8.1, 8.1, -1, 2 );
  TProfile effCvst1( "effCvst1", "effC vs time;trigger;eff C",
		     500, 0, 1E6, -1, 2 );
  TProfile effCvst5( "effCvst5", "effC vs time;trigger;eff C",
		     500, 0, 5E6, -1, 2 );
  TProfile effCvst40( "effCvst40", "effC vs time;trigger;eff C",
		      1000, 0, 40E6, -1, 2 );
  TProfile effCivst1( "effCivst1", "effC vs time;trigger;iso eff C",
		      500, 0, 1E6, -1, 2 );
  TProfile effCivst2( "effCivst2", "effC vs time;trigger;iso eff C",
		      400, 0, 2E6, -1, 2 );
  TProfile effCivst10( "effCivst10", "effC vs time;trigger;iso eff C",
		       500, 0, 10E6, -1, 2 );
  TProfile effCvsw( "effCvsw", "effC vs window;link window [mm];eff C",
		     40, 0.025, 2.025, -1, 2 );
  TProfile2D * effCmap1;
  effCmap1 = new TProfile2D( "effCmap1",
			     "C efficiency map;col;row;eff C",
			     8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1, -1, 2 );
  TProfile2D * effCmap4;
  effCmap4 = new TProfile2D( "effCmap4",
			     "C efficiency map;col;row;eff C",
			     4*54, -4*8.1, 4*8.1, 1*81, -8.1, 8.1, -1, 2 );

  TH1D hkxB( "kinkxB", "kink at B in x;kink x at B [mrad];tracks",
	     100, -10, 10 );
  TH1D hkyB( "kinkyB", "kink at B in y;kink y at B [mrad];tracks",
	     100, -10, 10 );
  TH1D hkxqB( "kinkxqB", "kink at B in x;kink x at B [mrad];Landau peak atracks",
	     100, -10, 10 );
  TH1D hkyqB( "kinkyqB", "kink at B in y;kink y at B [mrad];Landau peak tracks",
	     100, -10, 10 );

  TProfile2D * madkymapB =
    new TProfile2D( "madkymapB",
		    "B kink mad map;col;row;mad kink y [mrad]",
		    8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1, 0, 50 );
  TProfile2D * k2mapB =
    new TProfile2D( "k2mapB",
		    "B k_{xy}^{2} map;col;row;<k_{x}^{2}+k_{y}^{2} [mrad^{2}]",
		    8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1, 0, 1000 );

  TH1D hkxC( "kinkxC", "kink at C in x;kink x at C [mrad];tracks",
	     100, -10, 10 );
  TH1D hkyC( "kinkyC", "kink at C in y;kink y at C [mrad];tracks",
	     100, -10, 10 );
  TH1D hkxqC( "kinkxqC", "kink at C in x;kink x at C [mrad];Landau peak atracks",
	     100, -10, 10 );
  TH1D hkyqC( "kinkyqC", "kink at C in y;kink y at C [mrad];Landau peak tracks",
	     100, -10, 10 );

  TH1D hchi2( "chi2", "GBL chisq;chisq;track fits", 100, 0, 50 );
  TH1D hprob( "prob", "GBL fit prob;fit probability;track fits", 100, 0, 1 );

  TH1D hnADC( "nADC", "ADCplets;ADCplets;events", 21, -0.5, 20.5 );
  TH1D hnADB( "nADB", "ADBplets;ADBplets;events", 21, -0.5, 20.5 );
  TH1D hn4ev( "n4ev", "4plets;4plets;events", 21, -0.5, 20.5 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  int event_nr = 0;

  int n4 = 0;
  int nmille = 0;

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() )
      eudaq::PluginManager::Initialize(evt);

    int n4ev = 0;
    bool ldb = 0;
    if( event_nr == -1 )
      ldb = 1;

    if( ldb || event_nr%10000 == 0 )
      cout << "Quad processing event "<< event_nr << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    int xm = 0;
    int ym = 0;
    int adc = 0;
    double cal = 0;

    vector <cluster> cl[4];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      vector<double> pxl = plane.GetPixels<double>();

      if( ldb ) cout << "PLANE " << plane.ID() << ": ";

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int mod = 0; // QUAD
      if( plane.ID() == 6 ) mod = 1; // TRP
      if( plane.ID() == 7 ) mod = 2; // DUT
      if( plane.ID() == 8 ) mod = 3; // REF

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix) {

	if( ldb ) 
	  cout << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << " ";

	xm = plane.GetX(ipix); // global column 0..415
	ym = plane.GetY(ipix); // global row 0..159
	adc = plane.GetPixel(ipix); // ADC 0..255
	
	hcol[mod].Fill( xm );
	hrow[mod].Fill( ym );
	hmap[mod]->Fill( xm, ym );

	// leave space for big pixels:

	int roc = xm / 52; // 0..7
	int col = xm % 52; // 0..51
	int row = ym;
	int x = 1 + xm + 2*roc; // 1..52 per ROC with big pix
	int y = ym;
	if( ym > 79 ) y += 2;

	// flip for upper ROCs into local addresses:

	if( ym > 79 ) {
	  roc = 15 - roc; // 15..8
	  col = 51 - col; // 51..0
	  row = 159 - ym; // 79..0
	}

	cal = adc;
	if( xm < 0 || xm > 415 || ym < 0 || ym > 159 || adc < 0 || adc > 255 )
	  cout << "invalid pixel at event " << event_nr << endl;
	else if( haveGain[mod] ) {
	  double a0 = p0[mod][roc][col][row];
	  double a1 = p1[mod][roc][col][row];
	  double a2 = p2[mod][roc][col][row];
	  double a3 = p3[mod][roc][col][row];
	  double a4 = p4[mod][roc][col][row];
	  double a5 = p5[mod][roc][col][row];
	  cal = PHtoVcal( adc, a0, a1, a2, a3, a4, a5, mod ); // [Vcal]
	}
	
	hpxq[mod].Fill( cal*ke[mod] );

	// fill pixel block for clustering
	pb[npx].col = x;
	pb[npx].row = y;
	pb[npx].adc = adc;
	pb[npx].cal = cal;
	pb[npx].big = 0;
	++npx;

	// double big pixels:
	// 0+1
	// 2..51
	// 52+53

	col = xm % 52; // 0..51

	if( col == 0 ) {
	  pb[npx].col = x-1; // double
	  pb[npx].row = y;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

	if( col == 51 ) {
	  pb[npx].col = x+1; // double
	  pb[npx].row = y;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

	if( ym == 79 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 80;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

	if( ym == 80 ) {
	  pb[npx].col = x; // double
	  pb[npx].row = 81;
	  pb[npx-1].adc *= 0.5;
	  pb[npx-1].cal *= 0.5;
	  pb[npx].adc = 0.5*adc;
	  pb[npx].cal = 0.5*cal;
	  pb[npx].big = 1;
	  ++npx;
	}

      } // pix
      
      hnpx[mod].Fill(npx);

      if( ldb ) cout << endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[mod] = getClus();

      if( ldb ) cout << "A clusters " << cl[mod].size() << endl;

      hncl[mod].Fill( cl[mod].size() );

      for( vector<cluster>::iterator cA = cl[mod].begin(); cA != cl[mod].end(); ++cA ) {

	hsiz[mod].Fill( cA->size );
	hclq[mod].Fill( cA->charge*ke[mod] );
	hncol[mod].Fill( cA->ncol );
	hnrow[mod].Fill( cA->nrow );

      }

    } // planes = mod

    ++event_nr;

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // A-B cluster correlations:

    int mAB = 0;
    int mCB = 0;
    int mDB = 0;

    for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

      double xm = cB->col*0.15 - alignx[B] - 32.4;
      double ym = cB->row*0.10 - aligny[B] -  8.1;
      double xB = xm - ym*fx[B] - tx[B]*xm;
      double yB = ym + xm*fy[B] - ty[B]*ym;

      for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

	double xm = cA->col*0.15 - alignx[A] - 32.4;
	double ym = cA->row*0.10 - aligny[A] -  8.1;
	double xA = xm - ym*fx[A] - tx[A]*xm;
	double yA = ym + xm*fy[A] - ty[A]*ym;

	hxxAB.Fill( xB, xA );
	hyyAB.Fill( yB, yA );

	double dx = xA - xB;
	double dy = yA - yB;
	hdxAB.Fill( dx ); // includes angular spread
	hdyAB.Fill( dy );
	if( abs( dy ) < bicuty && cA->big == 0 && cB->big == 0 ) {
	  hdxcAB.Fill( dx );
	  dxvsxAB.Fill( xB, dx );
	  dxvsyAB.Fill( yB, dx );
	}
	if( abs( dx ) < bicutx && cA->big == 0 && cB->big == 0 ) {
	  hdycAB.Fill( dy );
	  dyvsxAB.Fill( xA, dy );
	  dyvsyAB.Fill( yA, dy );
	}
	if( abs( dx ) < bicutx && abs( dy ) < bicuty )
	  mAB = 1;

      } // clusters A

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // C-B cluster correlations:

      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	double xm = cC->col*0.15 - alignx[C] - 32.4;
	double ym = cC->row*0.10 - aligny[C] -  8.1;
	double xC = xm - ym*fx[C] - tx[C]*xm;
	double yC = ym + xm*fy[C] - ty[C]*ym;

	hxxCB.Fill( xB, xC );
	hyyCB.Fill( yB, yC );

	double dx = xC - xB;
	double dy = yC - yB;
	hdxCB.Fill( dx );
	hdyCB.Fill( dy );
	if( abs( dy ) < bicuty && cC->big == 0 && cB->big == 0 ) {
	  hdxcCB.Fill( dx );
	  dxvsxCB.Fill( xC, dx );
	  dxvsyCB.Fill( yC, dx );
	}
	if( abs( dx ) < bicutx && cC->big == 0 && cB->big == 0 ) {
	  hdycCB.Fill( dy );
	  dyvsxCB.Fill( xC, dy );
	  dyvsyCB.Fill( yC, dy );
	}

	if( abs( dx ) < bicutx && abs( dy ) < bicuty )
	  mCB = 1;

      } // clusters C

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // D-B correlations:

      for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	double xm = cD->col*0.15 - alignx[D] - 32.4;
	double ym = cD->row*0.10 - aligny[D] -  8.1;
	double xD = xm - ym*fx[D] - tx[D]*xm;
	double yD = ym + xm*fy[D] - ty[D]*ym;

	hxxDB.Fill( xB, xD );
	hyyDB.Fill( yB, yD );

	double dx = xD - xB;
	double dy = yD - yB;
	hdxDB.Fill( dx );
	hdyDB.Fill( dy );
	if( abs( dy ) < bicuty && cD->big == 0 && cB->big == 0 ) {
	  hdxcDB.Fill( dx );
	  dxvsxDB.Fill( xD, dx );
	  dxvsyDB.Fill( yD, dx );
	}
	if( abs( dx ) < bicutx && cD->big == 0 && cB->big == 0 ) {
	  hdycDB.Fill( dy );
	  dyvsxDB.Fill( xD, dy );
	  dyvsyDB.Fill( yD, dy );
	}

	if( abs( dx ) < bicutx && abs( dy ) < bicuty )
	  mDB = 1;

      } // cl D

    } // cl B

    mABvst.Fill( event_nr, mAB );
    mCBvst.Fill( event_nr, mCB );
    mDBvst.Fill( event_nr, mDB );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // unbiased triplet ACB:

    bool iso = cl[A].size() == 1 && cl[C].size() == 1;

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      double xm = cA->col*0.15 - alignx[A] - 32.4;
      double ym = cA->row*0.10 - aligny[A] -  8.1;
      double xA = xm - ym*fx[A] - tx[A]*xm;
      double yA = ym + xm*fy[A] - ty[A]*ym;

      double qA = cA->charge*ke[A];
      bool lqA = 1;
      if(      qA < qL ) lqA = 0;
      else if( qA > qR ) lqA = 0;

      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	double xm = cC->col*0.15 - alignx[C] - 32.4;
	double ym = cC->row*0.10 - aligny[C] -  8.1;
	double xC = xm - ym*fx[C] - tx[C]*xm;
	double yC = ym + xm*fy[C] - ty[C]*ym;

	double qC = cC->charge*ke[C];
	bool lqC = 1;
	if(      qC < 17 ) lqC = 0;
	else if( qC > 30 ) lqC = 0;

	// A-C track:

	double xavg2B = 0.5*(xA + xC); // interpolate
	double yavg2B = 0.5*(yA + yC); // equidistant

	double slpx = (xC - xA)/2/dz; // angle
	double slpy = (yC - yA)/2/dz; // angle

	for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	  double xm = cB->col*0.15 - alignx[B] - 32.4;
	  double ym = cB->row*0.10 - aligny[B] -  8.1;
	  double xB = xm - ym*fx[B] - tx[B]*xm;
	  double yB = ym + xm*fy[B] - ty[B]*ym;

	  double qB = cB->charge*ke[B];
	  bool lqB = 1;
	  if(      qB < 17 ) lqB = 0;
	  else if( qB > 30 ) lqB = 0;

	  // tri ACB:

	  double dx3 = xB - xavg2B;
	  double dy3 = yB - yavg2B;

	  hdxACB.Fill( dx3 );
	  hdyACB.Fill( dy3 );

	  if( abs( dy3 ) < tricuty
	      && cA->big == 0 && cC->big == 0 && cB->big == 0 ) {
	    hdxcACB.Fill( dx3*1E3 );
	    if( iso ) hdxciACB.Fill( dx3*1E3 );
	    dxvsxACB.Fill( xavg2B, dx3 );
	    dxvsyACB.Fill( yavg2B, dx3 );
	  }
	  if( abs( dx3 ) < tricutx
	      && cA->big == 0 && cC->big == 0 && cB->big == 0 ) {
	    hdycACB.Fill( dy3*1E3 );
	    if( iso ) hdyciACB.Fill( dy3*1E3 );
	    dyvsxACB.Fill( xavg2B, dy3 );
	    dyvsyACB.Fill( yavg2B, dy3 );
	    madyvsyACB.Fill( yavg2B, fabs(dy3)*1E3 );
	    if( yavg2B > -6 && yavg2B < 7 ) { // module handle cutout
	      hdycfACB.Fill( dy3*1E3 );
	      if( lqA && lqC && lqB )
		hdycfqACB.Fill( dy3*1E3 );
	    }
	  }

	  if( abs( dx3 ) > tricutx ) continue; // tight tri
	  if( abs( dy3 ) > tricuty ) continue;

	  hmapACB->Fill( xB, yB );

	  htxACB.Fill( slpx );
	  htyACB.Fill( slpy );

	  // efficiency of D vs BCA:

	  double xavg3D = (3*xC - xA)/2; // extrapolate
	  double yavg3D = (3*yC - yA)/2; // 

	  // Transform to local coordinate system
	  double yavg3Dlocal = (yavg3D-xavg3D*fy[D]/(1-ty[D]))/(1-ty[D]+fx[D]*fy[D]/(1-ty[D]));
	  double xavg3Dlocal = (xavg3D+yavg3Dlocal*fx[D])/(1-tx[D]);

	  bool fiducial = isFiducial(xavg3Dlocal, yavg3Dlocal);
	  if(!fiducial)continue;
	  
	  int nm[99] = {0};

	  for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	    double xm = cD->col*0.15 - alignx[D] - 32.4;
	    double ym = cD->row*0.10 - aligny[D] -  8.1;
	    double xD = xm - ym*fx[D] - tx[D]*xm;
	    double yD = ym + xm*fy[D] - ty[D]*ym;

	    double dx4 = xD - xavg3D;
	    double dy4 = yD - yavg3D;

	    for( int iw = 1; iw < 99; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) // for eff
		nm[iw] = 1;
	    
	  } // cl D

	  effDvst1.Fill( event_nr, nm[14] );
	  effDvst5.Fill( event_nr, nm[14] );
	  effDvst40.Fill( event_nr, nm[14] );
	  if( iso ) {
	    effDivst1.Fill( event_nr, nm[14] );
	    effDivst2.Fill( event_nr, nm[14] );
	    effDivst10.Fill( event_nr, nm[14] );
	  }
	  if( yavg3Dlocal < 0 )
	    effDvsx0.Fill( xavg3Dlocal, nm[14] );
	  else
	    effDvsx1.Fill( xavg3Dlocal, nm[14] );
	  effDvsy.Fill( yavg3Dlocal, nm[14] );
	  effDmap1->Fill( xavg3Dlocal, yavg3Dlocal, nm[14] );
	  effDmap4->Fill( xavg3Dlocal, yavg3Dlocal, nm[14] );
	  for( int iw = 1; iw < 99; ++ iw )
	    effDvsw.Fill( iw*0.050+0.005, nm[iw] );

	} // cl B

      } // cl C

    } // cl A

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // unbiased triplet BDC:

    iso = cl[B].size() == 1 && cl[D].size() == 1;

    for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

      double xm = cB->col*0.15 - alignx[B] - 32.4;
      double ym = cB->row*0.10 - aligny[B] -  8.1;
      double xB = xm - ym*fx[B] - tx[B]*xm;
      double yB = ym + xm*fy[B] - ty[B]*ym;

      double qB = cB->charge*ke[B];
      bool lqB = 1;
      if(      qB < qL ) lqB = 0;
      else if( qB > qR ) lqB = 0;

      for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	double xm = cD->col*0.15 - alignx[D] - 32.4;
	double ym = cD->row*0.10 - aligny[D] -  8.1;
	double xD = xm - ym*fx[D] - tx[D]*xm;
	double yD = ym + xm*fy[D] - ty[D]*ym;

	double qD = cD->charge*ke[D];
	bool lqD = 1;
	if(      qD < qL ) lqD = 0;
	else if( qD > qR ) lqD = 0;

	// B-D track:

	double xavg2C = 0.5*(xB + xD); // interpolate
	double yavg2C = 0.5*(yB + yD); // equidistant

	//double slpx = (xD - xB)/2/dz; // angle
	//double slpy = (yD - yB)/2/dz; // angle

	for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	  double xm = cC->col*0.15 - alignx[C] - 32.4;
	  double ym = cC->row*0.10 - aligny[C] -  8.1;
	  double xC = xm - ym*fx[C] - tx[C]*xm;
	  double yC = ym + xm*fy[C] - ty[C]*ym;

	  double qC = cC->charge*ke[C];
	  bool lqC = 1;
	  if(      qC < 17 ) lqC = 0;
	  else if( qC > 30 ) lqC = 0;

	  // tri BDC:

	  double dx3 = xC - xavg2C;
	  double dy3 = yC - yavg2C;

	  hdxBDC.Fill( dx3 );
	  hdyBDC.Fill( dy3 );

	  if( abs( dy3 ) < tricuty
	      && cB->big == 0 && cD->big == 0 && cC->big == 0 ) {
	    hdxcBDC.Fill( dx3*1E3 );
	    if( iso ) hdxciBDC.Fill( dx3*1E3 );
	    dxvsxBDC.Fill( xavg2C, dx3 );
	    dxvsyBDC.Fill( yavg2C, dx3 );
	  }
	  if( abs( dx3 ) < tricutx
	      && cB->big == 0 && cD->big == 0 && cC->big == 0 ) {
	    hdycBDC.Fill( dy3*1E3 );
	    if( iso ) hdyciBDC.Fill( dy3*1E3 );
	    dyvsxBDC.Fill( xavg2C, dy3 );
	    dyvsyBDC.Fill( yavg2C, dy3 );
	    madyvsyBDC.Fill( yavg2C, fabs(dy3)*1E3 );
	    if( yavg2C > -6 && yavg2C < 7 ) { // module handle cutout
	      madyvsxBDC.Fill( xavg2C, fabs(dy3)*1E3 );
	      hdycfBDC.Fill( dy3*1E3 );
	      if( lqB && lqD && lqC )
		hdycfqBDC.Fill( dy3*1E3 );
	    }
	  }

	  if( abs( dx3 ) > tricutx ) continue; // tight tri
	  if( abs( dy3 ) > tricuty ) continue;

	  hmapBDC->Fill( xavg2C, yavg2C );

	  // efficiency of A:

	  double xavg3A = (3*xB - xD)/2; // extrapolate
	  double yavg3A = (3*yB - yD)/2; // 

	  // Transform to local coordinate system
	  double yavg3Alocal = (yavg3A-xavg3A*fy[A]/(1-ty[A]))/(1-ty[A]+fx[A]*fy[A]/(1-ty[A]));
	  double xavg3Alocal = (xavg3A+yavg3Alocal*fx[A])/(1-tx[A]);

	  bool fiducial = isFiducial(xavg3Alocal, yavg3Alocal);
	  if(!fiducial) continue;

	  int nm[99] = {0};

	  for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

	    double xm = cA->col*0.15 - alignx[A] - 32.4;
	    double ym = cA->row*0.10 - aligny[A] -  8.1;
	    double xA = xm - ym*fx[A] - tx[A]*xm;
	    double yA = ym + xm*fy[A] - ty[A]*ym;

	    // tri ABC:

	    double dx4 = xA - xavg3A;
	    double dy4 = yA - yavg3A;

	    for( int iw = 1; iw < 99; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) // for eff
		nm[iw] = 1;

	  } // cl A

	  effAvst1.Fill( event_nr, nm[14] );
	  effAvst5.Fill( event_nr, nm[14] );
	  effAvst40.Fill( event_nr, nm[14] );

	  if( iso ) {
	    effAivst1.Fill( event_nr, nm[14] );
	    effAivst2.Fill( event_nr, nm[14] );
	    effAivst8.Fill( event_nr, nm[14] );
	    effAivst10.Fill( event_nr, nm[14] );
	  } // iso

	  if( yavg3Alocal < 0 )
	    effAvsx0.Fill( xavg3Alocal, nm[14] );
	  else
	    effAvsx1.Fill( xavg3Alocal, nm[14] );
	  effAvsy.Fill( yavg3Alocal, nm[14] );
	  effAmap1->Fill( xavg3Alocal, yavg3Alocal, nm[14] );
	  effAmap4->Fill( xavg3Alocal, yavg3Alocal, nm[14] );
	  for( int iw = 1; iw < 99; ++ iw )
	    effAvsw.Fill( iw*0.050+0.005, nm[iw] );

	} // cl C

      } // cl D

    } // cl B

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // eff of B and C vs A-D:

    int nADC = 0;
    int nADB = 0;

    iso = cl[A].size() == 1 && cl[D].size() == 1;

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      double xm = cA->col*0.15 - alignx[A] - 32.4;
      double ym = cA->row*0.10 - aligny[A] -  8.1;
      double xA = xm - ym*fx[A] - tx[A]*xm;
      double yA = ym + xm*fy[A] - ty[A]*ym;

      TMatrixD derivA( 2, 6 ); // -alignment derivatives x,y

      derivA[0][0] = 1.0; // -dresidx/dalignx
      derivA[1][0] = 0.0;

      derivA[0][1] = 0.0;
      derivA[1][1] = 1.0; // dresidy/daligny

      derivA[0][2] = ym; // dresidx/dfx
      derivA[1][2] = 0.0;

      derivA[0][3] = 0.0;
      derivA[1][3] =-xm; // dresidy/dfy

      derivA[0][4] = xm; // dresidx/dtx
      derivA[1][4] = 0.0;

      derivA[0][5] = 0.0;
      derivA[1][5] = ym; // dresidy/dty

      double qA = cA->charge*ke[A];
      bool lqA = 1;
      if(      qA < qL ) lqA = 0;
      else if( qA > qR ) lqA = 0;

      for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	double xm = cD->col*0.15 - alignx[D] - 32.4;
	double ym = cD->row*0.10 - aligny[D] -  8.1;
	double xD = xm - ym*fx[D] - tx[D]*xm;
	double yD = ym + xm*fy[D] - ty[D]*ym;

	TMatrixD derivD( 2, 6 ); // alignment derivatives x,y

	derivD[0][0] = 1.0; // dresidx/dalignx
	derivD[1][0] = 0.0;

	derivD[0][1] = 0.0;
	derivD[1][1] = 1.0; // dresidy/daligny

	derivD[0][2] = ym; // dresidx/dfx
	derivD[1][2] = 0.0;

	derivD[0][3] = 0.0;
	derivD[1][3] =-xm; // dresidy/dfy

	derivD[0][4] = xm; // dresidx/dtx
	derivD[1][4] = 0.0;

	derivD[0][5] = 0.0;
	derivD[1][5] = ym; // dresidy/dty

	double qD = cD->charge*ke[D];
	bool lqD = 1;
	if(      qD < qL ) lqD = 0;
	else if( qD > qR ) lqD = 0;

	hxxDA.Fill( xA, xD );
	hyyDA.Fill( yA, yD );

	double dx2 = xD - xA;
	double dy2 = yD - yA;
	hdxDA.Fill( dx2 );
	hdyDA.Fill( dy2 );
	if( abs( dy2 ) < bicuty && cD->big == 0 && cA->big == 0 ) {
	  hdxcDA.Fill( dx2 );
	  dxvsxDA.Fill( xD, dx2 );
	  dxvsyDA.Fill( yD, dx2 );
	}
	if( abs( dx2 ) < bicutx && cD->big == 0 && cA->big == 0 ) {
	  hdycDA.Fill( dy2 );
	  dyvsxDA.Fill( xD, dy2 );
	  dyvsyDA.Fill( yD, dy2 );
	}

	//if( abs( dx2 ) > bicutx ) continue; // angle cut
	//if( abs( dy2 ) > bicuty ) continue;

	double slpx = (xD - xA)/3/dz; // angle
	double slpy = (yD - yA)/3/dz; // angle

	double xavg3C = (xA + 2*xD)/3; // interpolate
	double yavg3C = (yA + 2*yD)/3; // A and D to C

	// tri ADC:

	for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	  double xm = cC->col*0.15 - alignx[C] - 32.4;
	  double ym = cC->row*0.10 - aligny[C] -  8.1;
	  double xC = xm - ym*fx[C] - tx[C]*xm;
	  double yC = ym + xm*fy[C] - ty[C]*ym;

	  double dx3 = xC - xavg3C;
	  double dy3 = yC - yavg3C;
	  hdxADC.Fill( dx3 );
	  hdyADC.Fill( dy3 );
	  if( abs( dy3 ) < tricuty && cA->big == 0 && cD->big == 0 && cC->big == 0 ) {
	    hdxcADC.Fill( dx3*1E3 );
	    if( iso ) hdxciADC.Fill( dx3*1E3 );
	    dxvsxADC.Fill( xavg3C, dx3 );
	    dxvsyADC.Fill( yavg3C, dx3 );
	  }
	  if( abs( dx3 ) < tricutx && cA->big == 0 && cD->big == 0 && cC->big == 0 ) {
	    hdycADC.Fill( dy3*1E3 );
	    if( iso ) hdyciADC.Fill( dy3*1E3 );
	    dyvsxADC.Fill( xavg3C, dy3 );
	    dyvsyADC.Fill( yavg3C, dy3 );
	  }
	  if( abs( dx3 ) < tricutx && abs( dy3 ) < tricuty ) {
	    hxADC.Fill( xavg3C );
	    hyADC.Fill( yavg3C );
	    hmapADC->Fill( xavg3C, yavg3C ); // D-C-A
	  }

	  if( abs( dx3 ) > tricutx ) continue; // tight tri
	  if( abs( dy3 ) > tricuty ) continue;

	  ++nADC;

	  htxADC.Fill( slpx*1E3 );
	  htyADC.Fill( slpy*1E3 );

	  // for tri ACB:

	  double xavg2B = 0.5*(xA + xC); // interpolate
	  double yavg2B = 0.5*(yA + yC); // equidistant

	  // Transform to local coordinate system
	  double yavg2Blocal = (yavg2B-xavg2B*fy[B]/(1-ty[B]))/(1-ty[B]+fx[B]*fy[B]/(1-ty[B]));
	  double xavg2Blocal = (xavg2B+yavg2Blocal*fx[B])/(1-tx[B]);

	  bool fiducial = isFiducial(xavg2Blocal, yavg2Blocal);
	  if(!fiducial) continue;

	  // efficiency of B:

	  int nm[99] = {0};

	  for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	    double xm = cB->col*0.15 - alignx[B] - 32.4;
	    double ym = cB->row*0.10 - aligny[B] -  8.1;
	    double xB = xm - ym*fx[B] - tx[B]*xm;
	    double yB = ym + xm*fy[B] - ty[B]*ym;

	    // tri ACB:

	    double dx4 = xB - xavg2B;
	    double dy4 = yB - yavg2B;

	    for( int iw = 1; iw < 41; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) // for eff
		nm[iw] = 1;

	    if( abs( dx4 ) < tricutx && abs( dy4 ) < tricuty  &&
		cA->big == 0 && cC->big == 0 && cB->big == 0 ) {
	      hsizB4.Fill( cB->size );
	      hclqB4.Fill( cB->charge*ke[B] );
	      hncolB4.Fill( cB->ncol );
	      hnrowB4.Fill( cB->nrow );
	    }

	  } // cl B

	  effBvst1.Fill( event_nr, nm[14] );
	  effBvst5.Fill( event_nr, nm[14] );
	  effBvst40.Fill( event_nr, nm[14] );

	  if( iso ) {
	    effBivst1.Fill( event_nr, nm[14] );
	    effBivst2.Fill( event_nr, nm[14] );
	    effBivst8.Fill( event_nr, nm[14] );
	    effBivst10.Fill( event_nr, nm[14] );
	  } // iso

	  if( yavg2Blocal < 0 )
	    effBvsx0.Fill( xavg2Blocal, nm[14] );
	  else
	    effBvsx1.Fill( xavg2Blocal, nm[14] );
	  effBvsy.Fill( yavg2Blocal, nm[14] );
	  effBmap1->Fill( xavg2Blocal, yavg2Blocal, nm[14] );
	  effBmap4->Fill( xavg2Blocal, yavg2Blocal, nm[14] );
	  for( int iw = 1; iw < 41; ++ iw )
	    effBvsw.Fill( iw*0.050+0.005, nm[iw] );

	} // cl C

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// tri ADB:

	double xavg3B = (2*xA + xD)/3; // interpolate
	double yavg3B = (2*yA + yD)/3; // A and D to B

	for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	  double xm = cB->col*0.15 - alignx[B] - 32.4;
	  double ym = cB->row*0.10 - aligny[B] -  8.1;
	  double xB = xm - ym*fx[B] - tx[B]*xm;
	  double yB = ym + xm*fy[B] - ty[B]*ym;

	  TMatrixD derivB( 2, 6 ); // alignment derivatives x,y

	  derivB[0][0] = 1.0; // dresidx/dalignx
	  derivB[1][0] = 0.0;

	  derivB[0][1] = 0.0;
	  derivB[1][1] = 1.0; // dresidy/daligny

	  derivB[0][2] = ym; // dresidx/dfx
	  derivB[1][2] = 0.0;

	  derivB[0][3] = 0.0;
	  derivB[1][3] =-xm; // dresidy/dfy

	  derivB[0][4] = xm; // dresidx/dtx
	  derivB[1][4] = 0.0;

	  derivB[0][5] = 0.0;
	  derivB[1][5] = ym; // dresidy/dty

	  double qB = cB->charge*ke[B];
	  bool lqB = 1;
	  if(      qB < qL ) lqB = 0;
	  else if( qB > qR ) lqB = 0;

	  double dx3 = xB - xavg3B;
	  double dy3 = yB - yavg3B;

	  hdxADB.Fill( dx3 );
	  hdyADB.Fill( dy3 );

	  if( abs( dy3 ) < tricuty && cA->big == 0 && cD->big == 0 && cB->big == 0 ) {
	    hdxcADB.Fill( dx3*1E3 );
	    if( iso ) hdxciADB.Fill( dx3*1E3 );
	    dxvsxADB.Fill( xavg3B, dx3 );
	    dxvsyADB.Fill( yavg3B, dx3 );
	  }
	  if( abs( dx3 ) < tricutx && cA->big == 0 && cD->big == 0 && cB->big == 0 ) {
	    hdycADB.Fill( dy3*1E3 );
	    if( iso ) hdyciADB.Fill( dy3*1E3 );
	    dyvsxADB.Fill( xavg3B, dy3 );
	    dyvsyADB.Fill( yavg3B, dy3 );
	  }
	  if( abs( dx3 ) < tricutx && abs( dy3 ) < tricuty ) {
	    hxADB.Fill( xavg3B );
	    hyADB.Fill( yavg3B );
	    hmapADB->Fill( xavg3B, yavg3B );
	  }

	  if( abs( dx3 ) > tricutx ) continue; // tight tri
	  if( abs( dy3 ) > tricuty ) continue;

	  ++nADB;

	  htxADB.Fill( slpx );
	  htyADB.Fill( slpy );

	  // B-D track:

	  double xavg2C = 0.5*(xB + xD); // equidistant
	  double yavg2C = 0.5*(yB + yD);

	  // Transform to local coordinate system
	  double yavg2Clocal = (yavg2C-xavg2C*fy[C]/(1-ty[C]))/(1-ty[C]+fx[C]*fy[C]/(1-ty[C]));
	  double xavg2Clocal = (xavg2C+yavg2Clocal*fx[C])/(1-tx[C]);
	  
	  bool fiducial = isFiducial(xavg2Clocal, yavg2Clocal);
	  if(!fiducial) continue;

	  // efficiency of C vs ADB:

	  int nm[99] = {0};

	  for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	    double xm = cC->col*0.15 - alignx[C] - 32.4;
	    double ym = cC->row*0.10 - aligny[C] -  8.1;
	    double xC = xm - ym*fx[C] - tx[C]*xm;
	    double yC = ym + xm*fy[C] - ty[C]*ym;

	    TMatrixD derivC( 2, 6 ); // alignment derivatives x,y

	    derivC[0][0] = 1.0; // dresidx/dalignx
	    derivC[1][0] = 0.0;

	    derivC[0][1] = 0.0;
	    derivC[1][1] = 1.0; // dresidy/daligny

	    derivC[0][2] = ym; // dresidx/dfx
	    derivC[1][2] = 0.0;

	    derivC[0][3] = 0.0;
	    derivC[1][3] =-xm; // dresidy/dfy

	    derivC[0][4] = xm; // dresidx/dtx
	    derivC[1][4] = 0.0;

	    derivC[0][5] = 0.0;
	    derivC[1][5] = ym; // dresidy/dty

	    double qC = cC->charge*ke[C];
	    bool lqC = 1;
	    if(      qC < qL ) lqC = 0;
	    else if( qC > qR ) lqC = 0;

	    double dx4 = xC - xavg2C;
	    double dy4 = yC - yavg2C;

	    for( int iw = 1; iw < 41; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) // for eff
		nm[iw] = 1;

	    if( abs( dx4 ) > tricutx ) continue; // quad
	    if( abs( dy4 ) > tricuty ) continue;

	    int minx = 999;
	    int maxx = 0;
	    int miny = 999;
	    int maxy = 0;

	    for( vector<pixel>::iterator p = cC->vpix.begin(); p != cC->vpix.end(); ++p ) {
	      if( p->col > maxx ) maxx = p->col;
	      if( p->col < minx ) minx = p->col;
	      if( p->row > maxy ) maxy = p->row;
	      if( p->row < miny ) miny = p->row;
	    }

	    if( cB->big == 0 && cD->big == 0 && cC->big == 0 ) {
	      hsizC4.Fill( cC->size );
	      hclqC4.Fill( cC->charge*ke[C] );
	      hncolC4.Fill( cC->ncol );
	      hnrowC4.Fill( cC->nrow );
	      hminxC4.Fill( (minx-1)%2 ); 
	      hmaxxC4.Fill( (maxx-1)%2 ); 
	    }

	    // we have a quad track !

	    ++n4ev;

	    if( cA->big > 0 ) continue;
	    if( cB->big > 0 ) continue;
	    if( cC->big > 0 ) continue;
	    if( cD->big > 0 ) continue;

	    ++n4;

	    // slopes and kinks:

	    double sxBA = (xB-xA) / dz; // track angle [rad]
	    double syBA = (yB-yA) / dz;
	    double sxCB = (xC-xB) / dz;
	    double syCB = (yC-yB) / dz;
	    double sxDC = (xD-xC) / dz;
	    double syDC = (yD-yC) / dz;

	    double kxB = sxCB - sxBA;
	    double kyB = syCB - syBA;
	    double k2B = kxB*kxB + kyB*kyB;
	    double kxC = sxDC - sxCB;
	    double kyC = syDC - syCB;

	    if( yB > -6 && yB < 7 ) { // Al handle cut out fiducial region
	      hkxB.Fill( kxB*1E3 );
	      hkyB.Fill( kyB*1E3 );
	      if( lqA && lqB && lqC ) {
		hkxqB.Fill( kxB*1E3 );
		hkyqB.Fill( kyB*1E3 );
	      }
	      k2mapB->Fill( xB, yB, k2B*1E6 );
	      madkymapB->Fill( xB, yB, abs(kyB)*1E3 );
	    }
	    if( yC > -6 && yC < 7 ) { // module handle cutout
	      hkxC.Fill( kxC*1E3 );
	      hkyC.Fill( kyC*1E3 );
	      if( lqB && lqD && lqC ) {
		hkxC.Fill( kxC*1E3 );
		hkyC.Fill( kyC*1E3 );
	      }
	    }

	    // GBL track fit:

	    vector<GblPoint> listOfPoints;
	    listOfPoints.reserve(4);
	    vector<double> sPoint;

	    // plane A:

	    TMatrixD jacPointToPoint(5, 5);
	    jacPointToPoint.UnitMatrix();
	    GblPoint *point = new GblPoint(jacPointToPoint);
	    meas[0] = 0;
	    meas[1] = 0;
	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsA, derivA ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // B:

	    jacPointToPoint = Jac5( dz );
	    point = new GblPoint(jacPointToPoint);
	    meas[0] = dx3;
	    meas[1] = dy3;
	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsB, derivB ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // C:

	    jacPointToPoint = Jac5( dz );
	    point = new GblPoint(jacPointToPoint);
	    meas[0] = xC - xavg3C;
	    meas[1] = yC - yavg3C;
	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsC, derivC ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // D:

	    jacPointToPoint = Jac5( dz );
	    point = new GblPoint(jacPointToPoint);
	    meas[0] = 0;
	    meas[1] = 0;
	    point->addMeasurement( proL2m, meas, measPrec );
	    //point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsD, derivD ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // track fit:

	    GblTrajectory traj( listOfPoints, 0 ); // 0 = no magnetic field
	    //traj.printPoints();

	    double Chi2;
	    int Ndf;
	    double lostWeight;

	    traj.fit( Chi2, Ndf, lostWeight );

	    double probchi = TMath::Prob( Chi2, Ndf );

	    hchi2.Fill( Chi2 );
	    hprob.Fill( probchi );

	    //cout << " Fit: " << Chi2 << ", " << Ndf << ", " << lostWeight << endl;

	    //traj.printTrajectory();
	    //traj.printPoints();

	    TVectorD aCorrection(5);
	    TMatrixDSym aCovariance(5);

	    // at plane C:

	    int ipos = 3; // starts at 1
	    traj.getResults( ipos, aCorrection, aCovariance );

	    //cout << " corrections: " << endl;
	    //aCorrection.Print();
	    //cout << " covariance: " << endl;
	    //aCovariance.Print();
	    //cout << "  sigma(x) = " << sqrt(aCovariance(3,3))*1E3 << " um";
	    //cout << ", sigma(y) = " << sqrt(aCovariance(4,4))*1E3 << " um";
	    //cout << endl;

	    // write to MP binary file

	    if( probchi > 0.01 ) { // bias with bad alignment ?
	      traj.milleOut( *mille );
	      ++nmille;
	    }
	    
	  } // cl C

	  effCvst1.Fill( event_nr, nm[14] );
	  effCvst5.Fill( event_nr, nm[14] );
	  effCvst40.Fill( event_nr, nm[14] );
	  if( iso ) {
	    effCivst1.Fill( event_nr, nm[14] );
	    effCivst2.Fill( event_nr, nm[14] );
	    effCivst10.Fill( event_nr, nm[14] );
	  }
	  if( yavg2Clocal < 0 )
	    effCvsx0.Fill( xavg2Clocal, nm[14] );
	  else
	    effCvsx1.Fill( xavg2Clocal, nm[14] );
	  effCvsy.Fill( yavg2Clocal, nm[14] );
	  effCmap1->Fill( xavg2Clocal, yavg2Clocal, nm[14] );
	  effCmap4->Fill( xavg2Clocal, yavg2Clocal, nm[14] );
	  for( int iw = 1; iw < 41; ++ iw )
	    effCvsw.Fill( iw*0.050+0.005, nm[iw] );

	} // cl B

      } // cl D

    } // cl A

    hnADC.Fill( nADC );
    hnADB.Fill( nADB );
    hn4ev.Fill( n4ev );

  } while( reader->NextEvent() && event_nr < lev );

  cout << endl << "events " << event_nr << endl;

  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignment fits:

  double nb, ne, nm;

  // A:

  nb = hdxAB.GetNbinsX();
  ne = hdxAB.GetSumOfWeights();
  nm = hdxAB.GetMaximum();
  cout << endl << hdxAB.GetTitle() << endl
       << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
       << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

  if( nm/(ne/nb) > 2 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, nm ); // amplitude
    fgp0->SetParameter( 1, hdxAB.GetBinCenter( hdxAB.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 0.05 ); // sigma
    fgp0->SetParameter( 3, hdxAB.GetBinContent(1) ); // BG
    hdxAB.Fit( "fgp0", "q" );
    cout << "  Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "  m " << fgp0->GetParameter(1)
	 << endl << "  s " << fgp0->GetParameter(2)
	 << endl << "  B " << fgp0->GetParameter(3)
	 << endl << "  alignx[0] = " << alignx[0] + fgp0->GetParameter(1)
	 << endl;
    alignx[0] += fgp0->GetParameter(1);

    // turn from profile:

    if( aligniteration ) {

      dxvsyAB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdxvsy = dxvsyAB.GetFunction( "pol1" );
      cout << endl << dxvsyAB.GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;
      fx[0] += fdxvsy->GetParameter(1); // same sign

      dxvsxAB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdxvsx = dxvsxAB.GetFunction( "pol1" );
      cout << endl << dxvsxAB.GetTitle() << " slope " << fdxvsx->GetParameter(1) << endl;
      tx[0] += fdxvsx->GetParameter(1); // same sign

    }

  }

  nb = hdyAB.GetNbinsX();
  ne = hdyAB.GetSumOfWeights();
  nm = hdyAB.GetMaximum();

  cout << endl << hdyAB.GetTitle() << endl
       << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
       << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

  if( nm/(ne/nb) > 2 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, nm ); // amplitude
    fgp0->SetParameter( 1, hdyAB.GetBinCenter( hdyAB.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 0.05 ); // sigma
    fgp0->SetParameter( 3, hdyAB.GetBinContent(1) ); // BG
    hdyAB.Fit( "fgp0", "q" );
    cout << "  Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "  m " << fgp0->GetParameter(1)
	 << endl << "  s " << fgp0->GetParameter(2)
	 << endl << "  B " << fgp0->GetParameter(3)
	 << endl << "  aligny[0] = " << aligny[0] + fgp0->GetParameter(1)
	 << endl;
    aligny[0] += fgp0->GetParameter(1);

    // x-y rotation from profile:

    if( aligniteration ) {

      dyvsxAB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdyvsx = dyvsxAB.GetFunction( "pol1" );
      cout << endl << dyvsxAB.GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;
      fy[0] -= fdyvsx->GetParameter(1); // opposite sign

      dyvsyAB.Fit( "pol1", "q", "", -7.5, 7.5 );
      TF1 * fdyvsy = dyvsyAB.GetFunction( "pol1" );
      cout << endl << dyvsyAB.GetTitle()
	   << " slope " << fdyvsy->GetParameter(1)
	   << endl;
      ty[0] += fdyvsy->GetParameter(1); // same sign

    }

  } // nm

  // C:

  nb = hdxCB.GetNbinsX();
  ne = hdxCB.GetSumOfWeights();
  nm = hdxCB.GetMaximum();

  cout << endl << hdxCB.GetTitle() << endl
       << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
       << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

  if( nm/(ne/nb) > 2 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, nm ); // amplitude
    fgp0->SetParameter( 1, hdxCB.GetBinCenter( hdxCB.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 0.05 ); // sigma
    fgp0->SetParameter( 3, hdxCB.GetBinContent(1) ); // BG
    hdxCB.Fit( "fgp0", "q" );
    cout << "  Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "  m " << fgp0->GetParameter(1)
	 << endl << "  s " << fgp0->GetParameter(2)
	 << endl << "  B " << fgp0->GetParameter(3)
	 << endl << "  alignx[2] = " << alignx[2] + fgp0->GetParameter(1)
	 << endl;
    alignx[2] += fgp0->GetParameter(1);

    // turn from profile:

    if( aligniteration ) {

      dxvsyCB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdxvsy = dxvsyCB.GetFunction( "pol1" );
      cout << endl << dxvsyCB.GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;
      fx[2] += fdxvsy->GetParameter(1); // same sign

      dxvsxCB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdxvsx = dxvsxCB.GetFunction( "pol1" );
      cout << endl << dxvsxCB.GetTitle() << " slope " << fdxvsx->GetParameter(1) << endl;
      tx[2] += fdxvsx->GetParameter(1); // same sign

    }

  } // nm

  nb = hdyCB.GetNbinsX();
  ne = hdyCB.GetSumOfWeights();
  nm = hdyCB.GetMaximum();

  cout << endl << hdyCB.GetTitle() << endl
       << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
       << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

  if( nm/(ne/nb) > 2 ) {

    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, nm ); // amplitude
    fgp0->SetParameter( 1, hdyCB.GetBinCenter( hdyCB.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 0.05 ); // sigma
    fgp0->SetParameter( 3, hdyCB.GetBinContent(1) ); // BG
    hdyCB.Fit( "fgp0", "q" );
    cout << "  Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "  m " << fgp0->GetParameter(1)
	 << endl << "  s " << fgp0->GetParameter(2)
	 << endl << "  B " << fgp0->GetParameter(3)
	 << endl << "  aligny[2] = " << aligny[2] + fgp0->GetParameter(1)
	 << endl;
    aligny[2] += fgp0->GetParameter(1);

    // x-y rotation from profile:

    if( aligniteration ) {

      dyvsxCB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdyvsx = dyvsxCB.GetFunction( "pol1" );
      cout << endl << dyvsxCB.GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;
      fy[2] -= fdyvsx->GetParameter(1); // opposite sign

      dyvsyCB.Fit( "pol1", "q", "", -7.5, 7.5 );
      TF1 * fdyvsy = dyvsyCB.GetFunction( "pol1" );
      cout << endl << dyvsyCB.GetTitle()
	   << " slope " << fdyvsy->GetParameter(1)
	   << endl;
      ty[2] += fdyvsy->GetParameter(1); // same sign

    }

  } // nm

  // D:

  nb = hdxDB.GetNbinsX();
  ne = hdxDB.GetSumOfWeights();
  nm = hdxDB.GetMaximum();

  cout << endl << hdxDB.GetTitle() << endl
       << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
       << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

  if( nm/(ne/nb) > 2 ) {
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, nm ); // amplitude
    fgp0->SetParameter( 1, hdxDB.GetBinCenter( hdxDB.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 0.05 ); // sigma
    fgp0->SetParameter( 3, hdxDB.GetBinContent(1) ); // BG
    hdxDB.Fit( "fgp0", "q" );
    cout << "  Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "  m " << fgp0->GetParameter(1)
	 << endl << "  s " << fgp0->GetParameter(2)
	 << endl << "  B " << fgp0->GetParameter(3)
	 << endl << "  alignx[3] = " << alignx[3] + fgp0->GetParameter(1)
	 << endl;
    alignx[3] += fgp0->GetParameter(1);

    // turn from profile:

    if( aligniteration ) {

      dxvsyDB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdxvsy = dxvsyDB.GetFunction( "pol1" );
      cout << endl << dxvsyDB.GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;
      fx[3] += fdxvsy->GetParameter(1); // same sign

      dxvsxDB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdxvsx = dxvsxDB.GetFunction( "pol1" );
      cout << endl << dxvsxDB.GetTitle() << " slope " << fdxvsx->GetParameter(1) << endl;
      tx[3] += fdxvsx->GetParameter(1); // same sign

    }

  } // nm

  nb = hdyDB.GetNbinsX();
  ne = hdyDB.GetSumOfWeights();
  nm = hdyDB.GetMaximum();

  cout << endl << hdyDB.GetTitle() << endl
       << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
       << "  Maximum " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

  if( nm/(ne/nb) > 2 ) {
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, nm ); // amplitude
    fgp0->SetParameter( 1, hdyDB.GetBinCenter( hdyDB.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 0.05 ); // sigma
    fgp0->SetParameter( 3, hdyDB.GetBinContent(1) ); // BG
    hdyDB.Fit( "fgp0", "q" );
    cout << "  Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "  m " << fgp0->GetParameter(1)
	 << endl << "  s " << fgp0->GetParameter(2)
	 << endl << "  B " << fgp0->GetParameter(3)
	 << endl << "  aligny[3] = " << aligny[3] + fgp0->GetParameter(1)
	 << endl;
    aligny[3] += fgp0->GetParameter(1);

    // x-y rotation from profile:

    if( aligniteration ) {

      dyvsxDB.Fit( "pol1", "q", "", -25, 25 );
      TF1 * fdyvsx = dyvsxDB.GetFunction( "pol1" );
      cout << endl << dyvsxDB.GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;
      fy[3] -= fdyvsx->GetParameter(1); // opposite sign

      dyvsyDB.Fit( "pol1", "q", "", -7.5, 7.5 );
      TF1 * fdyvsy = dyvsyDB.GetFunction( "pol1" );
      cout << endl << dyvsyDB.GetTitle()
	   << " slope " << fdyvsy->GetParameter(1)
	   << endl;
      ty[3] += fdyvsy->GetParameter(1); // same sign

    }

  } // nm

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write alignment to file:

  ofstream alignFile( alignFileName.str() );

  alignFile << "# modules alignment for run " << run << endl;
  ++aligniteration;
  alignFile << "iteration " << aligniteration << endl;

  cout << endl << "align iteration " << aligniteration << endl;

  for( int ipl = 0; ipl < 4; ++ipl ) {
    alignFile << endl;
    alignFile << "plane " << ipl << endl;
    alignFile << "alignx " << alignx[ipl] << endl;
    alignFile << "aligny " << aligny[ipl] << endl;
    alignFile << "fx " << fx[ipl] << endl;
    alignFile << "fy " << fy[ipl] << endl;
    alignFile << "tx " << tx[ipl] << endl;
    alignFile << "ty " << ty[ipl] << endl;

    cout << endl;
    cout << "plane " << ipl << endl;
    cout << "alignx " << alignx[ipl] << endl;
    cout << "aligny " << aligny[ipl] << endl;
    cout << "fx " << fx[ipl] << endl;
    cout << "fy " << fy[ipl] << endl;
    cout << "tx " << tx[ipl] << endl;
    cout << "ty " << ty[ipl] << endl;
  } // ipl

  alignFile.close();

  cout << endl
       << "written to " << alignFileName.str()
       << endl;

  if( aligniteration == 1 )
    cout << "need one more align iteration: please run again!" << endl;

  cout << endl;
  cout << "quad  tracks " << n4 << endl;
  cout << "mille tracks " << nmille << endl;


  cout << "Efficiencies:" << endl;
  cout << "Mod\tTotal\t\tupper\t\tlower" << setprecision(6) << endl;
  cout << "A\t" << effAvsw.GetBinContent(14) << "\t" << effAvsx1.GetMean(2) << "\t" << effAvsx0.GetMean(2) << endl;
  cout << "B\t" << effBvsw.GetBinContent(14) << "\t" << effBvsx1.GetMean(2) << "\t" << effBvsx0.GetMean(2) << endl;
  cout << "C\t" << effCvsw.GetBinContent(14) << "\t" << effCvsx1.GetMean(2) << "\t" << effCvsx0.GetMean(2) << endl;
  cout << "D\t" << effDvsw.GetBinContent(14) << "\t" << effDvsx1.GetMean(2) << "\t" << effDvsx0.GetMean(2) << endl;
  

  cout << endl << histoFile->GetName() << endl << endl;

  return 0;
}
