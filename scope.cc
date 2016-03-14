
// Daniel Pitzl, Feb 2016
// telescope analysis with eudaq and DUT

// make scope
// scope -l 9999 20833  # tilted DUT 504
// scope 19037

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TFile.h>
#include <TH1I.h> // counting
#include <TH1D.h> // weighted counts
#include <TH2I.h>
#include <TProfile.h>
#include <TProfile2D.h>
#include <TF1.h>

#include <sstream> // stringstream
#include <fstream> // filestream

using namespace std;
using namespace eudaq;

struct pixel {
  int col;
  int row;
  int adc;
  double q;
  int ord;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  double charge;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  bool lk;
};

// globals:

pixel pb[999]; // global declaration: vector of pixels with hit
int fNHit; // global

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
    do {
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

    c.size = c.vpix.size();
    c.col = 0;
    c.row = 0;
    double sumQ = 0;
    int minx = 999;
    int maxx = 0;
    int miny = 999;
    int maxy = 0;

    for( vector<pixel>::iterator p = c.vpix.begin();  p != c.vpix.end();  ++p ) {
      double Qpix = p->q; // calibrated [Vcal]
      if( Qpix < 0 ) Qpix = 1; // DP 1.7.2012
      c.charge += Qpix;
      sumQ += Qpix;
      c.col += (*p).col*Qpix;
      c.row += (*p).row*Qpix;
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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // further arguments:

  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

  } // argc

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double DUTtilt0 = 19.3;
  double pbeam = 5.6;
  int chip0 = 504;
  int refchip0 = 501;
  string gainFileName( "gain.dat" );
  string refgainFileName( "gain.dat" );
  int weib = 3;

  ifstream runsFile( "runs.dat" );

  if( runsFile.bad() || ! runsFile.is_open() ) {
    cout << "Error opening runs.dat" << endl;
    return 1;
  }
  // can there be instructions between if and else ? no

  else {

    cout << "read runs from runs.dat" << endl;

    string hash( "#" );
    string RUN( "run" );
    string GEO( "geo" );
    string GeV( "GeV" );
    string CHIP( "chip" );
    string WEIB( "weib" );
    string GAIN( "gain" );
    string REFCHIP( "refchip" );
    string REFGAIN( "refgain" );
    string TILT( "tilt" );
    bool found = 0;

    while( ! runsFile.eof() ) {

      string line;
      getline( runsFile, line );

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == RUN )  {
	int ival;
	tokenizer >> ival;
	if( ival == run ) {
	  found = 1;
	  break; // end file reading
	}
      }

      if( tag == TILT ) {
	tokenizer >> DUTtilt0;
	continue;
      }

      if( tag == GAIN ) {
	tokenizer >> gainFileName;
	continue;
      }

      if( tag == REFGAIN ) {
	tokenizer >> refgainFileName;
	continue;
      }

      if( tag == GEO ) {
	tokenizer >> geoFileName;
	continue;
      }

      if( tag == GeV ) {
	tokenizer >> pbeam;
	continue;
      }

      if( tag == CHIP ) {
	tokenizer >> chip0;
	continue;
      }

      if( tag == REFCHIP ) {
	tokenizer >> refchip0;
	continue;
      }

      if( tag == weib ) {
	tokenizer >> weib;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    if( found )
      cout 
	<< "settings for run " << run << ":" << endl
	<< "  beam " << pbeam << " GeV" << endl
	<< "  geo file " << geoFileName << endl
	<< "  nominal DUT tilt " << DUTtilt0 << " deg" << endl
	<< "  DUT chip " << chip0 << endl
	<< "  DUT gain file " << gainFileName << endl
	<< "  Weibull version " << weib << endl
	<< "  REF chip " << refchip0 << endl
	<< "  REF gain file " << refgainFileName << endl
	<< endl;
    else {
      cout << "run " << run << " not found in runs.dat" << endl;
      return 1;
    }

  } // alignFile

  runsFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // geometry:

  int nx[9]; // x-pixels per plane
  int ny[9]; // y-pixels per plane
  double sizex[9]; // x size per plane
  double sizey[9]; // y size per plane
  double ptchx[9]; // x-pixel size
  double ptchy[9]; // y-pixel size
  double midx[9]; // x mid
  double midy[9]; // y mid

  double zz[9];

  for( int ipl = 0; ipl < 9; ++ipl )
    nx[ipl] = 0; // missing plane flag

  ifstream geoFile( geoFileName );

  cout << endl;

  if( geoFile.bad() || ! geoFile.is_open() ) {
    cout << "Error opening " << geoFileName << endl;
    return 1;
  }

  cout << "read geometry from " << geoFileName << endl;

  { // open local scope

    string hash( "#" );
    string plane( "plane" );
    string type( "type" );
    string sizexs( "sizex" );
    string sizeys( "sizey" );
    string npixelx( "npixelx" );
    string npixely( "npixely" );
    string zpos( "zpos" );

    int ipl = 0;
    string chiptype;

    while( ! geoFile.eof() ) {

      string line;
      getline( geoFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane ) {
	tokenizer >> ipl;
	continue;
      }

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == type ) {
	tokenizer >> chiptype;
	continue;
      }

      if( tag == sizexs ) {
	double val;
	tokenizer >> val;
	sizex[ipl] = val;
	continue;
      }

      if( tag == sizeys ) {
	double val;
	tokenizer >> val;
	sizey[ipl] = val;
	continue;
      }

      if( tag == npixelx ) {
	int val;
	tokenizer >> val;
	nx[ipl] = val;
	continue;
      }

      if( tag == npixely ) {
	int val;
	tokenizer >> val;
	ny[ipl] = val;
	continue;
      }

      if( tag == zpos ) {
	double val;
	tokenizer >> val;
	zz[ipl] = val;
	continue;
      }

      // anything else on the line and in the file gets ignored

    } // while getline

    for( int ipl = 0; ipl < 9; ++ipl ) {
      if( nx[ipl] == 0 ) continue; // missing plane flag
      ptchx[ipl] = sizex[ipl] / nx[ipl]; // pixel size
      ptchy[ipl] = sizey[ipl] / ny[ipl];
      midx[ipl] = 0.5 * sizex[ipl]; // mid plane
      midy[ipl] = 0.5 * sizey[ipl]; // mid plane
    }

  } // geo scope

  geoFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // alignments:

  int aligniteration = 0;
  double alignx[9];
  double aligny[9];
  double rotx[9];
  double roty[9];

  ostringstream alignFileName; // output string stream

  alignFileName << "align_" << run << ".dat";

  ifstream ialignFile( alignFileName.str() );

  cout << endl;

  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    cout << "Error opening " << alignFileName.str() << endl
	 << "  please do: tele -g " << geoFileName << " " << run << endl
	 << endl;
    return 1;
  }
  // can there be instructions between if and else ?
  else {

    cout << "read alignment from " << alignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string plane( "plane" );
    string shiftx( "shiftx" );
    string shifty( "shifty" );
    string rotxvsy( "rotxvsy" );
    string rotyvsx( "rotyvsx" );

    int ipl = 0;

    while( ! ialignFile.eof() ) {

      string line;
      getline( ialignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> aligniteration;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 9 ) {
	cout << "wrong plane number " << ipl << endl;
	continue;
      }

      double val;
      tokenizer >> val;
      if(      tag == shiftx )
	alignx[ipl] = val;
      else if( tag == shifty )
	aligny[ipl] = val;
      else if( tag == rotxvsy )
	rotx[ipl] = val;
      else if( tag == rotyvsx )
	roty[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  double wt = atan(1.0) / 45.0; // pi/180 deg

  double thr = 1.5; // [ke]
  double qwid = 3.5; // [ke] for Moyal

  bool rot90 = 0; // 504
  //bool rot90 = 1; // 506

  int iDUT = 7;

  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0.0;
  double DUTturn = 0;
  double DUTtilt = 19.3; // [deg]
  double DUTz = 40 + zz[2];

  ostringstream DUTalignFileName; // output string stream

  DUTalignFileName << "alignDUT_" << run << ".dat";

  ifstream iDUTalignFile( DUTalignFileName.str() );

  cout << endl;

  if( iDUTalignFile.bad() || ! iDUTalignFile.is_open() ) {
    cout << "no " << DUTalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read DUTalignment from " << DUTalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iDUTalignFile.eof() ) {

      string line;
      getline( iDUTalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> DUTaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	DUTalignx = val;
      else if( tag == aligny )
	DUTaligny = val;
      else if( tag == rot )
	DUTrot = val;
      else if( tag == tilt )
	DUTtilt = val;
      else if( tag == turn )
	DUTturn = val;
      else if( tag == dz )
	DUTz = val + zz[2];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iDUTalignFile.close();

  if( DUTaligniteration == 0 )
    DUTtilt = DUTtilt0; // from runs.dat

  // normal vector on DUT surface:
  // N = ( 0, 0, -1 ) on DUT, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double co = cos( DUTturn*wt );
  const double so = sin( DUTturn*wt );
  const double ca = cos( DUTtilt*wt );
  const double sa = sin( DUTtilt*wt );
  const double cf = cos( DUTrot );
  const double sf = sin( DUTrot );

  const double Nx =-ca*so;
  const double Ny = sa;
  const double Nz =-ca*co;

  const double norm = cos( DUTturn*wt ) * cos( DUTtilt*wt ); // length of Nz

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // REF:

  int iREF = 8;

  int REFaligniteration = 0;
  double REFalignx = 0.0;
  double REFaligny = 0.0;
  double REFrot = 0.0;
  double REFz = 55 + zz[5];

  ostringstream REFalignFileName; // output string stream

  REFalignFileName << "alignREF_" << run << ".dat";

  ifstream iREFalignFile( REFalignFileName.str() );

  cout << endl;

  if( iREFalignFile.bad() || ! iREFalignFile.is_open() ) {
    cout << "no " << REFalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read REFalignment from " << REFalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string dz( "dz" );

    while( ! iREFalignFile.eof() ) {

      string line;
      getline( iREFalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> REFaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	REFalignx = val;
      else if( tag == aligny )
	REFaligny = val;
      else if( tag == rot )
	REFrot = val;
      else if( tag == dz )
	REFz = val + zz[5];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iREFalignFile.close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT gain:

  double p0[52][80]; // Weibull
  double p1[52][80];
  double p2[52][80];
  double p3[52][80];
  double p4[52][80];
  double p5[52][80];

  ifstream gainFile( gainFileName );

  if(! gainFile ) {
    cout << "gain file " << gainFileName << " not found" << endl;
    return 1;
  }
  else {
    cout << endl << "using DUT pixel gain file " << gainFileName << endl;
    char ih[99];
    int col;
    int row;
    while( gainFile >> ih ) {
      gainFile >> col;
      gainFile >> row;
      if( col < 0 || col > 51 || row < 0 || row > 79 ) {
	cout << "invalid pixel in gain file " << col << " " << row << endl;
	continue;
      }
      gainFile >> p0[col][row];
      gainFile >> p1[col][row];
      gainFile >> p2[col][row];
      gainFile >> p3[col][row];
      gainFile >> p4[col][row];
      gainFile >> p5[col][row]; // gain ratio
    } // while

  } // gainFile

  double ke = 0.250; // to get q0f peak at 22 ke

  if( chip0 == 500 )
    ke = 0.305; // 14393 to get q0f peak at 22 ke no eps in Q

  if( chip0 == 504 ) {
    ke = 0.254; // 14614 to get q0f peak at 22 ke
    if( run >= 19019 ) // Apr 2015
      ke = 0.235; // 19045 to get q0f peak at 22 ke
    if( run >= 20785 ) // Jul 2015
      ke = 0.235;
    if( run >= 20823 ) // 12.7.2015 chiller 17
      ke = 0.238; // to get q0f peak at 22 ke
  }

  if( chip0 == 506 ) {
    ke = 0.290; // 14654 to get q0f peak at 22 ke
    if( run >= 19441 ) // chiller off
      ke = 0.263;
    if( run >= 19582 ) // chiller off
      ke = 0.268;
    if( run >= 19702 ) // chiller on
      ke = 0.257;
    if( run >= 20744 ) // chiller off
      ke = 0.268;
  }

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // REF gain:

  double r0[52][80]; // Weibull
  double r1[52][80];
  double r2[52][80];
  double r3[52][80];
  double r4[52][80];
  double r5[52][80];

  ifstream refgainFile( refgainFileName );

  if(! refgainFile ) {
    cout << "REF gain file " << refgainFileName << " not found" << endl;
    return 1;
  }
  else {
    cout << endl << "using REF pixel gain file " << refgainFileName << endl;
    char ih[99];
    int col;
    int row;
    while( refgainFile >> ih ) {
      refgainFile >> col;
      refgainFile >> row;
      if( col < 0 || col > 51 || row < 0 || row > 79 ) {
	cout << "invalid pixel in gain file " << col << " " << row << endl;
	continue;
      }
      refgainFile >> r0[col][row];
      refgainFile >> r1[col][row];
      refgainFile >> r2[col][row];
      refgainFile >> r3[col][row];
      refgainFile >> r4[col][row];
      refgainFile >> r5[col][row]; // gain ratio
    } // while

  } // gainFile

  double refke = 0.250; // to get q0f peak at 22 ke

  if( refchip0 == 500 )
    refke = 0.305;

  else if( refchip0 == 501 )
    refke = 0.305;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // (re-)create root file:

  ostringstream rootFileName; // output string stream

  rootFileName << "scope" << run << ".root";

  TFile* histoFile = new TFile( rootFileName.str(  ).c_str(  ), "RECREATE" );

  // book histos:
 
  TH1I hdtus = TH1I( "dtus", "time between events;time between events [us];events", 100, 0, 1000 );
  TH1I hdtms = TH1I( "dtms", "time between events;time between events [ms];events", 100, 0, 1000 );

  TH1I hcol[9];
  TH1I hrow[9];
  TH1I hnpx[9];

  TH1I hncl[9];
  TH1I hsiz[9];
  TH1I hncol[9];
  TH1I hnrow[9];

  for( int ipl = 0; ipl < 9; ++ipl ) {

    hcol[ipl] = TH1I( Form( "col%i", ipl ),
		      Form( "%i col;col;%i pixels", ipl, ipl ), 
		      max( 52, nx[ipl]/4 ), 0, max( 52, nx[ipl]/4 ) );
    hrow[ipl] = TH1I( Form( "row%i", ipl ),
		      Form( "%i row;row;%i pixels", ipl, ipl ),
		      max( 80, ny[ipl]/2 ), 0, max( 80, ny[ipl]/2 ) );

    hnpx[ipl] = TH1I( Form( "npx%i", ipl ),
		      Form( "%i pixel per event;pixels;%i events", ipl, ipl ),
		      200, 0, 200 );

    hncl[ipl] = TH1I( Form( "ncl%i", ipl ),
		      Form( "plane %i cluster per event;cluster;plane %i events", ipl, ipl ),
		      51, -0.5, 50.5 );
    hsiz[ipl] = TH1I( Form( "clsz%i", ipl ),
		      Form( "%i cluster size;pixels/cluster;%i clusters", ipl, ipl ),
		      51, -0.5, 50.5 );
    hncol[ipl] = TH1I( Form( "ncol%i", ipl ), 
		       Form( "%i cluster size x;columns/cluster;%i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );
    hnrow[ipl] = TH1I( Form( "nrow%i", ipl ),
		       Form( "%i cluster size y;rows/cluster;%i clusters", ipl, ipl ),
		       21, -0.5, 20.5 );

  } // planes

  // driplets:

  TH1I hdx35 = TH1I( "dx35", "3-5 dx;3-5 dx [mm];cluster pairs", 100, -1, 1 );
  TH1I hdy35 = TH1I( "dy35", "3-5 dy;3-5 dy [mm];cluster pairs", 100, -1, 1 );

  TH1I hdxdri = TH1I( "dxdri", "driplet dx;driplet dx [mm];driplets", 100, -0.1, 0.1 );
  TH1I hdydri = TH1I( "dydri", "driplet dy;driplet dy [mm];driplets", 100, -0.1, 0.1 );

  TH1I hdxdric = TH1I( "dxdric", "driplet dx;driplet dx [mm];driplets", 100, -0.05, 0.05 );
  TH1I hdydric = TH1I( "dydric", "driplet dy;driplet dy [mm];driplets", 100, -0.05, 0.05 );

  TH1I ndriHisto = TH1I( "ndri", "driplets;driplets;events", 51, -0.5, 50.5 );

  // REF vs driplets:

  TH1I refsxaHisto = TH1I( "refsxa",
			   "REF + driplet x;REF cluster + driplet #Sigmax [mm];REF clusters",
			   440, -11, 11 );
  TH1I refdxaHisto = TH1I( "refdxa",
			   "REF - driplet x;REF cluster - driplet #Deltax [mm];REF clusters",
			   440, -11, 11 );

  TH1I refsyaHisto = TH1I( "refsya",
			   "REF + driplet y;REF cluster + driplet #Sigmay [mm];REF clusters",
			   220, -5.5, 5.5 );
  TH1I refdyaHisto = TH1I( "refdya",
			   "REF - driplet y;REF cluster - driplet #Deltay [mm];REF clusters",
			   220, -5.5, 5.5 );

  TH1I refdxHisto = TH1I( "refdx",
			   "REF - driplet x;REF cluster - driplet #Deltax [mm];REF clusters",
			   200, -1, 1 );
  TH1I refdyHisto = TH1I( "refdy",
			   "REF - driplet y;REF cluster - driplet #Deltay [mm];REF clusters",
			   200, -1, 1 );

  TProfile refdyvsx = TProfile( "refdyvsx",
			   "REF #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
				76, -3.8, 3.8, -0.5, 0.5 );
  TProfile refdyvsy = TProfile( "refdyvsy",
			   "REF #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
				76, -3.8, 3.8, -0.5, 0.5 );
  TProfile refdyvsty =
    TProfile( "refdyvsty",
	      "REF #Deltay vs #theta_{y};y track slope [rad];<cluster - triplet #Deltay> [mm]",
	      80, -0.002, 0.002, -0.5, 0.5 );

  TH1I refnpxHisto = TH1I( "refnpx",
			   "REF linked clusters;REF cluster size [pixels];linked REF cluster",
			   20, 0.5, 20.5 );

  TH1I refqHisto = TH1I( "refq",
			 "REF linked clusters;REF cluster charge [ke];linked REF cluster",
			 80, 0, 80 );

  TProfile2D refnpxvsxmym =
    TProfile2D( "refnpxvsxmym",
	      "REF cluster size vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];REF <cluster size> [pixels]",
		60, 0, 0.3, 40, 0, 0.2, 0, 20 );

  TProfile2D refqxvsxmym =
    TProfile2D( "refqxvsxmym",
	      "REF cluster charge vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];REF <cluster charge> [ke]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.1 );

  TH1I reflkxHisto = TH1I( "reflkx",
			   "linked driplet at REF x;driplet x at REF [mm];linked driplets",
			   220, -11, 11 );
  TH1I reflkyHisto = TH1I( "reflky",
			   "linked driplet at REF y;driplet y at REF [mm];linked driplets",
			   110, -5.5, 5.5 );

  TH1I reflkcolHisto = TH1I( "reflkcol",
			     "REF linked col;REF linked col;linked REF cluster",
			     52, 0, 52 );
  TH1I reflkrowHisto = TH1I( "reflkrow",
			     "REF linked row;REF linked row;linked REF cluster",
			     80, 0, 80 );

  // DUT:

  TH1I dutpxq1stHisto =
    TH1I( "dutpxq1st",
	  "DUT pixel charge 1st;1st pixel charge [ke];1st pixels",
	  100, 0, 25 );

  TH1I dutpxq2ndHisto =
    TH1I( "dutpxq2nd",
	  "DUT pixel charge 2nd;2nd pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd02Histo =
    TH1I( "dutpxq2nd02",
	  "DUT pixel charge 2nd - 0.02 q1;2nd - 0.02 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd03Histo =
    TH1I( "dutpxq2nd03",
	  "DUT pixel charge 2nd - 0.03 q1;2nd - 0.03 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd04Histo =
    TH1I( "dutpxq2nd04",
	  "DUT pixel charge 2nd - 0.04 q1;2nd - 0.04 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd05Histo =
    TH1I( "dutpxq2nd05",
	  "DUT pixel charge 2nd - 0.05 q1;2nd - 0.05 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd06Histo =
    TH1I( "dutpxq2nd06",
	  "DUT pixel charge 2nd - 0.06 q1;2nd - 0.06 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd07Histo =
    TH1I( "dutpxq2nd07",
	  "DUT pixel charge 2nd - 0.07 q1;2nd - 0.07 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd08Histo =
    TH1I( "dutpxq2nd08",
	  "DUT pixel charge 2nd - 0.08 q1;2nd - 0.08 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );
  TH1I dutpxq2nd09Histo =
    TH1I( "dutpxq2nd09",
	  "DUT pixel charge 2nd - 0.09 q1;2nd - 0.09 q1 pixel charge [ke];2nd pixels",
	  100, 0, 25 );

  TH1I dutcolHisto =
    TH1I( "dutcol",
	  "DUT pixel column;pixel column;pixels",
	  52, -0.5, 51.5 );
  TH1I dutrowHisto =
    TH1I( "dutrow",
	  "DUT pixel rowumn;pixel rowumn;pixels",
	  80, -0.5, 79.5 );

  TH1I dutq0fHisto =
    TH1I( "dutq0f",
	  "normal fiducial cluster charge;normal fiducial cluster charge [ke];fiducial clusters",
	  160, 0, 80 );

  TH1I dutncolHisto =
    TH1I( "dutncol",
	  "DUT cluster size;cluster size [columns];clusters",
	  52, 0.5, 52.5 );
  TH1I dutnrowHisto =
    TH1I( "dutnrow",
	  "DUT cluster size;cluster size [rows];clusters",
	  20, 0.5, 20.5 );
  TH1I dutcolminHisto =
    TH1I( "dutcolmin",
	  "DUT first cluster column;first cluster column;clusters",
	  52, -0.5, 51.5 );
  TH1I dutcolmaxHisto =
    TH1I( "dutcolmax",
	  "DUT last cluster column;last cluster column;clusters",
	  52, -0.5, 51.5 );
  TH1I dutcol0qHisto =
    TH1I( "dutcol0q",
	  "DUT first column charge;first column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol0oddqHisto =
    TH1I( "dutcol0oddq",
	  "DUT odd first column charge;odd first column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol0eveqHisto =
    TH1I( "dutcol0eveq",
	  "DUT eve first column charge;eve first column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol9qHisto =
    TH1I( "dutcol9q",
	  "DUT last column charge;last column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol1qHisto =
    TH1I( "dutcol1q",
	  "DUT 2nd column charge;2nd column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol2qHisto =
    TH1I( "dutcol2q",
	  "DUT 3rd column charge;3rd column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol3qHisto =
    TH1I( "dutcol3q",
	  "DUT 4th column charge;4th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol4qHisto =
    TH1I( "dutcol4q",
	  "DUT 5th column charge;5th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol5qHisto =
    TH1I( "dutcol5q",
	  "DUT 6th column charge;6th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol6qHisto =
    TH1I( "dutcol6q",
	  "DUT 7th column charge;7th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol7qHisto =
    TH1I( "dutcol7q",
	  "DUT 8th column charge;8th column charge [ke];clusters",
	  100, 0, 50 );
  TH1I dutcol8qHisto =
    TH1I( "dutcol8q",
	  "DUT 9th column charge;9th column charge [ke];clusters",
	  100, 0, 50 );

  // triplets:

  TH1I hdx02 = TH1I( "dx02", "0-2 dx;0-2 dx [mm];cluster pairs", 100, -1, 1 );
  TH1I hdy02 = TH1I( "dy02", "0-2 dy;0-2 dy [mm];cluster pairs", 100, -1, 1 );

  TH1I hdxtri = TH1I( "dxtri", "triplet dx;triplet dx [mm];triplets", 100, -0.1, 0.1 );
  TH1I hdytri = TH1I( "dytri", "triplet dy;triplet dy [mm];triplets", 100, -0.1, 0.1 );

  TH1I hdxtric = TH1I( "dxtric", "triplet dx;triplet dx [mm];triplets", 100, -0.05, 0.05 );
  TH1I hdytric = TH1I( "dytric", "triplet dy;triplet dy [mm];triplets", 100, -0.05, 0.05 );

  TH1I hdxtric1 = TH1I( "dxtric1", "triplet dx 1-col;1-col triplet dx [mm];1-col triplets",
			100, -0.05, 0.05 );
  TH1I hdxtric2 = TH1I( "dxtric2", "triplet dx 2-col;2-col triplet dx [mm];2-col triplets",
			100, -0.05, 0.05 );
  TH1I hdxtric3 = TH1I( "dxtric3", "triplet dx 3-col;3-col triplet dx [mm];3-col triplets",
			100, -0.05, 0.05 );
  TH1I hdxtric4 = TH1I( "dxtric4", "triplet dx 4-col;4-col triplet dx [mm];4-col triplets",
			100, -0.05, 0.05 );
  TH1I hdxtric5 = TH1I( "dxtric5", "triplet dx 5-col;5-col triplet dx [mm];5-col triplets",
			100, -0.05, 0.05 );

  TH1I hdxtris1 = TH1I( "dxtris1", "triplet dx 1-px;1-px triplet dx [mm];1-px triplets",
			100, -0.05, 0.05 );
  TH1I hdxtris2 = TH1I( "dxtris2", "triplet dx 2-px;2-px triplet dx [mm];2-px triplets",
			100, -0.05, 0.05 );
  TH1I hdxtris3 = TH1I( "dxtris3", "triplet dx 3-px;3-px triplet dx [mm];3-px triplets",
			100, -0.05, 0.05 );
  TH1I hdxtris4 = TH1I( "dxtris4", "triplet dx 4-px;4-px triplet dx [mm];4-px triplets",
			100, -0.05, 0.05 );
  TH1I hdxtris5 = TH1I( "dxtris5", "triplet dx 5-px;5-px triplet dx [mm];5-px triplets",
			100, -0.05, 0.05 );

  TH1I ntriHisto = TH1I( "ntri", "triplets;triplets;events", 51, -0.5, 50.5 );

  TH1I trixHisto = TH1I( "trix", "triplets at z DUT;x [mm];triplets",
			  240, -12, 12 );
  TH1I triyHisto = TH1I( "triy", "triplets at z DUT;y [mm];triplets",
			  120, -6, 6 );
  TH2I trixyHisto = TH2I( "trixy", "triplets at z DUT;x [mm];y [mm];triplets",
			  240, -12, 12, 120, -6, 6 );

  TH1I ttdxHisto = TH1I( "ttdx", "telescope triplets;triplets #Deltax [mm];triplet pairs",
			 100, -5, 5 );
  TH1I ttdx1Histo = TH1I( "ttdx1", "telescope triplets;triplets #Deltax [mm];triplet pairs",
			  100, -0.5, 0.5 );
  TH1I ttdminHisto = TH1I( "ttdmin",
			   "telescope triplets isolation;triplets #Delta [mm];triplet pairs",
			   100, 0, 1 );

  // dripets - triplets:

  TH1I hsixdx = TH1I( "sixdx", "six dx;dx [mm];triplet-driplet pairs", 100, -5, 5 );
  TH1I hsixdy = TH1I( "sixdy", "six dy;dy [mm];triplet-driplet pairs", 100, -5, 5 );
  TH1I hsixdxc = TH1I( "sixdxc", "six dx;dx [mm];triplet-driplet pairs", 100, -1, 1 );
  TH1I hsixdyc = TH1I( "sixdyc", "six dy;dy [mm];triplet-driplet pairs", 100, -1, 1 );

  TProfile sixdxvsy = TProfile( "sixdxvsy",
				"six #Deltax vs y;yB [mm];<driplet - triplet #Deltax> [mm]",
				100, -5, 5, -0.5, 0.5 );
  TProfile sixdyvsx = TProfile( "sixdyvsx",
				"six #Deltay vs x;xB [mm];<driplet - triplet #Deltay> [mm]",
				200, -10, 10, -0.5, 0.5 );

  TProfile sixdxvstx = TProfile( "sixdxvstx",
				"six #Deltax vs slope x;slope x [rad];<driplet - triplet #Deltax> [mm]",
				100, -0.002, 0.002, -0.5, 0.5 );
  TProfile sixdyvsty = TProfile( "sixdyvsty",
				"six #Deltay vs slope y;slope y [rad];<driplet - triplet #Deltay> [mm]",
				100, -0.002, 0.002, -0.5, 0.5 );

  TH1I hsixdslpx = TH1I( "sixdslpx",
			 "driplet slope x - triplet slope x;driplet slope x - triplet slope x;driplet-triplet pairs",
			 100, -0.005, 0.005 );     
  TH1I hsixdslpy = TH1I( "sixdslpy",
			 "driplet slope y - triplet slope y;driplet slope y - triplet slope y;driplet-triplet pairs",
			 100, -0.005, 0.005 );     

  TH2I sixxyHisto = TH2I( "sixxy", "sixplets at z DUT;x [mm];y [mm];sixplets",
			  240, -12, 12, 120, -6, 6 );

  // DUT pixel vs triplets:

  TH1I z3Histo = TH1I( "z3",
		       "z3 should be zero;z3 [mm];triplets",
		       100, -0.01, 0.01 );

  TH1I cmssxaHisto = TH1I( "cmssxa",
			   "Pixel + Telescope x;cluster + triplet #Sigmax [mm];clusters",
			   440, -11, 11 );
  TH1I cmsdxaHisto = TH1I( "cmsdxa",
			   "Pixel - Telescope x;cluster - triplet #Deltax [mm];clusters",
			   440, -11, 11 );

  TH1I cmssyaHisto = TH1I( "cmssya",
			   "Pixel + Telescope y;cluster + triplet #Sigmay [mm];clusters",
			   220, -5.5, 5.5 );
  TH1I cmsdyaHisto = TH1I( "cmsdya",
			   "Pixel - Telescope y;cluster - triplet #Deltay [mm];clusters",
			   220, -5.5, 5.5 );

  TH1I cmsdxHisto = TH1I( "cmsdx",
			   "Pixel - Telescope x;cluster - triplet #Deltax [mm];clusters",
			   200, -0.5, 0.5 );
  TH1I cmsdyHisto = TH1I( "cmsdy",
			   "Pixel - Telescope y;cluster - triplet #Deltay [mm];clusters",
			   500, -0.5, 0.5 );

  TH1I cmslkxHisto = TH1I( "cmslkx",
			   "linked triplet at DUT x;triplet x at DUT [mm];linked triplets",
			   220, -11, 11 );
  TH1I cmslkyHisto = TH1I( "cmslky",
			   "linked triplet at DUT y;triplet y at DUT [mm];linked triplets",
			   110, -5.5, 5.5 );

  TH1I cmscolHisto = TH1I( "cmscol",
			   "DUT linked columns;DUT linked cluster column;linked clusters",
			   52, 0, 52 );
  TH1I cmsrowHisto = TH1I( "cmsrow",
			   "DUT linked rows;DUT linked cluster row;linked clusters",
			   80, 0, 80 );

  TH1I cmsdxfcHisto =
    TH1I( "cmsdxfc",
	  "fiducial #Deltax cut y;cluster - triplet #Deltax [mm];fiducial clusters",
	  200, -0.5, 0.5 );

  TProfile cmsdxvsx = TProfile( "cmsdxvsx",
			   "#Deltax vs x;x track [mm];<cluster - triplet #Deltax> [mm]",
				80, -4, 4, -0.5, 0.5 );
  TProfile cmsdxvsy = TProfile( "cmsdxvsy",
			   "#Deltax vs y;y track [mm];<cluster - triplet #Deltax> [mm]",
				80, -4, 4, -0.5, 0.5 );
  TProfile cmsdxvstx =
    TProfile( "cmsdxvstx",
	      "#Deltax vs #theta_{x};x track slope [rad];<cluster - triplet #Deltax> [mm]",
	      80, -0.004, 0.004, -0.5, 0.5 );

  TH1I cmsdyfcHisto =
    TH1I( "cmsdyfc",
	  "fiducial #Deltay cut x;cluster - triplet #Deltay [mm];fiducial clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq0Histo =
    TH1I( "cmsdyfcq0",
	  "fiducial #Deltay cut x q0;cluster - triplet #Deltay [mm];fiducial q0 clusters",
	  200, -0.5, 0.5 );
  TH1I cmsdyfcq9Histo =
    TH1I( "cmsdyfcq9",
	  "fiducial #Deltay cut x q9;cluster - triplet #Deltay [mm];fiducial q9 clusters",
	  200, -0.5, 0.5 );
  TH1I cmsdyfcq1Histo =
    TH1I( "cmsdyfcq1",
	  "fiducial #Deltay cut x q1;cluster - triplet #Deltay [mm];fiducial q1 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq2Histo =
    TH1I( "cmsdyfcq2",
	  "fiducial #Deltay cut x q2;cluster - triplet #Deltay [mm];fiducial q2 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3Histo =
    TH1I( "cmsdyfcq3",
	  "fiducial #Deltay cut x q3;cluster - triplet #Deltay [mm];fiducial q3 clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3dotHisto =
    TH1I( "cmsdyfcq3dot",
	  "fiducial #Deltay cut x q3 dot;cluster - triplet #Deltay [mm];fiducial q3 dot clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3nodHisto =
    TH1I( "cmsdyfcq3nod",
	  "fiducial #Deltay cut x q3 no dot;cluster - triplet #Deltay [mm];fiducial q3 no dot clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3oddHisto =
    TH1I( "cmsdyfcq3odd",
	  "fiducial #Deltay cut x q3 odd col;cluster - triplet #Deltay [mm];fiducial q3 odd clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3eveHisto =
    TH1I( "cmsdyfcq3eve",
	  "fiducial #Deltay cut x q3 eve col;cluster - triplet #Deltay [mm];fiducial q3 eve clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3midHisto =
    TH1I( "cmsdyfcq3mid",
	  "fiducial #Deltay cut x q3 mid;cluster - triplet #Deltay [mm];fiducial q3 mid clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq3rimHisto =
    TH1I( "cmsdyfcq3rim",
	  "fiducial #Deltay cut x q3 rim;cluster - triplet #Deltay [mm];fiducial q3 rim clusters",
	  500, -0.5, 0.5 );

  TH1I cmsdyfcq4nodHisto =
    TH1I( "cmsdyfcq4nod",
	  "fiducial #Deltay cut x q4 no dot;cluster - triplet #Deltay [mm];fiducial q4 no dot clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq4oddHisto =
    TH1I( "cmsdyfcq4odd",
	  "fiducial #Deltay cut x q4 odd col;cluster - triplet #Deltay [mm];fiducial q4 odd clusters",
	  500, -0.5, 0.5 );
  TH1I cmsdyfcq4eveHisto =
    TH1I( "cmsdyfcq4eve",
	  "fiducial #Deltay cut x q4 eve col;cluster - triplet #Deltay [mm];fiducial q4 eve clusters",
	  500, -0.5, 0.5 );

  TProfile cmsrmsyvsq = TProfile( "cmsrmsyvsq",
			   "DUT #Deltay vs Q0;normal cluster charge [ke];MAD(#Deltay>) [mm]",
				150, 0, 150, 0, 0.2 );

  TProfile cmsdyvsx = TProfile( "cmsdyvsx",
			   "DUT #Deltay vs x;x track [mm];<cluster - triplet #Deltay> [mm]",
				76, -3.8, 3.8, -0.2, 0.2 );
  TProfile cmsdyvsy = TProfile( "cmsdyvsy",
			   "DUT #Deltay vs y;y track [mm];<cluster - triplet #Deltay> [mm]",
				76, -3.8, 3.8, -0.2, 0.2 );
  TProfile cmsdyvsty =
    TProfile( "cmsdyvsty",
	      "DUT #Deltay vs #theta_{y};y track slope [rad];<cluster - triplet #Deltay> [mm]",
	      80, -0.002, 0.002, -0.2, 0.2 );

  TH1I cmsq0fHisto =
    TH1I( "cmsq0f",
	  "normal fiducial cluster charge;normal cluster charge [ke];linked fiducial clusters",
	  160, 0, 80 );
  TH1I cmsq0fdotHisto =
    TH1I( "cmsq0fdot",
	  "normal fiducial cluster charge dot;normal cluster charge dot [ke];linked fiducial dot clusters",
	  160, 0, 80 );
  TH1I cmsq0fnodHisto =
    TH1I( "cmsq0fnod",
	  "normal fiducial cluster charge no dot;normal cluster charge no dot [ke];linked fiducial no dot clusters",
	  160, 0, 80 );
  TH1I cmsq0foddHisto =
    TH1I( "cmsq0fodd",
	  "normal fiducial cluster charge odd col;normal cluster charge odd col [ke];linked fiducial odd col clusters",
	  160, 0, 80 );
  TH1I cmsq0feveHisto =
    TH1I( "cmsq0feve",
	  "normal fiducial cluster charge eve col;normal cluster charge eve col [ke];linked fiducial eve col clusters",
	  160, 0, 80 );
  TH1I cmsq0fcHisto =
    TH1I( "cmsq0fc",
	  "normal fiducial cluster charge core;normal cluster charge core [ke];linked fiducial core clusters",
	  160, 0, 80 );
  TH1I cmsq0fpHisto =
    TH1I( "cmsq0fp",
	  "normal fiducial cluster charge peri;normal cluster charge peri [ke];linked fiducial peri clusters",
	  160, 0, 80 );

  TProfile cmsqxvsx = TProfile( "cmsqxvsx",
			   "DUT cluster charge vs x;x track [mm];<cluster charge> [ke]",
				76, -3.8, 3.8, 0, 0.1 );
  TProfile cmsqxvsy = TProfile( "cmsqxvsy",
			   "DUT cluster charge vs y;y track [mm];<cluster charge> [ke]",
				76, -3.8, 3.8, 0, 0.1 );
  TProfile cmsqxvsxm = TProfile( "cmsqxvsxm",
			   "DUT cluster charge vs xmod;x track mod 0.3 [mm];<cluster charge> [ke]",
				60, 0, 0.3, 0, 0.1 );
  TProfile cmsqxvsym = TProfile( "cmsqxvsym",
			   "DUT cluster charge vs ymod;y track mod 0.2 [mm];<cluster charge> [ke]",
				40, 0, 0.2, 0, 0.1 );
  TProfile cmsqxvsymn =
    TProfile( "cmsqxvsymn",
	      "DUT cluster charge vs ymod no dot;y track mod 0.2 [mm];<no dot cluster charge> [ke]",
	      40, 0, 0.2, 0, 0.1 );
  TProfile cmsqxvsymd =
    TProfile( "cmsqxvsymd",
	      "DUT cluster charge vs ymod dot;y track mod 0.2 [mm];<dot cluster charge> [ke]",
	      40, 0, 0.2, 0, 0.1 );
  TProfile cmsqxvsxmd =
    TProfile( "cmsqxvsxmd",
	      "DUT cluster charge vs xmod dot;x track mod 0.3 [mm];<dot cluster charge> [ke]",
	      60, 0, 0.3, 0, 0.1 );
  TProfile2D cmsqxvsxmym =
    TProfile2D( "cmsqxvsxmym",
	      "DUT cluster charge vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];<dot cluster charge> [ke]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.1 );

  TProfile cmsqxvst1 =
    TProfile( "cmsqxvst1",
	      "DUT cluster charge vs time;time [s];<cluster charge> [ke]",
	      70, 0, 350, 0, 0.1 );
  TProfile cmsqxvst2 =
    TProfile( "cmsqxvst2",
	      "DUT cluster charge vs time;time [s];<cluster charge> [ke]",
	      140, 0, 700, 0, 0.1 );
  TProfile cmsqxvst3 =
    TProfile( "cmsqxvst3",
	      "DUT cluster charge vs time;time [s];<cluster charge> [ke]",
	      200, 0, 2000, 0, 0.1 );
  TProfile cmsqxvst4 =
    TProfile( "cmsqxvst4",
	      "DUT cluster charge vs time;time [s];<cluster charge> [ke]",
	      400, 0, 4000, 0, 0.1 );

  TH1I cmsnpxHisto =
    TH1I( "cmsnpx",
	  "linked DUT cluster size;cluster size [pixels];linked fiducial clusters",
	  20, 0.5, 20.5 );
  TH1I cmsncolHisto =
    TH1I( "cmsncol",
	  "linked DUT cluster size;cluster size [columns];linked fiducial clusters",
	  52, 0.5, 52.5 );
  TH1I cmsnrowHisto =
    TH1I( "cmsnrow",
	  "linked DUT cluster size;cluster size [rows];linked fiducial clusters",
	  20, 0.5, 20.5 );
  TProfile cmsnpxvsq =
    TProfile( "cmsnpxvsq",
	      "DUT cluster size vs Q0;normal cluster charge [ke];<cluster size> [pixels]",
	      150, 0, 150, 0, 20 );

  TH1I cmspxqHisto =
    TH1I( "cmspxq",
	  "DUT pixel charge linked;pixel charge [ke];linked pixels",
	  100, 0, 25 );
  TH1I cmspxqoddHisto =
    TH1I( "cmspxqodd",
	  "DUT pixel charge linked odd col;pixel charge [ke];linked odd col pixels",
	  100, 0, 25 );
  TH1I cmspxqeveHisto =
    TH1I( "cmspxqeve",
	  "DUT pixel charge linked even col;pixel charge [ke];linked even col pixels",
	  100, 0, 25 );
  TH1I cmspxq1Histo =
    TH1I( "cmspxq1",
	  "DUT pixel charge linked 1-px;pixel charge [ke];linked 1-px pixels",
	  100, 0, 25 );
  TH1I cmspxq2Histo =
    TH1I( "cmspxq2",
	  "DUT pixel charge linked 2-px;pixel charge [ke];linked 2-px pixels",
	  100, 0, 25 );
  TH1I cmspxq3Histo =
    TH1I( "cmspxq3",
	  "DUT pixel charge linked 3-px;pixel charge [ke];linked 3-px pixels",
	  100, 0, 25 );
  TProfile cmspxqvsq =
    TProfile( "cmspxqvsq",
	      "DUT pixel charge vs Q0;normal cluster charge [ke];<pixel charge> [ke",
	      150, 0, 150, 0, 50 );
  TProfile cmspxqvsxm =
    TProfile( "cmspxqvsxm",
	      "DUT pixel charge vs xmod;x track mod 0.3 [mm];<pixel charge> [ke]",
	      60, 0, 0.3, 0, 50 );
  TProfile cmspxqvsym =
    TProfile( "cmspxqvsym",
	      "DUT pixel charge vs ymod;y track mod 0.2 [mm];<pixel charge> [ke]",
	      40, 0, 0.2, 0, 50 );

  TProfile cmsdxvsxm =
    TProfile( "cmsdxvsxm",
	      "DUT x residual vs xmod;x track mod 0.3 [mm];<x residual> [mm]",
	      60, 0, 0.3, -0.2, 0.2 );
  TProfile cmsdxvsym =
    TProfile( "cmsdxvsym",
	      "DUT x residual vs ymod;y track mod 0.2 [mm];<x residual> [mm]",
	      40, 0, 0.2, -0.2, 0.2 );
  TProfile cmsdyvsxm =
    TProfile( "cmsdyvsxm",
	      "DUT y residual vs xmod;x track mod 0.3 [mm];<y residual> [mm]",
	      60, 0, 0.3, -0.2, 0.2 );
  TProfile cmsdyvsym =
    TProfile( "cmsdyvsym",
	      "DUT y residual vs ymod;y track mod 0.2 [mm];<y residual> [mm]",
	      40, 0, 0.2, -0.2, 0.2 );
  TProfile cmsdyvsymodd =
    TProfile( "cmsdyvsymodd",
	      "DUT y residual vs ymod odd col;y track mod 0.2 [mm];<y residual> odd col [mm]",
	      40, 0, 0.2, -0.2, 0.2 );
  TProfile cmsdyvsymeve =
    TProfile( "cmsdyvsymeve",
	      "DUT y residual vs ymod even col;y track mod 0.2 [mm];<y residual> even col [mm]",
	      40, 0, 0.2, -0.2, 0.2 );

  TProfile cmsrmsxvsx =
    TProfile( "cmsrmsxvsx",
	      "DUT x resolution vs x;x track [mm];MAD(#Deltax) [mm]",
	      76, -3.8, 3.8, 0, 0.2 );
  TProfile cmsrmsyvsx =
    TProfile( "cmsrmsyvsx",
	      "DUT y resolution vs x;x track [mm];MAD(#Deltay) [mm]",
	      76, -3.8, 3.8, 0, 0.2 );
  TProfile cmsrmsxvsy =
    TProfile( "cmsrmsxvsy",
	      "DUT x resolution vs y;y track [mm];MAD(#Deltax) [mm]",
	      76, -3.8, 3.8, 0, 0.2 );
  TProfile cmsrmsyvsy =
    TProfile( "cmsrmsyvsy",
	      "DUT y resolution vs y;y track [mm];MAD(#Deltay) [mm]",
	      76, -3.8, 3.8, 0, 0.2 );

  TProfile cmsrmsxvsxm =
    TProfile( "cmsrmsxvsxm",
	      "DUT x resolution vs xmod;x track mod 0.3 [mm];MAD(#Deltax) [mm]",
	      60, 0, 0.3, 0, 0.2 );
  TProfile cmsrmsyvsxm =
    TProfile( "cmsrmsyvsxm",
	      "DUT y resolution vs xmod;x track mod 0.3 [mm];MAD(#Deltay) [mm]",
	      60, 0, 0.3, 0, 0.2 );
  TProfile cmsrmsxvsym =
    TProfile( "cmsrmsxvsym",
	      "DUT x resolution vs ymod;y track mod 0.2 [mm];MAD(#Deltax) [mm]",
	      40, 0, 0.2, 0, 0.2 );
  TProfile cmsrmsyvsym =
    TProfile( "cmsrmsyvsym",
	      "DUT y resolution vs ymod;y track mod 0.2 [mm];MAD(#Deltay) [mm]",
	      40, 0, 0.2, 0, 0.2 );

  TProfile2D cmsrmsxvsxmym =
    TProfile2D( "cmsrmsxvsxmym",
		"DUT x resolution vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];MAD(#Deltax) [mm]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.2 );
  TProfile2D cmsrmsyvsxmym =
    TProfile2D( "cmsrmsyvsxmym",
		"DUT y resolution vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];MAD(#Deltay) [mm]",
		60, 0, 0.3, 40, 0, 0.2, 0, 0.2 );

  TProfile cmsncolvsx =
    TProfile( "cmsncolvsx",
	      "DUT cluster size vs x;x track [mm];<cluster size> [columns]",
	      76, -3.8, 3.8, 0, 20 );
  TProfile cmsncolvsy =
    TProfile( "cmsncolvsy",
	      "DUT cluster size vs y;y track [mm];<cluster size> [columns]",
	      76, -3.8, 3.8, 0, 20 );

  TProfile cmsncolvsxm =
    TProfile( "cmsncolvsxm",
	      "DUT cluster size vs xmod;x track mod 0.3 [mm];<cluster size> [columns]",
	      60, 0, 0.3, 0, 20 );
  TProfile cmsnrowvsxm =
    TProfile( "cmsnrowvsxm",
	      "DUT cluster size vs xmod;x track mod 0.3 [mm];<cluster size> [rows]",
	      60, 0, 0.3, 0, 20 );

  TProfile cmsncolvsym =
    TProfile( "cmsncolvsym",
	      "DUT cluster size vs ymod;y track mod 0.2 [mm];<cluster size> [columns]",
	      40, 0, 0.2, 0, 20 );
  TProfile cmsnrowvsym =
    TProfile( "cmsnrowvsym",
	      "DUT cluster size vs ymod;y track mod 0.2 [mm];<cluster size> [rows]",
	      40, 0, 0.2, 0, 20 );
  TProfile2D cmsnpxvsxmym =
    TProfile2D( "cmsnpxvsxmym",
	      "DUT cluster size vs xmod ymod;x track mod 0.3 [mm];y track mod 0.2 [mm];<cluster size> [pixels]",
		60, 0, 0.3, 40, 0, 0.2, 0, 20 );

  TH2I trixylkHisto = TH2I( "trixylk",
			    "triplets with DUT;x [mm];y [mm];DUT-linked triplets",
			    240, -12, 12, 120, -6, 6 );
  TH2I sixxylkHisto = TH2I( "sixxylk",
			    "REF-linked sixplets at z DUT;x [mm];y [mm];REF-linked sixplets",
			    240, -12, 12, 120, -6, 6 );
  TH2I sixxyeffHisto = TH2I( "sixxyeff",
			    "REF-linked sixplets with DUT;x [mm];y [mm];DUT-REF-linked sixplets",
			    240, -12, 12, 120, -6, 6 );

  TProfile effvsx = TProfile( "effvsx",
			      "DUT efficiency vs x;x track at DUT [mm];efficiency",
			      54, -4.05, 4.05, -1, 2 ); // bin = col
  TProfile effvsy = TProfile( "effvsy",
			      "DUT efficiency vs y;y track at DUT [mm];efficiency",
			      80, -4, 4, -1, 2 ); // bin = row

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  cout << endl;

  FileReader * reader;
  if( run < 100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X");
  else if( run < 1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X");
  else if( run < 10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X");
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X");
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X");

  int event_nr = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock
  uint64_t prevTLU = 0;

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() )
      eudaq::PluginManager::Initialize(evt);

    bool ldbg = 0;

    if( event_nr <  0 )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( event_nr < 2  )
      evTLU0 = evTLU;
    double evsec = (evTLU - evTLU0) / fTLU;
    double evdt = (evTLU - prevTLU) / fTLU;
    hdtus.Fill( evdt * 1E6 ); // [us]
    hdtms.Fill( evdt * 1E3 ); // [ms]
    prevTLU = evTLU;

    if( event_nr < 10 )
      cout<<"Processing event " << event_nr << " time " << evsec << endl;
    else if( event_nr < 100 && event_nr%10 == 0 )
      cout<<"Processing event " << event_nr << " time " << evsec << endl;
    else if( event_nr < 1000 && event_nr%100 == 0 )
      cout<<"Processing event " << event_nr << " time " << evsec << endl;
    else if( event_nr%1000 == 0 )
      cout<<"Processing event " << event_nr << " time " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    int xm = 0;
    int ym = 0;
    int adc = 0;
    double q = 0;

    vector <cluster> cl[9];

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg ) std::cout << "PLANE " << plane.ID() << ": ";

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  std::cout << plane.GetX(ipix)
		    << " " << plane.GetY(ipix)
		    << " " << plane.GetPixel(ipix) << " ";

	xm = plane.GetX(ipix); // global column 0..415
	ym = plane.GetY(ipix); // global row 0..159
	adc = plane.GetPixel(ipix); // ADC 0..255

	q = adc;

	if( ipl == iDUT &&
	    adc > 0 &&
	    xm >= 0 && xm < 52 &&
	    ym >= 0 && ym < 80 ) {

	  double Ared = adc - p4[xm][ym]; // p4 is asymptotic maximum

	  if( Ared >= 0 )
	    Ared = -0.1; // avoid overflow

	  double a3 = p3[xm][ym]; // positive
	  if( weib == 3 )
	    q = p1[xm][ym] *
	      ( pow( -log( -Ared / a3 ), 1/p2[xm][ym] ) - p0[xm][ym] ) * ke;
	  // q = ( (-ln(-(A-p4)/p3))^1/p2 - p0 )*p1

	}  // DUT

	if( ipl == iREF &&
	    adc > 0 &&
	    xm >= 0 && xm < 52 &&
	    ym >= 0 && ym < 80 ) {

	  double Ared = adc - r4[xm][ym]; // r4 is asymptotic maximum

	  if( Ared >= 0 )
	    Ared = -0.1; // avoid overflow

	  double a3 = r3[xm][ym]; // positive
	  if( weib == 3 )
	    q = r1[xm][ym] *
	      ( pow( -log( -Ared / a3 ), 1/r2[xm][ym] ) - r0[xm][ym] ) * refke;
	  // q = ( (-ln(-(A-r4)/r3))^1/r2 - r0 )*r1

	}  // REF

	hcol[ipl].Fill( xm );
	hrow[ipl].Fill( ym );

	// fill pixel block for clustering:

	// mask telescope hot pixels: to do

	pb[npx].col = xm;
	pb[npx].row = ym;
	pb[npx].adc = adc;
	pb[npx].q = q;
	pb[npx].ord = npx; // readout order
	++npx;

	if( npx == 999 ) {
	  cout << "pixel buffer overflow in plane " << ipl
	       << ", event " << event_nr
	       << endl;
	  break;
	}

      } // pix

      hnpx[ipl].Fill(npx);

      if( ldbg ) std::cout << std::endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[ipl] = getClus();

      if( ldbg ) cout << "A clusters " << cl[ipl].size() << endl;

      hncl[ipl].Fill( cl[ipl].size() );

      for( vector<cluster>::iterator cA = cl[ipl].begin(); cA != cl[ipl].end(); ++cA ) {

	hsiz[ipl].Fill( cA->size );
	hncol[ipl].Fill( cA->ncol );
	hnrow[ipl].Fill( cA->nrow );

      } // cl A

    } // planes

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make driplets 3+5-4:

    vector <triplet> driplets;

    double driCut = 0.1; // [mm]

    for( vector<cluster>::iterator cA = cl[3].begin(); cA != cl[3].end(); ++cA ) {

      double xA = cA->col*ptchx[3] - alignx[3];
      double yA = cA->row*ptchy[3] - aligny[3];
      double xmid = xA - midx[3];
      double ymid = yA - midy[3];
      xA = xmid - ymid*rotx[3];
      yA = ymid + xmid*roty[3];

      for( vector<cluster>::iterator cC = cl[5].begin(); cC != cl[5].end(); ++cC ) {

	double xC = cC->col*ptchx[5] - alignx[5];
	double yC = cC->row*ptchy[5] - aligny[5];
	double xmid = xC - midx[5];
	double ymid = yC - midy[5];
	xC = xmid - ymid*rotx[5];
	yC = ymid + xmid*roty[5];

	double dx2 = xC - xA;
	double dy2 = yC - yA;
	double dz35 = zz[5] - zz[3]; // from 3 to 5 in z
	hdx35.Fill( dx2 );
	hdy35.Fill( dy2 );

	if( abs( dx2 ) > 0.005 * dz35 ) continue; // angle cut
	if( abs( dy2 ) > 0.005 * dz35 ) continue; // angle cut

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zz[3] + zz[5] ); // mid z
 
	double slpx = ( xC - xA ) / dz35; // slope x
	double slpy = ( yC - yA ) / dz35; // slope y

	// middle plane B = 4:

	for( vector<cluster>::iterator cB = cl[4].begin(); cB != cl[4].end(); ++cB ) {

	  double xB = cB->col*ptchx[4] - alignx[4];
	  double yB = cB->row*ptchy[4] - aligny[4];
	  double xmid = xB - midx[4];
	  double ymid = yB - midy[4];
	  xB = xmid - ymid*rotx[4];
	  yB = ymid + xmid*roty[4];

	  // interpolate track to B:

	  double dz = zz[4] - avz;
	  double xk = avx + slpx * dz; // driplet at k
	  double yk = avy + slpy * dz;

	  double dx3 = xB - xk;
	  double dy3 = yB - yk;
	  hdxdri.Fill( dx3 );
	  hdydri.Fill( dy3 );

	  if( abs( dy3 ) < 0.02 ) {
	    hdxdric.Fill( dx3 );
	  }
	  if( abs( dx3 ) < 0.02 )
	    hdydric.Fill( dy3 );

	  // telescope driplet cuts:

	  if( fabs(dx3) > driCut ) continue;
	  if( fabs(dy3) > driCut ) continue;

	  triplet dri;
	  dri.xm = avx;
	  dri.ym = avy;
	  dri.zm = avz;
	  dri.sx = slpx;
	  dri.sy = slpy;
	  dri.lk = 0;
	  driplets.push_back(dri);

	} // cl B

      } // cl C

    } // cl A

    ndriHisto.Fill( driplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // driplets vs REF:

    double xcutREF = 0.3;
    double ycutREF = 0.2;

    const double cr = cos( REFrot );
    const double sr = sin( REFrot );

    for( unsigned int iB = 0; iB < driplets.size(); ++iB ) { // iB = downstream

      double avx = driplets[iB].xm;
      double avy = driplets[iB].ym;
      double avz = driplets[iB].zm;
      double slx = driplets[iB].sx;
      double sly = driplets[iB].sy;

      double zB = REFz - avz; // z REF from mid of driplet
      double xB = avx + slx * zB; // driplet impact point on REF
      double yB = avy + sly * zB;

      // transform into REF system:

      double x3 = xB*cr + yB*sr; // rotate
      double y3 = yB*cr - xB*sr;

      double x4 = -x3 + REFalignx; // invert and shift to mid
      double y4 =  y3 + REFaligny; // shift to mid

      double xmod = fmod( 9.075 + x4, 0.3 ); // [0,0.3] mm, 2 pixel wide
      double ymod = fmod( 9.050 + y4, 0.2 ); // [0,0.2] mm

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // driplets vs REF clusters:

      for( vector<cluster>::iterator c = cl[iREF].begin(); c != cl[iREF].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;
	double refx = ( ccol - nx[iREF]/2 ) * ptchx[iREF]; // -3.9..3.9 mm
	double refy = ( crow - ny[iREF]/2 ) * ptchy[iREF]; // -4..4 mm
	double q = c->charge;
	double qx = exp( -q/qwid);

	int npx = c->size;

	// residuals for pre-alignment:

	refsxaHisto.Fill( refx + x3 ); // peak
	refdxaHisto.Fill( refx - x3 ); // 

	refsyaHisto.Fill( refy + y3 ); // 
	refdyaHisto.Fill( refy - y3 ); // peak

	double refdx = refx - x4;
	double refdy = refy - y4;

	if( abs( refdx ) < xcutREF ) {
	  refdyHisto.Fill( refdy );
	  refdyvsx.Fill( x4, refdy );
	  refdyvsy.Fill( y4, refdy );
	  refdyvsty.Fill( sly, refdy );
	}

	if( abs( refdy ) < ycutREF ) {
	  refdxHisto.Fill( refdx );
	}

	if( abs( refdx ) < xcutREF && abs( refdy ) < ycutREF ) {
	  refnpxHisto.Fill( npx );
	  refqHisto.Fill( q );
	  refnpxvsxmym.Fill( xmod, ymod, npx );
	  refqxvsxmym.Fill( xmod, ymod, qx );
	}

	if( fabs( refdy ) < 0.10 && fabs( refdx ) < 0.15 ) { // tight cuts for eff
	  driplets[iB].lk = 1;
	  reflkxHisto.Fill( xB );
	  reflkyHisto.Fill( yB );
	  reflkcolHisto.Fill( ccol );
	  reflkrowHisto.Fill( crow );
	}

      } // REF

    } // driplets

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

      if( c->vpix.size() < 2 ) continue; // skip 1-pix clusters

      if( c->vpix.size() > 2 ) continue; // want 2-pix clusters

      vector<pixel>::iterator pxa = c->vpix.begin();
      vector<pixel>::iterator pxb = pxa; // lower row (sorted along 80*col + row)
      ++pxb; // higher row in clustering

      if( pxa->ord > 0 && pxb->ord > 0 ) continue; // none is 1st

      double q1 = pxa->q; // read out first
      double q2 = pxb->q;
      if( pxb->ord == 0 ) {
	q1 = pxb->q; // read out first
	q2 = pxa->q;
      }

      dutpxq1stHisto.Fill( q1 ); // read out first on ROC
      dutpxq2ndHisto.Fill( q2 );

      dutpxq2nd02Histo.Fill( q2 - 0.02*q1 ); // 
      dutpxq2nd03Histo.Fill( q2 - 0.03*q1 ); // good for 506
      dutpxq2nd04Histo.Fill( q2 - 0.04*q1 ); // too little
      dutpxq2nd05Histo.Fill( q2 - 0.05*q1 ); // good
      dutpxq2nd06Histo.Fill( q2 - 0.06*q1 ); // best for digV2.1
      dutpxq2nd07Histo.Fill( q2 - 0.07*q1 ); // too much
      dutpxq2nd08Histo.Fill( q2 - 0.08*q1 ); // 
      dutpxq2nd09Histo.Fill( q2 - 0.09*q1 ); // 

    } // DUT cl

    double eps = 0.06;

    // apply:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

      // pix in clus:

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	if( px->ord == 0 ) continue;

	// look for previous, in all clusters:

	double qprv = 0;

	for( vector<cluster>::iterator d = cl[iDUT].begin(); d != cl[iDUT].end(); ++d )
	  for( vector<pixel>::iterator qx = d->vpix.begin(); qx != d->vpix.end(); ++qx )
	    if( qx->ord == px->ord - 1 ) 
	      qprv = qx->q;

	// apply tsunami correction:

	px->q -= eps*qprv; // overwrite!
	//px->q -= eps; // subtract constant, overwrite!

      } // pix

    } // DUT clusters

    // tsunami corrected clusters:

    for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

      int colmin = 99;
      int colmax = -1;
      int rowmin = 99;
      int rowmax = -1;

      double qsum = 0;
      double sumcol = 0;
      double sumrow = 0;

      double qcol[52];
      for( int icol = 0; icol < 52; ++icol ) qcol[icol] = 0;

      double qrow[80];
      for( int irow = 0; irow < 80; ++irow ) qrow[irow] = 0;

      for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	int icol = px->col;
	int irow = px->row;

	dutcolHisto.Fill( icol );
	dutrowHisto.Fill( irow );

	if( icol < colmin ) colmin = icol;
	if( icol > colmax ) colmax = icol;
	if( irow < rowmin ) rowmin = irow;
	if( irow > rowmax ) rowmax = irow;

	double q = px->q; // corrected
	if( q < 0 ) continue;

	qsum += q;
	qcol[icol] += q; // project cluster onto cols
	qrow[irow] += q; // project cluster onto rows

	sumcol += icol*q;
	sumrow += irow*q;

      } // pix

      int ncol = colmax - colmin + 1;
      int nrow = rowmax - rowmin + 1;

      // cluster charge after tsunami correction:

      //28.3.2015 c->charge = qsum; // overwrite creates shift in cmsqvsxmym
      c->col = sumcol / qsum; // overwrite !
      c->row = sumrow / qsum; // overwrite !

      if( colmin > 0 && colmax < 51 && rowmin > 0 && rowmax <79 ) {
	dutq0fHisto.Fill( qsum * norm );
	dutncolHisto.Fill( ncol );
	dutnrowHisto.Fill( nrow );
      }

      if( ncol > 2 ) {

	dutcolminHisto.Fill( colmin ); // strong even peaks at large turn 14848
	dutcolmaxHisto.Fill( colmax ); // weaker even peaks

	dutcol0qHisto.Fill(qcol[colmin]);
	if( colmin%2 )
	  dutcol0oddqHisto.Fill(qcol[colmin]);
	else
	  dutcol0eveqHisto.Fill(qcol[colmin]);
	dutcol9qHisto.Fill(qcol[colmax]);

	if( ncol > 1 ) dutcol1qHisto.Fill(qcol[colmin+1]);
	if( ncol > 2 ) dutcol2qHisto.Fill(qcol[colmin+2]);
	if( ncol > 3 ) dutcol3qHisto.Fill(qcol[colmin+3]);
	if( ncol > 4 ) dutcol4qHisto.Fill(qcol[colmin+4]);
	if( ncol > 5 ) dutcol5qHisto.Fill(qcol[colmin+5]);
	if( ncol > 6 ) dutcol6qHisto.Fill(qcol[colmin+6]);
	if( ncol > 7 ) dutcol7qHisto.Fill(qcol[colmin+7]);
	if( ncol > 8 ) dutcol8qHisto.Fill(qcol[colmin+8]);

      } // long

    } // DUT clusters

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // make triplets 2+0-1:

    vector <triplet> triplets;

    double triCut = 0.1; // [mm]

    for( vector<cluster>::iterator cA = cl[0].begin(); cA != cl[0].end(); ++cA ) {

      double xA = cA->col*ptchx[0] - alignx[0];
      double yA = cA->row*ptchy[0] - aligny[0];
      double xmid = xA - midx[0];
      double ymid = yA - midy[0];
      xA = xmid - ymid*rotx[0];
      yA = ymid + xmid*roty[0];

      for( vector<cluster>::iterator cC = cl[2].begin(); cC != cl[2].end(); ++cC ) {

	double xC = cC->col*ptchx[2] - alignx[2];
	double yC = cC->row*ptchy[2] - aligny[2];
	double xmid = xC - midx[2];
	double ymid = yC - midy[2];
	xC = xmid - ymid*rotx[2];
	yC = ymid + xmid*roty[2];

	double dx2 = xC - xA;
	double dy2 = yC - yA;
	double dz02 = zz[2] - zz[0]; // from 0 to 2 in z
	hdx02.Fill( dx2 );
	hdy02.Fill( dy2 );

	if( abs( dx2 ) > 0.005 * dz02 ) continue; // angle cut
	if( abs( dy2 ) > 0.005 * dz02 ) continue; // angle cut

	double avx = 0.5 * ( xA + xC ); // mid
	double avy = 0.5 * ( yA + yC );
	double avz = 0.5 * ( zz[0] + zz[2] ); // mid z
 
	double slpx = ( xC - xA ) / dz02; // slope x
	double slpy = ( yC - yA ) / dz02; // slope y

	// middle plane B = 1:

	for( vector<cluster>::iterator cB = cl[1].begin(); cB != cl[1].end(); ++cB ) {

	  double xB = cB->col*ptchx[1] - alignx[1];
	  double yB = cB->row*ptchy[1] - aligny[1];
	  double xmid = xB - midx[1];
	  double ymid = yB - midy[1];
	  xB = xmid - ymid*rotx[1];
	  yB = ymid + xmid*roty[1];

	  // interpolate track to B:

	  double dz = zz[1] - avz;
	  double xk = avx + slpx * dz; // triplet at k
	  double yk = avy + slpy * dz;

	  double dx3 = xB - xk;
	  double dy3 = yB - yk;
	  hdxtri.Fill( dx3 );
	  hdytri.Fill( dy3 );

	  if( abs( dy3 ) < 0.02 ) {

	    hdxtric.Fill( dx3 );

	    if(      cB->size == 1 )
	      hdxtris1.Fill( dx3 ); // 4.2 um
	    else if( cB->size == 2 )
	      hdxtris2.Fill( dx3 ); // 4.0 um
	    else if( cB->size == 3 )
	      hdxtris3.Fill( dx3 ); // 3.8 um
	    else if( cB->size == 4 )
	      hdxtris4.Fill( dx3 ); // 4.3 um
	    else
	      hdxtris5.Fill( dx3 ); // 3.6 um

	    if(      cB->ncol == 1 )
	      hdxtric1.Fill( dx3 ); // 4.0 um
	    else if( cB->ncol == 2 )
	      hdxtric2.Fill( dx3 ); // 4.1 um
	    else if( cB->ncol == 3 )
	      hdxtric3.Fill( dx3 ); // 3.6 um
	    else if( cB->ncol == 4 )
	      hdxtric4.Fill( dx3 ); // 3.5 um
	    else
	      hdxtric5.Fill( dx3 ); // 4.1 um
	  }
	  if( abs( dx3 ) < 0.02 )
	    hdytric.Fill( dy3 );

	  // telescope triplet cuts:

	  if( fabs(dx3) > triCut ) continue;
	  if( fabs(dy3) > triCut ) continue;

	  triplet tri;
	  tri.xm = avx;
	  tri.ym = avy;
	  tri.zm = avz;
	  tri.sx = slpx;
	  tri.sy = slpy;
	  triplets.push_back(tri);

	} // cl B

      } // cl C

    } // cl A

    ntriHisto.Fill( triplets.size() );

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets:

    double xcut = 0.4;
    double ycut = 0.2;

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double avx = triplets[iA].xm;
      double avy = triplets[iA].ym;
      double avz = triplets[iA].zm;
      double slx = triplets[iA].sx;
      double sly = triplets[iA].sy;

      double zA = DUTz - avz; // z DUT from mid of triplet
      double xA = avx + slx * zA; // triplet impact point on DUT
      double yA = avy + sly * zA;

      trixHisto.Fill( xA );
      triyHisto.Fill( yA );
      trixyHisto.Fill( xA, yA );

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // tri vs tri: isolation at DUT

      double dmin = 99.9;

      for( unsigned int jj = 0; jj < triplets.size(); ++jj ) {

	if( jj == iA ) continue;

	double zj = DUTz - triplets[jj].zm;
	double xj = triplets[jj].xm + triplets[jj].sx * zj; // triplet impact point on DUT
	double yj = triplets[jj].ym + triplets[jj].sy * zj;

	double dx = xA - xj;
	double dy = yA - yj;
	double dd = sqrt( dx*dx + dy*dy );
	if( dd < dmin ) dmin = dd;

	ttdxHisto.Fill( dx );
	ttdx1Histo.Fill( dx );

      } // jj

      ttdminHisto.Fill( dmin );

      bool liso = 0;
      if( dmin > 0.3 ) liso = 1;

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // match triplet and driplet for efficiency

      bool lsixlk = 0;
      double sixcut = 0.1; // [mm]

      for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // j = B = upstream

	double avxB = driplets[jB].xm;
	double avyB = driplets[jB].ym;
	double avzB = driplets[jB].zm;
	double slxB = driplets[jB].sx;
	double slyB = driplets[jB].sy;

	// driplet at mid:

	double zB = DUTz - avzB; // z from mid of triplet to mid
	double xB = avxB + slxB * zB; // triplet at mid
	double yB = avyB + slyB * zB;

	// driplet - triplet:

	double dx = xB - xA;
	double dy = yB - yA;

	hsixdx.Fill( dx ); // for align fit
	hsixdy.Fill( dy ); // for align fit
	if( abs(dy) < sixcut ) {
	  hsixdxc.Fill( dx );
	  sixdxvsy.Fill( yB, dx );
	  sixdxvstx.Fill( slx, dx );
	}
	if( abs(dx) < sixcut ) {
	  hsixdyc.Fill( dy );
	  sixdyvsx.Fill( xB, dy );
	  sixdyvsty.Fill( sly, dy );
	}

	// compare slopes:

	if( abs(dy) < 0.1 && abs(dx) < 0.1 ) {
	  hsixdslpx.Fill( slxB - slx );
	  hsixdslpy.Fill( slyB - sly ); // width: 0.3 mrad
	}

	if( fabs(dx) < sixcut && fabs(dy) < sixcut ) {
	  sixxyHisto.Fill( xA, yA );
	  if( driplets[jB].lk ) {
	    lsixlk = 1;
	  }
	}

      } // driplets

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // intersect inclined track with tilted DUT plane:

      double zc = (Nz*zA - Ny*avy - Nx*avx) / (Nx*slx + Ny*sly + Nz); // from avz
      double yc = avy + sly * zc;
      double xc = avx + slx * zc;

      double dzc = zc + avz - DUTz; // from DUT z0 [-8,8] mm

      // transform into DUT system: (passive).
      // large rotations don't commute: careful with order

      double x1 = co*xc - so*dzc; // turn o
      double y1 = yc;
      double z1 = so*xc + co*dzc;

      double x2 = x1;
      double y2 = ca*y1 + sa*z1; // tilt a
      double z2 =-sa*y1 + ca*z1; // should be zero (in DUT plane). is zero

      double x3 = cf*x2 + sf*y2; // rot
      double y3 =-sf*x2 + cf*y2;
      double z3 = z2; // should be zero (in DUT plane). is zero

      z3Histo.Fill( z3 ); // is zero

      double upsign = 1; // 504

      double upsignx = upsign;
      double upsigny = upsign;
    
      if( rot90 ) {
	upsignx =  1;
	upsigny = -1;
      }
      double x4 = upsignx*x3 + DUTalignx; // shift to mid
      double y4 =-upsigny*y3 + DUTaligny; // invert y, shift to mid

      bool fiducial = 1;
      if( abs( x4 ) > 3.9 ) fiducial = 0; // skip big col
      if( abs( y4 ) > 3.9 ) fiducial = 0; // skip 1st and last row

      // reduce to 2x2 pixel region:

      double xmod = fmod( 9.075 + x4, 0.3 ); // [0,0.3] mm, 2 pixel wide
      double ymod = fmod( 9.050 + y4, 0.2 ); // [0,0.2] mm
      if( rot90 ) { // x = col = yt, y = row = xt
	xmod = fmod( 9.075 + y4, 0.3 ); // [0,0.3] mm, 2 pixel wide
	ymod = fmod( 9.050 + x4, 0.2 ); // [0,0.2] mm
      }

      // bias dot, from cmsqvsxmym:

      bool ldot = 1;
      if( xmod < 0.105 ) ldot = 0; // dot at x = 125
      if( xmod > 0.195 ) ldot = 0; // and at x = 175
      if( DUTtilt < 6 ) {
	if( ymod < 0.055 ) ldot = 0; // dot at y =  75
	if( ymod > 0.195 ) ldot = 0; // dot at y = 175
	if( ymod > 0.095 && ymod < 0.155 ) ldot = 0; // band between dots
      }

      bool ymid = 0;
      if( ymod > 0.040 && ymod < 0.060 ) ymid = 1;
      if( ymod > 0.140 && ymod < 0.160 ) ymid = 1;

      // pixel core, 2x2 pixel region:

      bool lcore = 1;
      if( xmod < 0.020 ) lcore = 0; // outer edge, see cmsncolvsxm
      if( xmod > 0.280 ) lcore = 0; // outer edge
      if( ymod < 0.020 ) lcore = 0; // outer edge, see cmsnrowvsym
      if( ymod > 0.180 ) lcore = 0; // outer edge
      if( xmod > 0.130 && xmod < 0.170 ) lcore = 0; // inner edge
      if( ymod > 0.080 && ymod < 0.120 ) lcore = 0; // inner edge

      // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
      // DUT pixel clusters:

      int nm = 0; 

      for( vector<cluster>::iterator c = cl[iDUT].begin(); c != cl[iDUT].end(); ++c ) {

	double ccol = c->col;
	double crow = c->row;

	int oddcol = int(ccol) % 2;

	int colmin = 99;
	int colmax = -1;
	int rowmin = 99;
	int rowmax = -1;

	double qcol[52];
	for( int icol = 0; icol < 52; ++icol ) qcol[icol] = 0;

	double qrow[80];
	for( int irow = 0; irow < 80; ++irow ) qrow[irow] = 0;

	for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	  int icol = px->col;
	  if( icol <  0 ) continue;
	  if( icol > 51 ) continue;

	  int irow = px->row;
	  if( irow <  0 ) continue;
	  if( irow > 79 ) continue;

	  if( icol < colmin ) colmin = icol;
	  if( icol > colmax ) colmax = icol;
	  if( irow < rowmin ) rowmin = irow;
	  if( irow > rowmax ) rowmax = irow;

	  double q = px->q; // [ke] corrected
	  if( q < 0 ) continue;

	  qcol[icol] += q; // project cluster onto cols
	  qrow[irow] += q; // project cluster onto rows

	} // pix

	int ncol = colmax - colmin + 1;
	int nrow = rowmax - rowmin + 1;

	// eta-algo in rows:

	double q1 = 0; // highest charge
	double q2 = 0; // 2nd highest
	int i1 = 99;
	int i2 = 99;
	double sumq = 0;
	double sumrow = 0;
	double sumrow2 = 0;
	double sumrow3 = 0;

	for( int irow = rowmin; irow <= rowmax; ++irow ) {

	  double q = qrow[irow];
	  sumq += q;
	  sumrow += irow*q;
	  double drow = irow - crow;
	  sumrow2 += drow*drow*q;
	  sumrow3 += drow*drow*drow*q;

	  if( q > q1 ) {
	    q2 = q1;
	    q1 = q;
	    i2 = i1;
	    i1 = irow;
	  }
	  else if( q > q2 ) {
	    q2 = q;
	    i2 = irow;
	  }

	} // rows
	/*
	double rmsrow = 0;
	double skwrow = 0;
	if( nrow > 1 ) {
	  rmsrow = sqrt( sumrow2/sumq );
	  skwrow = sumrow3/sumq*8/nrow/nrow/nrow; // normalized 3rd moment
	  // from fitskw:
	  crow -= // overwrite!
	    ( skwoff + skwslp * (skwrow-skwmid) + // Sun 5.7.2015
	      skwrng * sin( (skwrow-skwmid) / skwwid ) ) * 1E-3 / pitchrow;
	}
	*/
	double q12 = q1 + q2;
	double eta = 0;
	if( q12 > 1 ) eta = ( q1 - q2 ) / q12;
	if( i2 > i1 ) eta = -eta;

	if( nrow == -2 ) { // inactive, need gain cal

	  // correct even/odd PH effect in readout direction:

	  double A1 = q1; // larger
	  double A2 = q2; // smaller

	  //double A12 = A1 + A2;

	  // cross talk for 2-row clusters:
	  // A1 = (1-cx)U1 + cx U2
	  // A2 = (1-cx)U2 + cx U1

	  // .x fitgp0.C("cmsdyfctq3d",-22,22) (no bias dot)
	  // minimize bias in cmsdyvsym.Draw()

	  //double cx =-0.02; // chip 205, run 10891, sy 7.07
	  double cx = 0.00; // chip 205, run 10891, sy 6.98+-0.05, run 6701 6.98
	  //double cx = 0.02; // chip 205, run 10891, sy 7.10, run 6701 7.06
	  //double cx = 0.04; // chip 205, run 10891, sy 7.21

	  //if( chip0 == 203 ) cx = 0.12; // cmsrmsyvseta worse
	  //if( chip0 == 203 ) cx =-0.12; // cmsrmsyvseta better but tilted sy 11.0
	  //if( chip0 == 203 ) cx =-0.25; // cmsrmsyvseta flat and hump at 0.5, sy 12.0

	  // analog chip 113, run 11553 .x fittp0.C("cmsdyfctq3d",-22,22)
	  // cx -0.02  sy 7.29
	  // cx  0.00  sy 7.22
	  // cx  0.02  sy 7.30 +- 0.035

	  //if( chip0 > 399 ) cx = 0.05; // digV2.1 sy 8.47 
	  //if( chip0 > 399 ) cx = 0.08; // digV2.1 sy 8.11 run 12521
	  if( chip0 > 399 ) cx = 0.10; // digV2.1 sy 8.07
	  //if( chip0 > 399 ) cx = 0.12; // digV2.1 sy 8.14

	  if( chip0 > 499 ) cx = 0.00; // chip 500  sy 8.28
	  //if( chip0 > 499 ) cx =  0.02; // chip 500  sy 8.52
	  //if( chip0 > 499 ) cx = -0.02; // chip 500  sy 8.32

	  // unfold = invert cross talk matrix:

	  double det = 1-2*cx; // Determinante

	  double U1 = ( (1-cx)*A1 - cx*A2 ) / det;
	  double U2 = ( (1-cx)*A2 - cx*A1 ) / det; // U1+U2 = A12

	  // threshold:

	  if( U1 < thr ) U1 = thr; // [ke] threshold
	  if( U2 < thr ) U2 = thr; // This alone improves resolution !

	  double U12 = U1 + U2;

	  if( fabs(cx) > 0.01 ) // for chips 400
	    crow = ( i1*U1 + i2*U2 ) / U12; // overwrite !

	  if( U12 > 1 ) eta = ( U1 - U2 ) / U12;
	  if( i2 > i1 ) eta = -eta;

	} // 2-row

	// column cluster:

	double qhtcol = qcol[colmin] + qcol[colmax];
	double asycol = ( qcol[colmax] - qcol[colmin] ) / qhtcol; // -1..1

	sumq = 0;
	double sumcol = 0;
	double sumcol2 = 0;
	double sumcol3 = 0;
	double p1 = 0; // highest charge
	double p2 = 0; // 2nd highest
	int j1 = 99;
	int j2 = 99;

	for( int icol = colmin; icol <= colmax; ++icol ) {

	  double q = qcol[icol];
	  if( q > p1 ) {
	    p2 = p1;
	    p1 = q;
	    j2 = j1;
	    j1 = icol;
	  }
	  else if( q > p2 ) {
	    p2 = q;
	    j2 = icol;
	  }

	  //Tue 21.7.2015 if( q > 17 ) q = 17; // truncate Landau tail [ke]
	  sumq += q;
	  sumcol += icol*q;
	  double dcol = icol - ccol; // distance from COG
	  sumcol2 += dcol*dcol*q; // 2nd central moment
	  sumcol3 += dcol*dcol*dcol*q; // 3rd central moment

	} // cols

	double rmscol = 0;
	double dm2col = 0;
	double sm2col = 0;
	double skwcol = 0;
	if( ncol > 1 ) {
	  rmscol = sqrt( sumcol2/sumq );
	  dm2col = ( rmscol - 0.288/ptchx[iDUT]*tan(DUTtilt*wt)/sqrt(12) ) / ncol;
	  sm2col = dm2col;
	  if( asycol < 0 )
	    sm2col = -sm2col;
	  //skwcol = sumcol3/sumq; // 3rd moment
	  //skwcol = sumcol3/sumq/rmscol/rmscol/rmscol; // normalized 3rd moment
	  skwcol = sumcol3/sumq*8/ncol/ncol/ncol; // normalized 3rd moment
	  // from fitskw:
	  /*
	  ccol -= // overwrite!
	    ( skwoff + skwslp * (skwcol-skwmid) + // Fri 3.7.2015
	      skwrng * sin( (skwcol-skwmid) / skwwid ) ) * 1E-3 / pitchcol;
	  */
	}

	double p12 = p1 + p2;
	double uta = 0;
	if( p12 > 1 ) uta = ( p1 - p2 ) / p12;
	if( j2 > j1 ) uta = -uta;

	if( rot90 )
	  eta = uta;

	double Q0 = c->charge * norm; // cluster charge normalized to vertical incidence
	double Qx = exp(-Q0/qwid);

	bool lq = 0;
	if( Q0 > 17 &&  Q0 < 30 ) lq = 1; // for Q0 Landau peak at 22

	// DUT - triplet:

	double cmsx = ( ccol - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	double cmsy = ( crow - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm

	if( rot90 ) {
	  cmsx = ( crow - ny[iDUT]/2 ) * ptchy[iDUT]; // -4..4 mm
	  cmsy = ( ccol - nx[iDUT]/2 ) * ptchx[iDUT]; // -3.9..3.9 mm
	}

	// residuals for pre-alignment:

	cmssxaHisto.Fill( cmsx + x3 ); // rot, tilt and turn but no shift
	cmsdxaHisto.Fill( cmsx - x3 ); // peak

	cmssyaHisto.Fill( cmsy + y3 ); // peak
	cmsdyaHisto.Fill( cmsy - y3 );

	// residuals:

	double cmsdx = cmsx - x4; // new style
	double cmsdy = cmsy - y4;

	cmsdxHisto.Fill( cmsdx );
	cmsdyHisto.Fill( cmsdy );

	if( chip0 == 500 ) { // 28.3.2015
	  if( oddcol )
	    cmsdy = cmsdy - 1.2E-3;
	  else
	    cmsdy = cmsdy + 1.2E-3;
	}
	if( chip0 == 504 ) {
	  if( oddcol )
	    cmsdy = cmsdy - 1.0E-3;
	  else
	    cmsdy = cmsdy + 1.0E-3;
	}

	if( abs(cmsdx) < xcut  && abs(cmsdy) < ycut ) {

	  cmslkxHisto.Fill( xc );
	  cmslkyHisto.Fill( yc );
	  cmscolHisto.Fill( ccol ); // map
	  cmsrowHisto.Fill( crow );

	}

	// cuts:

	if( fiducial && liso ) {

	  // for dx:

	  if( abs(cmsdy) < ycut ) {

	    cmsdxfcHisto.Fill( cmsdx );

	    // profiles for alignment:

	    if( lq ) {
	      cmsdxvsx.Fill( x4, cmsdx ); //cmsdxvsx->Fit("pol1","","",-3.5,3.5)
	      cmsdxvsy.Fill( y4, cmsdx );
	      cmsdxvstx.Fill( slx, cmsdx );
	    } // lq

	  } // dy

	  // for dy:

	  if( abs(cmsdx) < xcut ) {

	    cmsdyfcHisto.Fill( cmsdy );
	    cmsrmsyvsq.Fill( Q0, fabs(cmsdy) ); //resolution vs charge

	    if(      Q0 < 16 ) // broken cluster
	      cmsdyfcq0Histo.Fill( cmsdy );

	    else if( Q0 > 35 )  // Landau tail
	      cmsdyfcq9Histo.Fill( cmsdy );

	    else { // 16..35

	      cmsdyfcq1Histo.Fill( cmsdy );

	      if( Q0 < 30 ) {
		cmsdyfcq2Histo.Fill( cmsdy );
	      }

	      if( Q0 < 25 ) {

		cmsdyfcq3Histo.Fill( cmsdy );

		if( ldot )
		  cmsdyfcq3dotHisto.Fill( cmsdy );
		else {
		  cmsdyfcq3nodHisto.Fill( cmsdy );
		  if( oddcol )
		    cmsdyfcq3oddHisto.Fill( cmsdy );
		  else
		    cmsdyfcq3eveHisto.Fill( cmsdy );
		}

		if( ymid ) {
		  cmsdyfcq3midHisto.Fill( cmsdy );
		}
		else
		  cmsdyfcq3rimHisto.Fill( cmsdy );

	      }

	      if( Q0 > 17 && Q0 < 23 && !ldot) {

		cmsdyfcq4nodHisto.Fill( cmsdy );

		if( oddcol )
		  cmsdyfcq4oddHisto.Fill( cmsdy );
		else
		  cmsdyfcq4eveHisto.Fill( cmsdy );

	      }

	    } // Q0 > 16

	    if( lq ) {
	      cmsdyvsx.Fill( x4, cmsdy );
	      cmsdyvsy.Fill( y4, cmsdy );
	      cmsdyvsty.Fill( sly, cmsdy );
	    } // lq

	  } // dx

	  // cut in x and y:

	  if( abs(cmsdx) < xcut  && abs(cmsdy) < ycut ) {

	    cmsq0fHisto.Fill( Q0 ); // Landau

	    if( ldot ) { // sensor bias dot
	      cmsq0fdotHisto.Fill( Q0 );
	    }
	    else { // no dot

	      cmsq0fnodHisto.Fill( Q0 );

	      if( oddcol )
		cmsq0foddHisto.Fill( Q0 ); // more
	      else
		cmsq0feveHisto.Fill( Q0 ); // less

	      if( lcore ) { // pixel core
		cmsq0fcHisto.Fill( Q0 );
	      }
	      else {
		cmsq0fpHisto.Fill( Q0 );
	      } // not core

	    } // no dot

	    // cluster charge profiles, exponential weighting Qx:

	    cmsqxvsx.Fill( x4, Qx );
	    cmsqxvsy.Fill( y4, Qx );
	    cmsqxvsxm.Fill( xmod, Qx ); // Q within pixel
	    cmsqxvsym.Fill( ymod, Qx ); // Q within pixel

	    if( ! ldot )
	      cmsqxvsymn.Fill( ymod, Qx ); // q within pixel

	    if( ( xmod > 120 && xmod < 140 ) ||
		( xmod > 160 && xmod < 180 ) ) // dot
	      cmsqxvsymd.Fill( ymod, Qx ); // q within pixel

	    if( ( ymod >  60 && ymod <  80 ) ||
		( ymod > 160 && ymod < 180 ) ) // dot
	      cmsqxvsxmd.Fill( xmod, Qx ); // q within pixel

	    cmsqxvsxmym.Fill( xmod, ymod, Qx ); // cluster charge profile

	    cmsqxvst1.Fill( evsec, Qx ); // cluster charge vs time
	    cmsqxvst2.Fill( evsec, Qx ); // cluster charge vs time
	    cmsqxvst3.Fill( evsec, Qx ); // cluster charge vs time
	    cmsqxvst4.Fill( evsec, Qx ); // cluster charge vs time

	    if( lq ) {
	      cmsnpxHisto.Fill( c->size );
	      cmsncolHisto.Fill( ncol );
	      cmsnrowHisto.Fill( nrow );
	    }
	    cmsnpxvsq.Fill( Q0, c->size ); // rising

	    for( vector<pixel>::iterator px = c->vpix.begin(); px != c->vpix.end(); ++px ) {

	      cmspxqHisto.Fill( px->q );
	      if( oddcol )
		cmspxqoddHisto.Fill( px->q );
	      else
		cmspxqeveHisto.Fill( px->q );
	      /*
	      if( px->ord == 0 )
		cmspxq1stHisto.Fill( px->q );
	      else
		cmspxq2ndHisto.Fill( px->q );
	      */
	      if( c->size == 1 ) cmspxq1Histo.Fill( px->q );
	      if( c->size == 2 ) cmspxq2Histo.Fill( px->q ); // flat
	      if( c->size >= 3 ) cmspxq3Histo.Fill( px->q ); // hump at low q
	      cmspxqvsq.Fill( Q0, px->q );
	      cmspxqvsxm.Fill( xmod, px->q );
	      cmspxqvsym.Fill( ymod, px->q );

	    } // pix

	    if( lq ) { // Landau peak

	      cmsdxvsxm.Fill( xmod, cmsdx );
	      cmsdxvsym.Fill( ymod, cmsdx );
	      cmsdyvsxm.Fill( xmod, cmsdy );
	      if( ! ldot ) {
		cmsdyvsym.Fill( ymod, cmsdy );
		if( oddcol ) 
		  cmsdyvsymodd.Fill( ymod, cmsdy );
		else
		  cmsdyvsymeve.Fill( ymod, cmsdy );
	      }

	      cmsrmsxvsx.Fill( x4, fabs(cmsdx) ); //resolution across cols
	      cmsrmsyvsx.Fill( x4, fabs(cmsdy) ); //resolution across cols
	      cmsrmsxvsy.Fill( y4, fabs(cmsdx) ); //resolution across rows
	      cmsrmsyvsy.Fill( y4, fabs(cmsdy) ); //resolution across rows

	      cmsrmsxvsxm.Fill( xmod, fabs(cmsdx) ); //resolution within pixel
	      cmsrmsyvsxm.Fill( xmod, fabs(cmsdy) ); //resolution within pixel
	      if( !ldot ) {
		cmsrmsxvsym.Fill( ymod, fabs(cmsdx) ); //resolution within pixel
		cmsrmsyvsym.Fill( ymod, fabs(cmsdy) ); //resolution within pixel
	      }
	      cmsrmsyvsxmym.Fill( xmod, ymod, fabs(cmsdy) ); //resolution within pixel
	      cmsrmsxvsxmym.Fill( xmod, ymod, fabs(cmsdx) ); //resolution within pixel

	      cmsncolvsx.Fill( x4, ncol ); // no trend
	      cmsncolvsy.Fill( y4, ncol ); // 

	      cmsncolvsxm.Fill( xmod, ncol );
	      cmsnrowvsxm.Fill( xmod, nrow );
	      if( ( xmod > 0.020 && xmod < 0.130 ) ||
		  ( xmod > 0.170 && xmod < 0.280 ) ) {
		cmsncolvsym.Fill( ymod, ncol ); // within pixel
		cmsnrowvsym.Fill( ymod, nrow ); // within pixel
	      }
	      cmsnpxvsxmym.Fill( xmod, ymod, c->size ); // cluster size map

	    } // q Landau peak

	  } // dx and dy

	} // fiducial

	// wide linking cut in x and y for efficiency:

	if( abs(cmsdx) < xcut  && abs(cmsdy) < ycut ) {
	  nm = 1;
	  trixylkHisto.Fill( xA, yA );
	}

      } // loop DUT clusters

      // DUT efficiency vs isolated REF-linked fiducial tracks:

      if( liso && lsixlk ) {
	sixxylkHisto.Fill( xA, yA );
	if( nm ) sixxyeffHisto.Fill( xA, yA );
	if( abs( y4 ) < 3.9 )
	  effvsx.Fill( x4, nm );
	if( abs( x4 ) < 3.9 )
	  effvsy.Fill( y4, nm );
      }

    } // loop triplets iA

    ++event_nr;

  } while( reader->NextEvent() && event_nr < lev );

  cout << "done after " << event_nr << " events" << endl;
  histoFile->Write();
  histoFile->Close();

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT alignment:

  double newDUTalignx = DUTalignx;
  double newDUTaligny = DUTaligny;

  if( cmsdxaHisto.GetMaximum() > cmssxaHisto.GetMaximum() ) {
    cout << endl << cmsdxaHisto.GetTitle()
	 << " bin " << cmsdxaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, cmsdxaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdxaHisto.GetBinCenter( cmsdxaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, cmsdxaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdxaHisto.GetBinContent(1) ); // BG
    cmsdxaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTalignx = fgp0->GetParameter(1);
  }
  else {
    cout << endl << cmssxaHisto.GetTitle()
	 << " bin " << cmssxaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, cmssxaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmssxaHisto.GetBinCenter( cmssxaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, cmssxaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmssxaHisto.GetBinContent(1) ); // BG
    cmssxaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTalignx = fgp0->GetParameter(1);
  }

  if( cmsdyaHisto.GetMaximum() > cmssyaHisto.GetMaximum() ) {
    cout << endl << cmsdyaHisto.GetTitle()
	 << " bin " << cmsdyaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, cmsdyaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdyaHisto.GetBinCenter( cmsdyaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, cmsdyaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdyaHisto.GetBinContent(1) ); // BG
    cmsdyaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTaligny = fgp0->GetParameter(1);
  }
  else {
    cout << endl << cmssyaHisto.GetTitle()
	 << " bin " << cmssyaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, cmssyaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmssyaHisto.GetBinCenter( cmssyaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, cmssyaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmssyaHisto.GetBinContent(1) ); // BG
    cmssyaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTaligny = fgp0->GetParameter(1);
  }

  // finer alignment:

  if( abs( newDUTalignx - DUTalignx ) < 0.1 ) {

    cout << endl << cmsdxHisto.GetTitle()
	 << " bin " << cmsdxHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, cmsdxHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdxHisto.GetBinCenter( cmsdxHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 8*cmsdxHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdxHisto.GetBinContent(1) ); // BG
    cmsdxHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTalignx = DUTalignx + fgp0->GetParameter(1);

  }

  if( abs( newDUTaligny - DUTaligny ) < 0.1 ) {

    cout << endl << cmsdyHisto.GetTitle()
	 << " bin " << cmsdyHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, cmsdyHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, cmsdyHisto.GetBinCenter( cmsdyHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 5*cmsdyHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, cmsdyHisto.GetBinContent(1) ); // BG
    cmsdyHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newDUTaligny = DUTaligny + fgp0->GetParameter(1);

    // dyvsx -> rot

    cmsdyvsx.Fit( "pol1", "q", "", -midx[iDUT]+0.2, midx[iDUT]-0.2 );
    TF1 * fdyvsx = cmsdyvsx.GetFunction( "pol1" );
    cout << endl << cmsdyvsx.GetTitle()
	 << ": extra rot " << fdyvsx->GetParameter(1) << endl;
    DUTrot += fdyvsx->GetParameter(1);

    // dyvsy -> tilt:

    cmsdyvsy.Fit( "pol1", "q", "", -midy[iDUT]+0.2, midy[iDUT]-0.2 );
    TF1 * fdyvsy = cmsdyvsy.GetFunction( "pol1" );
    cout << endl << cmsdyvsy.GetTitle()
	 << ": slope " << fdyvsy->GetParameter(1)
	 << ", extra tilt " << fdyvsy->GetParameter(1)/wt/max(sa,0.01)
	 << " deg"
	 << endl;
    DUTtilt += fdyvsy->GetParameter(1)/wt/max(sa,0.01); // [deg]

    // dyvsty -> dz:

    cmsdyvsty.Fit( "pol1", "q", "", -0.002, 0.002 );
    TF1 * fdyvsty = cmsdyvsty.GetFunction( "pol1" );
    cout << endl << cmsdyvsty.GetTitle()
	 << ": z shift " << fdyvsty->GetParameter(1)
	 << " mm (subtract from DUTz)"
	 << endl;
    DUTz -= fdyvsty->GetParameter(1);
  }

  // write new DUT alignment:

  ofstream DUTalignFile( DUTalignFileName.str() );

  DUTalignFile << "# DUT alignment for run " << run << endl;
  ++DUTaligniteration;
  DUTalignFile << "iteration " << DUTaligniteration << endl;
  DUTalignFile << "alignx " << newDUTalignx << endl;
  DUTalignFile << "aligny " << newDUTaligny << endl;
  DUTalignFile << "rot " << DUTrot << endl;
  DUTalignFile << "tilt " << DUTtilt << endl;
  DUTalignFile << "turn " << DUTturn << endl;
  DUTalignFile << "dz " << DUTz - zz[2] << endl;

  DUTalignFile.close();

  cout << endl << "wrote DUT alignment iteration " << DUTaligniteration
       << " to " << DUTalignFileName.str()
       << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // REF alignment:

  double newREFalignx = REFalignx;
  double newREFaligny = REFaligny;

  if( refdxaHisto.GetMaximum() > refsxaHisto.GetMaximum() ) {
    cout << endl << refdxaHisto.GetTitle()
	 << " bin " << refdxaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdxaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdxaHisto.GetBinCenter( refdxaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refdxaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdxaHisto.GetBinContent(1) ); // BG
    refdxaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFalignx = fgp0->GetParameter(1);
  }
  else {
    cout << endl << refsxaHisto.GetTitle()
	 << " bin " << refsxaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refsxaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refsxaHisto.GetBinCenter( refsxaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refsxaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refsxaHisto.GetBinContent(1) ); // BG
    refsxaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFalignx = fgp0->GetParameter(1);
  }

  if( refdyaHisto.GetMaximum() > refsyaHisto.GetMaximum() ) {
    cout << endl << refdyaHisto.GetTitle()
	 << " bin " << refdyaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdyaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdyaHisto.GetBinCenter( refdyaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refdyaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdyaHisto.GetBinContent(1) ); // BG
    refdyaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFaligny = fgp0->GetParameter(1);
  }
  else {
    cout << endl << refsyaHisto.GetTitle()
	 << " bin " << refsyaHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refsyaHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refsyaHisto.GetBinCenter( refsyaHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, refsyaHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refsyaHisto.GetBinContent(1) ); // BG
    refsyaHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFaligny = fgp0->GetParameter(1);
  }

  // finer alignment:

  if( abs( newREFalignx - REFalignx ) < 0.1 ) {

    cout << endl << refdxHisto.GetTitle()
	 << " bin " << refdxHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdxHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdxHisto.GetBinCenter( refdxHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 8*refdxHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdxHisto.GetBinContent(1) ); // BG
    refdxHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFalignx = REFalignx + fgp0->GetParameter(1);

  }

  if( abs( newREFaligny - REFaligny ) < 0.1 ) {

    cout << endl << refdyHisto.GetTitle()
	 << " bin " << refdyHisto.GetBinWidth(1)
	 << endl;
    TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
    fgp0->SetParameter( 0, refdyHisto.GetMaximum() ); // amplitude
    fgp0->SetParameter( 1, refdyHisto.GetBinCenter( refdyHisto.GetMaximumBin() ) );
    fgp0->SetParameter( 2, 5*refdyHisto.GetBinWidth(1) ); // sigma
    fgp0->SetParameter( 3, refdyHisto.GetBinContent(1) ); // BG
    refdyHisto.Fit( "fgp0", "q" );
    cout << "Fit Gauss + BG:"
	 << endl << "  A " << fgp0->GetParameter(0)
	 << endl << "mid " << fgp0->GetParameter(1)
	 << endl << "sig " << fgp0->GetParameter(2)
	 << endl << " BG " << fgp0->GetParameter(3)
	 << endl;
    newREFaligny = REFaligny + fgp0->GetParameter(1);

    // dyvsx -> rot

    refdyvsx.Fit( "pol1", "q", "", -midx[iREF]+0.2, midx[iREF]-0.2 );
    TF1 * fdyvsx = refdyvsx.GetFunction( "pol1" );
    cout << endl << refdyvsx.GetTitle()
	 << ": extra rot " << fdyvsx->GetParameter(1) << endl;
    REFrot += fdyvsx->GetParameter(1);

    // dyvsty -> dz:

    refdyvsty.Fit( "pol1", "q", "", -0.002, 0.002 );
    TF1 * fdyvsty = refdyvsty.GetFunction( "pol1" );
    cout << endl << refdyvsty.GetTitle()
	 << ": z shift " << fdyvsty->GetParameter(1)
	 << " mm (add to REFz)"
	 << endl;
    REFz += fdyvsty->GetParameter(1);
  }

  // write new REF alignment:

  ofstream REFalignFile( REFalignFileName.str() );

  REFalignFile << "# REF alignment for run " << run << endl;
  ++REFaligniteration;
  REFalignFile << "iteration " << REFaligniteration << endl;
  REFalignFile << "alignx " << newREFalignx << endl;
  REFalignFile << "aligny " << newREFaligny << endl;
  REFalignFile << "rot " << REFrot << endl;
  REFalignFile << "dz " << REFz - zz[5] << endl;

  REFalignFile.close();

  cout << endl << "wrote REF alignment iteration " << REFaligniteration
       << " to " << REFalignFileName.str()
       << endl;

  cout << endl << histoFile->GetName() << endl;

  cout << endl;

  return 0;
}
