
// Daniel Pitzl, DESY, Apr 2016
// telescope event display using ROOT

// make evd
// evd -l 9 24500

#include "eudaq/FileReader.hh"
#include "eudaq/PluginManager.hh"

#include <TApplication.h>
#include <TGClient.h>
#include <TGFrame.h>
#include <TRootEmbeddedCanvas.h>
#include <TCanvas.h>
#include <TStyle.h>
#include <TGraph.h>
#include <TAxis.h>
#include <TFile.h>
#include <TH1I.h> // counting
#include <TH2I.h>
#include <TF1.h>
#include <TLine.h>

#include <sstream> // stringstream
#include <fstream> // filestream
#include <set>
#include <cmath> // fabs
#include <unistd.h> // usleep

using namespace std;
using namespace eudaq;

class MyMainFrame:public TGMainFrame
{
private:
  TGMainFrame * fMain;
  TRootEmbeddedCanvas * fEcanvas1;
  TRootEmbeddedCanvas * fEcanvas2;
public:
  MyMainFrame( const TGWindow * p, UInt_t w, UInt_t h );
  //virtual ~MyMainFrame(  );// crash
  ~MyMainFrame(  );
  TCanvas *GetCanvas1(  );
  TCanvas *GetCanvas2(  );
};

struct pixel {
  int col;
  int row;
  int adc;
  double q;
  int ord;
  bool big;
};

struct cluster {
  vector <pixel> vpix; // Armin Burgmeier: list
  int size;
  int ncol, nrow;
  double col, row;
  double charge;
  bool big;
};

struct triplet {
  double xm;
  double ym;
  double zm;
  double sx;
  double sy;
  bool lk;
  double ttdmin;
};

// globals:

pixel pb[999]; // global declaration: array of pixel hits
int fNHit; // global

//------------------------------------------------------------------------------
MyMainFrame::MyMainFrame( const TGWindow * p, UInt_t w, UInt_t h )
//  :TGMainFrame( p, w, h )
{
  cout << "MyMainFrame..." << endl;
  // Create a main frame:
  //fMain = new TGMainFrame( p, w, h );
  fMain = new TGMainFrame( p, w, h, kVerticalFrame );

  fMain->SetWMPosition( 99, 0 ); // no effect

  // Create canvas widget:
  fEcanvas1 = new TRootEmbeddedCanvas( "Ecanvas", fMain, w, h/2 );
  fEcanvas2 = new TRootEmbeddedCanvas( "Ecanvas", fMain, w, h/2 );

  fMain->AddFrame( fEcanvas1,
                   new TGLayoutHints( kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1 ) );
  fMain->AddFrame( fEcanvas2,
                   new TGLayoutHints( kLHintsExpandX | kLHintsExpandY, 1, 1, 1, 1 ) );

  // Set a name to the main frame:
  fMain->SetWindowName( "tele_scope event display" );

  // Map all subwindows of main frame:
  fMain->MapSubwindows(  );

  // Initialize the layout algorithm:
  fMain->Resize( fMain->GetDefaultSize(  ) );

  // Map main frame:
  fMain->MapWindow(  );
}

MyMainFrame::~MyMainFrame(  )
{
  // Clean up used widgets: frames, buttons, layouthints
  fMain->Cleanup(  );
  cout << "MyMainFrame: Cleanup" << endl;
  delete fMain;
  //delete fEcanvas; // crash
}

TCanvas *MyMainFrame::GetCanvas1(  )
{
  return ( fEcanvas1->GetCanvas(  ) );
}

TCanvas *MyMainFrame::GetCanvas2(  )
{
  return ( fEcanvas2->GetCanvas(  ) );
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
    c.big = 0;
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
      if( p->big ) c.big = 1;
      if( p->col > maxx ) maxx = p->col;
      if( p->col < minx ) minx = p->col;
      if( p->row > maxy ) maxy = p->row;
      if( p->row < miny ) miny = p->row;
    }

    //cout << "(cluster with " << c.vpix.size() << " pixels)" << endl;

    if( ! ( c.charge == 0 ) ) {
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

  int lev = 9; // last event displayed

  bool syncdut = 0; // re-sync required ?
  bool syncref = 0; // re-sync required ?
  bool syncmod = 0; // re-sync required ?

  double thr = 0; // offline pixel threshold [ke]

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] ); // last event

    if( !strcmp( argv[i], "-s" ) ) {
      syncdut = 1;
      syncref = 1;
      syncmod = 1;
    }

    if( !strcmp( argv[i], "-d" ) )
      syncdut = 1;

    if( !strcmp( argv[i], "-r" ) )
      syncref = 1;

  } // argc

  if( syncdut )
    cout << "re-sync DUT" << endl;
  if( syncref )
    cout << "re-sync REF" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // runs.dat:

  cout << endl;

  string geoFileName( "geo.dat" );
  double DUTtilt0 = 19.3;
  double MODtilt0 = 16.5;
  double MODturn0 = 27.8;
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

  } // runsFile

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
  // hot pixels:

  ostringstream hotFileName; // output string stream

  hotFileName << "hot_" << run << ".dat";

  ifstream ihotFile( hotFileName.str() );

  set <int> hotset[9];

  if( ihotFile.bad() || ! ihotFile.is_open() ) {
    cout << "no " << hotFileName.str() << " (created by tele)" << endl;
  }
  else {

    cout << "read hot pixel list from " << hotFileName.str() << endl;

    string hash( "#" );
    string plane( "plane" );
    string pix( "pix" );

    int ipl = 0;

    while( ! ihotFile.eof() ) {

      string line;
      getline( ihotFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == plane )
	tokenizer >> ipl;

      if( ipl < 0 || ipl >= 6 ) {
	//cout << "wrong plane number " << ipl << endl;
	continue;
      }

      if( tag == pix ) {
	int ix, iy;
	tokenizer >> ix;
	tokenizer >> iy;
	int ipx = ix*ny[ipl]+iy;
	hotset[ipl].insert(ipx);
      }

    } // while getline

  } // hotFile

  ihotFile.close();

  for( int ipl = 0; ipl < 6; ++ipl )
    cout << ipl << ": hot " << hotset[ipl].size() << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // DUT:

  double wt = atan(1.0) / 45.0; // pi/180 deg

  bool rot90 = 0; // 504
  //bool rot90 = 1; // 506

  double upsign = 1; // 504

  double upsignx = upsign;
  double upsigny = upsign;
  if( chip0 == 603 && run > 24290 ) // tilt 180 deg
    upsigny = -upsign;
  
  if( rot90 ) {
    upsignx =  1;
    upsigny = -1;
  }

  int iDUT = 7;

  int DUTaligniteration = 0;
  double DUTalignx = 0.0;
  double DUTaligny = 0.0;
  double DUTrot = 0.0;
  double DUTturn = 0;
  double DUTtilt = DUTtilt0; // [deg]
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

  zz[iDUT] = DUTz;
  alignx[iDUT] = DUTalignx;
  aligny[iDUT] = DUTaligny;
  rotx[iDUT] = DUTrot;
  roty[iDUT] = DUTrot;

  // normal vector on DUT surface:
  // N = ( 0, 0, -1 ) on DUT, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double co = cos( DUTturn*wt );
  const double so = sin( DUTturn*wt );
  const double ca = cos( DUTtilt*wt );
  const double sa = sin( DUTtilt*wt );

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

  zz[iREF] = REFz;
  alignx[iREF] = REFalignx;
  aligny[iREF] = REFaligny;
  rotx[iREF] = REFrot;
  roty[iREF] = REFrot;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // MOD:

  int iMOD = 6;

  int MODaligniteration = 0;
  double MODalignx = 0.0;
  double MODaligny = 0.0;
  double MODrot = 0.0;
  double MODtilt = MODtilt0; // [deg]
  double MODturn = MODturn0; // [deg]
  double MODz = 47 + zz[4];

  ostringstream MODalignFileName; // output string stream

  MODalignFileName << "alignMOD_" << run << ".dat";

  ifstream iMODalignFile( MODalignFileName.str() );

  cout << endl;

  if( iMODalignFile.bad() || ! iMODalignFile.is_open() ) {
    cout << "no " << MODalignFileName.str() << ", will bootstrap" << endl;
  }
  else {

    cout << "read MODalignment from " << MODalignFileName.str() << endl;

    string hash( "#" );
    string iteration( "iteration" );
    string alignx( "alignx" );
    string aligny( "aligny" );
    string rot( "rot" );
    string tilt( "tilt" );
    string turn( "turn" );
    string dz( "dz" );

    while( ! iMODalignFile.eof() ) {

      string line;
      getline( iMODalignFile, line );
      cout << line << endl;

      if( line.empty() ) continue;

      stringstream tokenizer( line );
      string tag;
      tokenizer >> tag; // leading white space is suppressed
      if( tag.substr(0,1) == hash ) // comments start with #
	continue;

      if( tag == iteration ) 
	tokenizer >> MODaligniteration;

      double val;
      tokenizer >> val;
      if(      tag == alignx )
	MODalignx = val;
      else if( tag == aligny )
	MODaligny = val;
      else if( tag == rot )
	MODrot = val;
      else if( tag == tilt )
	MODtilt = val;
      else if( tag == turn )
	MODturn = val;
      else if( tag == dz )
	MODz = val + zz[4];

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  iMODalignFile.close();

  zz[iMOD] = MODz;
  alignx[iMOD] = MODalignx;
  aligny[iMOD] = MODaligny;
  rotx[iMOD] = MODrot;
  roty[iMOD] = MODrot;

  // normal vector on MOD surface:
  // N = ( 0, 0, -1 ) on MOD, towards -z
  // transform into tele system:
  // tilt alpha around x
  // turn omega around y

  const double com = cos( MODturn*wt );
  const double som = sin( MODturn*wt );
  const double cam = cos( MODtilt*wt );
  const double sam = sin( MODtilt*wt );
  const double cfm = cos( MODrot );
  const double sfm = sin( MODrot );

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

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // ROOT:

  cout << "ROOT application..." << endl;

  TApplication theApp( "comet", &argc, argv );

  gStyle->SetTextFont( 62 ); // 62 = Helvetica bold
  gStyle->SetTextAlign( 11 );

  gStyle->SetTickLength( -0.02, "x" ); // tick marks outside
  gStyle->SetTickLength( -0.01, "y" );
  gStyle->SetTickLength( -0.01, "z" );

  gStyle->SetLabelOffset( 0.022, "x" );
  gStyle->SetLabelOffset( 0.013, "y" );
  gStyle->SetLabelOffset( 0.022, "z" );

  gStyle->SetTitleOffset( 1.6, "x" );
  gStyle->SetTitleOffset( 1.6, "y" );
  gStyle->SetTitleOffset( 1.7, "z" );

  gStyle->SetLabelFont( 62, "X" );
  gStyle->SetLabelFont( 62, "Y" );
  gStyle->SetLabelFont( 62, "z" );

  gStyle->SetTitleFont( 62, "X" );
  gStyle->SetTitleFont( 62, "Y" );
  gStyle->SetTitleFont( 62, "z" );

  gStyle->SetTitleBorderSize( 0 ); // no frame around global title
  gStyle->SetTitleAlign( 13 ); // 13 = left top align
  gStyle->SetTitleX( 0.12 ); // global title
  gStyle->SetTitleY( 0.98 ); // global title

  gStyle->SetLineWidth( 1 ); // frames
  gStyle->SetHistLineColor( 4 ); // 4=blau
  gStyle->SetHistLineWidth( 3 );
  gStyle->SetHistFillColor( 5 ); // 5=gelb
  //  gStyle->SetHistFillStyle(4050); // 4050 = half transparent
  gStyle->SetHistFillStyle( 1001 ); // 1001 = solid

  gStyle->SetFrameLineWidth( 2 );

  // statistics box:

  gStyle->SetOptStat( 0 );
  gStyle->SetStatFormat( "8.6g" ); // more digits, default is 6.4g
  gStyle->SetStatFont( 42 ); // 42 = Helvetica normal
  //  gStyle->SetStatFont(62); // 62 = Helvetica bold
  gStyle->SetStatBorderSize( 1 ); // no 'shadow'

  gStyle->SetStatX( 0.80 );
  gStyle->SetStatY( 0.95 );

  gStyle->SetPalette( 1 ); // rainbow colors

  gStyle->SetHistMinimumZero(  ); // no zero suppression

  gStyle->SetOptDate( 0 );

  cout << "open ROOT window..." << endl;
  MyMainFrame *myMF = new
    MyMainFrame( gClient->GetRoot(  ), 1400, 800 );

  cout << "open Canvas..." << endl;

  TCanvas * c1 = myMF->GetCanvas1(  );
  TCanvas * c2 = myMF->GetCanvas2(  );

  c1->SetBottomMargin( 0.08 );
  c1->SetLeftMargin( 0.03 );
  c1->SetRightMargin( 0.10 );

  c2->SetBottomMargin( 0.08 );
  c2->SetLeftMargin( 0.03 );
  c2->SetRightMargin( 0.10 );

  gPad->Update(  ); // required

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  cout << endl;

  FileReader * reader;
  if(      run <    100 )
    reader = new FileReader( runnum.c_str(), "data/run0000$2R$X" );
  else if( run <   1000 )
    reader = new FileReader( runnum.c_str(), "data/run000$3R$X" );
  else if( run <  10000 )
    reader = new FileReader( runnum.c_str(), "data/run00$4R$X" );
  else if( run < 100000 )
    reader = new FileReader( runnum.c_str(), "data/run0$5R$X" );
  else
    reader = new FileReader( runnum.c_str(), "data/run$6R$X" );

  int event_nr = 0;
  int nevd = 0;
  uint64_t evTLU0 = 0;
  const double fTLU = 384E6; // 384 MHz TLU clock

  vector <cluster> cl0[9]; // remember from previous event

  do {
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ) {
      eudaq::PluginManager::Initialize(evt);
      continue;
    }

    ++event_nr;

    bool ldbg = 0;

    if( event_nr == 1 )
      ldbg = 1;

    if( lev < 100 )
      ldbg = 1;

    uint64_t evTLU = evt.GetTimestamp(); // 384 MHz = 2.6 ns
    if( event_nr < 2  )
      evTLU0 = evTLU;
    double evsec = (evTLU - evTLU0) / fTLU;

    if( event_nr < 10 )
      cout << "scope processing  " << event_nr << "  taken " << evsec << endl;
    else if( event_nr < 100 && event_nr%10 == 0 )
      cout << "scope processing  " << event_nr << "  taken " << evsec << endl;
    else if( event_nr < 1000 && event_nr%100 == 0 )
      cout << "scope processing  " << event_nr << "  taken " << evsec << endl;
    else if( event_nr%1000 == 0 )
      cout << "scope processing  " << event_nr << "  taken " << evsec << endl;

    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);

    vector <cluster> cl[9];
    int ncl = 0;

    for( size_t iplane = 0; iplane < sevt.NumPlanes(); ++iplane ) {

      const eudaq::StandardPlane &plane = sevt.GetPlane(iplane);

      std::vector<double> pxl = plane.GetPixels<double>();

      if( ldbg ) cout << "plane " << plane.ID() << ": ";

      // /home/pitzl/eudaq/main/include/eudaq/CMSPixelHelper.hh

      int ipl = plane.ID();

      int npx = 0;

      for( size_t ipix = 0; ipix < pxl.size(); ++ipix ) {

	if( ldbg ) 
	  cout << plane.GetX(ipix)
	       << " " << plane.GetY(ipix)
	       << " " << plane.GetPixel(ipix) << " ";

	int ix = plane.GetX(ipix); // global column 0..415
	int iy = plane.GetY(ipix); // global row 0..159
	int adc = plane.GetPixel(ipix); // ADC 0..255

	// skip hot pixels:

	int ipx = ix*ny[ipl] + iy;
	if( hotset[ipl].count(ipx) ) continue;

	double q = adc;

	if( ipl == iDUT &&
	    adc > 0 &&
	    ix >= 0 && ix < 52 &&
	    iy >= 0 && iy < 80 ) {

	  double Ared = adc - p4[ix][iy]; // p4 is asymptotic maximum

	  if( Ared >= 0 )
	    Ared = -0.1; // avoid overflow

	  double a3 = p3[ix][iy]; // positive
	  if( weib == 3 )
	    q = p1[ix][iy] *
	      ( pow( -log( -Ared / a3 ), 1/p2[ix][iy] ) - p0[ix][iy] ) * ke;
	  // q = ( (-ln(-(A-p4)/p3))^1/p2 - p0 )*p1

	  if( q < thr ) continue; // offline threshold [ke]

	}  // DUT

	if( ipl == iREF &&
	    adc > 0 &&
	    ix >= 0 && ix < 52 &&
	    iy >= 0 && iy < 80 ) {

	  double Ared = adc - r4[ix][iy]; // r4 is asymptotic maximum

	  if( Ared >= 0 )
	    Ared = -0.1; // avoid overflow

	  double a3 = r3[ix][iy]; // positive
	  if( weib == 3 )
	    q = r1[ix][iy] *
	      ( pow( -log( -Ared / a3 ), 1/r2[ix][iy] ) - r0[ix][iy] ) * refke;
	  // q = ( (-ln(-(A-r4)/r3))^1/r2 - r0 )*r1

	}  // REF

	int xm = ix;
	int ym = iy;

	if( ipl == iMOD ) {

	  // leave space for big pixels:

	  int roc = ix / 52; // 0..7
	  int col = ix % 52; // 0..51
	  int row = iy;
	  int x = 1 + ix + 2*roc; // 1..52 per ROC with big pix
	  int y = iy;
	  if( iy > 79 ) y += 2;

	  // flip for upper ROCs into local addresses:

	  if( iy > 79 ) {
	    roc = 15 - roc; // 15..8
	    col = 51 - col; // 51..0
	    row = 159 - iy; // 79..0
	  }

	  ix = x;
	  iy = y;

	} // MOD

	// fill pixel block for clustering:

	pb[npx].col = ix;
	pb[npx].row = iy;
	pb[npx].adc = adc;
	pb[npx].q = q;
	pb[npx].ord = npx; // readout order
	pb[npx].big = 0;
	++npx;

	if( ipl == iMOD ) {

	  // double big pixels:
	  // 0+1
	  // 2..51
	  // 52+53

	  int col = xm % 52; // 0..51

	  if( col == 0 ) {
	    pb[npx].col = ix-1; // double
	    pb[npx].row = iy;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	  if( col == 51 ) {
	    pb[npx].col = ix+1; // double
	    pb[npx].row = iy;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	  if( ym == 79 ) {
	    pb[npx].col = ix; // double
	    pb[npx].row = 80;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	  if( ym == 80 ) {
	    pb[npx].col = ix; // double
	    pb[npx].row = 81;
	    pb[npx-1].adc *= 0.5;
	    pb[npx-1].q *= 0.5;
	    pb[npx].adc = 0.5*adc;
	    pb[npx].q = 0.5*q;
	    pb[npx].big = 1;
	    ++npx;
	  }

	} // MOD

	if( npx > 994 ) {
	  cout << "pixel buffer overflow in plane " << ipl
	       << ", event " << event_nr
	       << endl;
	  break;
	}

      } // pix

      if( ldbg ) cout << endl;

      // clustering:

      fNHit = npx; // for cluster search

      cl[ipl] = getClus();

      if( ldbg ) cout << "  z " << int(zz[ipl]+0.5)
		      << ", clusters " << cl[ipl].size()
		      << endl;

      ncl += cl[ipl].size();

    } // planes

    if( ! syncdut )
      cl0[iDUT] = cl[iDUT];
    if( ! syncref )
      cl0[iREF] = cl[iREF];
    if( ! syncmod )
      cl0[iMOD] = cl[iMOD];

    // event selection:

    if( ncl < 6 ) continue;
    if( cl0[iDUT].size() < 1 ) continue;
    if( cl0[iREF].size() < 1 ) continue;

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

	if( fabs( dx2 ) > 0.010 * dz35 ) continue; // angle cut
	if( fabs( dy2 ) > 0.010 * dz35 ) continue; // angle cut

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
	  dri.ttdmin = 99.9; // isolation [mm]
	  driplets.push_back(dri);

	} // cl B

      } // cl C

    } // cl A

    // event selection:
    //if( driplets.size() > 4 ) continue;
    //if( driplets.size() < 8 ) continue;

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

	if( fabs( dx2 ) > 0.010 * dz02 ) continue; // angle cut
	if( fabs( dy2 ) > 0.010 * dz02 ) continue; // angle cut

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

	  // telescope triplet cuts:

	  if( fabs(dx3) > triCut ) continue;
	  if( fabs(dy3) > triCut ) continue;

	  triplet tri;
	  tri.xm = avx;
	  tri.ym = avy;
	  tri.zm = avz;
	  tri.sx = slpx;
	  tri.sy = slpy;
	  tri.lk = 0;
	  tri.ttdmin = 99.9; // isolation [mm]
	  triplets.push_back(tri);

	} // cl B

      } // cl C

    } // cl A

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // book event histos:

    ++nevd;

    cout << "display " << nevd
	 << " with " << ncl << " clusters, "
	 << triplets.size() << " triplets, "
	 << driplets.size() << " driplets"
	 << endl;

    TH2D xzview( Form( "ev%ixz", nevd ),
		 Form( "display %i trigger %i x-z;z [mm];x [mm];hit", nevd, (int) event_nr ),
		 890, -10, 880, 220, -11, 11 );

    c1->cd();
    xzview.Draw("");

    TH2D yzview( Form( "ev%iyz", nevd ),
		 Form( "display %i trigger %i y-z;z [mm];y [mm];hit", nevd, (int) event_nr ),
		 890, -10, 880, 220, -11, 11 );

    c2->cd();
    yzview.Draw("");

    vector <double> vx(ncl);
    vector <double> vy(ncl);
    vector <double> vz(ncl);

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // Mimosa hits:

    int icl = 0;

    for( int ipl = 0; ipl < 6; ++ipl ) {

      for( vector<cluster>::iterator cA = cl[ipl].begin(); cA != cl[ipl].end(); ++cA ) {

	double xA = cA->col*ptchx[ipl] - alignx[ipl];
	double yA = cA->row*ptchy[ipl] - aligny[ipl];
	double xmid = xA - midx[ipl];
	double ymid = yA - midy[ipl];
	xA = xmid - ymid*rotx[ipl];
	yA = ymid + xmid*roty[ipl];
	double zA = zz[ipl];
	vx.at(icl) = xA;
	vy.at(icl) = yA;
	vz.at(icl) = zA;
	++icl;

      } // cl

    } // planes

    c1->cd();

    TGraph gxz( icl, &vz[0], &vx[0] );
    gxz.SetMarkerColor(4);
    gxz.SetMarkerStyle(20);
    gxz.SetMarkerSize(0.8);
    gxz.Draw("P"); // without axis option: overlay

    c2->cd();

    TGraph gyz( icl, &vz[0], &vy[0] );
    gyz.SetMarkerColor(4);
    gyz.SetMarkerStyle(20);
    gyz.SetMarkerSize(0.8);
    gyz.Draw("P"); // without axis option: overlay

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // DUT hits:

    icl = 0;

    for( vector<cluster>::iterator cA = cl0[iDUT].begin(); cA != cl0[iDUT].end(); ++cA ) {

      double xA = cA->col*ptchx[iDUT] - alignx[iDUT];
      double yA = cA->row*ptchy[iDUT] - aligny[iDUT];
      double xmid = xA - midx[iDUT];
      double ymid = yA - midy[iDUT];
      xmid = upsignx*xmid;
      ymid =-upsigny*ymid;
      xA = xmid - ymid*rotx[iDUT];
      yA = ymid + xmid*roty[iDUT];
      double zA = zz[iDUT];

      double x2 = xA;
      double y2 = ca*yA; // tilt -a
      double z2 = sa*yA;

      xA = co*x2 + so*z2; // turn -o
      yA = y2;
      zA+=-so*x2 + co*z2;

      vx.at(icl) = xA;
      vy.at(icl) = yA;
      vz.at(icl) = zA;
      ++icl;

    } // cl

    TGraph * DUTgxz;
    TGraph * DUTgyz;

    if( icl > 0 ) {
      c1->cd();
      DUTgxz = new TGraph( icl, &vz[0], &vx[0] );
      DUTgxz->SetMarkerColor(2);
      DUTgxz->SetMarkerStyle(21);
      DUTgxz->SetMarkerSize(0.8);
      DUTgxz->Draw("P"); // without axis option: overlay

      c2->cd();
      DUTgyz = new TGraph( icl, &vz[0], &vy[0] );
      DUTgyz->SetMarkerColor(2);
      DUTgyz->SetMarkerStyle(21);
      DUTgyz->SetMarkerSize(0.8);
      DUTgyz->Draw("P"); // without axis option: overlay
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // MOD hits:

    icl = 0;

    for( vector<cluster>::iterator cA = cl0[iMOD].begin(); cA != cl0[iMOD].end(); ++cA ) {

      double xA = cA->col*ptchx[iMOD] - alignx[iMOD];
      double yA = cA->row*ptchy[iMOD] - aligny[iMOD];

      double xmid = xA - midx[iMOD];
      double ymid = yA - midy[iMOD];

      xmid =-xmid;
      ymid = ymid;

      xA = cfm*xmid - sfm*ymid; // -rot
      yA = sfm*xmid + cfm*ymid;
      double zA = zz[iMOD];

      double x2 = xA;
      double y2 = cam*yA; // tilt -a
      double z2 = sam*yA;

      xA = com*x2 + som*z2; // turn -o
      yA = y2;
      zA+=-som*x2 + com*z2;

      vx.at(icl) = xA;
      vy.at(icl) = yA;
      vz.at(icl) = zA;
      ++icl;

    } // cl

    TGraph * MODgxz;
    TGraph * MODgyz;

    if( icl > 0 ) {

      c1->cd();
      MODgxz = new TGraph( icl, &vz[0], &vx[0] );
      MODgxz->SetMarkerColor(7);
      MODgxz->SetMarkerStyle(21);
      MODgxz->SetMarkerSize(0.8);
      MODgxz->Draw("P"); // without axis option: overlay

      c2->cd();
      MODgyz = new TGraph( icl, &vz[0], &vy[0] );
      MODgyz->SetMarkerColor(7);
      MODgyz->SetMarkerStyle(21);
      MODgyz->SetMarkerSize(0.8);
      MODgyz->Draw("P"); // without axis option: overlay
    }
    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // REF hits:

    icl = 0;

    for( vector<cluster>::iterator cA = cl0[iREF].begin(); cA != cl0[iREF].end(); ++cA ) {

      double xA = cA->col*ptchx[iREF] - alignx[iREF];
      double yA = cA->row*ptchy[iREF] - aligny[iREF];
      double xmid = xA - midx[iREF];
      double ymid = yA - midy[iREF];
      xmid =-xmid;
      ymid = ymid;
      xA = xmid - ymid*rotx[iREF];
      yA = ymid + xmid*roty[iREF];
      double zA = zz[iREF];

      vx.at(icl) = xA;
      vy.at(icl) = yA;
      vz.at(icl) = zA;
      ++icl;

    } // cl

    TGraph * REFgxz;
    TGraph * REFgyz;
    if( icl > 0 ) {
      c1->cd();
      REFgxz = new TGraph( icl, &vz[0], &vx[0] );
      REFgxz->SetMarkerColor(8);
      REFgxz->SetMarkerStyle(21);
      REFgxz->SetMarkerSize(0.8);
      REFgxz->Draw("P"); // without axis option: overlay

      c2->cd();
      REFgyz = new TGraph( icl, &vz[0], &vy[0] );
      REFgyz->SetMarkerColor(8);
      REFgyz->SetMarkerStyle(21);
      REFgyz->SetMarkerSize(0.8);
      REFgyz->Draw("P"); // without axis option: overlay
    }

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // triplets:

    vector <TLine> xlines; // in main scope
    vector <TLine> ylines; // in main scope

    for( unsigned int iA = 0; iA < triplets.size(); ++iA ) { // iA = upstream

      double avx = triplets[iA].xm;
      double avy = triplets[iA].ym;
      double avz = triplets[iA].zm;
      double slx = triplets[iA].sx;
      double sly = triplets[iA].sy;

      double d0 = zz[0] - avz; // z DUT from mid of triplet
      double x0 = avx + slx * d0; // triplet impact point on DUT
      double y0 = avy + sly * d0;

      double dA = DUTz - avz; // z DUT from mid of triplet
      double xA = avx + slx * dA; // triplet impact point on DUT
      double yA = avy + sly * dA;

      TLine lx( zz[0], x0, DUTz, xA );
      lx.SetLineColor(1);
      xlines.push_back(lx);

      TLine ly( zz[0], y0, DUTz, yA );
      ly.SetLineColor(1);
      ylines.push_back(ly);

    } // loop triplets iA

    // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    // driplets:

    for( unsigned int jB = 0; jB < driplets.size(); ++jB ) { // jB = downstream

      double avx = driplets[jB].xm;
      double avy = driplets[jB].ym;
      double avz = driplets[jB].zm;
      double slx = driplets[jB].sx;
      double sly = driplets[jB].sy;

      double dA = DUTz - avz; // z DUT from mid of driplet
      double xA = avx + slx * dA; // driplet at DUT
      double yA = avy + sly * dA;

      double dB = REFz - avz; // z REF from mid of driplet
      double xB = avx + slx * dB; // driplet at REF
      double yB = avy + sly * dB;

      TLine lx( DUTz, xA, REFz, xB );
      lx.SetLineColor(6);
      xlines.push_back(lx);

      TLine ly( DUTz, yA, REFz, yB );
      ly.SetLineColor(6);
      ylines.push_back(ly);

    } // loop driplets jB

    c1->cd();
    for( size_t ii = 0; ii < xlines.size(); ++ii )
      xlines[ii].Draw("same");
    c1->Update();

    c2->cd();
    for( size_t ii = 0; ii < ylines.size(); ++ii )
      ylines[ii].Draw("same");
    c2->Update();

    if( syncdut )
      cl0[iDUT] = cl[iDUT]; // remember for re-sync
    if( syncref )
      cl0[iREF] = cl[iREF];
    if( syncmod )
      cl0[iMOD] = cl[iMOD];

    // keyboard interaction:
    /*
    string input;
    cout << "hit enter for next...";
    getline( cin, input, '\n' );
    */

    usleep( 500000 ); // sleep some micro-seconds

  } while( reader->NextEvent() && nevd < lev );

  cout << "done after " << event_nr << " events" << endl;

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // done

  delete myMF;

  cout << endl;

  return 0;
}
