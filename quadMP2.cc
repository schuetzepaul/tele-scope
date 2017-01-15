
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

#include "lcio.h"
#include "IO/LCReader.h"
#include "EVENT/LCRunHeader.h"
#include "EVENT/LCCollection.h"
#include "EVENT/LCEvent.h"
#include "EVENT/TrackerData.h"

#include "stdio.h"
#include "stdlib.h"
#include <sys/types.h>
#include <sys/stat.h>


using namespace std;
using namespace gbl;
using namespace eudaq;
using namespace lcio;

struct pixel {
  int col;
  int row;
  int adc;
  double cal;
  bool big;
  int roc;
};

struct cluster {
  vector <pixel> vpix;
  int roc;
  int size;
  int sumA; // DP
  double charge;
  double col,row;
  bool big;
  int ncol, nrow;
  double x5,y5,z5;
  double xg,yg,zg;
};

// globals:

pixel pb[66560]; // global declaration: vector of pixels with hit
int fNHit; // global

//------------------------------------------------------------------------------
// inverse decorrelated Weibull PH -> large Vcal DAC
double PHtoVcal( double ph, double a0, double a1, double a2, double a3, double a4, double a5)
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
    bool sameRoc = true;
    int prevRoc = 0;
    int npx = 0;

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
      if( npx > 0 && prevRoc != p->roc) sameRoc = false;
      prevRoc = p->roc;
      npx++;
    }

    if(sameRoc){
      c.roc = prevRoc;
    }else{
      c.roc = -1;
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
TMatrixD Jac5( double ds, double bfield ) // for GBL
{
  /*
    B-field
    track = 
    q/p, x', y', x, y
    0,   1,  2,  3, 4
  */

  // angle change should be DeltaAlpha = dz*B*cos(incAngle)*(e/p) (see manual draft and calculation via a isoscele triangle)
  
  bfield = bfield * 0.3;

  TMatrixD jac(5, 5);
  jac.UnitMatrix();
  jac[3][1] = ds; // x = xp * ds
  jac[4][2] = ds; // y = yp * ds
  jac[1][0] = -1.*bfield*ds;
  jac[3][0] = -0.5*bfield*ds*ds;

  return jac;
}

//------------------------------------------------------------------------------
Double_t landau_gauss_peak(TH1* h);


//------------------------------------------------------------------------------
bool isFiducial( double x, double y)
{
  bool ffiducial = true;
  
  double addBorder = 0.5;
  
  if(y < -(8.1-0.3-0.06-addBorder) || y > (8.1-0.3-0.06-addBorder)
     || x < -(32.4-0.45-0.06-addBorder) || x > (32.4-0.45-0.06-addBorder)) ffiducial = false;
  return ffiducial;
}


//------------------------------------------------------------------------------
int searchBScanParameters(int runnr, double &Bfield){

  ifstream parameterFile( "BScanParameters.dat" );

  cout << endl;
  if( parameterFile.bad() || !parameterFile.is_open() ) {
    cout << "File could not be found." << endl;
    return -1;
  }
  else {

    //    cout << "Reading parameters." << endl;

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
        Bfield = currentBcurrent*0.00095; // convert current to field.
        parameterFile.close();

	cout << "Found run " << runnr << ". B = " << Bfield << " T" << endl;

        return 1;
      }
      
    }
    cout << "Run " << runnr << " not found in list. Please add it." << endl;
    parameterFile.close();
    return -2;
  }
}


//------------------------------------------------------------------------------
int searchRunlist(int runnr, double &momentum, int *modName, bool &CCSuppressed, int &alignmentRun, double &turn){

  ifstream runlistFile( "runlist-quad.dat" );

  cout << endl;
  if( runlistFile.bad() || ! runlistFile.is_open() ) {
    cout << "runlist-quad.dat could not be found." << endl;
    return -1;
  }
  else {

    cout << "Reading runlist." << endl;

    int currentRunnr;
    double currentMomentum;
    int currentModNames[4];
    int currentCCSuppressed = 0;
    double currentTurn;
    alignmentRun = 0;
    
    while( ! runlistFile.eof() ) {

      currentCCSuppressed = 0;

      string line;
      getline( runlistFile, line );
      
      if( line.empty() ) continue;
      if( line.at(0) == '#' ) continue;

      stringstream thisline(line);
      
      thisline >> currentRunnr >> currentMomentum;
      for(int mod = 0; mod < 4; mod++){
	thisline >> currentModNames[mod];
      }
      thisline >> currentCCSuppressed >> alignmentRun >> currentTurn;
      
      if(!(currentRunnr && currentModNames[0] && currentModNames[1] && currentModNames[2] && currentModNames[3])){
	continue; // No correct data in runlist
      }
      
      if(currentRunnr == runnr){
	cout << "Found entry in runlist:" << endl;
	cout << line << endl;
	momentum = currentMomentum;
	for(int mod = 0; mod < 4; mod++){
	  modName[mod] = currentModNames[mod];
	}
	if(alignmentRun == 0) alignmentRun = runnr;
	CCSuppressed = currentCCSuppressed;
	turn = currentTurn;
	runlistFile.close();
	return 1;
      }
      
    }
    cout << "Run " << runnr << " not found in runlist-quad.dat. Please add it." << endl;
    runlistFile.close();
    return -2;
  }
}

//------------------------------------------------------------------------------
int getShifts(int runnr, int * shifts){
  
  ifstream shiftFile( "shiftParameters.dat" );

  if(shiftFile.bad() || !shiftFile.is_open() ){
    cout << "Shift file could not be found!" << endl;
    return -1;
  }
  else{
    
    int currentRunnr;

    while( !shiftFile.eof() ) {
      
      string line;
      getline(shiftFile, line);

      if(line.empty()) continue;
      if(line.at(0) == '#') continue;

      stringstream thisline(line);

      thisline >> currentRunnr;
      for(int mod = 0; mod < 4; mod++){
	thisline >> shifts[mod];
      }
      
      if(currentRunnr == runnr){
	cout << "Found following shifts for this run: " << endl << line <<  endl;
	cout << "However - I will just use different ones." << endl;
	if(runnr >= 2187 && runnr <= 2307){
	  shifts[0] = 0;
	  shifts[1] = 13;
	  shifts[2] = 13;
	  shifts[3] = -1.;
	}
	return 1;
      }

    }
    cout << "Run " << runnr << " not found in shiftFile." << endl;
    shiftFile.close();
    return -2;
  }

}


//------------------------------------------------------------------------------
int linefit(double x0, double y0,double x1, double y1, double x2, double y2, double &resultM, double &resultB){

  double sumxi = x0 + x1 + x2;
  double sumyi = y0 + y1 + y2;
  double sumxiyi = x0*y0 + x1*y1 + x2*y2;
  double sumxixi = x0*x0 + x1*x1 + x2*x2;;
  double sumxixj = sumxi*sumxi;
  double sumxiyj = sumxi*sumyi;

  try{
    resultM = (sumxiyi - sumxiyj/3.) / (sumxixi - sumxixj/3.);
    resultB = (sumyi - resultM*sumxi)/3.;
  }catch(...){
    return 0;
  }

  return 1;
}

//----------------------------------

double residSquaredAt(double m, double b, double x, double y){

  return (m*x+b-y)*(m*x+b-y);

}


int linefit4(double x0, double y0,double x1, double y1, double x2, double y2, double x3, double y3, double &resultM, double &resultB, double &sigmay){

  double sumxi = x0 + x1 + x2 + x3;
  double sumyi = y0 + y1 + y2 + y3;
  double sumxiyi = x0*y0 + x1*y1 + x2*y2 + x3*y3;
  double sumxixi = x0*x0 + x1*x1 + x2*x2 + x3*x3;
  double sumxixj = sumxi*sumxi;
  double sumxiyj = sumxi*sumyi;

  try{
    resultM = (sumxiyi - sumxiyj/4.) / (sumxixi - sumxixj/4.);
    resultB = (sumyi - resultM*sumxi)/4.;
  }catch(...){
    return 0;
  }

  sigmay = (residSquaredAt(resultM,resultB,x0,y0) + residSquaredAt(resultM,resultB,x1,y1) + residSquaredAt(resultM,resultB,x2,y2) + residSquaredAt(resultM,resultB,x3,y3))/(4-2);
  sigmay = TMath::Sqrt(sigmay);

  return 1;
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

  int alignmentRun;

  // Searching for entry in runlist
  
  double p;
  int modName[4];
  bool CCSupressed = false;
  int conversionRun = 0;
  bool writeEfficiency = false;
  bool fitSupressed = false;
  double turnTable = 0;
  
  bool useLCIO = false;

  // further arguments:

  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] );

    // Suppress conversion factor correction
    if( !strcmp( argv[i], "-c" ) )
      {
	if( strchr( argv[++i], '-') == NULL){
	  conversionRun = atoi( argv[i] );
	  cout << "Conversion factors taken from run " << conversionRun << endl;
	}else{
	  --i;
	}
	CCSupressed = true;
      }

    if( !strcmp( argv[i], "-e" ) )
      writeEfficiency = true;

    if( !strcmp( argv[i], "-f" ) )
      fitSupressed = true;

    if( !strcmp( argv[i], "-a" ) )
      useLCIO = true;


  } // argc

  if(searchRunlist(run, p, modName, CCSupressed, alignmentRun, turnTable) < 0){
    exit(0);
  }

  double bfield;
  searchBScanParameters(run, bfield);


  cout << "run " << run << endl;
  FileReader * reader;
  LCReader* lcReader = LCFactory::getInstance()->createLCReader();

  if(useLCIO){
    char lciofilename[200];
    sprintf(lciofilename, "lcio/run%06d-converter.slcio", run);
    lcReader->open(lciofilename);

    LCRunHeader *runHdr = lcReader->readNextRunHeader();

    cout << "  Run : " << runHdr->getRunNumber() << " - "      << runHdr->getDetectorName() << ":  "      << runHdr->getDescription() << endl ;    

  }else{
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
  }


  double dz = 22.5; // [mm] projected z spacing
  if( run > 150 )
    dz = 40.0;
  if( !(run > 42000 && run < 42700) )
    dz = 32.0;


  int aligniteration = 0;
  double alignx[4];
  double aligny[4];
  double tiltA[4];
  double turnA[4];
  double rotA[4];
  double alignz[4];

  for( int ipl = 0; ipl < 4; ++ipl ) {
    alignx[ipl] = 0;
    aligny[ipl] = 0;
    tiltA[ipl] = 0;
    turnA[ipl] = 0;
    rotA[ipl] = 0;
    alignz[ipl] = 0;
  }

  ostringstream alignFileName; // output string stream

  if(!alignmentRun){
    alignFileName << "alignmentMP/align_" << run << ".dat";
  }else{
    alignFileName << "alignmentMP/align_" << alignmentRun << ".dat";
  }

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
    string Rotation( "rot" );
    string Tilt( "tilt" );
    string Turn( "turn" );
    string Alignz( "alignz" );

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
      else if( tag == Rotation )
	rotA[ipl] = val;
      else if( tag == Tilt )
	tiltA[ipl] = val;
      else if( tag == Turn )
	turnA[ipl] = val;
      else if( tag == Alignz )
	alignz[ipl] = val;

      // anything else on the line and in the file gets ignored

    } // while getline

  } // alignFile

  ialignFile.close();


  // Read Millepede file:

  // create directory if necessary:
  stringstream dirname;
  dirname << "run" << run;
  string rundir = dirname.str();

  struct stat sb;

  if (stat(rundir.c_str(), &sb) == 0 && S_ISDIR(sb.st_mode)){
    cout << "Directory exists." << endl;
  }else{
    cout << "Directory does not exist. Creating it." << endl;
    stringstream createdir;
    createdir << "mkdir run" << run;
    system(createdir.str().c_str());
  }

  dirname << "/millepede.res";

  ifstream milleRes( dirname.str() );

  const int ndimMP = 5;

  cout << endl;
  if( milleRes.bad() || ! milleRes.is_open() ) {
    cout << "No " << dirname.str() << " found. Skipping this." << endl;
    cout << endl;
  }
  else {

    cout << "Read Millepede results from " << dirname.str() << endl;

    std::vector<double > tokens;
    std::stringstream tokenizer;
    std::string line;

    getline(milleRes, line);

    unsigned int numpars = 4*ndimMP; // 2* dx, dy, drot, dtilt, dturn, dz
    std::map< unsigned int, double > alpar; // map = associative array

    for( unsigned int ipar = 0; ipar < numpars; ++ipar ) {

      if( milleRes.eof() ) break;
      getline( milleRes, line );
      if( line.empty() ) continue;

      tokens.clear();
      tokenizer.clear();
      tokenizer.str( line );

      double buffer;

      while( tokenizer >> buffer ) tokens.push_back( buffer );
      int lpar = static_cast< int >( tokens[0] + 0.5 ); // par label

      alpar[lpar] = tokens[1];

    }

    bool showAlPars[4] = {true, false, true, false};

    for(size_t mod = 0; mod < 4; mod++){
      if(showAlPars[mod]){
	cout << setprecision(6) << "MP alignment corrections plane " << mod << ":" << endl;
	cout <<  "\tdx    =  " << alpar[mod*ndimMP+1]*1E3 << " um" << endl;
	cout <<  "\tdy    =  " << alpar[mod*ndimMP+2]*1E3 << " um" << endl;
	cout <<  "\tdrot  =  " << alpar[mod*ndimMP+3]*1E3 << " mrad" << endl;
	cout <<  "\tdtilt =  " << alpar[mod*ndimMP+4]*1E3 << " mrad" << endl;
	cout <<  "\tdturn =  " << alpar[mod*ndimMP+5]*1E3 << " mrad" << endl;
      }
    }

    alignx[0] += alpar[1];
    aligny[0] += alpar[2];
    rotA[0]   += alpar[3];
    tiltA[0]  += alpar[4];
    turnA[0]  += alpar[5];

    alignx[2] += alpar[2*ndimMP+1];
    aligny[2] += alpar[2*ndimMP+2];
    rotA[2]   += alpar[2*ndimMP+3];
    tiltA[2]  += alpar[2*ndimMP+4];
    turnA[2]  += alpar[2*ndimMP+5];

  }

  // calculate sine / cosine of alignment parameters:

  double ctu[4] ;
  double stu[4] ;
  double cti[4] ;
  double sti[4] ;
  double cro[4] ;
  double sro[4] ;

  double wt = 180./TMath::Pi();


  for(size_t mod = 0; mod < 4; mod++){
    ctu[mod] = cos( turnA[mod] );
    stu[mod] = sin( turnA[mod] );
    cti[mod] = cos( tiltA[mod] );
    sti[mod] = sin( tiltA[mod] );
    cro[mod] = cos( rotA[mod] );
    sro[mod] = sin( rotA[mod] );
  }



  // pair matching cuts:

  //  double bicutx = 5E-3*dz; // [mm]
  //  double bicuty = 5E-3*dz; // [mm]

  double bicutx = 150E-3*dz; // [mm]
  double bicuty = 20E-3*dz; // [mm]

  // triplet linking cuts:

  //  double tricutx = 0.06; // [mm]
  //  double tricuty = 0.06; // [mm]

  double tricutx = 0.7; // [mm]
  double tricuty = 0.35; // [mm]

  /*
  if( tricutx < 3*dz*tetSi ) {
    tricutx = 3*dz*tetSi;
    tricuty = 3*dz*tetSi;
    }*/

  cout << "\nlinking cuts " << tricutx << ", " << tricuty << " mm\n" << endl;


  double chCutLow = 14;
  double chCutHigh = 40;
  double chCutResLow = 18;
  double chCutResHigh = 25;
  //double chCutLow = 16;
  //double chCutHigh = 25;

  // global labels for Pede:


  vector<int> labelsA( ndimMP );
  labelsA[0] =  1; // dx
  labelsA[1] =  2; // dy
  labelsA[2] =  3; // drot
  labelsA[3] =  4; // dtilt
  labelsA[4] =  5; // dturn
  labelsA[5] =  6; // dz

  vector<int> labelsB( ndimMP );
  labelsB[0] =  1*ndimMP+1; // dx
  labelsB[1] =  1*ndimMP+2; // dy
  labelsB[2] =  1*ndimMP+3; // drot
  labelsB[3] =  1*ndimMP+4; // dtilt
  labelsB[4] =  1*ndimMP+5; // dturn
  labelsB[5] =  1*ndimMP+6; // dz

  vector<int> labelsC( ndimMP );
  labelsC[0] =  2*ndimMP+1; // dx
  labelsC[1] =  2*ndimMP+2; // dy
  labelsC[2] =  2*ndimMP+3; // drot
  labelsC[3] =  2*ndimMP+4; // dtilt
  labelsC[4] =  2*ndimMP+5; // dturn
  labelsC[5] =  2*ndimMP+6; // dz

  vector<int> labelsD( ndimMP );
  labelsD[0] =  3*ndimMP+1; // dx
  labelsD[1] =  3*ndimMP+2; // dy
  labelsD[2] =  3*ndimMP+3; // drot
  labelsD[3] =  3*ndimMP+4; // dtilt
  labelsD[4] =  3*ndimMP+5; // dturn
  labelsD[5] =  3*ndimMP+6; // dz

  // Prepare the mille binary file:                                                                                                                                                                                 
  std::stringstream name;
  name << "run" << run << "/mille.bin";
  string m_millefilename = name.str();
  MilleBinary * mille;
  mille = new gbl::MilleBinary( m_millefilename.c_str() );

  // Landau peak cuts: Mon 27.7.2015


  string gainFileName[4];
  double ke[4][16];
  for(int mod = 0; mod < 4; mod++){
    for(int roc = 0; roc < 16; roc++){
      ke[mod][roc] = -1.;
    }
  }


  const int A = 0;
  const int B = 1;
  const int C = 2;
  const int D = 3;

  // Get gaincal files

  if( run >= 435 ) { // 2016 May
    
    stringstream gainstream;
    for(int mod = 0; mod < 4; mod++){
      gainstream << "gaincal/D" << modName[mod] << "-tb24-gaincal.dat";
      gainFileName[mod] = gainstream.str();
      gainstream.str("");
    }
  }
  
  

  // Get ke (conversion small Vcal -> e-)

  double turn = turnTable;

  if(useLCIO) turn = -1.*turn;

  double turnOnTable = 0.;
  if(run > 2180 && run < 2310) turnOnTable = -27.8;
  if(run >= 44000) turnOnTable = 27.8;

  double tilt = 19.3;

  double costilt = cos(tilt/wt);
  double sintilt = sin(tilt/wt);
  double costurn = cos(turn/wt);
  double sinturn = sin(turn/wt);
  double costurnon = cos(turnOnTable/wt);
  double sinturnon = sin(turnOnTable/wt);


  double norm = 1*TMath::Cos(TMath::Pi()*tilt/180.)*TMath::Cos(TMath::Pi()*(turn+turnOnTable)/180.);


  // for GBL:

  double resx = 9.9E-3; // [mm] col hit resolution
  double resy = 7.7E-3; // [mm] row hit resolution

  resx = 43.E-3;
  resy = 29.E-3;
  if(fabs(turn+turnOnTable) > 13) resx = 9.E-3;
  if(fabs(tilt) > 15) resy = 6.5E-3;


  // X0 Si = 21.82/2.33 = 9.365 cm
  // X0 Al = 24.01/2.70 = 8.89 cm
  // X0 Cu = 12.86/8.96 = 1.435 cm
  // X0 air = 36.62/1.204E-3 = 304 m

  double X0Si = (( 0.3 + 0.175 ) / 94.) / norm; // Sensor + ROC + HDI(appr.)

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

  TVectorD wscatSi(2);
  wscatSi[0] = 1.0 / ( tetSi * tetSi ); // weight
  wscatSi[1] = 1.0 / ( tetSi * tetSi );



  ostringstream conversionFileName; // output string stream

  conversionRun = alignmentRun;
  

  if(alignmentRun == run){
    conversionFileName << "conversions/conversion_" << run << ".dat";
  }else{
    conversionFileName << "conversions/conversion_" << conversionRun << ".dat";    
  }

  ifstream conversionFile( conversionFileName.str() );

  cout << endl;
  if( conversionFile.bad() || ! conversionFile.is_open() ) {
    cout << "no " << conversionFileName.str() << ", will bootstrap" << endl;
    cout << endl;
  }else{
    int modNr;
    int rocNr;
    double keVal;
    while(conversionFile >> modNr >> rocNr >> keVal){
      ke[modNr][rocNr] = keVal;
    }
  }
  for(int mod = 0; mod < 4; mod++){
    for(int roc = 0; roc < 16; roc++){
      if(ke[mod][roc] < 0.) ke[mod][roc] = 0.045;
    }
  }
  cout << "Roc-wise conversion factors:";
  for(int mod = 0; mod < 4; mod++){
    cout << endl << modName[mod] << ":\t";
    for(int roc = 0; roc < 16; roc++){
      cout << setprecision(4) << ke[mod][roc] << " ";
    }
  }
  cout << endl;

  conversionFile.close();
  
  // Shifts for pixel charge:

  int shifts[4] = {0,0,0,0};

  if(getShifts(run, shifts) == -2){
    for(int mod=0; mod < 4; mod++){
      shifts[mod] = 0;
    }
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

  fname << "histogramsMP/quad-" << run << ".root";

  TFile* histoFile = new TFile( fname.str(  ).c_str(  ), "RECREATE" );

  // book histos:
 
  TH1D hcol[4];
  TH1D hrow[4];
  TH1D hclx[4];
  TH1D hcly[4];
  TH1D hpxdig[4];
  TH1D hpxq[4];
  TH2D * hmap[4];
  TH2D * hxyGlobal[4];
  TH2D * hxzGlobal[4];
  TH1D hnpx[4];
  TH1D hsiz[4];
  TH1D hclq[4];
  TH1D hclq0[4];
  TH1D hclq0r[4][16];
  TH1D hclq0g[4];
  TH1D hncol[4];
  TH1D hncolq[4];
  TH1D hncolqf4[4];
  TH1D hnrow[4];
  TH1D heffRoc[4];
  TProfile ncolvsclq[4];
  TProfile effvsRoc[4];

  //  TDirectory * gbldir = histoFile->mkdir("GBL");  


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
    hclx[mod] = TH1D( Form("clx%c", modtos),
		      Form("%c Cluster x;x;%c clusters", modtos, modtos), 
		      416, -31.2, 31.2 );
    hcly[mod] = TH1D( Form("cly%c", modtos),
		      Form("%c Cluster y;y;%c clusters", modtos, modtos), 
		      160, -8.0, 8.0 );
    hpxq[mod] = TH1D( Form("pxq%c",modtos),
		      Form("%c pixel charge;pixel q [ke];%c pixels",modtos,modtos),
		      100, 0, 25 );
    hpxdig[mod] = TH1D( Form("pxdig%c",modtos),
		      Form("%c pixel pulse height;pixel PH [ADC counts];%c pixels",modtos,modtos),
		      256, 0, 255 );
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
    hclq0[mod] = TH1D( Form("clq0%c",modtos),
		       Form("%c normalized cluster charge;norm. cluster charge [ke];%c clusters",modtos,modtos),
		      100, 0, 100 );
    for( int roc = 0; roc < 16; ++roc ) {    
      hclq0r[mod][roc] = TH1D( Form("clq0%c%d",modtos,roc),
			       Form("%c, ROC %d, normalized cluster charge;norm. cluster charge [ke];%c clusters",modtos,roc,modtos),
			       100, 0, 100 );
    }
    hclq0g[mod] = TH1D( Form("clq0%cg",modtos),
		       Form("%c, between ROCs, normalized cluster charge;norm. cluster charge [ke];%c clusters",modtos,modtos),
		      100, 0, 100 );
    hncol[mod]= TH1D( Form("ncol%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hncolq[mod]= TH1D( Form("ncolq%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hncolqf4[mod]= TH1D( Form("ncolqf4%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hnrow[mod]= TH1D( Form("nrow%c",modtos),
		      Form("%c cluster size;rows/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    heffRoc[mod]= TH1D( Form("effRoc%c",modtos),
			Form("eff%c per ROC;efficiency;ROCs",modtos),
			25, 0.995, 1.0);
    ncolvsclq[mod]= TProfile( Form("ncolvsclq%c",modtos),
			     Form("ncol%c vs norm cluster charge;cluster charge [ke];ncol %c",modtos,modtos),
			     100,0,100);
    effvsRoc[mod]= TProfile( Form("eff%cvsRoc",modtos),
			     Form("eff%c vs ROC;ROC;eff %c",modtos,modtos),
			     16,-0.5,15.5, -1, 2);
    hxyGlobal[mod] = new  TH2D( Form("Globalxy%c",modtos),
			   Form("%c cluster map, global var.;x [mm];y [mm];%c clusters",modtos,modtos),
			   500, -40., 40., 160, -10., 10. );
    hxzGlobal[mod] = new  TH2D( Form("Globalxz%c",modtos),
			   Form("%c cluster map, global var.;x [mm];z [mm];%c clusters",modtos,modtos),
			   500, -40., 40., 160, -20., 20. );

  } // module planes

  TH2D * hxzAllGlobal = new  TH2D( "GlobalxzAll",
				 "cluster map, global var.;x [mm];z [mm];clusters",
				 500, -40., 40., 300, -15., 110. );


  TH2D hxxAB( "xxAB", "A vs B;col B;col A;clusters",
	      432, -32.4, 32.4, 432, -32.4, 32.4 );
  TH2D hyyAB( "yyAB", "A vs B;row B;row A;clusters",
	      162, -8.1, 8.1, 162, -8.1, 8.1 );
  TH1D hdxAB( "dxAB", "Ax-Bx;x-x [mm];cluster pairs", 150, -2, 2 );
  TH1D hdyAB( "dyAB", "Ay-By;y-y [mm];cluster pairs", 200, -2, 2 );
  TH1D hdxcAB( "dxcAB", "Ax-Bx;x-x [mm];cluster pairs", 200, -2, 2 );
  TH1D hdycAB( "dycAB", "Ay-By;y-y [mm];cluster pairs", 200, -2, 2 );
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
  TH1D hdxCB( "dxCB", "Cx-Bx;x-x [mm];cluster pairs", 200, -2, 2 );
  TH1D hdyCB( "dyCB", "Cy-By;y-y [mm];cluster pairs", 200, -2, 2 );
  TH1D hdxcCB( "dxcCB", "Cx-Bx;x-x [mm];cluster pairs", 200, -2, 2 );
  TH1D hdycCB( "dycCB", "Cy-By;y-y [mm];cluster pairs", 200, -2, 2 );
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
  TH1D hdxDB( "dxDB", "Dx-Bx;x-x [mm];cluster pairs", 200, -4, 4 );
  TH1D hdyDB( "dyDB", "Dy-By;y-y [mm];cluster pairs", 200, -2, 2 );
  TH1D hdxcDB( "dxcDB", "Dx-Bx;x-x [mm];cluster pairs", 200, -4, 4 );
  TH1D hdycDB( "dycDB", "Dy-By;y-y [mm];cluster pairs", 200, -4, 4 );
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
  TH1D hdxcACB( "dxcACB", "ACB dx;x-x [um];ACBs", 200, -700, 700 );
  TH1D hdycACB( "dycACB", "ACB dy;y-y [um];ACBs", 200, -700, 700 );
  TH1D hdxciACB( "dxciACB", "ACB dx;x-x [um];isolated ACBs",
		 200, -700, 700 );
  TH1D hdyciACB( "dyciACB", "ACB dy;y-y [um];isolated ACBs",
		 200, -700, 700 );
  TH1D hdycfACB( "dycfACB", "ACB dy;y-y [um];inner ACBs",
		 200, -700, 700 );
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
  TH1D hclq0D3( "clq0D3", "normalized D cluster charge;norm. cluster charge [ke];D3 clusters",
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
  TH1D hdxcfqBDC( "dxcfqBDC", "BDC dx;x-x [um];Landau peak inner BDCs",
		  200, -200, 200 );
  TH1D hdxcfqtBDC( "dxcfqtBDC", "BDC dx;x-x [um];Landau peak inner BDCs",
		  200, -200, 200 );
  TH1D hdxcfqtiBDC( "dxcfqtiBDC", "BDC dx;x-x [um];Landau peak inner BDCs",
		  200, -200, 200 );
  TProfile rmsxBDCvsq0("rmsxBDCvsq0",
		      "BDC dx vs q0; q0 [ke], BDC <abs(dx)>",
		      80, 0, 80, 0, 50);
  TH1D hdyciBDC( "dyciBDC", "BDC dy;y-y [um];isolated BDCs",
		 200, -200, 200 );
  TH1D hdycfBDC( "dycfBDC", "BDC dy;y-y [um];inner BDCs",
		 200, -200, 200 );
  TH1D hdycfqBDC( "dycfqBDC", "BDC dy;y-y [um];Landau peak inner BDCs",
		  200, -200, 200 );
  TH1D hdycfqtBDC( "dycfqtBDC", "BDC dy;y-y [um];Landau peak inner BDCs",
		  200, -200, 200 );
  TH1D hdycfqtiBDC( "dycfqtiBDC", "BDC dy;y-y [um];Landau peak inner BDCs",
		  200, -200, 200 );
  TProfile rmsyBDCvsq0("rmsyBDCvsq0",
		      "BDC dy vs q0; q0 [ke], BDC <abs(dy)>",
		      80, 0, 80, 0, 50);


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

  TH1D htxBDC( "txBDC", "tri angle BDC x;tri angle BDC x [mrad];BDC triplets", 100, -10, 10 );
  TH1D htyBDC( "tyBDC", "tri angle BDC y;tri angle BDC y [mrad];BDC triplets", 100, -10, 10 );


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
  TH1D hclq0A3( "clq0A3", "normalized A cluster charge;norm. cluster charge [ke];A3 clusters",
	       100, 0, 100 );
  TH1D hncolA3( "ncolA3", "A cluster size;columns/cluster;A3 clusters",
		21, -0.5, 20.5 );
  TH1D hncolA3q1( "ncolA3q1", "A cluster size;columns/cluster;A3 clusters",
		21, -0.5, 20.5 );
  TH1D hncolA3q2( "ncolA3q2", "A cluster size;columns/cluster;A3 clusters",
		21, -0.5, 20.5 );
  TH1D hncolA3q3( "ncolA3q3", "A cluster size;columns/cluster;A3 clusters",
		21, -0.5, 20.5 );
  TH1D hnrowA3( "nrowA3", "A cluster size;rows/cluster;A3 clusters",
		21, -0.5, 20.5 );

  // A-D tracks:

  TH2D hxxDA( "xxDA", "D vs A;col A;col D;clusters",
	      432, -32.4, 32.4, 432, -32.4, 32.4 );
  TH2D hyyDA( "yyDA", "D vs A;row A;row D;clusters",
	      162, -8.1, 8.1, 162, -8.1, 8.1 );
  TH1D hdxDA( "dxDA", "Dx-xA;x-x [mm];cluster pairs", 200, -5, 5 );
  TH1D hdyDA( "dyDA", "Dy-Ay;y-y [mm];cluster pairs", 200, -5, 5 );
  TH1D hdxDAwoa( "dxDAwoa", "Dx-Ax w/o alignment;x-x [mm];cluster pairs", 200, -20, 20 );
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
  TH1D hclq0B4( "clq0B4", "normalized B cluster charge;norm. cluster charge [ke];B4 clusters",
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
  TH1D hclq0C4( "clq0C4", "normalized C cluster charge;norm. cluster charge [ke];C4 clusters",
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


  TH1D hdx4C( "dx4C", "4C dx;x-x [mm];4plets", 200, -1, 1 );
  TH1D hdy4C( "dy4C", "4C dy;y-y [mm];4plets", 200, -1, 1 );
  TH1D hdxc4C( "dxc4C", "4C dx;x-x [um];4plets", 200, -200, 200 );
  TH1D hdyc4C( "dyc4C", "4C dy;y-y [um];4plets", 200, -200, 200 );

  TH1D hdxci4C( "dxci4C", "4C dx;x-x [um];isolated 4plets",
		 200, -200, 200 );
  TH1D hdxcfq4C( "dxcfq4C", "4C dx;x-x [um];Landau peak inner 4plets",
		  200, -200, 200 );
  TH1D hdxcfqt4C( "dxcfqt4C", "4C dx;x-x [um];Landau peak inner 4plets",
		  200, -200, 200 );
  TH1D hdxcfqti4C( "dxcfqti4C", "4C dx;x-x [um];Landau peak inner 4plets",
		  200, -200, 200 );
  TProfile rmsx4Cvsq0("rmsx4Cvsq0",
		      "4C dx vs q0; q0 [ke], 4C <abs(dx)>",
		      80, 0, 80, 0, 50);
  TH1D hdyci4C( "dyci4C", "4C dy;y-y [um];isolated 4plets",
		 200, -200, 200 );
  TH1D hdycf4C( "dycf4C", "4C dy;y-y [um];inner 4plets",
		 200, -200, 200 );
  TH1D hdycfq4C( "dycfq4C", "4C dy;y-y [um];Landau peak inner 4plets",
		  200, -200, 200 );
  TH1D hdycfqt4C( "dycfqt4C", "4C dy;y-y [um];Landau peak inner 4plets",
		  200, -200, 200 );
  TH1D hdycfqti4C( "dycfqti4C", "4C dy;y-y [um];Landau peak inner 4plets",
		  200, -200, 200 );
  TProfile rmsy4Cvsq0("rmsy4Cvsq0",
		      "4C dy vs q0; q0 [ke], 4plet <abs(dy)>",
		      80, 0, 80, 0, 50);
  TH1D hmx4C( "mx4C", "4C fit slope x;slope [mrad];4plets", 200, -5, 5 );
  TH1D hbx4C( "bx4C", "4C fit x0;fitted x0 [mm];4plets", 200, -32, 32 );
  TH1D hmy4C( "my4C", "4C fit slope y;slope [mrad];4plets", 200, -5, 5 );
  TH1D hby4C( "by4C", "4C fit y0;fitted y0 [mm];4plets", 200, -8, 8 );

  TH1D hmx4fit( "mx4fit", "4fit fit slope x;slope [mrad];4plets", 200, -5, 5 );
  TH1D hbx4fit( "bx4fit", "4fit fit x0;fitted x0 [mm];4plets", 200, -32, 32 );
  TH1D hmy4fit( "my4fit", "4fit fit slope y;slope [mrad];4plets", 200, -5, 5 );
  TH1D hby4fit( "by4fit", "4fit fit y0;fitted y0 [mm];4plets", 200, -8, 8 );

  TH1D hsigmax4fit( "sx4fit", "4fit fit sigma x;sigma x [um];4plets", 100, 0, 100 );
  TH1D hsigmay4fit( "sy4fit", "4fit fit sigma y;sigma y [um];4plets", 100, 0, 100 );

  TH1D hcircdistXBmC( "circdistXBmC", "x distances from ADline, C-B;Delta dx [um];4plets", 100, -100, 100 );
  TH1D hcircdistYBmC( "circdistYBmC", "y distances from ADline, C-B;Delta dy [um];4plets", 100, -100, 100 );

  TH1D hcircdistXB( "circdistXB", "x distance from ADline, B ; dx [um];4plets", 100, -200, 200 );
  TH1D hcircdistYB( "circdistYB", "y distance from ADline, B ; dy [um];4plets", 100, -100, 100 );
  TH1D hcircdistXC( "circdistXC", "x distance from ADline, C ; dx [um];4plets", 100, -200, 200 );
  TH1D hcircdistYC( "circdistYC", "y distance from ADline, C ; dy [um];4plets", 100, -100, 100 );


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


  TH1D hchi2( "GBL/chi2", "GBL chisq;chisq;track fits", 100, 0, 50 );
  TH1D hprob( "GBL/prob", "GBL fit prob;fit probability;track fits", 100, 0, 1 );

  TH1D fitrxHisto[4];
  TH1D fitryHisto[4];
  TH1D gblaxHisto[4];
  TH1D gblayHisto[4];
  TH1D gbldxHisto[4];
  TH1D gbldyHisto[4];
  TH1D gblrxHisto[4];
  TH1D gblryHisto[4];
  TH1D gblrxqHisto[4];
  TH1D gblryqHisto[4];
  TH1D gblpxHisto[4];
  TH1D gblpyHisto[4];
  TH1D gblqxHisto[4];
  TH1D gblqyHisto[4];
  TH1D gbltxHisto[4];
  TH1D gbltyHisto[4];
  TH1D gblqpHisto[4];
  
  for( int mod = 0; mod < 4; ++mod ) {
    char modtos;
    switch( mod ) {
    case 0: modtos = 'A'; break;
    case 1: modtos = 'B'; break;
    case 2: modtos = 'C'; break;
    case 3: modtos = 'D'; break;
    default:modtos = 'X';
    }
    
    fitrxHisto[mod] = TH1D(Form("fitrx%c",modtos), Form("Linear fit residual x %c;resid x [um];track fits",modtos), 200, -150, 150);
    fitryHisto[mod] = TH1D(Form("fitry%c",modtos), Form("Linear fit residual y %c;resid y [um];track fits",modtos), 200, -150, 150);

    gblaxHisto[mod] = TH1D(Form("gblax%c",modtos), Form("GBL angle x %c;angle x [mrad];track fits",modtos), 500, -100, 100);
    gblayHisto[mod] = TH1D(Form("gblay%c",modtos), Form("GBL angle y %c;angle y [mrad];track fits",modtos), 100, -5, 5);
    gbldxHisto[mod] = TH1D(Form("gbldx%c",modtos), Form("GBL shift x %c;shift x [um];track fits",modtos), 100, -200, 200);
    gbldyHisto[mod] = TH1D(Form("gbldy%c",modtos), Form("GBL shift y %c;shift y [um];track fits",modtos), 100, -100, 100);
    gblrxHisto[mod] = TH1D(Form("gblrx%c",modtos), Form("GBL residual x %c;resid x [um];track fits",modtos), 100, -50, 50);
    gblryHisto[mod] = TH1D(Form("gblry%c",modtos), Form("GBL residual y %c;resid y [um];track fits",modtos), 100, -50, 50);
    gblrxqHisto[mod] = TH1D(Form("gblrxq%c",modtos), Form("GBL residual x, q cut %c;resid x [um];track fits",modtos), 100, -50, 50);
    gblryqHisto[mod] = TH1D(Form("gblryq%c",modtos), Form("GBL residual y, q cut %c;resid y [um];track fits",modtos), 100, -50, 50);
    gblpxHisto[mod] = TH1D(Form("gblpx%c",modtos), Form("GBL pull x %c;pull x [sigma];track fits",modtos), 100, -10, 10);
    gblpyHisto[mod] = TH1D(Form("gblpy%c",modtos), Form("GBL pull y %c;pull y [sigma];track fits",modtos), 100, -10, 10);
    gblqxHisto[mod] = TH1D(Form("gblqx%c",modtos), Form("GBL kink x %c;kink x [mrad];track fits",modtos), 100, -1, 1);
    gblqyHisto[mod] = TH1D(Form("gblqy%c",modtos), Form("GBL kink y %c;kink y [mrad];track fits",modtos), 100, -1, 1);
    gbltxHisto[mod] = TH1D(Form("gbltx%c",modtos), Form("GBL kink pull x %c;kink pull;track fits",modtos), 100, -10, 10);
    gbltyHisto[mod] = TH1D(Form("gblty%c",modtos), Form("GBL kink pull y %c;kink pull;track fits",modtos), 100, -10, 10);
    gblqpHisto[mod] = TH1D(Form("gblqp%c",modtos), Form("GBL p %c;p;track fits",modtos), 200, 0., 10.);
    
  
    
  } // module planes
  

  TH1D hnADC( "nADC", "ADCplets;ADCplets;events", 21, -0.5, 20.5 );
  TH1D hnADB( "nADB", "ADBplets;ADBplets;events", 21, -0.5, 20.5 );
  TH1D hn4ev( "n4ev", "4plets;4plets;events", 21, -0.5, 20.5 );

  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  int event_nr = 0;

  int n4 = 0;
  int nmille = 0;
  int ntry = 0;
  do {

    vector <cluster> cl[4];
    int n4ev = 0;


  
    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();
  
    if( evt.IsBORE() ){
      eudaq::PluginManager::Initialize(evt);
      cout << "I'm BOREing" << endl;
    }

    bool ldb = false;
    if( event_nr == -1 )
      ldb = 1;
  
    if( ldb || event_nr%10000 == 0 )
      cout << "Quad processing event (raw) "<< event_nr << endl;
  
    StandardEvent sevt = eudaq::PluginManager::ConvertToStandard(evt);
    
    int xm = 0;
    int ym = 0;
    int adc = 0;
    double cal = 0;
  
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

	// Shift adc to known threshold

	adc += shifts[mod];

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
	  if(modName[mod]>=4000 && modName[mod]<4005){ // Even out the gain calibration for each pixel for the simulation
	    a0 = 0.999738;
	    a1 = 322248.1;
	    a2 = 4894.3;
	    a3 = 206.7;
	    a4 = 184.2;
	    a5 = 6.846;
	  }
	  cal = PHtoVcal( adc, a0, a1, a2, a3, a4, a5 ) * ke[mod][roc]; // [ke]
	}
	
	hpxdig[mod].Fill( adc );
	hpxq[mod].Fill( cal );

	// fill pixel block for clustering
	pb[npx].col = x;
	pb[npx].row = y;
	pb[npx].roc = roc;
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
	hclq[mod].Fill( cA->charge );
	hclq0[mod].Fill( cA->charge*norm );
	hncol[mod].Fill( cA->ncol );
	hnrow[mod].Fill( cA->nrow );

	double x1 = cA->col*0.15 - 32.4;
	double y1 = -1.*(cA->row*0.10 - 8.1);
	double z1 = 0.;

	// tilt around x:

	double x2 = x1;
	double y2 = costilt*y1 + sintilt*z1;
	double z2 =-sintilt*y1 + costilt*z1;

	// turn around y (turning the modules on the table):

	double x3 = costurnon*x2 + sinturnon*z2;
	double y3 = y2;
	double z3 =-sinturnon*x2 + costurnon*z2;

	// spread along z:

	z3 += dz*mod; // 0 32 64 96
	  

	// turn around y (turning the whole table):

	double x4 = costurn*x3 + sinturn*z3;
	double y4 = y3;
	double z4 =-sinturn*x3 + costurn*z3;

	// Prealignment & Mille:

	double x5 = x4 - alignx[mod];
	double y5 = y4 - aligny[mod];
	double z5 = z4 - alignz[mod];

	cA->x5 = x5;
	cA->y5 = y5;
	cA->z5 = z5;

	double x6 = cro[mod]*x5 - sro[mod]*y5;
	double y6 = sro[mod]*x5 + cro[mod]*y5;
	double z6 = z5;

	double x7 = x6;
	double y7 = cti[mod]*y6 - sti[mod]*z6;
	double z7 = sti[mod]*y6 + cti[mod]*z6;

	double x8 = ctu[mod]*x7 + stu[mod]*z7;
	double y8 = y7;
	double z8 =-stu[mod]*x7 + ctu[mod]*z7;

	cA->xg = x8;
	cA->yg = y8;
	cA->zg = z8;


	hclx[mod].Fill( x5 );
	hcly[mod].Fill( y5 );
	  
	ncolvsclq[mod].Fill( cA->charge*norm, cA->ncol );
	if(cA->roc == -1){
	  hclq0g[mod].Fill( cA->charge*norm );
	}else{
	  hclq0r[mod][cA->roc].Fill( cA->charge*norm );
	}
	if(cA->charge*norm > chCutLow && cA->charge*norm < chCutHigh){
	  hncolq[mod].Fill( cA->ncol );
	}

	hxyGlobal[mod]->Fill( x5, y5 ); // front view
	hxzGlobal[mod]->Fill( x5, z5-32.*mod ); // top view

	hxzAllGlobal->Fill( x5, z5 ); // top view, all

      }
    }
    if(!(reader->NextEvent())) break;
    
        
    ++event_nr;

    
    // Prealignment:

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      double xA = cA->xg;
      double yA = cA->yg;

      for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

      double xB = cB->xg;
      double yB = cB->yg;
     
      double dx = xA-xB;
      double dy = yA-yB;

      hdxAB.Fill(dx);
      hdyAB.Fill(dy);
 
      }
    }

    for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

      double xC = cC->xg;
      double yC = cC->yg;

      for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

      double xB = cB->xg;
      double yB = cB->yg;
     
      double dx = xC-xB;
      double dy = yC-yB;

      hdxCB.Fill(dx);
      hdyCB.Fill(dy);
 
      }
    }

    for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

      double xD = cD->xg;
      double yD = cD->yg;

      for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

      double xB = cB->xg;
      double yB = cB->yg;
     
      double dx = xD-xB;
      double dy = yD-yB;

      hdxDB.Fill(dx);
      hdyDB.Fill(dy);
 
      }
    }

  
    int nADC = 0;
    int nADB = 0;

    bool iso = cl[A].size() == 1 && cl[D].size() == 1;

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {


      double xA = cA->xg;
      double yA = cA->yg;
      double x5A = cA->x5;
      double y5A = cA->y5;
      double z5A = cA->z5;

      TMatrixD derivA( 2, ndimMP ); // -alignment derivatives x,y

      // First order -> Terms with only cosine used. If !avail, use terms with one sine.

      derivA[0][0] = 1.0; // dx/dx
      derivA[1][0] = 0.0;

      derivA[0][1] = 0.0;
      derivA[1][1] = 1.0; // dy/dy

      derivA[0][2] = ctu[A]*cro[A]*y5A; // dx/drot
      derivA[1][2] =-cti[A]*cro[A]*x5A;

      derivA[0][3] =-stu[A]*cro[A]*cti[A]*y5A; // dx/dtilt
      derivA[1][3] = cti[A]*z5A;

      derivA[0][4] =-ctu[A]*cti[A]*z5A; // dx/dturn
      derivA[1][4] = 0.0;


      double qA = cA->charge*norm;
      bool lqA = 1;
      if(      qA < chCutResLow ) lqA = 0;
      else if( qA > chCutResHigh ) lqA = 0;

      for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	double xD = cD->xg;
	double yD = cD->yg;

	TMatrixD derivD( 2, ndimMP ); // alignment derivatives x,y

	// First order -> Terms with only cosine used. If !avail, use terms with one sine.
	
	double x5D = cD->x5;
	double y5D = cD->y5;
	double z5D = cD->z5;

	derivD[0][0] = 1.0; // dx/dx
	derivD[1][0] = 0.0;

	derivD[0][1] = 0.0;
	derivD[1][1] = 1.0; // dy/dy

	derivD[0][2] = ctu[D]*cro[D]*y5D; // dx/drot
	derivD[1][2] =-cti[D]*cro[D]*x5D;

	derivD[0][3] =-stu[D]*cro[D]*cti[D]*y5D; // dx/dtilt
	derivD[1][3] = cti[D]*z5D;

	derivD[0][4] =-ctu[D]*cti[D]*z5D; // dx/dturn
	derivD[1][4] = 0.0;


	double qD = cD->charge*norm;
	bool lqD = 1;
	if(      qD < chCutResLow ) lqD = 0;
	else if( qD > chCutResHigh ) lqD = 0;

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

	  double xC = cC->xg;
	  double yC = cC->yg;


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

	  //if( abs( dx3 ) > tricutx ) continue; // tight tri
	  //if( abs( dy3 ) > tricuty ) continue;

	  ++nADC;

	  htxADC.Fill( slpx*1E3 );
	  htyADC.Fill( slpy*1E3 );

	  // for tri ACB:

	  double xavg2B = 0.5*(xA + xC); // interpolate
	  double yavg2B = 0.5*(yA + yC); // equidistant

	
	  // efficiency of B:

	  int nm[99] = {0};

	  for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	    double xB = cB->xg;
	    double yB = cB->yg;

	    // tri ACB:

	    double dx4 = xB - xavg2B;
	    double dy4 = yB - yavg2B;

	    for( int iw = 1; iw < 41; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) // for eff
		nm[iw] = 1;

	    if( abs( dx4 ) < tricutx && abs( dy4 ) < tricuty  &&
		cA->big == 0 && cC->big == 0 && cB->big == 0 ) {
	      hsizB4.Fill( cB->size );
	      hclqB4.Fill( cB->charge );
	      hclq0B4.Fill( cB->charge*norm );
	      hncolB4.Fill( cB->ncol );
	      hnrowB4.Fill( cB->nrow );
	      if(cB->charge*norm > chCutLow && cB->charge*norm < chCutHigh){
		hncolqf4[1].Fill(cB->ncol);
	      }


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

	  for( int iw = 1; iw < 41; ++ iw )
	    effBvsw.Fill( iw*0.050+0.005, nm[iw] );


	} // cl C

	// - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
	// tri ADB:

	double xavg3B = (2*xA + xD)/3; // interpolate
	double yavg3B = (2*yA + yD)/3; // A and D to B

	for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	  double xB = cB->xg;
	  double yB = cB->yg;

	  TMatrixD derivB( 2, ndimMP ); // alignment derivatives x,y

	  double x5B = cB->x5;
	  double y5B = cB->y5;
	  double z5B = cB->z5;

	  derivB[0][0] = 1.0; // dx/dx
	  derivB[1][0] = 0.0;

	  derivB[0][1] = 0.0;
	  derivB[1][1] = 1.0; // dy/dy

	  derivB[0][2] = ctu[B]*cro[B]*y5B; // dx/drot
	  derivB[1][2] =-cti[B]*cro[B]*x5B;
	  
	  derivB[0][3] =-stu[B]*cro[B]*cti[B]*y5B; // dx/dtilt
	  derivB[1][3] = cti[B]*z5B;
	  
	  derivB[0][4] =-ctu[B]*cti[B]*z5B; // dx/dturn
	  derivB[1][4] = 0.0;

	  double qB = cB->charge*norm;
	  bool lqB = 1;
	  if(      qB < chCutResLow ) lqB = 0;
	  else if( qB > chCutResHigh ) lqB = 0;

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

	  //if( abs( dx3 ) > tricutx ) continue; // tight tri
	  //if( abs( dy3 ) > tricuty ) continue;

	  ++nADB;

	  htxADB.Fill( slpx );
	  htyADB.Fill( slpy );

	  // B-D track:

	  double xavg2C = 0.5*(xB + xD); // equidistant
	  double yavg2C = 0.5*(yB + yD);

	  // efficiency of C vs ADB:

	  int nm[99] = {0};

	  for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	    double xC = cC->xg;
	    double yC = cC->yg;

	    TMatrixD derivC( 2, ndimMP ); // alignment derivatives x,y

	    double x5C = cC->x5;
	    double y5C = cC->y5;
	    double z5C = cC->z5;

	    derivC[0][0] = 1.0; // dx/dx
	    derivC[1][0] = 0.0;

	    derivC[0][1] = 0.0;
	    derivC[1][1] = 1.0; // dy/dy

	    derivC[0][2] = ctu[C]*cro[C]*y5C; // dx/drot
	    derivC[1][2] =-cti[C]*cro[C]*x5C;
	    
	    derivC[0][3] =-stu[C]*cro[C]*cti[C]*y5C; // dx/dtilt
	    derivC[1][3] = cti[C]*z5C;
	    
	    derivC[0][4] =-ctu[C]*cti[C]*z5C; // dx/dturn
	    derivC[1][4] = 0.0;

	    double qC = cC->charge*norm;
	    bool lqC = 1;
	    if(      qC < chCutResLow ) lqC = 0;
	    else if( qC > chCutResHigh ) lqC = 0;

	    double dx4 = xC - xavg2C;
	    double dy4 = yC - yavg2C;

	    for( int iw = 1; iw < 41; ++ iw )
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ) // for eff
		nm[iw] = 1;

	    //if( abs( dx4 ) > tricutx ) continue; // quad
	    //if( abs( dy4 ) > tricuty ) continue;

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
	      hclqC4.Fill( cC->charge );
	      hclq0C4.Fill( cC->charge*norm );
	      hncolC4.Fill( cC->ncol );
	      hnrowC4.Fill( cC->nrow );
	      hminxC4.Fill( (minx-1)%2 ); 
	      hmaxxC4.Fill( (maxx-1)%2 ); 
	      if(cC->charge*norm > chCutLow && cC->charge*norm < chCutHigh){
		hncolqf4[2].Fill(cC->ncol);
	      }
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


	    // Only use right charge cuts:

	    // if(!(lqA && lqB && lqC && lqD)) continue;

	    double rowdown = -50.;
	    double rowup = 60.;

	    if(!useLCIO){
	      if(cA->row < rowdown || cA->row > rowup ) continue;
	      if(cB->row < rowdown || cB->row > rowup ) continue;
	      if(cC->row < rowdown || cC->row > rowup ) continue;
	    }
	    
	    // Initial fit:
	    double za=cA->zg;
	    double zb=cB->zg;
	    double zc=cC->zg;
	    double zd=cD->zg;
	    
	    double mx,bx,my,by,sigmax,sigmay;

	    linefit4(za,xA,zb,xB,zc,xC,zd,xD,mx,bx,sigmax);
	    linefit4(za,yA,zb,yB,zc,yC,zd,yD,my,by,sigmay);

	    fitrxHisto[A].Fill((xA-(mx*za+bx))*1.E3);
	    fitryHisto[A].Fill((yA-(my*za+by))*1.E3);
	    fitrxHisto[B].Fill((xB-(mx*zb+bx))*1.E3);
	    fitryHisto[B].Fill((yB-(my*zb+by))*1.E3);
	    fitrxHisto[C].Fill((xC-(mx*zc+bx))*1.E3);
	    fitryHisto[C].Fill((yC-(my*zc+by))*1.E3);
	    fitrxHisto[D].Fill((xD-(mx*zd+bx))*1.E3);
	    fitryHisto[D].Fill((yD-(my*zd+by))*1.E3);

	    vector<GblPoint> listOfPoints;
	    listOfPoints.reserve(4);
	    vector<double> sPoint;

	    ntry++;

	    // plane A:

	    TMatrixD jacPointToPoint(5, 5);
	    jacPointToPoint.UnitMatrix();
	    GblPoint *point = new GblPoint(jacPointToPoint);
	    meas[0] = xA-(mx*za+bx);
	    meas[1] = yA-(my*za+by);
	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsA, derivA ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // B:

	    jacPointToPoint = Jac5( zb-za, bfield  );
	    point = new GblPoint(jacPointToPoint);
	    meas[0] = xB-(mx*zb+bx);
	    meas[1] = yB-(my*zb+by);
	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsB, derivB ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // C:

	    jacPointToPoint = Jac5( zc-zb, bfield );
	    point = new GblPoint(jacPointToPoint);
	    meas[0] = xC-(mx*zc+bx);
	    meas[1] = yC-(my*zc+by);
	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsC, derivC ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // D:

	    jacPointToPoint = Jac5( zd-zc, bfield );
	    point = new GblPoint(jacPointToPoint);
	    meas[0] = xD-(mx*zd+bx);
	    meas[1] = yD-(my*zd+by);
	    point->addMeasurement( proL2m, meas, measPrec );
	    point->addScatterer( scat, wscatSi );
	    point->addGlobals( labelsD, derivD ); // for MillePede alignment
	    listOfPoints.push_back(*point);
	    delete point;

	    // track fit:

	    bool curved = true;
	    if(bfield < 0.001) curved = false;
	    GblTrajectory traj( listOfPoints, curved ); // 0 = no magnetic field
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



	    // Let's see the results:

	    TVectorD aCorrection(5);
	    TMatrixDSym aCovariance(5);

	    unsigned int ndim = 2;
	    TVectorD aResiduals(ndim);
	    TVectorD aMeasErrors(ndim);
	    TVectorD aResErrors(ndim);
	    TVectorD aDownWeights(ndim);

	    TVectorD kKinks(ndim);
	    TVectorD kKinkErrors(ndim);
	    TVectorD kResErrors(ndim);
	    TVectorD kDownWeights(ndim);
	    
	    //track = q/p, x', y', x, y
	    //        0,   1,  2,  3, 4


	    if(probchi > 0.01){
	      
	      for(size_t iplane = 1; iplane < 5; iplane++){
		
		traj.getResults( iplane, aCorrection, aCovariance );
		traj.getMeasResults( static_cast<unsigned int>(iplane), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
		traj.getScatResults( static_cast<unsigned int>(iplane), ndim, kKinks, kKinkErrors, kResErrors, kDownWeights );
		
		gblaxHisto[iplane-1].Fill( (mx-aCorrection[1])*1E3 ); // angle x [mrad]
		gblayHisto[iplane-1].Fill( (my-aCorrection[2])*1E3 ); // angle y [mrad]
		gbldxHisto[iplane-1].Fill( aCorrection[3]*1E3 ); // shift x [um]
		gbldyHisto[iplane-1].Fill( aCorrection[4]*1E3 ); // shift y [um]
		gblrxHisto[iplane-1].Fill( aResiduals[0] * 1E3 ); // residual x [um]
		gblryHisto[iplane-1].Fill( aResiduals[1] * 1E3 ); // residual y [um]
		gblpxHisto[iplane-1].Fill( aResiduals[0] / aResErrors[0] ); // pull x
		gblpyHisto[iplane-1].Fill( aResiduals[1] / aResErrors[1] ); // pull y
		gblqxHisto[iplane-1].Fill( kKinks[0]*1E3 ); // kink x 
		gblqyHisto[iplane-1].Fill( kKinks[1]*1E3 ); // kink y
		gbltxHisto[iplane-1].Fill( kKinks[0]/kResErrors[0] ); // x kink pull ?
		gbltyHisto[iplane-1].Fill( kKinks[1]/kResErrors[1] ); // y kink pull ?
		gblqpHisto[iplane-1].Fill( 1/aCorrection[0]*1E-3 );

		if(lqA && lqB && lqC && lqD){
		  gblrxqHisto[iplane-1].Fill( aResiduals[0] * 1E3 ); // residual x [um]
		  gblryqHisto[iplane-1].Fill( aResiduals[1] * 1E3 ); // residual y [um]
		}
		
	      }
	    }



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
	  for( int iw = 1; iw < 41; ++ iw )
	    effCvsw.Fill( iw*0.050+0.005, nm[iw] );
	

	} // cl B

      } // cl D

    } // cl A

    hnADC.Fill( nADC );
    hnADB.Fill( nADB );
    hn4ev.Fill( n4ev );


  } while( event_nr < lev );

  cout << endl << "events " << event_nr << endl;

  // Analyze efficiency plots:
  

  histoFile->Write();
  histoFile->Close();

  cout << nmille << " GBL tracks." << endl;
  cout << ntry << " attempts for GBL fits." << endl;

  cout << endl << histoFile->GetName() << endl << endl;

  delete mille;

  // MILLEPEDE:


  ofstream steerFile;
  std::stringstream steername;
  steername << "run" << run << "/steerPede.txt";
  steerFile.open( steername.str().c_str() );

  if( steerFile.is_open() ) {



    steerFile << "! generated by quadMP" << std::endl;
    steerFile << "Cfiles" << std::endl;
    steerFile << "mille.bin" << std::endl;
    steerFile << std::endl;

    steerFile << "Parameter" << std::endl;

    bool fixed[4] = {false, true, false, true};

    for(size_t mod = 0; mod < 4 ; mod++){
      
      if(fixed[mod]){
	steerFile << mod*ndimMP+1 << "  0.0  -1.0" << std::endl; // dx
	steerFile << mod*ndimMP+2 << "  0.0  -1.0" << std::endl; // dy
	steerFile << mod*ndimMP+3 << "  0.0  -1.0" << std::endl; // drot
	steerFile << mod*ndimMP+4 << "  0.0  -1.0" << std::endl; // dtilt
	steerFile << mod*ndimMP+5 << "  0.0  -1.0" << std::endl; // dturn
	steerFile << std::endl;
    
      }else{
	steerFile << mod*ndimMP+1 << "  0.0  0.0" << std::endl; // dx
	steerFile << mod*ndimMP+2 << "  0.0  0.0" << std::endl; // dy
	steerFile << mod*ndimMP+3 << "  0.0  0.0" << std::endl; // drot
	steerFile << mod*ndimMP+4 << "  0.0  0.0" << std::endl; // dtilt
	steerFile << mod*ndimMP+5 << "  0.0  0.0" << std::endl; // dturn
	steerFile << std::endl;
      }
    }

    steerFile << "! chiscut 5.0 2.5" << std::endl;
    steerFile << "outlierdownweighting 4" << std::endl;
    steerFile << std::endl;
    steerFile << "method inversion 10  0.1" << std::endl;
    steerFile << "threads 10 1" << std::endl;
    steerFile << std::endl;
    steerFile << "! histprint" << std::endl;
    steerFile << std::endl;
    steerFile << "end" << std::endl;

    steerFile.close();



  }



  // Prealignment

  double goodtracks = (double)nmille/(double)ntry;

  if(goodtracks < 0.25 && alignmentRun == run){
    
    cout << "Not enough tracks for MP. Do Prealignment." << endl;

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
	   << endl << "  m " << fgp0->GetParameter(1)
	   << endl << "  alignx[0] = " << alignx[0] + fgp0->GetParameter(1)
	   << endl;
      alignx[0] += fgp0->GetParameter(1);

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
	   << endl << "  m " << fgp0->GetParameter(1)
	   << endl << "  aligny[0] = " << aligny[0] + fgp0->GetParameter(1)
	   << endl;
      aligny[0] += fgp0->GetParameter(1);
    }

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
	   << endl << "  m " << fgp0->GetParameter(1)
	   << endl << "  alignx[2] = " << alignx[2] + fgp0->GetParameter(1)
	   << endl;
      alignx[2] += fgp0->GetParameter(1);

    }

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
	   << endl << "  m " << fgp0->GetParameter(1)
	   << endl << "  aligny[2] = " << aligny[2] + fgp0->GetParameter(1)
	   << endl;
      aligny[2] += fgp0->GetParameter(1);
    }

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
	   << endl << "  m " << fgp0->GetParameter(1)
	   << endl << "  alignx[3] = " << alignx[3] + fgp0->GetParameter(1)
	   << endl;
      alignx[3] += fgp0->GetParameter(1);
    
    }

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
	   << endl << "  m " << fgp0->GetParameter(1)
	   << endl << "  aligny[3] = " << aligny[3] + fgp0->GetParameter(1)
	   << endl;
      aligny[3] += fgp0->GetParameter(1);

    }

  }else if(alignmentRun == run){ 
    // if prealignment was not performed, execute millepede.

    if(chdir(rundir.c_str()) != 0){
      cout << "Could not enter run directory." << endl;
    }else{
      
      stringstream pedecmd;
      pedecmd << "./../millepede/pede steerPede.txt";
      
      system(pedecmd.str().c_str());
      
      if(chdir("..") != 0){
	cout << "Could not go back to main directory." << endl;
      }
    }

  }


  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // write alignment to file:
  
  if(alignmentRun == run){

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
      alignFile << "rot " << rotA[ipl] << endl;
      alignFile << "tilt " << tiltA[ipl] << endl;
      alignFile << "turn " << turnA[ipl] << endl;

      cout << endl;
      cout << "plane " << ipl << endl;
      cout << "alignx " << alignx[ipl] << endl;
      cout << "aligny " << aligny[ipl] << endl;
      cout << "rot " << rotA[ipl] << endl;
      cout << "tilt " << tiltA[ipl] << endl;
      cout << "turn " << turnA[ipl] << endl;
    } // ipl

    alignFile.close();

  
  
    cout << endl
	 << "written to " << alignFileName.str()
	 << endl << endl;

  }else{
    cout << endl
	 << "Not writing an alignment file."
	 << endl << endl;    
  }


  cout << "Just to sum it up: " << endl;
  cout << nmille << " GBL tracks." << endl;
  cout << ntry << " attempts for GBL fits." << endl;
  cout << "Ratio: " << goodtracks << endl << endl;


  // Correct conversion factors

  if(fitSupressed){
    cout << "Skipping landau fits." << endl;
  }else{

    ofstream conversionFileOut;

    if(!CCSupressed || conversionRun == run){
      conversionFileOut.open(conversionFileName.str(),std::ofstream::out);

      if(conversionFileOut.is_open()){
	cout << "Conversion file opened. Let's put the right factors there." << endl;
      }else{
	cout << "Conversion file NOT opened. I repeat: NOT OPENED!" << endl;
      }
    }

    double landau_peak[4][16];
    double correction[4][16];

    //cout << "conv run: " << conversionRun << endl;

    cout << "\nLandau peaks [ke]:";
    for(int mod = 0; mod < 4; mod++){
      cout << endl << modName[mod] << ":\t";
      if(haveGain[mod]){
	for(int roc = 0; roc < 16; roc++){
	  landau_peak[mod][roc] = landau_gauss_peak(&hclq0r[mod][roc]);
	  correction[mod][roc] = 22./landau_peak[mod][roc];
	  if(!CCSupressed || conversionRun == run){
	    if(correction[mod][roc] > 0.0001)ke[mod][roc] *= correction[mod][roc];
	    conversionFileOut << mod << "\t" << roc << "\t" << ke[mod][roc] << endl;
	  }
	  cout << landau_peak[mod][roc] << " ";
	}
      }else{
	cout << "No gain file available.";
      }
    }
    conversionFileOut.close();
  }
  if(CCSupressed){
    cout << "\nNo conversion correction wanted." << endl;
  }






  return 0;
}


Double_t fitLandauGauss( Double_t *x, Double_t *par ) {

  static int nn=0;
  nn++;
  static double xbin = 1;
  static double b1 = 0;
  if( nn == 1 ) { b1 = x[0]; }
  if( nn == 2 ) { xbin = x[0] - b1; } // bin width needed for normalization

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

Double_t landau_gauss_peak(TH1* h) {

  double aa = h->GetEntries();//normalization

  if(h->GetEntries() < 100) return 22.;

  // find peak:
  h->GetXaxis()->SetRange(7, 50);
  int ipk = h->GetMaximumBin();
  h->GetXaxis()->UnZoom();
  double xpk = h->GetBinCenter(ipk);
  double sm = xpk / 9; // sigma
  double ns = sm; // noise

  // fit range:
  int ib0 = h->FindBin(10);
  int ib9 = h->FindBin(40);
  double x0 = h->GetBinLowEdge(ib0);
  double x9 = h->GetBinLowEdge(ib9) + h->GetBinWidth(ib9);

  // create a TF1 with the range from x0 to x9 and 4 parameters
  TF1 *fitFcn = new TF1( "fitFcn", fitLandauGauss, x0, x9, 4 );

  // set start values:
  fitFcn->SetParameter( 0, xpk ); // peak position, defined above
  fitFcn->SetParameter( 1, sm ); // width
  fitFcn->SetParameter( 2, aa ); // area
  fitFcn->SetParameter( 3, ns ); // noise

  h->Fit("fitFcn", "R Q", "ep" );// R = range from fitFcn
  TF1 *fit = h->GetFunction("fitFcn");
  return fit->GetParameter(0);
}
