// Paul Schuetze, DESY, 2016-2018
// based on code of Daniel Pitzl, DESY


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
  double xi,yi,zi;
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
  
  if( y < -5 || y > 5 ) { // Al handle cut out fiducial region
    ffiducial=false;
  }
  
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

	//cout << "Found run " << runnr << ". B = " << Bfield << " T" << endl;

        return 1;
      }
      
    }
    cout << "Run " << runnr << " not found in list. Please add it." << endl;
    parameterFile.close();
    return -2;
  }
}


//------------------------------------------------------------------------------
int searchRunlist(int runnr, double &momentum, int *modName, bool &CCSuppressed, int &alignmentRun, double &turn, double &bfield){

  ifstream runlistFile( "runlist-quad.dat" );

  cout << endl;
  if( runlistFile.bad() || ! runlistFile.is_open() ) {
    cout << "runlist-quad.dat could not be found." << endl;
    return -1;
  }
  else {

    int currentRunnr;
    double currentMomentum;
    int currentModNames[4];
    int currentCCSuppressed = 0;
    double currentTurn;
    double currentBfield;
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
      thisline >> currentCCSuppressed >> alignmentRun >> currentBfield >> currentTurn;
      
      if(!(currentRunnr && currentModNames[0] && currentModNames[1] && currentModNames[2] && currentModNames[3])){
	continue; // No correct data in runlist
      }
      
      if(currentRunnr == runnr){
	//cout << "Found entry in runlist:" << endl;
	//cout << line << endl;
	momentum = currentMomentum;
	for(int mod = 0; mod < 4; mod++){
	  modName[mod] = currentModNames[mod];
	}
	if(alignmentRun == 0) alignmentRun = runnr;
	CCSuppressed = currentCCSuppressed;
	turn = currentTurn;
	bfield = 0.00095*currentBfield;
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


//----------------------------------

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


//----------------------------------


int main( int argc, char* argv[] )
{
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
  int angleRun = 0;
  bool fitSupressed = false;
  double turnTable = 0;

  bool doAngleRun = false;

  bool isDirrad = false;
  
  bool useLCIO = false;

  int skipEvents = -1;
  if(run==3207){
    skipEvents = 40000;
  }else if(run==3223){
    skipEvents = 80000;
  }else if(run==3193){
    skipEvents = 120000;
  }else if(run==3210){
    skipEvents = 30000;//120000;
  }

  // further arguments:

  int lev = 999222111; // last event

  for( int i = 1; i < argc; ++i ) {

    if( !strcmp( argv[i], "-l" ) )
      lev = atoi( argv[++i] );

    if( !strcmp( argv[i], "-s" ) )
      skipEvents = atoi( argv[++i] );

    // Suppress conversion factor correction
    if( !strcmp( argv[i], "-c" ) ){
      if( strchr( argv[++i], '-') == NULL){
	conversionRun = atoi( argv[i] );
	cout << "Conversion factors taken from run " << conversionRun << endl;
      }else{
	--i;
      }
      CCSupressed = true;
    }

    if( !strcmp( argv[i], "-f" ) )
      fitSupressed = true;

    if( !strcmp( argv[i], "-a" ) )
      useLCIO = true;

    if( !strcmp( argv[i], "-x" ) ){
      if( strchr( argv[++i], '-') == NULL){
	angleRun = atoi( argv[i] );
	cout << "Alignment taken from run " << angleRun << endl;
      }else{
	--i;
      }
      doAngleRun = true;
    }

  } // argc

  double bfield;
  if(searchRunlist(run, p, modName, CCSupressed, alignmentRun, turnTable, bfield) < 0){
    exit(0);
  }

  if(modName[3]==4582){
    isDirrad = true;
  }

  if(bfield < 0.01){
    searchBScanParameters(run, bfield);
  }

  std::cout <<
    "Run number: " << run << std::endl <<
    "  Modules: " << modName[0] << " " << modName[1] << " " << modName[2] << " " << modName[3] << std::endl <<
    "  Irradiated: " << isDirrad << std::endl <<
    "  Table turn angle: " << turnTable << std::endl <<
    "  BField: " << bfield << std::endl << 
    "  Alignment run: " << alignmentRun << std::endl << std::endl;

  
  char filename [100];
  sprintf(filename, "data/run%06d.raw", run);

  FileReader * reader = new FileReader (runnum.c_str(), filename);
  
  double dz = 32.;

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

  if(doAngleRun){
    alignmentRun = angleRun;
  }

  if(!alignmentRun){
    alignFileName << "alignmentRespin/align_" << run << ".dat";
  }else{
    alignFileName << "alignmentRespin/align_" << alignmentRun << ".dat";
  }


  ifstream ialignFile( alignFileName.str() );

  bool isFirstAlignment = false;
  cout << endl;
  if( ialignFile.bad() || ! ialignFile.is_open() ) {
    isFirstAlignment = true;
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


  // calculate sine / cosine of alignment parameters:

  double ctu[4] ;
  double stu[4] ;
  double cti[4] ;
  double sti[4] ;
  double cro[4] ;
  double sro[4] ;

  double wt = 180./TMath::Pi();


  for(size_t mod = 0; mod < 4; mod++){
    if(isFirstAlignment){
      tiltA[mod] = 19.3/wt;
      if (mod==99)
	tiltA[mod] = 13./wt;
    }
    ctu[mod] = cos( turnA[mod] );
    stu[mod] = sin( turnA[mod] );
    cti[mod] = cos( tiltA[mod] );
    sti[mod] = sin( tiltA[mod] );
    cro[mod] = cos( rotA[mod] );
    sro[mod] = sin( rotA[mod] );
  }

  double bicutx = 150E-3*dz; // [mm]
  double bicuty = 20E-3*dz; // [mm]

  double tricutx = 0.7; // [mm]
  double tricuty = 0.35; // [mm]

  double chCutLow = 18; // [ke]
  double chCutHigh = 25; // [ke]

  double chCutLooseLow = 14; // [ke]
  double chCutLooseHigh = 40; // [ke]

  double angleCut = 0.15; // [rad]

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
  string milleResName = dirname.str();

  
  const int ndimMP = 5;


  vector<int> labelsA( ndimMP );
  labelsA[0] =  1; // dx
  labelsA[1] =  2; // dy
  labelsA[2] =  3; // drot
  labelsA[3] =  4; // dtilt
  labelsA[4] =  5; // dturn

  vector<int> labelsB( ndimMP );
  labelsB[0] =  1*ndimMP+1; // dx
  labelsB[1] =  1*ndimMP+2; // dy
  labelsB[2] =  1*ndimMP+3; // drot
  labelsB[3] =  1*ndimMP+4; // dtilt
  labelsB[4] =  1*ndimMP+5; // dturn

  vector<int> labelsC( ndimMP );
  labelsC[0] =  2*ndimMP+1; // dx
  labelsC[1] =  2*ndimMP+2; // dy
  labelsC[2] =  2*ndimMP+3; // drot
  labelsC[3] =  2*ndimMP+4; // dtilt
  labelsC[4] =  2*ndimMP+5; // dturn

  vector<int> labelsD( ndimMP );
  labelsD[0] =  3*ndimMP+1; // dx
  labelsD[1] =  3*ndimMP+2; // dy
  labelsD[2] =  3*ndimMP+3; // drot
  labelsD[3] =  3*ndimMP+4; // dtilt
  labelsD[4] =  3*ndimMP+5; // dturn

  // Prepare the mille binary file:                                                                                                                                                                                 
  std::stringstream milleBinName;
  milleBinName << "run" << run << "/mille2.bin";
  string m_millefilename = milleBinName.str();
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

  stringstream gainstream;
  for(int mod = 0; mod < 4; mod++){
    gainstream << "gaincal/D" << modName[mod] << "-tb24-gaincal.dat";
    gainFileName[mod] = gainstream.str();
    gainstream.str("");
  }
  
  double turn = turnTable;

  if(doAngleRun){
    //turnTable = 0.;
  }


  // for GBL:

  double norm = 1.*TMath::Cos(tiltA[0])*TMath::Cos(turnA[0]+turnTable/wt);

  double resx = 9.9E-3; // [mm] col hit resolution
  double resy = 7.7E-3; // [mm] row hit resolution

  //resx = 43.E-3;
  resy = 29.E-3;
  //if(fabs(turn+turnOnTable) > 13) resx = 11.5E-3;
  //resx = (37.-1.0*fabs(turn))*1.E-3;
  resx = (37.-1.0*fabs(turn))*1.E-3;
  if(fabs(turn)>28.){
    resx = (-33.88+1.48*fabs(turn))*1.E-3;
  }

  if(fabs(tiltA[0]) > 15/wt) resy = 7.7E-3;
  //if(fabs(tilt) > 15) resy = 5.0E-3;


  // X0 Si = 21.82/2.33 = 9.365 cm
  // X0 Al = 24.01/2.70 = 8.89 cm
  // X0 Cu = 12.86/8.96 = 1.435 cm
  // X0 air = 36.62/1.204E-3 = 304 m

  // measurement = residual
  TVectorD meas(2);
  meas.Zero(); // ideal

  TVectorD measPrec(2); // precision = 1/resolution^2
  measPrec[0] = 1.0 / resx / resx;
  measPrec[1] = 1.0 / resy / resy;

  TMatrixD proL2m( 2, 2 ); // measurement projection matrix
  proL2m.UnitMatrix();
  
  /*
  proL2m[0][0] = costurn*costurnon - sinturn*sinturnon;
  proL2m[0][1] =-costurn*sinturnon*sintilt - sinturn*costurnon*sintilt;
  proL2m[1][0] = 0;
  proL2m[1][1] = costilt;
  */

  // scatter:
  TVectorD scat(2);
  scat.Zero(); // mean is zero


  double X0Si = 0.75*(( 0.3 + 0.175 + 1.0 ) / 94.) / norm; // Sensor + ROC + HDI(appr. - gives a good kink pull)

  double tetSi = 0.0136 * sqrt(X0Si) / p * ( 1 + 0.038*log(X0Si) );

  TVectorD wscatSi(2);
  wscatSi[0] = 1.0 / ( tetSi * tetSi ); // weight
  wscatSi[1] = 1.0 / ( tetSi * tetSi );


  ostringstream conversionFileName; // output string stream

  //conversionRun = alignmentRun;
  

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

  if(doAngleRun){
    fname << "histogramsAngle/quadAngle-" << run << ".root";
  }else{
    fname << "histogramsRespin/quad-" << run << ".root";
  }

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
  TH2D * hxyLocal[4];
  TH2D * hxzGlobal[4];
  TH1D hnpx[4];
  TH1D hsiz[4];
  TH1D hclq[4];
  TH1D hclq0[4];
  TH1D hclq0r[4][16];
  TH1D hclq0g[4];
  TH1D hclq0t[4];
  TH1D hclq0tstart[4];
  TH1D hclq0tend[4];
  TProfile hclq0tvst[4];
  TH1D hncol[4];
  TH1D hncolq[4];
  TH1D hncolqf4[4];
  TH1D hnrowqf4[4];
  TH1D hncol3qf4[4];
  TH1D hncol4qf4[4];
  TH1D hnrow[4];
  TH1D hnrow4[4];
  TProfile ncolvsclq[4];
  TProfile nrowvsclq[4];

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

    hncl[mod] = TH1D( Form("ncl%c", modtos),
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
    hclq0t[mod] = TH1D( Form("clq0%ct",modtos),
		       Form("%c, On track, normalized cluster charge;norm. cluster charge [ke];%c clusters",modtos,modtos),
		      100, 0, 100 );
    hclq0tstart[mod] = TH1D( Form("clq0%ctstart",modtos),
		       Form("%c, On track, normalized cluster charge, first 100k events;norm. cluster charge [ke];%c clusters",modtos,modtos),
		      100, 0, 100 );
    hclq0tend[mod] = TH1D( Form("clq0%ctend",modtos),
			Form("%c, On track, normalized cluster charge, >900k events;norm. cluster charge [ke];%c clusters",modtos,modtos),
		      100, 0, 100 );
    hclq0tvst[mod] = TProfile( Form("clq0%ctvst",modtos),
			       Form( "clq%cvst", "clq %c vs time;trigger;<charge> %c", modtos, modtos, modtos),
			       500, 0, 1E6, 0, 100 );
    hncol[mod]= TH1D( Form("ncol%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hncolq[mod]= TH1D( Form("ncolq%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hncolqf4[mod]= TH1D( Form("ncolqf4%c",modtos), 
		      Form("%c cluster size;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hncol3qf4[mod]= TH1D( Form("ncol3qf4%c",modtos), 
		      Form("%c cluster size w/ charge cuts in ABC;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hncol4qf4[mod]= TH1D( Form("ncol4qf4%c",modtos), 
		      Form("%c cluster size w/ charge cuts in ABCD;columns/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hnrowqf4[mod]= TH1D( Form("nrowqf4%c",modtos), 
		      Form("%c cluster size w/ charge cuts;rows/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hnrow[mod]= TH1D( Form("nrow%c",modtos),
		      Form("%c cluster size;rows/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    hnrow4[mod]= TH1D( Form("nrow4%c",modtos),
		      Form("%c cluster size, good fits;rows/cluster;%c clusters",modtos,modtos),
		      21, -0.5, 20.5 );
    ncolvsclq[mod]= TProfile( Form("ncolvsclq%c",modtos),
			     Form("ncol%c vs norm cluster charge;cluster charge [ke];ncol %c",modtos,modtos),
			     100,0,100);
    nrowvsclq[mod]= TProfile( Form("nrowvsclq%c",modtos),
			     Form("nrow%c vs norm cluster charge;cluster charge [ke];nrow %c",modtos,modtos),
			     100,0,100);
    hxyGlobal[mod] = new  TH2D( Form("Globalxy%c",modtos),
			   Form("%c cluster map, global var.;x [mm];y [mm];%c clusters",modtos,modtos),
			   500, -40., 40., 160, -10., 10. );
    hxyLocal[mod] = new  TH2D( Form("Localxy%c",modtos),
			   Form("%c cluster map, local var.;x [mm];y [mm];%c clusters",modtos,modtos),
			   432, -32.4, 32.4, 162, -8.1, 8.1 );
    hxzGlobal[mod] = new  TH2D( Form("Globalxz%c",modtos),
			   Form("%c cluster map, global var.;x [mm];z [mm];%c clusters",modtos,modtos),
			   500, -40., 40., 750, -60., 60. );

    
  } // module planes


  histoFile->cd("");
  TH2D * hxzAllGlobal = new  TH2D( "GlobalxzAll",
				 "cluster map, global var.;x [mm];z [mm];clusters",
				 500, -40., 40., 300, -60., 60. );


  // Prealignment:

  //TDirectory * preDir = histoFile->mkdir("Prealignment");
  //preDir->cd();

  TH1D hdxAC( "dxAC", "Ax-Cx;x-x [mm];cluster pairs", 150, -2, 2 );
  TH1D hdyAC( "dyAC", "Ay-Cy;y-y [mm];cluster pairs", 200, -2, 2 );

  TH1D hdxBC( "dxBC", "Bx-Cx;x-x [mm];cluster pairs", 150, -2, 2 );
  TH1D hdyBC( "dyBC", "By-Cy;y-y [mm];cluster pairs", 200, -2, 2 );

  TH1D hdxDC( "dxDC", "dx-Cx;x-x [mm];cluster pairs", 150, -2, 2 );
  TH1D hdyDC( "dyDC", "dy-Cy;y-y [mm];cluster pairs", 200, -2, 2 );

  TProfile dxvsxAC( "dxvsxAC", "A-C dx vs x;x [mm];A-C <dx>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyAC( "dxvsyAC", "A-C dx vs y;y [mm];A-C <dx>",
		    81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxAC( "dyvsxAC", "A-C dy vs x;x [mm];A-C <dy>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyAC( "dyvsyAC", "A-C dy vs y;y [mm];A-C <dy>",
		    81, -8.1, 8.1, -1, 1 );

  TProfile dxvsxBC( "dxvsxBC", "B-C dx vs x;x [mm];B-C <dx>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyBC( "dxvsyBC", "B-C dx vs y;y [mm];B-C <dx>",
		    81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxBC( "dyvsxBC", "B-C dy vs x;x [mm];B-C <dy>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyBC( "dyvsyBC", "B-C dy vs y;y [mm];B-C <dy>",
		    81, -8.1, 8.1, -1, 1 );

  TProfile dxvsxDC( "dxvsxDC", "d-C dx vs x;x [mm];D-C <dx>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dxvsyDC( "dxvsyDC", "d-C dx vs y;y [mm];D-C <dx>",
		    81, -8.1, 8.1, -1, 1 );
  TProfile dyvsxDC( "dyvsxDC", "D-C dy vs x;x [mm];D-C <dy>",
		    216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyDC( "dyvsyDC", "D-C dy vs y;y [mm];D-C <dy>",
		    81, -8.1, 8.1, -1, 1 );

  histoFile->cd("");


  // ACB Histos

  //TDirectory * ACBdir = histoFile->mkdir("ACB");
  //gDirectory->cd("ACB");
  
  TH1D hdxACB( "dxACB", "ACB dx;x-x [mm];ACBs", 200, -1, 1 );
  TH1D hdyACB( "dyACB", "ACB dy;y-y [mm];ACBs", 200, -1, 1 );

  TH1D hdxbACB( "dxbACB", "ACB both tricuts dx;x-x [mm];ACBs", 200, -1, 1 );
  TH1D hdybACB( "dybACB", "ACB both tricuts dy;y-y [mm];ACBs", 200, -1, 1 );

  TH1D hdxcACB( "dxcACB", "ACB dx;x-x [um];ACBs", 200, -700, 700 );
  TH1D hdycACB( "dycACB", "ACB dy;y-y [um];ACBs", 200, -700, 700 );
  TH1D hdxciACB( "dxciACB", "ACB dx;x-x [um];isolated ACBs",
		 200, -700, 700 );
  TH1D hdxcfACB( "dxcfACB", "ACB dx;x-x [um];Landau peak inner ACBs",
		  200, -200, 200 );
  TH1D hdxcfqACB( "dxcfqACB", "ACB dx;x-x [um];Landau peak inner ACBs",
		  200, -200, 200 );
  TH1D hdxcfqtACB( "dxcfqtACB", "ACB dx;x-x [um];Landau peak inner ACBs",
		  200, -200, 200 );
  TH1D hdxcfqtiACB( "dxcfqtiACB", "ACB dx;x-x [um];Landau peak inner ACBs",
		  200, -200, 200 );
  TProfile rmsxACBvsq0("rmsxACBvsq0",
		      "ACB dx vs q0; q0 [ke], ACB <abs(dx)>",
		      80, 0, 80, 0, 50);
  TH1D hdyciACB( "dyciACB", "ACB dy;y-y [um];isolated ACBs",
		 200, -700, 700 );
  TH1D hdycfACB( "dycfACB", "ACB dy;y-y [um];inner ACBs",
		 200, -700, 700 );
  TH1D hdycfqACB( "dycfqACB", "ACB dy;y-y [um];Landau peak inner ACBs",
		  200, -200, 200 );
  TH1D hdycfqtACB( "dycfqtACB", "ACB dy;y-y [um];Landau peak inner ACBs",
		  200, -200, 200 );
  TH1D hdycfqtiACB( "dycfqtiACB", "ACB dy;y-y [um];Landau peak inner ACBs",
		  200, -200, 200 );
  TProfile rmsyACBvsq0("rmsyACBvsq0",
		      "ACB dy vs q0; q0 [ke], ACB <abs(dy)>",
		      80, 0, 80, 0, 50);
  TProfile rmsyACBvst("rmsyACBvst",
		      "ACB dy vs t; Trigger, ACB <abs(dy)>",
		      80, 0, 1e6, 0, 50);
  TProfile dyvsxACB( "dyvsxACB",
		     "ACB dy vs x;x [mm];ACB <dy>",
		     216, -32.4, 32.4, -1, 1 );
  TProfile dyvsyACB( "dyvsyACB",
		     "ACB dy vs y;y [mm];ACB <dy>",
		     81, -8.1, 8.1, -1, 1 );

  TH2D * hmapACB;
  hmapACB = new TH2D( "mapACB",
		      "ACBplet map;ACBplet col;ACBplet row;ACBplets",
		      8*54, -4*8.1, 4*8.1, 2*81, -8.1, 8.1 );


  // Efficiencies B:

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
  TProfile effBvsw( "effBvsw", "effB vs window;link window [mm];eff B",
		     40, 0.025, 2.025, -1, 2 );
  TProfile effBinvvsw( "effBinvvsw", "ineffB vs window;link window [mm];ineff B",
		     40, 0.025, 2.025, -1, 2 );

  // Efficiencies C:

  TProfile effCvsx0( "effCvsx0", "effC vs lower x;lower ADCplet x [mm];eff C",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effCvsx1( "effCvsx1", "effC vs upper x;upper ADCplet x [mm];eff C",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effCvsy( "effCvsy", "effC vs y;ADCplet y [mm];eff C",
		    81, -8.1, 8.1, -1, 2 );
  TProfile effCvst1( "effCvst1", "effC vs time;trigger;eff C",
		     500, 0, 1E6, -1, 2 );
  TProfile effCvst5( "effCvst5", "effC vs time;trigger;eff C",
		     500, 0, 5E6, -1, 2 );
  TProfile effCvst40( "effCvst40", "effC vs time;trigger;eff C",
		      1000, 0, 40E6, -1, 2 );
  TProfile effCvsw( "effCvsw", "effC vs window;link window [mm];eff C",
		     40, 0.025, 2.025, -1, 2 );
  TProfile effCinvvsw( "effCinvvsw", "ineffC vs window;link window [mm];ineff C",
		     40, 0.025, 2.025, -1, 2 );


  // Efficiencies D:

  TProfile effDvsx0( "effDvsx0", "effD vs lower x;lower ADCplet x [mm];eff D",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effDvsx1( "effDvsx1", "effD vs upper x;upper ADCplet x [mm];eff D",
		     216, -32.4, 32.4, -1, 2 );
  TProfile effDvsy( "effDvsy", "effD vs y;ADCplet y [mm];eff D",
		    81, -8.1, 8.1, -1, 2 );
  TProfile effDvst1( "effDvst1", "effD vs time;trigger;eff D",
		     500, 0, 1E6, -1, 2 );
  TProfile effDvst5( "effDvst5", "effD vs time;trigger;eff D",
		     500, 0, 5E6, -1, 2 );
  TProfile effDvst40( "effDvst40", "effD vs time;trigger;eff D",
		      1000, 0, 40E6, -1, 2 );
  TProfile effDvsw( "effDvsw", "effD vs window;link window [mm];eff D",
		     40, 0.025, 2.025, -1, 2 );
  TProfile effDinvvsw( "effDinvvsw", "ineffD vs window;link window [mm];ineff D",
		     40, 0.025, 2.025, -1, 2 );


  TH1D hdxADC( "dxADC", "ADC dx;x-x [mm];ADCplets", 200, -1, 1 );
  TH1D hdyADC( "dyADC", "ADC dy;y-y [mm];ADCplets", 200, -1, 1 );

  TH1D hdxADC4B( "dxADC4B", "ADC4B dx;x-x [mm];ADCplets", 200, -1, 1 );
  TH1D hdyADC4B( "dyADC4B", "ADC4B dy;y-y [mm];ADCplets", 200, -1, 1 );



  histoFile->cd("");


  //TDirectory * gblDir = histoFile->mkdir("GBL");
  //gblDir->cd();

  TH1D fit3ry("fit3ry", "Linear fit residual y;resid y [um];track fits", 200, -150, 150);


  TH1D hchi2( "GBLchi2", "GBL chisq;chisq;track fits", 100, 0, 50 );
  TH1D hndf( "GBLndf", "GBL ndf;ndf;track fits", 50, 0, 50 );
  TH1D hprob( "GBLprob", "GBL fit prob;fit probability;track fits", 100, 0, 1 );

  TH1D derxrotHisto( "derxrot", "dx/drot;dx/drot [mm/mrad];tries", 200, -10.,10.);
  TH1D deryrotHisto( "deryrot", "dy/drot;dy/drot [mm/mrad];tries", 200, -10.,10.);

  TH1D derxtiltHisto( "derxtilt", "dx/dtilt;dx/dtilt [mm/mrad];tries", 200, -10.,10.);
  TH1D derytiltHisto( "derytilt", "dy/dtilt;dy/dtilt [mm/mrad];tries", 200, -10.,10.);

  TH1D derxturnHisto( "derxturn", "dx/dturn;dx/dturn [mm/mrad];tries", 200, -10.,10.);
  TH1D deryturnHisto( "deryturn", "dy/dturn;dy/dturn [mm/mrad];tries", 200, -10.,10.);
  
  TH1D fitrxHisto[4];
  TH1D fitryHisto[4];
  TH1D fitrxAllHisto[4];
  TH1D fitryAllHisto[4];
  TH1D fitrxqHisto[4];
  TH1D fitryqHisto[4];
  TH1D fit2rxHisto[4];
  TH1D fit2ryHisto[4];

  TH1D gblaxSHisto[4];
  TH1D gblaxHisto[4];
  TH1D gblayHisto[4];
  TH1D gbldxHisto[4];
  TH1D gbldyHisto[4];
  TH1D gblrxHisto[4];
  TH1D gblryHisto[4];
  TH1D gblrxqHisto[4];
  TH1D gblryqHisto[4];
  TH1D gblrxqAllHisto[4];
  TH1D gblryqAllHisto[4];
  TH1D gblpxHisto[4];
  TH1D gblpyHisto[4];
  TH1D gblqxHisto[4];
  TH1D gblqyHisto[4];
  TH1D gbltxHisto[4];
  TH1D gbltyHisto[4];
  TH1D gblqpHisto[4];
  TH1D gblpqHisto[4];
  TProfile gblrxvsq0[4];
  TProfile gblryvsq0[4];

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
    fitrxqHisto[mod] = TH1D(Form("fitrxq%c",modtos), Form("Linear fit residual x, charge cut %c;resid x [um];track fits",modtos), 200, -150, 150);
    fitryqHisto[mod] = TH1D(Form("fitryq%c",modtos), Form("Linear fit residual y, charge cut %c;resid y [um];track fits",modtos), 200, -150, 150);
    fitrxAllHisto[mod] = TH1D(Form("fitrxAll%c",modtos), Form("Linear fit residual x %c, not only GBL-tracks;resid x [um];track fits",modtos), 200, -150, 150);
    fitryAllHisto[mod] = TH1D(Form("fitryAll%c",modtos), Form("Linear fit residual y %c, not only GBL-tracks;resid y [um];track fits",modtos), 200, -150, 150);


    gblaxSHisto[mod] = TH1D(Form("gblaxS%c",modtos), Form("GBL angle x %c w/o mx;angle x [mrad];track fits",modtos), 500, -100, 100);
    gblaxHisto[mod] = TH1D(Form("gblax%c",modtos), Form("GBL angle x %c;angle x [mrad];track fits",modtos), 500, -100, 100);
    gblayHisto[mod] = TH1D(Form("gblay%c",modtos), Form("GBL angle y %c;angle y [mrad];track fits",modtos), 100, -5, 5);
    gbldxHisto[mod] = TH1D(Form("gbldx%c",modtos), Form("GBL shift x %c;shift x [um];track fits",modtos), 100, -200, 200);
    gbldyHisto[mod] = TH1D(Form("gbldy%c",modtos), Form("GBL shift y %c;shift y [um];track fits",modtos), 100, -100, 100);
    gblrxHisto[mod] = TH1D(Form("gblrx%c",modtos), Form("GBL residual x %c;resid x [um];track fits",modtos), 100, -50, 50);
    gblryHisto[mod] = TH1D(Form("gblry%c",modtos), Form("GBL residual y %c;resid y [um];track fits",modtos), 100, -50, 50);
    gblrxqHisto[mod] = TH1D(Form("gblrxq%c",modtos), Form("GBL residual x, q cut %c;resid x [um];track fits",modtos), 100, -50, 50);
    gblryqHisto[mod] = TH1D(Form("gblryq%c",modtos), Form("GBL residual y, q cut %c;resid y [um];track fits",modtos), 100, -50, 50);
    gblrxqAllHisto[mod] = TH1D(Form("gblrxqAll%c",modtos), Form("GBL residual x, all mod q cut %c;resid x [um];track fits",modtos), 100, -50, 50);
    gblryqAllHisto[mod] = TH1D(Form("gblryqAll%c",modtos), Form("GBL residual y, all mod q cut %c;resid y [um];track fits",modtos), 100, -50, 50);
    gblpxHisto[mod] = TH1D(Form("gblpx%c",modtos), Form("GBL pull x %c;pull x [sigma];track fits",modtos), 100, -10, 10);
    gblpyHisto[mod] = TH1D(Form("gblpy%c",modtos), Form("GBL pull y %c;pull y [sigma];track fits",modtos), 100, -10, 10);
    gblqxHisto[mod] = TH1D(Form("gblqx%c",modtos), Form("GBL kink x %c;kink x [mrad];track fits",modtos), 100, -1, 1);
    gblqyHisto[mod] = TH1D(Form("gblqy%c",modtos), Form("GBL kink y %c;kink y [mrad];track fits",modtos), 100, -1, 1);
    gbltxHisto[mod] = TH1D(Form("gbltx%c",modtos), Form("GBL kink pull x %c;kink pull;track fits",modtos), 100, -10, 10);
    gbltyHisto[mod] = TH1D(Form("gblty%c",modtos), Form("GBL kink pull y %c;kink pull;track fits",modtos), 100, -10, 10);

    
    gblqpHisto[mod] = TH1D(Form("gblqp%c",modtos), Form("GBL p %c;p;track fits",modtos), 600, 0., 30.);
    gblpqHisto[mod] = TH1D(Form("gblpq%c",modtos), Form("GBL 1/p %c;p;track fits",modtos), 600, 0., 2.);
    
    gblrxvsq0[mod] = TProfile(Form("gblrx%cvsq0",modtos), Form("GBL residual x vs q0, %c; q0 [ke]; resid x [um]",modtos), 80, 0, 80, 0, 50);
    gblryvsq0[mod] = TProfile(Form("gblry%cvsq0",modtos), Form("GBL residual y vs q0, %c; q0 [ke]; resid y [um]",modtos), 80, 0, 80, 0, 50);
      
    
  } // module planes


  //histoFile->cd("");




  //TDirectory * fitDir = histoFile->mkdir("Linefit");
  //fitDir->cd();

  TH1D hmx4fit( "mx4fit", "4fit fit slope x;slope [mrad];4plets", 200, -50, 50 );
  TH1D hmx4fitall( "mx4fitall", "4fit fit slope x;slope [mrad];4plets", 1200, -300, 300 );
  TH1D hmx4fitchi( "mx4fitchi", "4fit fit slope x;slope [mrad];4plets", 1200, -300, 300 );
  TH1D hmx4fitchiinv( "mx4fitchiinv", "4fit fit slope x;slope [mrad];4plets", 1200, -300, 300 );
  TH1D hmx4fitl( "mx4fitl", "4fit fit slope x;slope [mrad];4plets", 1200, -300, 300 );
  TH1D hmx4fitlused( "mx4fitlused", "4fit fit slope x, used tracks;slope [mrad];4plets", 3000, -300, 300 );
  TH1D hmx4fitlusedinv( "mx4fitlusedinv", "4fit fit slope x, used tracks, inverted;slope [mrad];4plets", 3000, -300, 300 );
  TH1D hbx4fit( "bx4fit", "4fit fit x0;fitted x0 [mm];4plets", 200, -32, 32 );
  TH1D hmy4fit( "my4fit", "4fit fit slope y;slope [mrad];4plets", 200, -200, 200 );
  TH1D hby4fit( "by4fit", "4fit fit y0;fitted y0 [mm];4plets", 200, -8, 8 );

  TH1D hsigmax4fit( "sx4fit", "4fit fit sigma x;sigma x [um];4plets", 200, 0, 100 );
  TH1D hsigmay4fit( "sy4fit", "4fit fit sigma y;sigma y [um];4plets", 200, 0, 100 );

  TH1D hchi2fitUsed( "chi2fitUsed", "chi2 linear fit;chi2 [um*2];4plets", 500, 0, 25000 );
  TH1D hchi2fitPassed( "chi2fitPassed", "chi2 linear fit;chi2 [um*2];4plets", 500, 0, 25000 );
  TH1D hchi2fitFailed( "chi2fitFailed", "chi2 linear fit;chi2 [um*2];4plets", 1000, 0, 100000 );

  TProfile mx4vsxfit( "mx4vsxfit", "Track slope vs position;x0 [mm];<slope x> [mrad]",
		      200, 32, -32, -300,300 );

  //histoFile->cd("");

  // Events

  TH1D hn4candidates( "hn4candidates", "Track candidates; Nr. of track candidates;Events", 25, -0.5, 24.5 );
  TH1D hn4good( "hn4good", "Fitted tracks; Nr. of good track fits;Events", 25, -0.5, 24.5 );




  // - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
  // event loop:

  int event_nr = 0;

  int n4 = 0;
  int nmille = 0;
  int ntry = 0;
  
  int n4ev = 0;
  
  do {

    vector <cluster> cl[4];


    // Get next event:
    DetectorEvent evt = reader->GetDetectorEvent();

    if( evt.IsBORE() ){
      eudaq::PluginManager::Initialize(evt);
      cout << "I'm BOREing" << endl;
    }

    if( event_nr < skipEvents ){
      if( event_nr%10000 == 0 ){
	cout << "Skipping event (raw) "<< event_nr << endl;
      }
      event_nr++;
      reader->NextEvent();
      continue;
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

	if(modName[mod]==4582 && ((223 <= xm && xm <= 224 && 109 <= ym && ym <= 121) || (225 <= xm && xm <= 226 && 114 <= ym && ym <= 117)) ){
	  continue;
	} // Noisy spot

	  // Shift adc to known threshold

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
      
      for( vector<cluster>::iterator clus = cl[mod].begin(); clus != cl[mod].end(); ++clus ) {

	hsiz[mod].Fill( clus->size );
	hclq[mod].Fill( clus->charge );
	hclq0[mod].Fill( clus->charge*norm );
	hncol[mod].Fill( clus->ncol );
	hnrow[mod].Fill( clus->nrow );

	double x1 = clus->col*0.15 - 32.4;
	double y1 = -1.*(clus->row*0.10 - 8.1);
	double z1 = 0.;
	
	hxyLocal[mod]->Fill( x1, y1 ); // Local, front view

	
	// Rot around z
	
	
	double x2 = cro[mod]*x1 - sro[mod]*y1;
	double y2 = sro[mod]*x1 + cro[mod]*y1;
	double z2 = z1;

	// Turn around y

	double x3 = ctu[mod]*x2 + stu[mod]*z2;
	double y3 = y2;
	double z3 =-stu[mod]*x2 + ctu[mod]*z2;

	// Tilt around x
	
	double x4 = x3;
	double y4 = cti[mod]*y3 - sti[mod]*z3;
	double z4 = sti[mod]*y3 + cti[mod]*z3;


	// Prealignment & Mille:
	
	double x5 = x4 - alignx[mod];
	double y5 = y4 - aligny[mod];
	double z5 = z4 + dz*mod - 1.5*dz;


	if(mod==3)
	  z5 -= 2.;
	

	clus->xi = x1;
	clus->yi = y1;
	clus->zi = z1;

	clus->xg = x5;
	clus->yg = y5;
	clus->zg = z5;

	hclx[mod].Fill( x5 );
	hcly[mod].Fill( y5 );
	  
	ncolvsclq[mod].Fill( clus->charge*norm, clus->ncol );
	nrowvsclq[mod].Fill( clus->charge*norm, clus->nrow );
	if(clus->roc == -1){
	  hclq0g[mod].Fill( clus->charge*norm );
	}else{
	  hclq0r[mod][clus->roc].Fill( clus->charge*norm );
	}
	if(clus->charge*norm > chCutLow && clus->charge*norm < chCutHigh){
	  hncolq[mod].Fill( clus->ncol );
	}

	hxyGlobal[mod]->Fill( x5, y5 ); // front view
	hxzGlobal[mod]->Fill( x5, z5 ); // top view

	hxzAllGlobal->Fill( x5, z5 ); // top view, all


      } // clusters


    } // modules
    if(!(reader->NextEvent())) break;
    ++event_nr;

    // Now let's do sth with this :-)


    // Prealignment:

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {
      double xA = cA->xg;
      double yA = cA->yg;
      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {
	double xC = cC->xg;
	double yC = cC->yg;
	double dx = xA-xC;
	double dy = yA-yC;
	hdxAC.Fill(dx);
	hdyAC.Fill(dy);
	dxvsxAC.Fill(xA,dx);
	dxvsyAC.Fill(yA,dx);
	dyvsxAC.Fill(xA,dy);
	dyvsyAC.Fill(yA,dy);
	
      }
    }

    for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {
      double xB = cB->xg;
      double yB = cB->yg;
      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {
	double xC = cC->xg;
	double yC = cC->yg;
	double dx = xB-xC;
	double dy = yB-yC;
	hdxBC.Fill(dx);
	hdyBC.Fill(dy);
	dxvsxBC.Fill(xB,dx);
	dxvsyBC.Fill(yB,dx);
	dyvsxBC.Fill(xB,dy);
	dyvsyBC.Fill(yB,dy);
	
      }
    }

    for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {
      double xD = cD->xg;
      double yD = cD->yg;
      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {
	double xC = cC->xg;
	double yC = cC->yg;
	double dx = xD-xC;
	double dy = yD-yC;
	hdxDC.Fill(dx);
	hdyDC.Fill(dy);
	dxvsxDC.Fill(xD,dx);
	dxvsyDC.Fill(yD,dx);
	dyvsxDC.Fill(xD,dy);
	dyvsyDC.Fill(yD,dy);
	
      }
    }

    // ACB triplet
    bool iso = cl[A].size() == 1 && cl[C].size() == 1;

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      double xA = cA->xg;
      double yA = cA->yg;
      double zA = cA->zg;

      double qA = cA->charge*norm;
      bool lqA = 1;
      if(      qA < chCutLow ) lqA = 0;
      else if( qA > chCutHigh ) lqA = 0;

      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	double xC = cC->xg;
	double yC = cC->yg;
	double zC = cC->zg;

	double qC = cC->charge*norm;
	bool lqC = 1;
	if(      qC < chCutLow ) lqC = 0;
	else if( qC > chCutHigh ) lqC = 0;

	// A-C track:

	// Where does the track hit B along z?

	double zBi = -0.5*dz;
	double nxB = stu[B];
	double nyB = -sti[B]*ctu[B];
	double nzB = cti[B]*ctu[B];
	double lx = xC-xA;
	double ly = yC-yA;
	double lz = zC-zA;

	double dIntersect = ((-xA)*nxB + (-yA)*nyB + (zBi-zA)*nzB)/(lx*nxB + ly*nyB + lz*nzB);
	// dIntersect is vecA + dIntersect*(vecC-vecA)
	double xpB = xA+dIntersect*lx;
	double ypB = yA+dIntersect*ly;
	double zpB = zA+dIntersect*lz;

	double xavg2B = xA + (zpB-zA)*(xC-xA)/(zC-zA);
	double yavg2B = yA + (zpB-zA)*(yC-yA)/(zC-zA);

	//double xavg2B = 0.5*(xA + xC); // interpolate
	//double yavg2B = 0.5*(yA + yC); // equidistant


	double slpx = (xC - xA)/2/dz + turnTable/wt; // angle
	double slpy = (yC - yA)/2/dz; // angle

	for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	  double xB = cB->xg;
	  double yB = cB->yg;
	  double zB = cB->zg;

	  double qB = cB->charge*norm;
	  bool lqB = 1;
	  if(      qB < chCutLow ) lqB = 0;
	  else if( qB > chCutHigh ) lqB = 0;

	  // tri ACB:

	  //double dx3 = xB - xavg2B;
	  //double dy3 = yB - yavg2B;

	  //double dx3 = (xB - xavg2B)/ctu[B];
	  //double dy3 = (yB - yavg2B)/cti[B];

	  double dx3 = (xB - xpB);
	  double dy3 = (yB - ypB);

	  hdxACB.Fill( dx3 );
	  hdyACB.Fill( dy3 );
	  
	  // Resolution studies
	  
	  if( abs( dy3 ) < tricuty
	      && cA->big == 0 && cC->big == 0 && cB->big == 0 ) {
	    hdxcACB.Fill( dx3*1E3 );

	    if(fabs(slpx) < angleCut && fabs(slpy) < angleCut && iso){
	      rmsxACBvsq0.Fill( qB , fabs(dx3)*1E3 );
	    }
	    if(isFiducial(xpB,ypB)){
	      hdxcfACB.Fill( dx3*1E3 );

	      if( lqA && lqC && lqB){
		hdxcfqACB.Fill( dx3*1E3 );
		if(fabs(slpx) < angleCut && fabs(slpy) < angleCut ){
		  hdxcfqtACB.Fill( dx3*1E3 );
		  if(iso){
		    //if(cA->ncol > 1 && cB->ncol > 1 && cC->ncol > 1)
		    hdxcfqtiACB.Fill( dx3*1E3 );
		  }
		}
	      }
	    }
	  }

	  if( abs( dx3 ) < tricutx
	      && cA->big == 0 && cC->big == 0 && cB->big == 0 ) {
	    hdycACB.Fill( dy3*1E3 );
	    
	    if(fabs(slpx) < angleCut && fabs(slpy) < angleCut && iso){
	      rmsyACBvsq0.Fill( qB , fabs(dy3)*1E3 );
	    }

	    if(isFiducial(xpB,ypB)){
	      hdycfACB.Fill( dy3*1E3 );
	      if( lqA && lqC && lqB){
		hdycfqACB.Fill( dy3*1E3 );
		if(fabs(slpx) < angleCut && fabs(slpy) < angleCut ){
		  hdycfqtACB.Fill( dy3*1E3 );
		  if(iso){
		    //if(cA->ncol > 1 && cB->ncol > 1 && cC->ncol > 1)
		    hdycfqtiACB.Fill( dy3*1E3 );
		    rmsyACBvst.Fill( event_nr , fabs(dy3)*1E3 );
		  }
		}
	      }
	    }
	  }

	  if( abs( dx3 ) > tricutx || abs( dy3 ) > tricuty  ){
	    continue; // tight tri
	  }

	  dyvsxACB.Fill( xpB, fabs(dy3)*1E3 );
	  dyvsyACB.Fill( ypB, fabs(dy3)*1E3 );
	  hmapACB->Fill( xavg2B, yavg2B );


	}
      }
    }


    // ACD -> B for efficiency
    iso = cl[A].size() == 1 && cl[C].size() == 1 && cl[D].size() == 1;
    
    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      bool lq[4] = {true,true,true,true};

      double xA = cA->xg;
      double yA = cA->yg;

      double qA = cA->charge*norm;
      bool lqA = 1;
      if(      qA < chCutLow ) lqA = 0;
      else if( qA > chCutHigh ) lqA = 0;
      lq[A] = lqA;

      for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	double xD = cD->xg;
	double yD = cD->yg;

	double qD = cD->charge*norm;
	bool lqD = 1;
	if(      qD < chCutLow ) lqD = 0;
	else if( qD > chCutHigh ) lqD = 0;
	lq[D] = lqD;

	double slpx = (xD - xA)/3/dz  + turnTable/wt; // angle
	double slpy = (yD - yA)/3/dz; // angle


	if(fabs(slpx) > angleCut || fabs(slpy) > angleCut){
	  continue;
	}	


	double xavg2C = (xA + 2.*xD)/3.; // interpolate
	double yavg2C = (yA + 2.*yD)/3.; // equidistant


	for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	  double xC = cC->xg;
	  double yC = cC->yg;

	  double qC = cC->charge*norm;
	  bool lqC = 1;
	  if(      qC < chCutLow ) lqC = 0;
	  else if( qC > chCutHigh ) lqC = 0;
	  lq[C] = lqC;
	  
	  double dx3 = xC-xavg2C;
	  double dy3 = yC-yavg2C;

	  hdxADC.Fill(dx3);
	  hdyADC.Fill(dy3);
	  
	  if(fabs(dx3)>tricutx || fabs(dy3)>tricuty){
	    continue;
	  }

	  if(!iso){
	    continue;
	  }

	  // This is a valid triplet ADC - is it also there in B?
	  
	  double xavg2B = (2.*xA + xD)/3.; // interpolate
	  double yavg2B = (2.*yA + yD)/3.; // equidistant

	  if(!isFiducial(xavg2B,yavg2B)){
	    continue;
	  }

	  int nm[99] = {0};

	  for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	    double xB = cB->xg;
	    double yB = cB->yg;

	    double qB = cB->charge*norm;
	    bool lqB = 1;
	    if(      qB < chCutLow ) lqB = 0;
	    else if( qB > chCutHigh ) lqB = 0;
	    lq[B] = lqB;

	    double dx4 = xB-xavg2B;
	    double dy4 = yB-yavg2B;

	    hdxADC4B.Fill(dx4);
	    hdyADC4B.Fill(dy4);

	    for( int iw = 1; iw < 41; ++ iw ){
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ){ // for eff
		nm[iw] = 1;
	      }
	    }
	  }

	  int effCut = 10;

	  effBvsx0.Fill( xA, nm[effCut] );
	  effBvsx1.Fill( xD, nm[effCut] );
	  effBvsy.Fill( yA, nm[effCut] );

	  effBvst1.Fill( event_nr, nm[effCut] );
	  effBvst5.Fill( event_nr, nm[effCut] );
	  effBvst40.Fill( event_nr, nm[effCut] );

	  for( int iw = 1; iw < 41; ++ iw ){
	    effBvsw.Fill( iw*0.050+0.005, nm[iw] );
	    effBinvvsw.Fill( iw*0.050+0.005, 1-nm[iw] );
	  }
	}
      }
    }


    // Same efficiency for C
    
    iso = cl[A].size() == 1 && cl[B].size() == 1 && cl[D].size() == 1;
    
    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      bool lq[4] = {true,true,true,true};

      double xA = cA->xg;
      double yA = cA->yg;

      double qA = cA->charge*norm;
      bool lqA = 1;
      if(      qA < chCutLow ) lqA = 0;
      else if( qA > chCutHigh ) lqA = 0;
      lq[A] = lqA;

      for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	double xD = cD->xg;
	double yD = cD->yg;

	double qD = cD->charge*norm;
	bool lqD = 1;
	if(      qD < chCutLow ) lqD = 0;
	else if( qD > chCutHigh ) lqD = 0;
	lq[D] = lqD;

	double slpx = (xD - xA)/3/dz  + turnTable/wt; // angle
	double slpy = (yD - yA)/3/dz; // angle


	if(fabs(slpx) > angleCut || fabs(slpy) > angleCut){
	  continue;
	}	


	double xavg2B = (2.*xA + xD)/3.; // interpolate
	double yavg2B = (2.*yA + yD)/3.; // equidistant


	for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	  double xB = cB->xg;
	  double yB = cB->yg;

	  double qB = cB->charge*norm;
	  bool lqB = 1;
	  if(      qB < chCutLow ) lqB = 0;
	  else if( qB > chCutHigh ) lqB = 0;
	  lq[B] = lqB;
	  
	  double dx3 = xB-xavg2B;
	  double dy3 = yB-yavg2B;
	  
	  if(fabs(dx3)>tricutx || fabs(dy3)>tricuty){
	    continue;
	  }

	  if(!iso){
	    continue;
	  }

	  // This is a valid triplet ADC - is it also there in B?
	  
	  double xavg2C = (2.*xA + xD)/3.; // interpolate
	  double yavg2C = (2.*yA + yD)/3.; // equidistant

	  if(!isFiducial(xavg2C,yavg2C)){
	    continue;
	  }

	  int nm[99] = {0};

	  for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	    double xC = cC->xg;
	    double yC = cC->yg;

	    double qC = cC->charge*norm;
	    bool lqC = 1;
	    if(      qC < chCutLow ) lqC = 0;
	    else if( qC > chCutHigh ) lqC = 0;
	    lq[C] = lqC;

	    double dx4 = xC-xavg2C;
	    double dy4 = yC-yavg2C;

	    for( int iw = 1; iw < 41; ++ iw ){
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ){ // for eff
		nm[iw] = 1;
	      }
	    }
	  }

	  int effCut = 10;

	  effCvsx0.Fill( xA, nm[effCut] );
	  effCvsx1.Fill( xD, nm[effCut] );
	  effCvsy.Fill( yA, nm[effCut] );

	  effCvst1.Fill( event_nr, nm[effCut] );
	  effCvst5.Fill( event_nr, nm[effCut] );
	  effCvst40.Fill( event_nr, nm[effCut] );

	  for( int iw = 1; iw < 41; ++ iw ){
	    effCvsw.Fill( iw*0.050+0.005, nm[iw] );
	    effCinvvsw.Fill( iw*0.050+0.005, 1-nm[iw] );
	  }
	}
      }
    }


    // Same efficiency for D
    
    iso = cl[A].size() == 1 && cl[B].size() == 1 && cl[C].size() == 1;
    
    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      bool lq[4] = {true,true,true,true};

      double xA = cA->xg;
      double yA = cA->yg;
      double zA = cA->zg;

      double qA = cA->charge*norm;
      bool lqA = 1;
      if(      qA < chCutLow ) lqA = 0;
      else if( qA > chCutHigh ) lqA = 0;
      lq[A] = lqA;

      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	double xC = cC->xg;
	double yC = cC->yg;
	double zC = cC->zg;

	double qC = cC->charge*norm;
	bool lqC = 1;
	if(      qC < chCutLow ) lqC = 0;
	else if( qC > chCutHigh ) lqC = 0;
	lq[C] = lqC;

	double slpx = (xC - xA)/2/dz  + turnTable/wt; // angle
	double slpy = (yC - yA)/2/dz; // angle


	if(fabs(slpx) > angleCut || fabs(slpy) > angleCut){
	  continue;
	}	


	double xavg2B = (xA + xC)/2.; // interpolate
	double yavg2B = (yA + yC)/2.; // equidistant


	for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	  double xB = cB->xg;
	  double yB = cB->yg;
	  double zB = cB->zg;

	  double qB = cB->charge*norm;
	  bool lqB = 1;
	  if(      qB < chCutLow ) lqB = 0;
	  else if( qB > chCutHigh ) lqB = 0;
	  lq[B] = lqB;
	  
	  double dx3 = xB-xavg2B;
	  double dy3 = yB-yavg2B;
	  
	  if(fabs(dx3)>tricutx || fabs(dy3)>tricuty){
	    continue;
	  }

	  if(!iso){
	    continue;
	  }

	  // This is a valid triplet ADC - is it also there in B?
	  
	  double xavg2D = (-1.*xA + 3.*xC)/2.; // interpolate
	  double yavg2D = (-1.*yA + 3.*yC)/2.; // equidistant

	  if(!isFiducial(xavg2D,yavg2D)){
	    continue;
	  }

	  int nm[99] = {0};

	  double my3,by3;
	  linefit(zA,yA,zB,yB,zC,yC,my3,by3);

	  for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	    double xD = cD->xg;
	    double yD = cD->yg;
	    double zD = cD->zg;

	    double qD = cD->charge*norm;
	    bool lqD = 1;
	    if(      qD < chCutLow ) lqD = 0;
	    else if( qD > chCutHigh ) lqD = 0;
	    lq[D] = lqD;

	    double dx4 = xD-xavg2D;
	    double dy4 = yD-yavg2D;

	    //Fit residual:
	    fit3ry.Fill((yD-(my3*zD+by3))*1.E3);

	    for( int iw = 1; iw < 41; ++ iw ){
	      if( abs( dx4 ) < iw*0.050 && abs( dy4 ) < iw*0.050 ){ // for eff
		nm[iw] = 1;
	      }
	    }
	  }

	  int effCut = 10;

	  effDvsx0.Fill( xavg2D, nm[effCut] );
	  effDvsx1.Fill( xavg2D, nm[effCut] );
	  effDvsy.Fill( yavg2D, nm[effCut] );

	  effDvst1.Fill( event_nr, nm[effCut] );
	  effDvst5.Fill( event_nr, nm[effCut] );
	  effDvst40.Fill( event_nr, nm[effCut] );

	  for( int iw = 1; iw < 41; ++ iw ){
	    effDvsw.Fill( iw*0.050+0.005, nm[iw] );
	    effDinvvsw.Fill( iw*0.050+0.005, 1-nm[iw] );
	  }
	}
      }
    }


    
    // Make a 4plet

    int nCandidatesPerEvent = 0;
    int nGoodPerEvent = 0;

    bool lq[4] = {true,true,true,true};

    for( vector<cluster>::iterator cA = cl[A].begin(); cA != cl[A].end(); ++cA ) {

      double xA = cA->xg;
      double yA = cA->yg;

      double qA = cA->charge*norm;
      bool lqA = 1;
      if(      qA < chCutLow ) lqA = 0;
      else if( qA > chCutHigh ) lqA = 0;
      lq[A] = lqA;

      bool lqAl = 1;
      if(      qA < chCutLooseLow ) lqAl = 0;
      else if( qA > chCutLooseHigh ) lqAl = 0;

      double xA1 = cA->xi;
      double yA1 = cA->yi;
      double zA1 = cA->zi;
      
      TMatrixD derivA( 2, ndimMP ); // -alignment derivatives x,y

      //d/dx
      derivA[0][0] = 1.;
      derivA[1][0] = 0.;

      //d/dy
      derivA[0][1] = 0.;
      derivA[1][1] = 1.;

      //d/drot
      derivA[0][2] = -1.*(-ctu[A]*sro[A]*xA1 - ctu[A]*cro[A]*yA1);
      derivA[1][2] = -1.*(cti[A]*cro[A]*xA1 - cti[A]*sro[A]*yA1 - sti[A]*stu[A]*sro[A]*xA1 - sti[A]*stu[A]*cro[A]*yA1);

      //d/dtilt
      derivA[0][3] = 0.;
      derivA[1][3] = -1.*(-sti[A]*sro[A]*xA1 - sti[A]*cro[A]*yA1 + cti[A]*stu[A]*cro[A]*xA1 - cti[A]*stu[A]*sro[A]*yA1 - cti[A]*ctu[A]*zA1);

      //d/dturn
      derivA[0][4] = -1.*(-stu[A]*cro[A]*xA1 + stu[A]*sro[A]*yA1 + ctu[A]*zA1);
      derivA[1][4] = -1.*(sti[A]*ctu[A]*cro[A]*xA1 - sti[A]*ctu[A]*sro[A]*yA1 + sti[A]*stu[A]*zA1);


      for( vector<cluster>::iterator cC = cl[C].begin(); cC != cl[C].end(); ++cC ) {

	double xC = cC->xg;
	double yC = cC->yg;

	double qC = cC->charge*norm;
	bool lqC = 1;
	if(      qC < chCutLow ) lqC = 0;
	else if( qC > chCutHigh ) lqC = 0;
	lq[C] = lqC;


	bool lqCl = 1;
	if(      qC < chCutLooseLow ) lqCl = 0;
	else if( qC > chCutLooseHigh ) lqCl = 0;

	double slpx = (xC - xA)/2/dz  + turnTable/wt; // angle
	double slpy = (yC - yA)/2/dz; // angle

	if(fabs(slpx) > angleCut || fabs(slpy) > angleCut){
	  continue;
	}

	double xavg2B = 0.5*(xA + xC); // interpolate
	double yavg2B = 0.5*(yA + yC); // equidistant

	double xC1 = cC->xi;
	double yC1 = cC->yi;
	double zC1 = cC->zi;
      
	TMatrixD derivC( 2, ndimMP ); // -alignment derivatives x,y

	//d/dx
	derivC[0][0] = 1.;
	derivC[1][0] = 0.;

	//d/dy
	derivC[0][1] = 0.;
	derivC[1][1] = 1.;

	//d/drot
	derivC[0][2] = -1.*(-ctu[C]*sro[C]*xC1 - ctu[C]*cro[C]*yC1);
	derivC[1][2] = -1.*(cti[C]*cro[C]*xC1 - cti[C]*sro[C]*yC1 - sti[C]*stu[C]*sro[C]*xC1 - sti[C]*stu[C]*cro[C]*yC1);

	//d/dtilt
	derivC[0][3] = 0.;
	derivC[1][3] = -1.*(-sti[C]*sro[C]*xC1 - sti[C]*cro[C]*yC1 + cti[C]*stu[C]*cro[C]*xC1 - cti[C]*stu[C]*sro[C]*yC1 - cti[C]*ctu[C]*zC1);

	//d/dturn
	derivC[0][4] = -1.*(-stu[C]*cro[C]*xC1 + stu[C]*sro[C]*yC1 + ctu[C]*zC1);
	derivC[1][4] = -1.*(sti[C]*ctu[C]*cro[C]*xC1 - sti[C]*ctu[C]*sro[C]*yC1 + sti[C]*stu[C]*zC1);
	

	for( vector<cluster>::iterator cB = cl[B].begin(); cB != cl[B].end(); ++cB ) {

	  double xB = cB->xg;
	  double yB = cB->yg;

	  double qB = cB->charge*norm;
	  bool lqB = 1;
	  if(      qB < chCutLow ) lqB = 0;
	  else if( qB > chCutHigh ) lqB = 0;
	  lq[B] = lqB;

	  bool lqBl = 1;
	  if(      qB < chCutLooseLow ) lqBl = 0;
	  else if( qB > chCutLooseHigh ) lqBl = 0;

	  double dx3 = xB - xavg2B;
	  double dy3 = yB - yavg2B;

	  if(fabs(dx3) > tricutx || fabs(dy3) > tricuty){
	    continue;
	  }

	  hdxbACB.Fill( dx3 );
	  hdybACB.Fill( dy3 );	  

	  double xB1 = cB->xi;
	  double yB1 = cB->yi;
	  double zB1 = cB->zi;
      
	  TMatrixD derivB( 2, ndimMP ); // -alignment derivatives x,y

	  //d/dx
	  derivB[0][0] = 1.;
	  derivB[1][0] = 0.;

	  //d/dy
	  derivB[0][1] = 0.;
	  derivB[1][1] = 1.;

	  //d/drot
	  derivB[0][2] = -1.*(-ctu[B]*sro[B]*xB1 - ctu[B]*cro[B]*yB1);
	  derivB[1][2] = -1.*(cti[B]*cro[B]*xB1 - cti[B]*sro[B]*yB1 - sti[B]*stu[B]*sro[B]*xB1 - sti[B]*stu[B]*cro[B]*yB1);

	  //d/dtilt
	  derivB[0][3] = 0.;
	  derivB[1][3] = -1.*(-sti[B]*sro[B]*xB1 - sti[B]*cro[B]*yB1 + cti[B]*stu[B]*cro[B]*xB1 - cti[B]*stu[B]*sro[B]*yB1 - cti[B]*ctu[B]*zB1);

	  //d/dturn
	  derivB[0][4] = -1.*(-stu[B]*cro[B]*xB1 + stu[B]*sro[B]*yB1 + ctu[B]*zB1);
	  derivB[1][4] = -1.*(sti[B]*ctu[B]*cro[B]*xB1 - sti[B]*ctu[B]*sro[B]*yB1 + sti[B]*stu[B]*zB1);

	  derxrotHisto.Fill(derivB[0][2]);
	  deryrotHisto.Fill(derivB[1][2]);
	  derxtiltHisto.Fill(derivB[0][3]);
	  derytiltHisto.Fill(derivB[1][3]);
	  derxturnHisto.Fill(derivB[0][4]);
	  deryturnHisto.Fill(derivB[1][4]);

	  
	  for( vector<cluster>::iterator cD = cl[D].begin(); cD != cl[D].end(); ++cD ) {

	    double xD = cD->xg;
	    double yD = cD->yg;

	    double qD = cD->charge*norm;
	    bool lqD = 1;
	    if(      qD < chCutLow ) lqD = 0;
	    else if( qD > chCutHigh ) lqD = 0;
	    lq[D] = lqD;

	    bool lqDl = 1;
	    if(      qD < chCutLooseLow ) lqDl = 0;
	    else if( qD > chCutLooseHigh ) lqDl = 0;

	    double xD1 = cD->xi;
	    double yD1 = cD->yi;
	    double zD1 = cD->zi;
      
	    TMatrixD derivD( 2, ndimMP ); // -alignment derivatives x,y

	    //d/dx
	    derivD[0][0] = 1.;
	    derivD[1][0] = 0.;

	    //d/dy
	    derivD[0][1] = 0.;
	    derivD[1][1] = 1.;

	    //d/drot
	    derivD[0][2] = -1.*(-ctu[D]*sro[D]*xD1 - ctu[D]*cro[D]*yD1);
	    derivD[1][2] = -1.*(cti[D]*cro[D]*xD1 - cti[D]*sro[D]*yD1 - sti[D]*stu[D]*sro[D]*xD1 - sti[D]*stu[D]*cro[D]*yD1);

	    //d/dtilt
	    derivD[0][3] = 0.;
	    derivD[1][3] = -1.*(-sti[D]*sro[D]*xD1 - sti[D]*cro[D]*yD1 + cti[D]*stu[D]*cro[D]*xD1 - cti[D]*stu[D]*sro[D]*yD1 - cti[D]*ctu[D]*zD1);

	    //d/dturn
	    derivD[0][4] = -1.*(-stu[D]*cro[D]*xD1 + stu[D]*sro[D]*yD1 + ctu[D]*zD1);
	    derivD[1][4] = -1.*(sti[D]*ctu[D]*cro[D]*xD1 - sti[D]*ctu[D]*sro[D]*yD1 + sti[D]*stu[D]*zD1);



	    if(!isFiducial(xB,yB) || !isFiducial(xC,yC)){
	      continue;
	    }
	    
	    if( cA->big > 0 ) continue;
	    if( cB->big > 0 ) continue;
	    if( cC->big > 0 ) continue;
	    if( cD->big > 0 ) continue;

	    double za=cA->zg;
	    double zb=cB->zg;
	    double zc=cC->zg;
	    double zd=cD->zg;

	    bool curved = true;
	    if(bfield < 0.001) curved = false;

	    double mx,bx,my,by,sigmax,sigmay;

	    linefit4(za,xA,zb,xB,zc,xC,zd,xD,mx,bx,sigmax);
	    linefit4(za,yA,zb,yB,zc,yC,zd,yD,my,by,sigmay);


	    n4ev++;

	    vector<GblPoint> listOfPoints;
	    listOfPoints.reserve(4);
	    vector<double> sPoint;

	    ntry++;
	    nCandidatesPerEvent++;

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

	    GblTrajectory traj( listOfPoints, curved ); // 0 = no magnetic field
	    //traj.printPoints();

	    double Chi2;
	    int Ndf;
	    double lostWeight;

	    traj.fit( Chi2, Ndf, lostWeight );

	    double probchi = TMath::Prob( Chi2, Ndf );

	    hchi2.Fill( Chi2 );
	    hndf.Fill( Ndf );
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

	    fitrxAllHisto[A].Fill((xA-(mx*za+bx))*1.E3);
	    fitryAllHisto[A].Fill((yA-(my*za+by))*1.E3);
	    fitrxAllHisto[B].Fill((xB-(mx*zb+bx))*1.E3);
	    fitryAllHisto[B].Fill((yB-(my*zb+by))*1.E3);
	    fitrxAllHisto[C].Fill((xC-(mx*zc+bx))*1.E3);
	    fitryAllHisto[C].Fill((yC-(my*zc+by))*1.E3);
	    fitrxAllHisto[D].Fill((xD-(mx*zd+bx))*1.E3);
	    fitryAllHisto[D].Fill((yD-(my*zd+by))*1.E3);

	    hmx4fitall.Fill(mx*1000.);
	    if(sigmax*sigmax<0.01){
	      hmx4fitchi.Fill(mx*1000.);
	      hmx4fitchiinv.Fill(-mx*1000.);
	    }

	    if(probchi > 0.01){

	      nGoodPerEvent++;

	      hmx4fit.Fill(mx*1000.);
	      hmx4fitl.Fill(mx*1000.);
	      mx4vsxfit.Fill(bx,mx*1000.);
	      hbx4fit.Fill(bx);
	      hmy4fit.Fill(my*1000.);
	      hby4fit.Fill(by);
	      hsigmax4fit.Fill(sigmax*1000.);
	      hsigmay4fit.Fill(sigmay*1000.);
	      hchi2fitPassed.Fill(sigmax*sigmax*1e6);

	      fitrxHisto[A].Fill((xA-(mx*za+bx))*1.E3);
	      fitryHisto[A].Fill((yA-(my*za+by))*1.E3);
	      fitrxHisto[B].Fill((xB-(mx*zb+bx))*1.E3);
	      fitryHisto[B].Fill((yB-(my*zb+by))*1.E3);
	      fitrxHisto[C].Fill((xC-(mx*zc+bx))*1.E3);
	      fitryHisto[C].Fill((yC-(my*zc+by))*1.E3);
	      fitrxHisto[D].Fill((xD-(mx*zd+bx))*1.E3);
	      fitryHisto[D].Fill((yD-(my*zd+by))*1.E3);

	      if(lq[A]){
		fitrxqHisto[A].Fill((xA-(mx*za+bx))*1.E3);
		fitryqHisto[A].Fill((yA-(my*za+by))*1.E3);
	      }
	      if(lq[B]){
		fitrxqHisto[B].Fill((xB-(mx*zb+bx))*1.E3);
		fitryqHisto[B].Fill((yB-(my*zb+by))*1.E3);
	      }
	      if(lq[C]){	      
		fitrxqHisto[C].Fill((xC-(mx*zc+bx))*1.E3);
		fitryqHisto[C].Fill((yC-(my*zc+by))*1.E3);
	      }
	      if(lq[D]){
		fitrxqHisto[D].Fill((xD-(mx*zd+bx))*1.E3);
		fitryqHisto[D].Fill((yD-(my*zd+by))*1.E3);
	      }


	      for(size_t iplane = 1; iplane < 5; iplane++){
		
		traj.getResults( iplane, aCorrection, aCovariance );
		traj.getMeasResults( static_cast<unsigned int>(iplane), ndim, aResiduals, aMeasErrors, aResErrors, aDownWeights );
		traj.getScatResults( static_cast<unsigned int>(iplane), ndim, kKinks, kKinkErrors, kResErrors, kDownWeights );
		
		gblaxSHisto[iplane-1].Fill( (-aCorrection[1])*1E3 ); // angle x [mrad]
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
		gblpqHisto[iplane-1].Fill( aCorrection[0]*1E3 );

		if(lqA && lqB && lqC && lqD){
		  gblrxqAllHisto[iplane-1].Fill( aResiduals[0] * 1E3 );
		  gblryqAllHisto[iplane-1].Fill( aResiduals[1] * 1E3 );
		}
		if(lq[iplane-1]){
		  gblrxqHisto[iplane-1].Fill( aResiduals[0] * 1E3 );
		  gblryqHisto[iplane-1].Fill( aResiduals[1] * 1E3 );
		}

		double currentCharge = 0.;
		if(iplane-1==0){
		  currentCharge=cA->charge*norm;
		}else if(iplane-1==1){
		  currentCharge=cB->charge*norm;
		}else if(iplane-1==2){
		  currentCharge=cC->charge*norm;
		}else if(iplane-1==3){
		  currentCharge=cD->charge*norm;
		}

		gblrxvsq0[iplane-1].Fill(currentCharge, aResiduals[0] * 1E3 );
		gblryvsq0[iplane-1].Fill(currentCharge, aResiduals[1] * 1E3 );

	      
		hclq0t[iplane-1].Fill( currentCharge );
		if(event_nr<100000){
		  hclq0tstart[iplane-1].Fill( currentCharge );
		}
		if(event_nr>900000 && event_nr<1000000){
		  hclq0tend[iplane-1].Fill( currentCharge );
		}
		hclq0tvst[iplane-1].Fill(event_nr, currentCharge);

		double currentNrow = 0.;
		if(iplane-1==0){
		  currentNrow=cA->nrow;
		}else if(iplane-1==1){
		  currentNrow=cB->nrow;
		}else if(iplane-1==2){
		  currentNrow=cC->nrow;
		}else if(iplane-1==3){
		  currentNrow=cD->nrow;
		}
		hnrow4[iplane-1].Fill(currentNrow);

		double currentNcol = 0.;
		if(iplane-1==0){
		  currentNcol=cA->ncol;
		}else if(iplane-1==1){
		  currentNcol=cB->ncol;
		}else if(iplane-1==2){
		  currentNcol=cC->ncol;
		}else if(iplane-1==3){
		  currentNcol=cD->ncol;
		}

		if(lqAl && lqBl && lqCl && lqDl){
		  hncol4qf4[iplane-1].Fill(currentNcol);
		}
		if(lqAl && lqBl && lqCl){
		  hncol3qf4[iplane-1].Fill(currentNcol);
		  hmx4fitlused.Fill(mx*1000.);
		  hmx4fitlusedinv.Fill(-mx*1000.);
		  hchi2fitUsed.Fill(sigmax*sigmax*1e6);
		}


		if(iplane-1==0){
		  if(lqA){
		    hncolqf4[0].Fill(cA->ncol);
		    hnrowqf4[0].Fill(cA->nrow);
		  }
		}else if(iplane-1==1){
		  if(lqB){
		    hncolqf4[1].Fill(cB->ncol);
		    hnrowqf4[1].Fill(cB->nrow);
		  }
		}else if(iplane-1==2){
		  if(lqC){
		    hncolqf4[2].Fill(cC->ncol);
		    hnrowqf4[2].Fill(cC->nrow);
		  }
		}else if(iplane-1==3){
		  if(isDirrad){
		    if(cD->charge*norm > 5. && cD->charge*norm < chCutHigh){
		      hncolqf4[3].Fill(cD->ncol);
		      hnrowqf4[3].Fill(cD->nrow);
		    }
		  }else{
		    if(lqD){
		      hncolqf4[3].Fill(cD->ncol);
		      hnrowqf4[3].Fill(cD->nrow);
		    }
		  }
		}


	      }

	    }else{
	      hchi2fitFailed.Fill(sigmax*sigmax*1e6);
	    }

	    // write to MP binary file

	    if( probchi > 0.01 ) { // bias with bad alignment ?
	      traj.milleOut( *mille );
	      ++nmille;
	    }



	  }	  
	}
      }
    }

    hn4candidates.Fill(nCandidatesPerEvent);
    hn4good.Fill(nGoodPerEvent);


  } while( event_nr < lev ); // <-----  Event loop


  delete mille;

  // Prealignment

  if((nmille < 900 && alignmentRun == run && !doAngleRun) || !aligniteration){

    cout << "Not enough tracks for MP. Do Prealignment." << endl;

    double nb, ne, nm;

    nb = hdxAC.GetNbinsX();
    ne = hdxAC.GetSumOfWeights();
    nm = hdxAC.GetMaximum();
    cout << endl << hdxAC.GetTitle() << endl
	 << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
	 << "  Maximum = " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

    if( nm/(ne/nb) > 2 ) {

      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, nm ); // amplitude
      fgp0->SetParameter( 1, hdxAC.GetBinCenter( hdxAC.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 0.05 ); // sigma
      fgp0->SetParameter( 3, hdxAC.GetBinContent(1) ); // BG
      hdxAC.Fit( "fgp0", "q" );
      cout << "  Fit Gauss + BG:"
	   << endl << "  m = " << fgp0->GetParameter(1)
	   << endl << "  alignx[0] = " << alignx[0] + fgp0->GetParameter(1)
	   << endl;
      alignx[0] += fgp0->GetParameter(1);

    }

    if( aligniteration ) {

      dxvsyAC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdxvsy = dxvsyAC.GetFunction( "pol1" );
      cout << endl << dxvsyAC.GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;
      //rotA[0] += fdxvsy->GetParameter(1); // same sign

      dxvsxAC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdxvsx = dxvsxAC.GetFunction( "pol1" );
      cout << endl << dxvsxAC.GetTitle() << " slope " << fdxvsx->GetParameter(1) << endl;
      //turnA[0] -= fdxvsx->GetParameter(1); // same sign

    }

    nb = hdyAC.GetNbinsX();
    ne = hdyAC.GetSumOfWeights();
    nm = hdyAC.GetMaximum();

    cout << endl << hdyAC.GetTitle() << endl
	 << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
	 << "  Maximum = " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

    if( nm/(ne/nb) > 2 ) {

      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, nm ); // amplitude
      fgp0->SetParameter( 1, hdyAC.GetBinCenter( hdyAC.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 0.05 ); // sigma
      fgp0->SetParameter( 3, hdyAC.GetBinContent(1) ); // BG
      hdyAC.Fit( "fgp0", "q" );
      cout << "  Fit Gauss + BG:"
	   << endl << "  m = " << fgp0->GetParameter(1)
	   << endl << "  aligny[0] = " << aligny[0] + fgp0->GetParameter(1)
	   << endl;
      aligny[0] += fgp0->GetParameter(1);
    }

    if( aligniteration ) {

      dyvsxAC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdyvsx = dyvsxAC.GetFunction( "pol1" );
      cout << endl << dyvsxAC.GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;
      //rotA[0] -= fdyvsx->GetParameter(1); // opposite sign

      dyvsyAC.Fit( "pol1", "q", "", -7.5, 7.5 );
      TF1 * fdyvsy = dyvsyAC.GetFunction( "pol1" );
      cout << endl << dyvsyAC.GetTitle()
	   << " slope " << fdyvsy->GetParameter(1)
	   << endl;
      //tiltA[0] += fdyvsy->GetParameter(1); // same sign

    }

    // B:

    nb = hdxBC.GetNbinsX();
    ne = hdxBC.GetSumOfWeights();
    nm = hdxBC.GetMaximum();

    cout << endl << hdxBC.GetTitle() << endl
	 << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
	 << "  Maximum = " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

    if( nm/(ne/nb) > 2 ) {

      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, nm ); // amplitude
      fgp0->SetParameter( 1, hdxBC.GetBinCenter( hdxBC.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 0.05 ); // sigma
      fgp0->SetParameter( 3, hdxBC.GetBinContent(1) ); // BG
      hdxBC.Fit( "fgp0", "q" );
      cout << "  Fit Gauss + BG:"
	   << endl << "  m = " << fgp0->GetParameter(1)
	   << endl << "  alignx[2] = " << alignx[2] + fgp0->GetParameter(1)
	   << endl;
      alignx[1] += fgp0->GetParameter(1);

    }

    if( aligniteration ) {

      dxvsyBC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdxvsy = dxvsyBC.GetFunction( "pol1" );
      cout << endl << dxvsyBC.GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;
      //rotA[1] += fdxvsy->GetParameter(1); // same sign

      dxvsxBC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdxvsx = dxvsxBC.GetFunction( "pol1" );
      cout << endl << dxvsxBC.GetTitle() << " slope " << fdxvsx->GetParameter(1) << endl;
      //turnA[1] -= fdxvsx->GetParameter(1); // same sign

    }

    nb = hdyBC.GetNbinsX();
    ne = hdyBC.GetSumOfWeights();
    nm = hdyBC.GetMaximum();

    cout << endl << hdyBC.GetTitle() << endl
	 << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
	 << "  Maximum = " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

    if( nm/(ne/nb) > 2 ) {

      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, nm ); // amplitude
      fgp0->SetParameter( 1, hdyBC.GetBinCenter( hdyBC.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 0.05 ); // sigma
      fgp0->SetParameter( 3, hdyBC.GetBinContent(1) ); // BG
      hdyBC.Fit( "fgp0", "q" );
      cout << "  Fit Gauss + BG:"
	   << endl << "  m = " << fgp0->GetParameter(1)
	   << endl << "  aligny[2] = " << aligny[2] + fgp0->GetParameter(1)
	   << endl;
      aligny[1] += fgp0->GetParameter(1);
    }

    if( aligniteration ) {

      dyvsxBC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdyvsx = dyvsxBC.GetFunction( "pol1" );
      cout << endl << dyvsxBC.GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;
      //rotA[1] -= fdyvsx->GetParameter(1); // opposite sign

      dyvsyBC.Fit( "pol1", "q", "", -7.5, 7.5 );
      TF1 * fdyvsy = dyvsyBC.GetFunction( "pol1" );
      cout << endl << dyvsyBC.GetTitle()
	   << " slope " << fdyvsy->GetParameter(1)
	   << endl;
      //tiltA[1] += fdyvsy->GetParameter(1); // same sign
    }

    // D:
    
    nb = hdxDC.GetNbinsX();
    ne = hdxDC.GetSumOfWeights();
    nm = hdxDC.GetMaximum();

    cout << endl << hdxDC.GetTitle() << endl
	 << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
	 << "  Maximum = " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

    if( nm/(ne/nb) > 2 ) {
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, nm ); // amplitude
      fgp0->SetParameter( 1, hdxDC.GetBinCenter( hdxDC.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 0.05 ); // sigma
      fgp0->SetParameter( 3, hdxDC.GetBinContent(1) ); // BG
      hdxDC.Fit( "fgp0", "q" );
      cout << "  Fit Gauss + BG:"
	   << endl << "  m = " << fgp0->GetParameter(1)
	   << endl << "  alignx[3] = " << alignx[3] + fgp0->GetParameter(1)
	   << endl;
      alignx[3] += fgp0->GetParameter(1);
    
    }

    if( aligniteration ) {

      dxvsyDC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdxvsy = dxvsyDC.GetFunction( "pol1" );
      cout << endl << dxvsyDC.GetTitle() << " slope " << fdxvsy->GetParameter(1) << endl;
      //rotA[3] += fdxvsy->GetParameter(1); // same sign

      dxvsxDC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdxvsx = dxvsxDC.GetFunction( "pol1" );
      cout << endl << dxvsxDC.GetTitle() << " slope " << fdxvsx->GetParameter(1) << endl;
      //turnA[3] -= fdxvsx->GetParameter(1); // same sign

    }

    nb = hdyDC.GetNbinsX();
    ne = hdyDC.GetSumOfWeights();
    nm = hdyDC.GetMaximum();

    cout << endl << hdyDC.GetTitle() << endl
	 << "  Inside  " << ne << " (" << ne/nb << " per bin)" << endl
	 << "  Maximum = " << nm << " (factor " << nm/ne*nb << " above mean)" << endl;

    if( nm/(ne/nb) > 2 ) {
      TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)+[3]", -1, 1 );
      fgp0->SetParameter( 0, nm ); // amplitude
      fgp0->SetParameter( 1, hdyDC.GetBinCenter( hdyDC.GetMaximumBin() ) );
      fgp0->SetParameter( 2, 0.05 ); // sigma
      fgp0->SetParameter( 3, hdyDC.GetBinContent(1) ); // BG
      hdyDC.Fit( "fgp0", "q" );
      cout << "  Fit Gauss + BG:"
	   << endl << "  m = " << fgp0->GetParameter(1)
	   << endl << "  aligny[3] = " << aligny[3] + fgp0->GetParameter(1)
	   << endl;
      aligny[3] += fgp0->GetParameter(1);

    }

    if( aligniteration ) {

      dyvsxDC.Fit( "pol1", "q", "", 0, 15 );
      TF1 * fdyvsx = dyvsxDC.GetFunction( "pol1" );
      cout << endl << dyvsxDC.GetTitle() << " slope " << fdyvsx->GetParameter(1) << endl;
      //rotA[3] -= fdyvsx->GetParameter(1); // opposite sign

      dyvsyDC.Fit( "pol1", "q", "", -7.5, 7.5 );
      TF1 * fdyvsy = dyvsyDC.GetFunction( "pol1" );
      cout << endl << dyvsyDC.GetTitle()
	   << " slope " << fdyvsy->GetParameter(1)
	   << endl;
      //tiltA[3] += fdyvsy->GetParameter(1); // same sign

    }
  }else if(alignmentRun == run && !doAngleRun){ 
    // if prealignment was not performed, execute millepede.

    // MILLEPEDE:

    ofstream steerFile;
    std::stringstream steername;
    steername << "run" << run << "/steerPede.txt";
    steerFile.open( steername.str().c_str() );

    if( steerFile.is_open() ) {



      steerFile << "! generated by quadMP" << std::endl;
      steerFile << "Cfiles" << std::endl;
      steerFile << "mille2.bin" << std::endl;
      steerFile << std::endl;

      steerFile << "Parameter" << std::endl;

      bool fixed[4] = {true, false, true, false};

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


    if(chdir(rundir.c_str()) != 0){
      cout << "Could not enter run directory." << endl;
    }else{
      remove("millepede.res");

      stringstream pedecmd;
      pedecmd << "./../millepede/pede steerPede.txt";
      
      std::cout << "Starting mille" << std::endl;

      system(pedecmd.str().c_str());

      std::cout << "Mille ran." << std::endl;
      
      if(chdir("..") != 0){
	cout << "Could not go back to main directory." << endl;
      }
    }

    //Read results
    
    ifstream milleRes( milleResName );

    cout << endl;
    if( milleRes.bad() || ! milleRes.is_open() ) {
      cout << "No " << milleResName << " found. Skipping this." << endl;
      cout << endl;
    }
    else {

      cout << "Read Millepede results from " << milleResName << endl;

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

      bool showAlPars[4] = {false, false, true, true};

      for(size_t mod = 0; mod < 4; mod++){
	cout << setprecision(6) << "MP alignment corrections plane " << mod << ":" << endl;
	cout <<  "\tdx    =  " << alpar[mod*ndimMP+1]*1E3 << " um" << endl;
	cout <<  "\tdy    =  " << alpar[mod*ndimMP+2]*1E3 << " um" << endl;
	cout <<  "\tdrot  =  " << alpar[mod*ndimMP+3]*1E3 << " mrad" << endl;
	cout <<  "\tdtilt =  " << alpar[mod*ndimMP+4]*1E3 << " mrad" << endl;
	cout <<  "\tdturn =  " << alpar[mod*ndimMP+5]*1E3 << " mrad" << endl;
      }

      size_t mod = 0;
      alignx[mod] += alpar[mod*ndimMP+1];
      aligny[mod] += alpar[mod*ndimMP+2];
      rotA[mod]   += alpar[mod*ndimMP+3];
      tiltA[mod]  += alpar[mod*ndimMP+4];
      turnA[mod]  += alpar[mod*ndimMP+5];

      mod = 1;
      alignx[mod] += alpar[mod*ndimMP+1];
      aligny[mod] += alpar[mod*ndimMP+2];
      rotA[mod]   += alpar[mod*ndimMP+3];
      tiltA[mod]  += alpar[mod*ndimMP+4];
      turnA[mod]  += alpar[mod*ndimMP+5];

      mod = 2;
      alignx[mod] += alpar[mod*ndimMP+1];
      aligny[mod] += alpar[mod*ndimMP+2];
      rotA[mod]   += alpar[mod*ndimMP+3];
      tiltA[mod]  += alpar[mod*ndimMP+4];
      turnA[mod]  += alpar[mod*ndimMP+5];

      mod = 3;
      alignx[mod] += alpar[mod*ndimMP+1];
      aligny[mod] += alpar[mod*ndimMP+2];
      rotA[mod]   += alpar[mod*ndimMP+3];
      tiltA[mod]  += alpar[mod*ndimMP+4];
      turnA[mod]  += alpar[mod*ndimMP+5];

    }

  }



  if(alignmentRun == run && !doAngleRun){

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

  }


  /*
  TF1 * fgp0 = new TF1( "fgp0", "[0]*exp(-0.5*((x-[1])/[2])^2)", -1, 1 );
  fgp0->SetParameter( 0, ntry ); // amplitude
  fgp0->SetParameter( 1, hmx4fitl.GetBinCenter( hmx4fitl.GetMaximumBin() ) );
  fgp0->SetParameter( 2, 2. ); // sigma
  hmx4fitl.Fit( "fgp0", "q" );
  */
  cout << "Mean track slope: " << hmx4fitl.GetMean() << " mrad" << endl << endl;


  std::cout << std::endl <<
    "Track fits: " << n4ev << std::endl <<
    "GBL fits: " << ntry << std::endl <<
    "GBL good: " << nmille << std::endl <<
    "Ratio GBL: " << double(nmille)/double(ntry) << std::endl << std::endl;

  

  
  if(alignmentRun != run) fitSupressed = true;
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

    
    cout << setprecision(4) << "\nLandau peaks [ke]:";
    for(int mod = 0; mod < 4; mod++){
      cout << endl << modName[mod] << ":\t";
      if(haveGain[mod]){
	for(int roc = 0; roc < 16; roc++){
	  if(modName[mod]==4582 && roc==8){
	    continue;
	  }
	  landau_peak[mod][roc] = landau_gauss_peak(&hclq0r[mod][roc]);
	  correction[mod][roc] = 22./landau_peak[mod][roc];
	  if(!CCSupressed || conversionRun == run){
	    if(modName[mod]!=4582){
	      if(correction[mod][roc] > 0.0001)ke[mod][roc] *= correction[mod][roc];
	      conversionFileOut << mod << "\t" << roc << "\t" << ke[mod][roc] << endl;
	    }
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



  histoFile->Write();
  histoFile->Close();
  
  

  return 1;
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
