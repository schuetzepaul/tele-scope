#include "TFile.h"
#include "TChain.h"
#include "THStack.h"
#include "TF1.h"
#include "TH1.h"
#include "TH3.h"
#include "TH2F.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TStyle.h"
#include "TPaveStats.h"
#include "TGaxis.h"
#include "TAxis.h"
#include "TList.h"
#include "TLatex.h"
#include "TLine.h"
#include "TObject.h"
#include "TDirectory.h"
#include "TGraphAsymmErrors.h"
#include "TEfficiency.h"
#include "TKey.h"
#include <iostream>
#include <algorithm>
#include <vector>
#include <exception>
#include <cmath> 
#include <iomanip>
#include <fstream>
#include <string>
#include <sys/stat.h>
#include <sstream>
#include <map>
#include "TMath.h"
#include <TPRegexp.h>


//*****************************************************************************

//*****************************************************************************

/// Macro produces a consolidated pdf with selected plots for the purpose of checking the results of the data

/// Example run by executing,
///      root -l -n 'macro/checkDataPlots.C("scope23168.root")'

///@ param: input file name
void  checkDataPlots(std::string file = "scope23168.root") {

  /// Suppress class constructor print outs (gErrorIgnoreLevel = kPrint, kInfo, kWarning, kError, kBreak, kSysError, kFatal)
  gROOT->ProcessLine("gErrorIgnoreLevel = 1001;");

  /// Extract file name
  TObjArray* token  = TString(file).Tokenize(".");
  std::string  scopeNumber = ((TObjString*)token->At(0))->GetString().Data();

  /// Create directory to store pdf file
  gSystem->MakeDirectory(TString(scopeNumber));

  /// Create maps for TH1I, TProfile, and TProfile2D
  std::map<std::string, TH1I*> mapOfTH1I;
  std::map<std::string, TProfile*> mapOfTProfile;
  std::map<std::string, TProfile2D*> mapOfTProfile2D;
  std::map<std::string, TCanvas*> mapOfTCanvas;

  /// Create input file handler
  TFile* sampleFile = new TFile(file.c_str());
  std::cout << "reading file " << sampleFile->GetName() << std::endl;

  /// Check input file is not a zombine (corrupt or does not exit)
  if(sampleFile->IsZombie()) {
    std::cout<<"\nFile"+TString(sampleFile->GetName())+" not found!!\n";
  }  
  else {
    //// Create list from list of histogram from root file
    TList* list = sampleFile->GetListOfKeys();

    /// For debugging purposes (print list of keys)
    bool debug = false;
    if(debug)list->Print();

    /// Check to see if key is not found in list
    if(!list) {
      std::cout << TString::Format("<Error> No keys found in file\n");
      exit(1);
    }

    /// Create iterator from list and store to map based on histogram type
    TIter next((TList*)list);
    TKey* key;

    /// Iterate and extract histograms from the root file and store in separate types of maps
    while((key=(TKey*)next())){

      if (std::strstr(key->GetClassName(),"TH1I")) {
	mapOfTH1I[key->GetName()] = (TH1I*)sampleFile->Get(key->GetName());
      }
      else if (std::strstr(key->GetClassName(),"TProfile2D")) {
        mapOfTProfile2D[key->GetName()] = (TProfile2D*)sampleFile->Get(key->GetName());
      }
      else if (std::strstr(key->GetClassName(),"TProfile")) {
	mapOfTProfile[key->GetName()] = (TProfile*)sampleFile->Get(key->GetName());
      }
    }
  }

  /// Declare output file name  
  std::string outputfilename = "plot_"+scopeNumber+".pdf";

  /// Create canvanses based on number of row and column for plots
  mapOfTCanvas["canvas2x1"] = new TCanvas("canvas2x1","Data monitoring",10,10,700,700);
  mapOfTCanvas["canvas2x1"]->Divide(2,1);

  mapOfTCanvas["canvas2x2"] = new TCanvas("canvas2x2","Data monitoring",10,10,700,700);
  mapOfTCanvas["canvas2x2"]->Divide(2,2);

  mapOfTCanvas["canvas3x1"] = new TCanvas("canvas3x1","Data monitoring",10,10,900,600);
  mapOfTCanvas["canvas3x1"]->Divide(3,1);

  mapOfTCanvas["canvas3x2"] = new TCanvas("canvas3x2","Data monitoring",10,10,900,600);
  mapOfTCanvas["canvas3x2"]->Divide(3,2);
  
  /// Open file to print and save canvases on different pages
  mapOfTCanvas["canvas3x2"]->Print(TString(scopeNumber+"/"+outputfilename+"["));

  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTH1I["trix"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTH1I["triy"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas2x1"]->cd(1);
  mapOfTProfile["effvsx"]->Draw();
  mapOfTCanvas["canvas2x1"]->cd(2);
  mapOfTProfile["effvsy"]->Draw();
  mapOfTCanvas["canvas2x1"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas2x1"]->cd(1);
  mapOfTH1I["sixdxc"]->Draw();
  mapOfTCanvas["canvas2x1"]->cd(2);
  mapOfTH1I["sixdyc"]->Draw();
  mapOfTCanvas["canvas2x1"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas3x1"]->cd(1);
  mapOfTProfile["refdyvsx"]->Draw();
  mapOfTCanvas["canvas3x1"]->cd(2);
  mapOfTProfile["refdyvsy"]->Draw();
  mapOfTCanvas["canvas3x1"]->cd(3);
  mapOfTProfile["refdyvsty"]->Draw();
  mapOfTCanvas["canvas3x1"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas3x1"]->cd(1);
  mapOfTProfile["cmsdyvsx"]->Draw();
  mapOfTCanvas["canvas3x1"]->cd(2);
  mapOfTProfile["cmsdyvsy"]->Draw();
  mapOfTCanvas["canvas3x1"]->cd(3);
  mapOfTProfile["cmsdyvsty"]->Draw();
  mapOfTCanvas["canvas3x1"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas3x1"]->cd(1);
  mapOfTProfile["cmspxqvsq"]->Draw();
  mapOfTCanvas["canvas3x1"]->cd(2);
  mapOfTProfile["cmspxqvsxm"]->Draw();
  mapOfTCanvas["canvas3x1"]->cd(3);
  mapOfTProfile["cmspxqvsym"]->Draw();
  mapOfTCanvas["canvas3x1"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTH1I["dutpxq2nd02"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTH1I["dutpxq2nd03"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(3);
  mapOfTH1I["dutpxq2nd04"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(4);
  mapOfTH1I["dutpxq2nd05"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));
 
  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTProfile["cmsrmsxvsx"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTProfile["cmsrmsyvsx"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(3);
  mapOfTProfile["cmsrmsxvsy"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(4);
  mapOfTProfile["cmsrmsyvsy"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));
  
  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTProfile["cmsdxvsxm"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTProfile["cmsdxvsym"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(3);
  mapOfTProfile["cmsdyvsxm"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(4);
  mapOfTProfile["cmsdyvsym"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTProfile["cmsrmsxvsx"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTProfile["cmsrmsyvsx"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(3);
  mapOfTProfile["cmsrmsxvsy"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(4);
  mapOfTProfile["cmsrmsyvsy"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTProfile["cmsrmsxvsxm"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTProfile["cmsrmsyvsxm"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(3);
  mapOfTProfile["cmsrmsxvsym"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(4);
  mapOfTProfile["cmsrmsyvsym"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTProfile["cmsqxvst1"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTProfile["cmsqxvst2"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(3);
  mapOfTProfile["cmsqxvst3"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas3x2"]->cd(1);
  mapOfTH1I["cmsdyfcq0"]->Draw();
  mapOfTCanvas["canvas3x2"]->cd(2);
  mapOfTH1I["cmsdyfcq1"]->Draw();
  mapOfTCanvas["canvas3x2"]->cd(3);
  mapOfTH1I["cmsdyfcq2"]->Draw();
  mapOfTCanvas["canvas3x2"]->cd(4);
  mapOfTH1I["cmsdyfcq3"]->Draw();
  mapOfTCanvas["canvas3x2"]->cd(5);
  mapOfTH1I["cmsdyfcq9"]->Draw();
  mapOfTCanvas["canvas3x2"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas3x2"]->cd(1);
  mapOfTH1I["cmsq0f"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(2);
  mapOfTH1I["cmsq0fdot"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(3);
  mapOfTH1I["cmsq0fnod"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(4);
  mapOfTH1I["cmsq0fc"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(5);
  mapOfTH1I["cmsq0fp"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->Print(TString(scopeNumber+"/"+outputfilename));

  mapOfTCanvas["canvas3x2"]->cd(1);
  mapOfTProfile["cmsqxvsx"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(2);
  mapOfTProfile["cmsqxvsy"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(3);
  mapOfTProfile["cmsqxvsxm"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(4);
  mapOfTProfile["cmsqxvsym"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(5);
  mapOfTProfile["cmsqxvsymn"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->Print(TString(scopeNumber+"/"+outputfilename));  

  mapOfTCanvas["canvas2x2"]->cd(1);
  mapOfTH1I["cmsdx"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(2);
  mapOfTH1I["cmssxa"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(3);
  mapOfTH1I["cmsdxfc"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(4);
  mapOfTH1I["cmsdy"]->Draw();
  mapOfTCanvas["canvas2x2"]->cd(5);
  mapOfTH1I["cmssya"]->Draw();
  mapOfTCanvas["canvas2x2"]->Print(TString(scopeNumber+"/"+outputfilename));
  mapOfTCanvas["canvas2x2"]->cd(6);
  mapOfTH1I["cmsdyfc"]->Draw();

  mapOfTCanvas["canvas3x2"]->cd(1);
  mapOfTProfile2D["cmsqxvsxmym"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(2);
  mapOfTProfile2D["cmsnpxvsxmym"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(3);
  mapOfTProfile2D["cmsrmsxvsxmym"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(4);
  mapOfTProfile2D["cmsrmsyvsxmym"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(5);
  mapOfTProfile2D["refnpxvsxmym"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->cd(6);
  mapOfTProfile2D["refqxvsxmym"]->Draw("colz");
  mapOfTCanvas["canvas3x2"]->Print(TString(scopeNumber+"/"+outputfilename));


  /// Close file used to print and save canvases on different pages   
  mapOfTCanvas["canvas3x2"]->Print(TString(scopeNumber+"/"+outputfilename+"]"));

  /// Print the location and name of output file
  std::cout << "Output file name " << scopeNumber << "/" << outputfilename << std::endl;

  /// Exit ROOT prompt
  gROOT->ProcessLine(".q");

  return;
}
