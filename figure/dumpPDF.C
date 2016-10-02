#include <TROOT.h>
#include "TFile.h"
#include "TH1.h"
#include "TTree.h"
#include "TProfile.h"
#include "TKey.h"
#include <string>
#include <iostream>
#include <TSystem.h>
#include <TCanvas.h>

using namespace std;

void dumpPDF(std::string inputFile_)
{

  TString endfix=gSystem->GetFromPipe(Form("file=%s; echo \"${file%%.root}\"",inputFile_.data()));

  cout << endfix << endl;
  std::string dirName = "fig";
  gSystem->mkdir(dirName.data());

  TCanvas* c1 = new TCanvas("c1","",500,500);

  TFile *_file = TFile::Open(inputFile_.data());
  _file->cd();
  TDirectory *current_sourcedir = gDirectory;
  int nPage=0;
  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key;
  while ((key = (TKey*)nextkey()) ) {

    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( "TH1" ) ) {
      // descendant of TH1 -> scale it

      cout << "outputing histogram " << obj->GetName() << endl;
      TH1 *h1 = (TH1*)obj;
      h1->SetMarkerStyle(8);
      h1->SetMarkerSize(1);
      h1->SetLineWidth(3);
      h1->SetLineColor(4);
      h1->Draw();
      if(nPage==0)
	c1->Print(Form("%s/%s.pdf(",dirName.data(),endfix.Data()),"pdf");
      else 
	c1->Print(Form("%s/%s.pdf",dirName.data(),endfix.Data()),"pdf");
      nPage++;
    } // if the object is a histogram


  } // loop over keys

  c1->Print(Form("%s/%s.pdf)",dirName.data(),endfix.Data()),"pdf");
    


}
