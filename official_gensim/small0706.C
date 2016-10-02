#include <TLegend.h>
#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TRandom.h>
#include <TLorentzVector.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH1F.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TROOT.h>
#include "TImage.h"
#include "TSystem.h"
#include "TStyle.h"
#include "untuplizer.h"
#include <TClonesArray.h>
#include <fstream>
#include <cmath>
#include <TSystem.h>
#include <string>
#include <sstream>
#include "setNCUStyle.C"
#include<TH2.h>
#include "TLine.h"
#include "TF1.h"
#include"TGraphAsymmErrors.h"
#include "TLatex.h"

TH2D* readTxt(string inputDir[2],string outputName,int option=0){
	TCanvas* c1,*c2;
	setNCUStyle();
	c1 = new TCanvas("c1","",889*1.5,768);
	int massZ[15]={10,15,20,50,95,100,200,295,300,500,995,1000,1995,2000,10000};
	//int inputZ[8]={2,4,6,8,10,13,16,21};
	int massA[6]={1,10,50,150,500,1000};
	
	TH2D* th2[5];
	th2[0]=new TH2D("eff","eff",15,0,15,6,0,6);
	for(int i=0;i<15;i++){
		for(int j=0;j<6;j++){
			fstream file1(Form("%s_MZp-%d_MChi-%d_13TeV_%s.txt",inputDir[0].data(),massZ[i],massA[j],inputDir[1].data() ));
			cout<<massA[j]*2-massZ[i]<<endl;
			double db1=0;
			file1>>db1;
			if(massA[j]*2==massZ[i]){
				fstream file2(Form("%s_MZp-%d_MChi-%d_13TeV_%s.txt",inputDir[0].data(),massZ[i]-5,massA[j]),inputDir[1].data());
				cout<<"Y"<<endl;
				file2>>db1;
			}
			
			//cout<<Form("%s_MZp%d_MChi%d_hbb_%s.txt",inputDir[0].data(),massZ[i],massA[j],inputDir[1].data())<<endl;
			th2[0]->Fill(i,j,db1);
		}
	}
	
	for(int i=0;i<15;i++)th2[0]->GetXaxis()->SetBinLabel(i+1,Form("%d",massZ[i]));
	for(int i=0;i<6;i++)th2[0]->GetYaxis()->SetBinLabel(i+1,Form("%d",massA[i]));
	th2[0]->SetXTitle("m_{Z'}[GeV]");
	th2[0]->SetYTitle("m_{DM}[GeV]");
	th2[0]->SetMarkerSize(2);
	th2[0]->SetTitle(Form("%s",outputName.data()));
	gStyle->SetPaintTextFormat(" .2g ");
	th2[0]->Draw("colztext");
	c1->Print(Form("plot/%s.pdf",outputName.data()));
	c1->SaveAs(Form("plot/%s.png",outputName.data()));
	return th2[0];
}

void small0706(){
	
	
	string st[2]={
		"crab_MonoHbb_ZpBaryonic",
		"resolved"
	};
		
	readTxt(st,"MonoHbb_ZpBaryonic_resolve");
	st[1] = "boosted";
       	readTxt(st,"MonoHbb_ZpBaryonic_booested");

}




