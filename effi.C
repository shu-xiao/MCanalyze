#include "TFile.h"
#include "iostream"
#include "TH1F.h"
#include "fstream"
#include "TFile.h"


using namespace std;
void effi() {
    
    int ma0[]={300,400,500,600,700,800};
    int mzp[]={600,800,1000,1200,1400,1700,2000,2500};
    TFile *f_effroot;
    TH1F *h_met, *h_num;
    float NUM, DEN;
    ofstream myfile;
    for (int i = 0; i < 6; i++){
        for (int j = 0; j < 8 ; j++){
            f_effroot = TFile::Open(Form("HistogramsAllRegion_sr/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-%d_MA0-%d_13TeV-madgraph-SkimTree.root",mzp[j],ma0[i]));
	    if (!f_effroot || !f_effroot->IsOpen()) continue;
            h_met = (TH1F*)f_effroot->Get("h_met_0");
            h_num = (TH1F*)f_effroot->Get("h_total_weight");
            NUM = h_met->Integral();
            DEN = h_num->Integral();
            cout << NUM/DEN << "i=" << i <<", j="<< j << endl;
            myfile.open(Form("SReffitxt/ZptoA0h_MZp-%d_MA0-%d.txt",mzp[j],ma0[i]));
            myfile << NUM/DEN << endl;
            myfile.close();

        }
    }
    //f_effroot = TFile::Open("HistogramsAllRegion_sr/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-1200_MA0-300_13TeV-madgraph-SkimTree.root");
    //h_met = (TH1F*)f_effroot->Get("h_met_0");
    //h_num = (TH1F*)f_effroot->Get("h_total_weight");
    //NUM = h_met->Integral();
    //DEN = h_num->Integral();
    //cout << NUM/DEN <<endl;
}
