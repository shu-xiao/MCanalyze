#include "iostream"
#include "fstream"
#include "TH1F.h"
#include "TH2.h"
#include "TGraph.h"
#include "TLegend.h"
using namespace std;

void readtxt() {
c1 = new TCanvas("c1","",1200,800);
Float_t mzp[] = {600,800,1000,1200,1400,1700,2000,2500};
int ma0[] = {300,400,500,600,700,800};
int line[6] = {8,8,7,6,6,6};
Int_t linecolor[6] = {51,61,71,81,91,100};
TGraph *effplot[6];
int skip;
int a;
// fig 0, ma0=300
cout << "ma0=300" << endl;
j=0;
Float_t eff0[8];
skip=0;
for (int i=0;i<8;i++) {
	//cout << "mzp"<< mzp[i] << endl;
	fstream file(Form("efftxtfile/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-%i_MA0-%d_13TeV-madgraph-SkimTree.root_HistogramsSignalRegion.txt",(int)mzp[i],ma0[j]));
	if (!file.is_open()) {
		skip += 1;
		continue;}
	a=i-skip;
	file >> eff0[a];
	cout << eff0[a] << endl;
}
//cout << eff0[7] <<endl;
effplot[j] = new TGraph(8,mzp,eff0);
for (int i=0;i<8;i++) {cout << mzp[i] << "  " <<eff0[i]<<endl;}
effplot[j]->SetLineWidth(2);
effplot[j]->SetLineColor(linecolor[j]);
effplot[j]->SetMarkerStyle(8);
effplot[j]->SetMarkerColor(linecolor[j]);
//effplot[j]->Draw("LP");
  //end of loop


// fig 1, ma0=400
cout << "ma0=400" << endl;
j=1;
Float_t eff1[8];
skip=0;
for (int i=0;i<8;i++) {
	//cout << "mzp"<< mzp[i] << endl;
	fstream file(Form("efftxtfile/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-%i_MA0-%d_13TeV-madgraph-SkimTree.root_HistogramsSignalRegion.txt",(int)mzp[i],ma0[j]));
	if (!file.is_open()) {
		skip+=1;
		continue;}
	a=i-skip;
	file >> eff1[a];
	cout << eff1[a] << endl;
}
effplot[j] = new TGraph(line[j],mzp,eff1);
for (int i=0;i<8;i++) {cout << mzp[i] << "  " <<eff1[i]<<endl;}
effplot[j]->SetLineWidth(2);
effplot[j]->SetLineColor(linecolor[j]);
effplot[j]->SetMarkerStyle(8);
effplot[j]->SetMarkerColor(linecolor[j]);
  //end of loop


// fig 2, ma0=500
cout << "ma0=500" << endl;
j=2;
Float_t eff2[7];
skip=0;
for (int i=0;i<8;i++) {
	//cout << "mzp"<< mzp[i] << endl;
	fstream file(Form("efftxtfile/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-%i_MA0-%d_13TeV-madgraph-SkimTree.root_HistogramsSignalRegion.txt",(int)mzp[i],ma0[j]));
	if (!file.is_open()) {
		skip+=1;
		continue;}
	a=i-skip;
	file >> eff2[a];
	cout << eff2[a] << endl;
}
Float_t mzp1[7] = {800,1000,1200,1400,1700,2000,2500};
effplot[j] = new TGraph(line[j],mzp1,eff2);
for (int i=0;i<8;i++) {cout << mzp1[i] << "  " <<eff2[i]<<endl;}
effplot[j]->SetLineWidth(2);
effplot[j]->SetLineColor(linecolor[j]);
effplot[j]->SetMarkerStyle(8);
effplot[j]->SetMarkerColor(linecolor[j]);
  //end of loop

// fig 3, ma0=600
cout << "ma0=600" << endl;
j=3;
Float_t eff3[6];
skip=0;
for (int i=0;i<8;i++) {
	//cout << "mzp"<< mzp[i] << endl;
	fstream file(Form("efftxtfile/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-%i_MA0-%d_13TeV-madgraph-SkimTree.root_HistogramsSignalRegion.txt",(int)mzp[i],ma0[j]));
	if (!file.is_open()) {
		skip+=1;
		continue;}
	a=i-skip;
	file >> eff3[a];
	cout << eff3[a] << endl;
}
Float_t mzp2[6] = {800,1000,1200,1400,1700,2500};
effplot[j] = new TGraph(line[j],mzp2,eff3);
for (int i=0;i<8;i++) {cout << mzp2[i] << "  " <<eff3[i]<<endl;}
effplot[j]->SetLineWidth(2);
effplot[j]->SetLineColor(linecolor[j]);
effplot[j]->SetMarkerStyle(8);
effplot[j]->SetMarkerColor(linecolor[j]);
  //end of loop

// fig 4, ma0=700
cout << "ma0=700" << endl;
j=4;
Float_t eff4[6];
skip=0;
for (int i=0;i<8;i++) {
	//cout << "mzp"<< mzp[i] << endl;
	fstream file(Form("efftxtfile/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-%i_MA0-%d_13TeV-madgraph-SkimTree.root_HistogramsSignalRegion.txt",(int)mzp[i],ma0[j]));
	if (!file.is_open()) {
		skip+=1;
		continue;}
	a=i-skip;
	file >> eff4[a];
	cout << eff4[a] << endl;
}
Float_t mzp3[6] = {1000,1200,1400,1700,2000,2500};
effplot[j] = new TGraph(line[j],mzp3,eff4);
for (int i=0;i<8;i++) {cout << mzp3[i] << "  " <<eff4[i]<<endl;}
effplot[j]->SetLineWidth(2);
effplot[j]->SetLineColor(linecolor[j]);
effplot[j]->SetMarkerStyle(8);
effplot[j]->SetMarkerColor(linecolor[j]);
  //end of loop

// fig 5, ma0=800
cout << "ma0=800" << endl;
j=5;
Float_t eff5[6];
skip=0;
for (int i=0;i<8;i++) {
	//cout << "mzp"<< mzp[i] << endl;
	fstream file(Form("efftxtfile/Merged_ZprimeToA0hToA0chichihbb_2HDM_MZp-%i_MA0-%d_13TeV-madgraph-SkimTree.root_HistogramsSignalRegion.txt",(int)mzp[i],ma0[j]));
	if (!file.is_open()) {
		skip+=1;
		continue;}
	a=i-skip;
	file >> eff5[a];
	cout << eff5[a] << endl;
}
effplot[j] = new TGraph(line[j],mzp3,eff5);
for (int i=0;i<8;i++) {cout << mzp3[i] << "  " <<eff5[i]<<endl;}
effplot[j]->SetLineWidth(2);
effplot[j]->SetLineColor(linecolor[j]);
effplot[j]->SetMarkerStyle(8);
effplot[j]->SetMarkerColor(linecolor[j]);
  //end of loop

effplot[5]->SetTitle("");
effplot[5]->GetXaxis()->SetTitle("M_{Z'} [GeV]");
effplot[5]->GetYaxis()->SetTitle("efficiency");

for (int i=0;i<6;i++) {
	effplot[i]->GetXaxis()->SetLimits(500.,2600.);
	effplot[i]->SetMinimum(0.);
	effplot[i]->SetMaximum(0.5);
	effplot[i]->Draw();
}

TLegend *leg = new TLegend(0.7, 0.63, 0.85, 0.9);
leg->SetBorderSize(0);
leg->SetFillColor(0);
leg->SetFillStyle(0);
leg->SetTextSize(0.04);
for (int i=0;i<6;i++){
	leg->AddEntry(effplot[i],Form("m_{A0} = %d GeV",ma0[i]),"LP");
	leg->Draw("same");
}


effplot[0]->Draw("LP");
for (int i=1;i<6;i++) {
	effplot[i]->Draw("LPsame");
}
c1->Update();
c1->Print("efficiencyPlot.pdf");
c1->SaveAs("efficiencyPlot.png");
}
