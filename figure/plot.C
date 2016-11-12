#include "TFile.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "cstdlib"
#include "iostream"
#include "TStyle.h"
#include "TLegend.h"
#include "THStack.h"
#include "TLatex.h"
#include "TStyle.h"
#include "TAxis.h"
#include "TStyle.h"
#include "TROOT.h"
#include "setNCUStyle.C"
#include "algorithm"


void plot()
{
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetFrameLineWidth(3);
    setNCUStyle();
	
    TCanvas *c1 = new TCanvas("c1","",900,700);
	TLegend* leg = new TLegend(0.7,0.7,0.93,0.88);
    
    TFile *f_mzp1 = TFile::Open("ZpBaryonic_ZpBaryonic_MZp50_MChi50_hbb_boosted.root");
    TFile *f_mzp2 = TFile::Open("ZpBaryonic_ZpBaryonic_MZp50_MChi10_hbb_boosted.root"); //old
    TFile *f_mzp3 = TFile::Open("ZpBaryonic_ZpBaryonic_MZp50_MChi1_hbb_boosted.root");
    //TFile *f_mzp4 = TFile::Open("ZpBaryonic_ZpBaryonic_MZp10_MChi10_hbb_boosted.root"); //dm=100
    //TFile *f_mzp5 = TFile::Open("ZpBaryonic_ZpBaryonic_MZp10_MChi1_hbb_boosted.root");
    //TFile *f_mzp6 = TFile::Open("MZp600Ma0_filedm100_25.root");
    string outputname = "ZpBaryonic_MZp50";

    string hname[10] = {
        "h_pfMet" ,             // 0
        "h_higgsPt" ,           // 1
        "h_higgsEta" ,         // 2
        "h_deltaR_subjet" , // 3
        "h_extraEle",           // 4
        "h_extraMuo",         // 5
        "h_extraTau" ,          // 6
        "h_extraBJet",          // 7
        "h_extraAK4Jet",      // 8
        "h_higgsJetMass"};  // 9
    string hXaxis[10]={
        "E_{T}^{mass} (GeV)",
        "P_{T} (GeV)",
        "#eta",
        "#DeltaR",
        "",
        "",
        "",
        "",
        "",
        "M (GeV)"
    };
        int hnum;


    for (hnum = 0; hnum<10; hnum++) {
//        if (hnum >= 4 and hnum <= 8) {continue;}
    TH1F* h_A0m1 = (TH1F*)f_mzp1->Get(hname[hnum].data());
    TH1F* h_A0m2 = (TH1F*)f_mzp2->Get(hname[hnum].data());
    TH1F* h_A0m3 = (TH1F*)f_mzp3->Get(hname[hnum].data());    
    //TH1F* h_A0m4 = (TH1F*)f_mzp4->Get(hname[hnum].data());
    //TH1F* h_A0m5 = (TH1F*)f_mzp5->Get(hname[hnum].data());
    //TH1F* h_A0m6 = (TH1F*)f_mzp1->Get("monoHbbM1200_800");
//dR(bb):h_D_dR0
//PT(A0): h_Bpt1
//PT(H): h_Bpt0
//TM(x):h_XmT
//h_Xm; h_Xy;
//h_Xpz;h_Bpz0;h_Bpz1;h_Dpz0;h_Dpz1;
//h_cosThetaStar;h_cosPhi;h_Bm0;h_BmT0;h_BmT1;h_By0;h_By1

    float scale1 = 1.0/h_A0m1->Integral();
    float scale2 = 1.0/h_A0m2->Integral();
    float scale3 = 1.0/h_A0m3->Integral();
    //float scale4 = 1.0/h_A0m4->Integral();
    //float scale5 = 1.0/h_A0m5->Integral();
 //   float scale6 = 1.0/h_A0m6->Integral();

    h_A0m1->Scale(scale1);
    h_A0m2->Scale(scale2);
    h_A0m3->Scale(scale3);
    //h_A0m4->Scale(scale4);
    //h_A0m5->Scale(scale5);
    //h_A0m6->Scale(scale6);
    // h_zpM700_hA0_bbxx->Scale(scale7);

    h_A0m1->Sumw2();
    h_A0m2->Sumw2();
    h_A0m3->Sumw2();
    //h_A0m4->Sumw2();
    //h_A0m5->Sumw2();
    //h_A0m6->Sumw2();


    h_A0m1->SetLineWidth(2);
    h_A0m2->SetLineWidth(2);
    h_A0m3->SetLineWidth(2);
    //h_A0m4->SetLineWidth(2);
    //h_A0m5->SetLineWidth(2);
   // h_A0m6->SetLineWidth(2);
    // h->SetLineColor(1);
    double hmax[] = {h_A0m1->GetMaximum(),h_A0m2->GetMaximum(),h_A0m3->GetMaximum()};//,h_A0m4->GetMaximum(),h_A0m5->GetMaximum()};
    h_A0m1->SetMaximum(*max_element(hmax,hmax+3)*1.1);

    // h_A0m1->GetYaxis()->SetRange(0,700);
    // h_A0m2->GetYaxis()->SetRange(0,700);
    // h_A0m3->GetYaxis()->SetRange(0,700);
    h_A0m1->GetYaxis()->SetTitle("Normalized to 1");
    h_A0m1->GetYaxis()->SetTitleSize(0.04);
    h_A0m1->GetYaxis()->CenterTitle();
    //h_A0m1->GetXaxis()->SetTitle("P_{T} (GeV)");
    h_A0m1->GetXaxis()->SetTitle(hXaxis[hnum].data());
    h_A0m1->GetXaxis()->SetTitleSize(0.04);
    //h_A0m3->GetXaxis()->CenterTitle();
   
    // h_A0m5->GetYaxis()->SetRange(0,700);
    // h_A0m6->GetYaxis()->SetRange(0,700); 

    h_A0m1->Draw("hist");        
    h_A0m2->Draw("histsame");
    h_A0m3->Draw("histsame"); 
    //h_A0m4->Draw("histsame");
    //h_A0m5->Draw("histsame");
    //h_A0m6->Draw("histsame");

    //h_A0m5->SetLineColor(11);
    //h_A0m5->SetFillColor(98);
  //  h_A0m5->SetFillStyle(3022);
    h_A0m1->SetLineColor(79);
    //h_A0m2->SetFillStyle(3004);
    h_A0m2->SetLineColor(2);
    h_A0m3->SetLineColor(92);
    //h_A0m4->SetLineColor(1);
    //h_A0m5->SetLineColor(51);


    /*leg->AddEntry(h_A0m1,"M_{Z'} = 295 GeV, M_{#chi} = 150 GeV");
    leg->AddEntry(h_A0m2,"M_{Z'} = 300 GeV, M_{#chi} = 50 GeV");
    leg->AddEntry(h_A0m3,"M_{Z'} = 300 GeV, M_{#chi} = 1 GeV");*/
    leg->AddEntry(h_A0m1,"M_{#chi} = 150 GeV");
    leg->AddEntry(h_A0m2,"M_{#chi} = 50 GeV");
    leg->AddEntry(h_A0m3,"M_{#chi} = 1 GeV");
    //leg->AddEntry(h_A0m4,"M_{#chi} = 10 GeV");
    //leg->AddEntry(h_A0m5,"M_{#chi} = 1 GeV");
    //leg->AddEntry(h_A0m6,"M_{A_{0}} = 800 GeV");

    leg->Draw();


    c1->Update();
    // Latex
    TString latexCMSname= "CMS Simulation Preliminary #sqrt{s} = 13 TeV";
    //TString latexCMSname2= "Z' #rightarrow A_{0} + H";
    TString latexCMSname3= "MZp50";
    //TString latexCMSname4= "tan#beta = 1, M_{#chi} = 100 GeV";
    TLatex Tl; Tl.SetTextFont(72); Tl.SetTextSize(0.035); 
    Tl.SetNDC(kTRUE); 
    Tl.SetTextAlign(22);
    Tl.DrawLatex(0.5,0.96,latexCMSname);
    Tl.DrawLatex(0.85,0.6,"ZpBaryonic" );
    Tl.DrawLatex(0.85, 0.55,latexCMSname3);
    Tl.DrawLatex(0.85, 0.5, hname[hnum].substr(2).data());
    string pdfname = "plot/" + outputname + ".pdf";
    if (hnum == 0) {pdfname += "(";}
    if (hnum == 9) {pdfname += ")";}
    string pngname = "plot/" + outputname + "_" + hname[hnum].substr(2) + ".png";
    c1->Print(pdfname.data());
    c1->SaveAs(pngname.data());
    c1->Clear();
    c1->cd();
    leg->Clear();}
}
