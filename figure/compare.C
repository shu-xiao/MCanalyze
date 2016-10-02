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

using namespace std;
void comparebase(string MZp, string MChi, string macro) {
    gStyle->SetOptStat(0);
    gStyle->SetOptTitle(0);
    gStyle->SetFrameLineWidth(3);
    setNCUStyle();
	
    TCanvas *c1 = new TCanvas("c1","hmzp",900,700);
    TLegend* leg = new TLegend(0.7,0.75,0.9,0.8761);
    string model1 = "ZpBaryonic_ZpBaryonic_MZp" + MZp + "_MChi" + MChi + "_hbb_" + macro + ".root";
    string model2 = "ZpHS_ZpHS_MZp" + MZp + "_MChi" + MChi + "_hbb_" + macro + ".root";
    string pdfName = "MZp" + MZp + "_MChi" + MChi + "_" + macro + ".pdf";
    //string model1 = "ZpBaryonic_ZpBaryonic_MZp10000_MChi1000_hbb_resolved.root";
    //string model2 = "ZpHS_ZpHS_MZp10000_MChi1000_hbb_resolved.root"; 
    //string pdfName = "MZp10000_MChi1000_resolved.pdf";
    TFile *f_ZpBY = TFile::Open(model1.data());
    TFile *f_ZpHS = TFile::Open(model2.data());
 

 // pfMET
    TH1F *h_pfMet1 = (TH1F*)f_ZpBY->Get("h_pfMet");
    TH1F *h_pfMet2 = (TH1F*)f_ZpHS->Get("h_pfMet");
    h_pfMet1->SetLineColor(2);
    h_pfMet2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);

    float scale1 = 1.0/h_pfMet1->Integral();
    float scale2 = 1.0/h_pfMet2->Integral();

    h_pfMet1->Scale(scale1);
    h_pfMet2->Scale(scale2);

    h_pfMet1->Sumw2();
    h_pfMet2->Sumw2();

    h_pfMet2->GetYaxis()->SetTitle("Normalized to 1");
    h_pfMet2->GetYaxis()->SetTitleSize(0.04);
    h_pfMet2->GetYaxis()->CenterTitle();
    h_pfMet2->GetXaxis()->SetTitle("E (GeV)");
    h_pfMet2->GetXaxis()->SetTitleSize(0.04);

    h_pfMet2->Draw("hist");
    h_pfMet1->Draw("histsame");

    leg->AddEntry(h_pfMet1,"ZpBaryonic");
    leg->AddEntry(h_pfMet2,"ZpHS");
    leg->Draw();    

    TLatex Tl; Tl.SetTextFont(72); Tl.SetTextSize(0.04); 
    Tl.SetNDC(kTRUE); 
    Tl.SetTextAlign(22);
    Tl.DrawLatex(0.5,0.96,"MET");
    pdfName += "(";
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();
    pdfName.erase(pdfName.end()-1, pdfName.end()); // remove '('
    
    // higgs Pt
    TH1F *h_higgsPt1 = (TH1F*)f_ZpBY->Get("h_higgsPt");
    TH1F *h_higgsPt2 = (TH1F*)f_ZpHS->Get("h_higgsPt");
    h_higgsPt1->SetLineColor(2);
    h_higgsPt2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_higgsPt1->Integral();
    scale2 = 1.0/h_higgsPt2->Integral();

    h_higgsPt1->Scale(scale1);
    h_higgsPt2->Scale(scale2);

    h_higgsPt1->Sumw2();
    h_higgsPt2->Sumw2();

    h_higgsPt2->GetYaxis()->SetTitle("Normalized to 1");
    h_higgsPt2->GetYaxis()->SetTitleSize(0.04);
    h_higgsPt2->GetYaxis()->CenterTitle();
    h_higgsPt2->GetXaxis()->SetTitle("Pt (GeV)");
    h_higgsPt2->GetXaxis()->SetTitleSize(0.04);

    h_higgsPt2->Draw("hist");
    h_higgsPt1->Draw("histsame");


    leg->AddEntry(h_higgsPt1,"ZpBaryonic");
    leg->AddEntry(h_higgsPt2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"higgs Pt");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();

    // higgs Eta
    TH1F *h_higgsEta1 = (TH1F*)f_ZpBY->Get("h_higgsEta");
    TH1F *h_higgsEta2 = (TH1F*)f_ZpHS->Get("h_higgsEta");
    h_higgsEta1->SetLineColor(2);
    h_higgsEta2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_higgsEta1->Integral();
    scale2 = 1.0/h_higgsEta2->Integral();

    h_higgsEta1->Scale(scale1);
    h_higgsEta2->Scale(scale2);

    h_higgsEta1->Sumw2();
    h_higgsEta2->Sumw2();

    h_higgsEta2->GetYaxis()->SetTitle("Normalized to 1");
    h_higgsEta2->GetYaxis()->SetTitleSize(0.04);
    h_higgsEta2->GetYaxis()->CenterTitle();
    h_higgsEta2->GetXaxis()->SetTitle("#eta");
    h_higgsEta2->GetXaxis()->SetTitleSize(0.04);

    h_higgsEta2->Draw("hist");
    h_higgsEta1->Draw("histsame");


    leg->AddEntry(h_higgsEta1,"ZpBaryonic");
    leg->AddEntry(h_higgsEta2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"higgs Eta");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();


    // delta R of 2 subJet
    TH1F *h_deltaR1 = (TH1F*)f_ZpBY->Get("h_deltaR_subjet");
    TH1F *h_deltaR2 = (TH1F*)f_ZpHS->Get("h_deltaR_subjet");
    h_deltaR1->SetLineColor(2);
    h_deltaR2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_deltaR1->Integral();
    scale2 = 1.0/h_deltaR2->Integral();

    h_deltaR1->Scale(scale1);
    h_deltaR2->Scale(scale2);

    h_deltaR1->Sumw2();
    h_deltaR2->Sumw2();

    h_deltaR1->GetYaxis()->SetTitle("Normalized to 1");
    h_deltaR1->GetYaxis()->SetTitleSize(0.04);
    h_deltaR1->GetYaxis()->CenterTitle();
    h_deltaR1->GetXaxis()->SetTitle("deltaR");
    h_deltaR1->GetXaxis()->SetTitleSize(0.04);

    h_deltaR1->Draw("hist");
    h_deltaR2->Draw("histsame");


    leg->AddEntry(h_deltaR1,"ZpBaryonic");
    leg->AddEntry(h_deltaR2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"delta R of subjet");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();

        // extra electrons
    TH1F *h_extraEle1 = (TH1F*)f_ZpBY->Get("h_extraEle");
    TH1F *h_extraEle2 = (TH1F*)f_ZpHS->Get("h_extraEle");
    h_extraEle1->SetLineColor(2);
    h_extraEle2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_extraEle1->Integral();
    scale2 = 1.0/h_extraEle2->Integral();

    h_extraEle1->Scale(scale1);
    h_extraEle2->Scale(scale2);

    h_extraEle1->Sumw2();
    h_extraEle2->Sumw2();

    h_extraEle2->GetYaxis()->SetTitle("Normalized to 1");
    h_extraEle2->GetYaxis()->SetTitleSize(0.04);
    h_extraEle2->GetYaxis()->CenterTitle();
    h_extraEle2->GetXaxis()->SetTitle("extra Electrons");
    h_extraEle2->GetXaxis()->SetTitleSize(0.04);

    h_extraEle2->Draw("hist");
    h_extraEle1->Draw("histsame");


    leg->AddEntry(h_extraEle1,"ZpBaryonic");
    leg->AddEntry(h_extraEle2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"extra Electrons");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();

    // Extra Muons
    TH1F *h_extraMuo1 = (TH1F*)f_ZpBY->Get("h_extraMuo");
    TH1F *h_extraMuo2 = (TH1F*)f_ZpHS->Get("h_extraMuo");
    h_extraMuo1->SetLineColor(2);
    h_extraMuo2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_extraMuo1->Integral();
    scale2 = 1.0/h_extraMuo2->Integral();

    h_extraMuo1->Scale(scale1);
    h_extraMuo2->Scale(scale2);

    h_extraMuo1->Sumw2();
    h_extraMuo2->Sumw2();

    h_extraMuo2->GetYaxis()->SetTitle("Normalized to 1");
    h_extraMuo2->GetYaxis()->SetTitleSize(0.04);
    h_extraMuo2->GetYaxis()->CenterTitle();
    h_extraMuo2->GetXaxis()->SetTitle("extra muons");
    h_extraMuo2->GetXaxis()->SetTitleSize(0.04);

    h_extraMuo2->Draw("hist");
    h_extraMuo1->Draw("histsame");


    leg->AddEntry(h_extraMuo1,"ZpBaryonic");
    leg->AddEntry(h_extraMuo2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"extra Muons");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();

    // Extra Taus
    TH1F *h_extraTau1 = (TH1F*)f_ZpBY->Get("h_extraTau");
    TH1F *h_extraTau2 = (TH1F*)f_ZpHS->Get("h_extraTau");
    h_extraTau1->SetLineColor(2);
    h_extraTau2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_extraTau1->Integral();
    scale2 = 1.0/h_extraTau2->Integral();

    h_extraTau1->Scale(scale1);
    h_extraTau2->Scale(scale2);

    h_extraTau1->Sumw2();
    h_extraTau2->Sumw2();

    h_extraTau2->GetYaxis()->SetTitle("Normalized to 1");
    h_extraTau2->GetYaxis()->SetTitleSize(0.04);
    h_extraTau2->GetYaxis()->CenterTitle();
    h_extraTau2->GetXaxis()->SetTitle("extra Taus");
    h_extraTau2->GetXaxis()->SetTitleSize(0.04);

    h_extraTau2->Draw("hist");
    h_extraTau1->Draw("histsame");


    leg->AddEntry(h_extraTau1,"ZpBaryonic");
    leg->AddEntry(h_extraTau2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"extra Taus");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();

        // Extra b jets
    TH1F *h_extraBJet1 = (TH1F*)f_ZpBY->Get("h_extraBJet");
    TH1F *h_extraBJet2 = (TH1F*)f_ZpHS->Get("h_extraBJet");
    h_extraBJet1->SetLineColor(2);
    h_extraBJet2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_extraBJet1->Integral();
    scale2 = 1.0/h_extraBJet2->Integral();

    h_extraBJet1->Scale(scale1);
    h_extraBJet2->Scale(scale2);

    h_extraBJet1->Sumw2();
    h_extraBJet2->Sumw2();

    h_extraBJet2->GetYaxis()->SetTitle("Normalized to 1");
    h_extraBJet2->GetYaxis()->SetTitleSize(0.04);
    h_extraBJet2->GetYaxis()->CenterTitle();
    h_extraBJet2->GetXaxis()->SetTitle("extra b Jets");
    h_extraBJet2->GetXaxis()->SetTitleSize(0.04);

    h_extraBJet2->Draw("hist");
    h_extraBJet1->Draw("histsame");


    leg->AddEntry(h_extraBJet1,"ZpBaryonic");
    leg->AddEntry(h_extraBJet2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"extra b Jets");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();

      // Extra AK4 jets
    TH1F *h_extraAK4Jet1 = (TH1F*)f_ZpBY->Get("h_extraAK4Jet");
    TH1F *h_extraAK4Jet2 = (TH1F*)f_ZpHS->Get("h_extraAK4Jet");
    h_extraAK4Jet1->SetLineColor(2);
    h_extraAK4Jet2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_extraAK4Jet1->Integral();
    scale2 = 1.0/h_extraAK4Jet2->Integral();

    h_extraAK4Jet1->Scale(scale1);
    h_extraAK4Jet2->Scale(scale2);

    h_extraAK4Jet1->Sumw2();
    h_extraAK4Jet2->Sumw2();

    h_extraAK4Jet1->GetYaxis()->SetTitle("Normalized to 1");
    h_extraAK4Jet1->GetYaxis()->SetTitleSize(0.04);
    h_extraAK4Jet1->GetYaxis()->CenterTitle();
    h_extraAK4Jet1->GetXaxis()->SetTitle("extra AK4 Jets");
    h_extraAK4Jet1->GetXaxis()->SetTitleSize(0.04);

    h_extraAK4Jet1->Draw("hist");
    h_extraAK4Jet2->Draw("histsame");


    leg->AddEntry(h_extraAK4Jet1,"ZpBaryonic");
    leg->AddEntry(h_extraAK4Jet2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"extra AK4 Jets");
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();

     // higgs jet mass
    TH1F *h_higgsJetMass1 = (TH1F*)f_ZpBY->Get("h_higgsJetMass");
    TH1F *h_higgsJetMass2 = (TH1F*)f_ZpHS->Get("h_higgsJetMass");
    h_higgsJetMass1->SetLineColor(2);
    h_higgsJetMass2->SetLineColor(4);
    //h_pfMet1->GetXaxis()->SetTitleFont(22);
    //h_pfMet1->GetYaxis()->SetTitleFont(22);
   
    scale1 = 1.0/h_higgsJetMass1->Integral();
    scale2 = 1.0/h_higgsJetMass2->Integral();

    h_higgsJetMass1->Scale(scale1);
    h_higgsJetMass2->Scale(scale2);

    h_higgsJetMass1->Sumw2();
    h_higgsJetMass2->Sumw2();

    h_higgsJetMass1->GetYaxis()->SetTitle("Normalized to 1");
    h_higgsJetMass1->GetYaxis()->SetTitleSize(0.04);
    h_higgsJetMass1->GetYaxis()->CenterTitle();
    h_higgsJetMass1->GetXaxis()->SetTitle("M (GeV)");
    h_higgsJetMass1->GetXaxis()->SetTitleSize(0.04);
    h_higgsJetMass1->Draw("hist");
    h_higgsJetMass2->Draw("histsame");


    leg->AddEntry(h_extraAK4Jet1,"ZpBaryonic");
    leg->AddEntry(h_extraAK4Jet2,"ZpHS");
  
    leg->Draw();    
    Tl.DrawLatex(0.5,0.96,"higgs Jet Mass");
    pdfName += ")";
    c1->Print(pdfName.data());
    c1->Clear();
    c1->cd();
    leg->Clear();
}

void compare() {
    string macrolist[2] = {"resolved", "boosted"};
    for (int i=0; i<2; i++){
        comparebase("10000", "1000", macrolist[i].data());
        comparebase("10000", "150", macrolist[i].data());
        comparebase("10000", "500", macrolist[i].data());
        comparebase("10000", "50", macrolist[i].data());
        comparebase("1000", "1000", macrolist[i].data());
        comparebase("1000", "150", macrolist[i].data());
        comparebase("100", "10", macrolist[i].data());
        comparebase("10", "1000", macrolist[i].data());
        comparebase("10", "10", macrolist[i].data());
        comparebase("10", "1", macrolist[i].data());
        comparebase("10", "500", macrolist[i].data());
        comparebase("10", "50", macrolist[i].data());
        comparebase("15", "10", macrolist[i].data());
        comparebase("2000", "1", macrolist[i].data());
        comparebase("2000", "500", macrolist[i].data());
        comparebase("200", "150", macrolist[i].data());
        comparebase("200", "1", macrolist[i].data());
        comparebase("200", "50", macrolist[i].data());
        comparebase("20", "1", macrolist[i].data());
        comparebase("295", "150", macrolist[i].data());
        comparebase("300", "1", macrolist[i].data());
        comparebase("300", "50", macrolist[i].data());
        comparebase("500", "150", macrolist[i].data());
        comparebase("500", "1", macrolist[i].data());
        comparebase("500", "500", macrolist[i].data());
        comparebase("50", "10", macrolist[i].data());
        comparebase("50", "1", macrolist[i].data());
        comparebase("50", "50", macrolist[i].data());
    }
}