// example code to run 2015 mono-Higgs boosted selections on signal (EXO-16-012)

#include <vector>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <TH1D.h>
#include <TFile.h>
#include "untuplizer.h"
#include <TClonesArray.h>
#include <TLorentzVector.h>
#include <fstream>
#include <TCanvas.h>
#include <TTreeReader.h>

string outputFile, outputRootFile;

void efferr(float nsig,float ntotal,float factor=1)
{
  float eff = nsig/ntotal;
  float err = sqrt( (1-eff)*eff/ntotal);
  cout << "efficiency = " << eff*factor << " +- " << err*factor << endl;
  ofstream myfile;
  myfile.open (outputFile.data());
  myfile<<eff*factor<<endl;
  myfile<<err*factor;
}


using namespace std;
void boosted_xAna_monoHiggsBase(std::string inputFile){

  //get TTree from file ...
//  TreeReader data(Form("/data7/khurana/MonoHPrivateSamples/PrivateSamplesForEfficiency/MonoH-%s/0000/NCUGlobalTuples_1.root",inputFile.data()));
  TreeReader data(Form("MonoH-%s/0000/NCUGlobalTuples_1.root",inputFile.data()));


  TString endfix;
  endfix=gSystem->GetFromPipe(Form("file=%s; test=${file%%/crab*}; echo \"${test}\"",inputFile.data()));
  outputFile=Form("%s_boosted.txt",endfix.Data());

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  TCanvas* c1 = new TCanvas("c1","",889*1.5,768);
 // TCanvas* c2 = new TCanvas("c2","",889*1.5,768);
  //TH1F* h_pfMet_fin = new TH1F("h_pfMet_fin", "phMet_fin", 25,0,2500);
  TH1F* h_pfMet = new TH1F("h_pfMet", "phMet", 25,0,2500);
  TH1F* h_higgsPt = new TH1F("h_higgsPt", "higgs Pt", 20,0,2000);
  TH1F* h_higgsEta = new TH1F("h_higgsEta", "higgs Eta", 30,-3,3);

  //TH1F* h_deltaR_0 = new TH1F("h_deltaR_0", "deltaR_0", 15,0,1.5);
  //TH1F* h_deltaR_1 = new TH1F("h_deltaR_1", "deltaR_1", 15,0,1.5);
  TH1F* h_deltaR_subjet = new TH1F("h_deltaR_subjet", "deltaR_subjet", 20,0,2);
  TH1F* h_extraEle = new TH1F("h_extraEle", "extra electrons", 6,-0.5,5.5);
  TH1F* h_extraMuo = new TH1F("h_extraMuo", "extra muons", 6,-0.5,5.5);
  TH1F* h_extraTau = new TH1F("h_extraTau", "extra tau", 6,-0.5,5.5);
  TH1F* h_extrabJ = new TH1F("h_extraBJet", "extra b jet", 8,-0.5,7.5);
  TH1F* h_extraAK4j = new TH1F("h_extraAK4Jet", "extra AK4 jet", 8,-0.5,7.5);
  TH1F* h_mindphi = new TH1F("h_mindphi", "Minimum delta phi", 17,0,3.4);

  TH1F* h_higgsJetM = new TH1F("h_higgsJetMass", "higgs jet Mass", 24,50,170);
  int Hindex[2]={-1,-1};

  for(Long64_t jEntry=0; jEntry<data.GetEntriesFast() ;jEntry++){

    if (jEntry % 50000 == 0)
      fprintf(stderr, "Processing event %lli of %lli\n", jEntry + 1, data.GetEntriesFast());

    data.GetEntry(jEntry);
    nTotal++;
    //0. has a good vertex
    int nVtx        = data.GetInt("nVtx");
    if(nVtx<1)continue;
    nPass[0]++;

    //1. trigger 
    std::string* trigName = data.GetPtrString("hlt_trigName");
    vector<bool> &trigResult = *((vector<bool>*) data.GetPtr("hlt_trigResult"));

    bool passTrigger=false;
    for(unsigned int it=0; it< trigResult.size(); it++)
      {
  std::string thisTrig= trigName[it];
  bool results = trigResult[it];

  if( (thisTrig.find("HLT_PFMET90_PFMHT90_")!= 
       std::string::npos && results==1) || 
      (thisTrig.find("HLT_PFMET170_NoiseCleaned")!= 
       std::string::npos && results==1))
    {
      //      cout << thisTrig << endl;
      passTrigger=true;
      break;
    }


      }

    //if(!passTrigger)continue;
    nPass[1]++;

    // apply noise filters if it is data
    bool isData = data.GetBool("isData");
    std::string* filterName = data.GetPtrString("hlt_filterName");
    vector<bool> &filterResult = *((vector<bool>*) data.GetPtr("hlt_filterResult"));
    bool passFilter=false;
    for(unsigned int it=0; it< filterResult.size(); it++)
      {
  std::string thisFilter= filterName[it];
  bool results = filterResult[it];

  if( (thisFilter.find("Flag_CSCTightHaloFilter")!= 
       std::string::npos && results==1) &&
      (thisFilter.find("Flag_eeBadScFilter")!= 
       std::string::npos && results==1) &&
      (thisFilter.find("Flag_HBHENoiseFilter")!= 
       std::string::npos && results==1) &&
      (thisFilter.find("Flag_HBHENoiseIsoFilter")!= 
       std::string::npos && results==1) )
    {
      passFilter=true;
      break;
    } 
      }
    //if( isData && !passFilter )continue;
    nPass[2]++;
  



    float pfMet = data.GetFloat("pfMetCorrPt");
    float pfMetPhi = data.GetFloat("pfMetCorrPhi");
    h_pfMet->Fill(pfMet);
    if(pfMet<200.)continue;
    nPass[3]++;



    //veto extra electrons
    int    nEle       = data.GetInt("nEle");
    TClonesArray* eleP4 = (TClonesArray*) data.GetPtrTObject("eleP4");
    vector<bool>& eleIsPassLoose = *((vector<bool>*) data.GetPtr("eleIsPassLoose"));
    vector<int> myEles;
    myEles.clear();
    for(int ie = 0; ie < nEle; ie++){    
      TLorentzVector* myEle = (TLorentzVector*)eleP4->At(ie);
      if( myEle->Pt()<10 )continue;
      if( fabs(myEle->Eta())>2.5 )continue;
      if( !eleIsPassLoose[ie] )continue;
      myEles.push_back(ie);
    } 
    //if(myEles.size()>0)continue;
    nPass[4]++;

    //veto extra muons
    int    nMu       = data.GetInt("nMu");
    TClonesArray* muP4 = (TClonesArray*) data.GetPtrTObject("muP4");
    vector<bool>& isLooseMuon = *((vector<bool>*) data.GetPtr("isLooseMuon"));
    float* muChHadIso = data.GetPtrFloat("muChHadIso");
    float* muNeHadIso = data.GetPtrFloat("muNeHadIso");
    float* muGamIso   = data.GetPtrFloat("muGamIso");
    float* muPUPt     = data.GetPtrFloat("muPUPt");
       
    vector<int> myMuos;
    myMuos.clear();
    for(int im = 0; im < nMu; im++){
      TLorentzVector* myMu = (TLorentzVector*)muP4->At(im);
      if( myMu->Pt()<10 )continue;
      if( fabs(myMu->Eta())>2.4 ) continue;
      if( !isLooseMuon[im] )continue;
      
      float relPFIso = (muChHadIso[im]+ 
      TMath::Max(0., muNeHadIso[im] + muGamIso[im] - 0.5*muPUPt[im]))/myMu->Pt();
      if(relPFIso>0.4)continue;
      myMuos.push_back(im);
    }
    //if(myMuos.size()>0)continue;
    nPass[5]++;

    //veto extra taus
    int    nTau      = data.GetInt("HPSTau_n");
    TClonesArray* tauP4 = (TClonesArray*) data.GetPtrTObject("HPSTau_4Momentum");
    vector<bool>& isDecayModeFinding = *((vector<bool>*) data.GetPtr("disc_decayModeFinding"));
    vector<bool>& passLooseTauIso = *((vector<bool>*) data.GetPtr("disc_byLooseIsolationMVA3oldDMwLT"));
   
    vector<int> myTaus;
    for(int it=0; it < nTau; it++)
      {
  TLorentzVector* myTau = (TLorentzVector*)tauP4->At(it);
  if( myTau->Pt()<20 )continue;
  if( fabs(myTau->Eta())>2.3 )continue;
  if( !isDecayModeFinding[it] )continue;
  if( !passLooseTauIso[it] )continue;
  myTaus.push_back(it);
      }

    h_extraEle->Fill(myEles.size());
    h_extraMuo->Fill(myMuos.size());
    h_extraTau->Fill(myTaus.size());
    if(myTaus.size()>0)continue;
    if(myEles.size()>0)continue;
    if(myMuos.size()>0)continue;
    nPass[6]++;
              

    int nFATJets     = data.GetInt("FATnJet");
    TClonesArray* fatjetP4 = (TClonesArray*) data.GetPtrTObject("FATjetP4");
    float*  fatjetPRmassL2L3Corr = data.GetPtrFloat("FATjetPRmassL2L3Corr");
    int*   nSubSoftDropJet = data.GetPtrInt("FATnSubSDJet");
    vector<float>   *subjetSDCSV =  data.GetPtrVectorFloat("FATsubjetSDCSV", nFATJets);
    vector<float>   *subjetSDPx  =  data.GetPtrVectorFloat("FATsubjetSDPx", nFATJets);
    vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJets); // integral number need fix
    vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJets);
    vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJets);  // subjet 4 vector
    vector<bool>    &passFatJetTightID = *((vector<bool>*) data.GetPtr("FATjetPassIDTight"));
    

    int HIndex=-1;
    for(int ij=0; ij<nFATJets; ij++)
      {
      
      TLorentzVector* thisJet = (TLorentzVector*)fatjetP4->At(ij);
      if(thisJet->Pt()<200)continue;
  if(fabs(thisJet->Eta())>2.4)continue;
  if(!passFatJetTightID[ij])continue;
  HIndex=ij;
  break; // only take the first one also the leading-pt one
      }
    
    if(HIndex<0)continue;
    nPass[7]++;
    TLorentzVector* subjetP4[2];
    //int i=HIndex;
    //for(int i=0;i<2;i++){
      for(int j=0;j<2;j++){

            subjetP4[j]=new TLorentzVector;
            subjetP4[j]->SetPxPyPzE(subjetSDPx[HIndex][j],subjetSDPy[HIndex][j],subjetSDPz[HIndex][j],subjetSDE[HIndex][j]);
           }
       h_deltaR_subjet->Fill(subjetP4[0]->DeltaR(*subjetP4[1]));
      //}
    TLorentzVector* higgsJet = (TLorentzVector*)fatjetP4->At(HIndex);
    float higgsPt, higgsEta;
    higgsPt = higgsJet->Pt();
    higgsEta = higgsJet->Eta();

    //veto n>=1 AK4 jets and n>=0 AK4 b jet
    const int nTHINJets     = data.GetInt("THINnJet");
    TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
    float* thinJetCSV =  data.GetPtrFloat("THINjetCISVV2");
    vector<bool>& passThinJetLooseID = *((vector<bool>*) data.GetPtr("THINjetPassIDLoose"));
    vector<bool>& passThinJetPUID = *((vector<bool>*) data.GetPtr("THINisPUJetID"));    
    
    unsigned int nGoodTHINBJets=0;
    unsigned int nGoodTHINJets=0;
    int jetIndex=-1;
    std::vector<int> indexForDPhi;
    indexForDPhi.clear();

    /*TLorentzVector  bjet[2];
    for(int ib=0; ib<2;ib++)bjet[ib] = 
            *((TLorentzVector*)thinjetP4->At(Hindex[ib]));
    TLorentzVector  higgsJetTV = bjet[0]+bjet[1];
    float higgsPt, higgsEta, higgsM;
    higgsPt = higgsJetTV->Pt();
    higgsEta = higgsJetTV->Eta();
    higgsM = higgsJetTV->M();*/

    for(int ij=0; ij < nTHINJets; ij++){
      TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij); //AK4 jet
      if(thisJet->Pt()<30)continue;
      if(fabs(thisJet->Eta())>4.5)continue;
      if(!passThinJetLooseID[ij])continue;
      if(!passThinJetPUID[ij])continue;
      indexForDPhi.push_back(ij);
      if(thisJet->DeltaR(*higgsJet)<0.8)continue;
      
      nGoodTHINJets++;
      jetIndex=ij;

      // for b-jet
      if(fabs(thisJet->Eta())>2.4)continue;
      if(thinJetCSV[ij]<0.46)continue; // loose b-tag
      nGoodTHINBJets++;

    } // end of loop

    h_extraAK4j->Fill(nGoodTHINJets); // extra AK4 jet
    h_extrabJ->Fill(nGoodTHINBJets); // extra b jet
    if(nGoodTHINJets>1)continue;
    if(nGoodTHINBJets>0)continue;
    nPass[8]++;
    nPass[9]++;

    bool passDphi=true;
    float mindphi = 99;
    for(unsigned int i=0; i<indexForDPhi.size(); i++)
      {
  int jetIndex=indexForDPhi[i];
  TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(jetIndex);
  double dphi=TVector2::Phi_mpi_pi(pfMetPhi-thisJet->Phi());
      if (mindphi > fabs(dphi)) mindphi = fabs(dphi);
  if(fabs(dphi)<0.4)
    {
      passDphi=false;
      //break;
    }
      }
    h_mindphi->Fill(mindphi);
    if(!passDphi)continue;
    nPass[10]++;
    
    h_higgsJetM->Fill(fatjetPRmassL2L3Corr[HIndex]); // wrong
    // additional mass cut and subjet cuts on Higgs
    if(fatjetPRmassL2L3Corr[HIndex]<100)continue;
    if(fatjetPRmassL2L3Corr[HIndex]>150)continue; 
    // require subjets b-tagged
    int nSubBJet=0;
    for(int is=0; is < nSubSoftDropJet[HIndex]; is++){
      if(subjetSDCSV[HIndex][is] < 0.46)continue; // loose b-tag
      nSubBJet++;
    }
    if(nSubBJet<2)continue;

    nPass[11]++;


    //h_pfMet_fin->Fill(pfMet);
    h_higgsPt->Fill(higgsPt);
    h_higgsEta->Fill(higgsEta);
  //  deltaR_0->Fill(higgsJet->DeltaR(bjet[0]));hh
    //deltaR_1->Fill(higgsJet->DeltaR(bjet[1]));
  } // end of loop over entries


  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;


  efferr(nPass[11],nTotal);
  //std::cout << "Minimum deltaPhi between the AK4 jets and MET = " << mindphi << endl;
  /*h_pfMet_fin->Draw("hist");
  c1->SaveAs("h_phMet_fin_boosted.png");
  h_pfMet->Draw("hist");
  c1->SaveAs("h_phMet_resolved.png");
  h_higgsPt->Draw("hist");
  c1->SaveAs("higgsPt_boosted.png");
  h_higgsEta->Draw("hist");
  c1->SaveAs("higgsEta_boosted.png");*/

  /*h_deltaR_0->Draw("hist");
  c1->SaveAs("deltaR_0_boosted.png");
  h_deltaR_1->Draw("hist");
  c1->SaveAs("deltaR_1_boosted.png");*/

  /*h_deltaR_subjet->Draw("hist");
  c1->SaveAs("deltaR_subjet_boosted.png");

  h_extraEle->Draw("hist");
  c1->SaveAs("extraElectron_boosted.png");
  h_extraMuo->Draw("hist");
  c1->SaveAs("extraMuon_boosted.png");
  h_extraTau->Draw("hist");
  c1->SaveAs("extraTau_boosted.png");

  h_extrabJ->Draw("hist");
  c1->SaveAs("extraBJet_boosted.png");
  h_extraAK4j->Draw("hist");
  c1->SaveAs("extraAK4jet_boosted.png");

  h_higgsJetM->Draw("hist");
  c1->SaveAs("h_higgsJetM_boosted.png");
  h_mindphi->Draw("hist");
  c1->SaveAs("mindphi_boosted.png");
  */
  outputRootFile=Form("%s_boosted.root",endfix.Data());
  TFile* outFile = new TFile(outputRootFile.data(),"recreate");
  //h_pfMet_fin->Write();
  h_pfMet->Write();
  h_higgsPt->Write();
  h_higgsEta->Write();
  //h_deltaR_0->Write();
  //h_deltaR_1->Write();
  h_deltaR_subjet->Write();
  h_extraEle->Write();
  h_extraMuo->Write();
  h_extraTau->Write();
  h_extrabJ->Write();
  h_extraAK4j->Write();
  h_higgsJetM->Write();
  h_mindphi->Write();
  outFile->Close();

  // TFile* outFile = new TFile("test.root","recreate");

  // h_hh->Write();
  // h_hh2->Write();
  // h_SD1->Write();
  // h_SD2->Write();

  // outFile->Close();



}

void boosted_xAna_monoHiggs(){

  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi1000_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi1000_hbb/160904_134335");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi10_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi10_hbb/160904_134511");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi150_hbb/160904_134422");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi1_hbb/160904_134649");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi500_hbb/160904_134558");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi50_hbb/160904_134825");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp1000_MChi1000_hbb/crab_MonoH-Scalar_Scalar_MZp1000_MChi1000_hbb/160904_134737");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp1000_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp1000_MChi150_hbb/160904_135001");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp1000_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp1000_MChi1_hbb/160904_134913");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi1000_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi1000_hbb/160904_135049");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi150_hbb/160904_135316");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi1_hbb/160904_135226");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi500_hbb/160904_135457");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi50_hbb/160904_135404");

  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp15_MChi10_hbb/crab_MonoH-Scalar_Scalar_MZp15_MChi10_hbb/160904_135633");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp200_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp200_MChi150_hbb/160904_135544");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp200_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp200_MChi1_hbb/160904_135813");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp200_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp200_MChi50_hbb/160904_135721");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp20_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp20_MChi1_hbb/160904_135947");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp295_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp295_MChi150_hbb/160904_135859");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp300_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp300_MChi1_hbb/160904_140124");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp500_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp500_MChi150_hbb/160904_140035");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp500_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp500_MChi500_hbb/160904_140302");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp50_MChi10_hbb/crab_MonoH-Scalar_Scalar_MZp50_MChi10_hbb/160904_140211");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp50_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp50_MChi1_hbb/160904_140438");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp50_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp50_MChi50_hbb/160904_140348");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp95_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp95_MChi50_hbb/160904_140624");
  boosted_xAna_monoHiggsBase("Scalar_Scalar_MZp995_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp995_MChi500_hbb/160904_140534");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi1000_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi1000_hbb/160904_140712");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi150_hbb/160904_140804");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi500_hbb/160904_140852");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi50_hbb/160904_140940");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp1000_MChi1000_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp1000_MChi1000_hbb/160904_141029");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp1000_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp1000_MChi150_hbb/160904_141118");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp100_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp100_MChi10_hbb/160904_141206");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi1000_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi1000_hbb/160904_141255");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi10_hbb/160904_141341");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi1_hbb/160904_141431");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi500_hbb/160904_141524");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi50_hbb/160904_141638");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp15_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp15_MChi10_hbb/160904_141729");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp2000_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp2000_MChi1_hbb/160904_141817");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp2000_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp2000_MChi500_hbb/160904_141908");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp200_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp200_MChi150_hbb/160904_141956");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp200_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp200_MChi1_hbb/160904_142058");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp200_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp200_MChi50_hbb/160904_142150");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp20_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp20_MChi1_hbb/160904_142237");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp295_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp295_MChi150_hbb/160904_142327");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp300_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp300_MChi1_hbb/160904_142415");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp300_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp300_MChi50_hbb/160904_142504");
  
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp500_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp500_MChi150_hbb/160904_142555");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp500_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp500_MChi1_hbb/160904_142645");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp500_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp500_MChi500_hbb/160904_142732");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp50_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp50_MChi10_hbb/160904_142822");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp50_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp50_MChi1_hbb/160904_142909");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp50_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp50_MChi50_hbb/160904_143002");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp95_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp95_MChi50_hbb/160904_143102");
  boosted_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp995_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp995_MChi500_hbb/160904_143203");
  
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi1000_hbb/160905_192611");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi10_hbb/160905_192728");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi150_hbb/160905_192843");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi1_hbb/160905_193000");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi500_hbb/160905_193119");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi50_hbb/160905_193235");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1000_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp1000_MChi1000_hbb/160905_193352");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1000_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp1000_MChi150_hbb/160905_193509");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1000_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp1000_MChi1_hbb/160905_193625");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp100_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp100_MChi10_hbb/160905_193741");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp100_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp100_MChi1_hbb/160905_193859");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi1000_hbb/160905_194016");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi10_hbb/160905_194139");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi150_hbb/160905_194254");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi1_hbb/160905_194410");
  
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi500_hbb/160905_194527");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi50_hbb/160905_194644");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp15_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp15_MChi10_hbb/160905_194805");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1995_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp1995_MChi1000_hbb/160905_194927");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp2000_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp2000_MChi1_hbb/160905_195046");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp2000_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp2000_MChi500_hbb/160905_195208");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp200_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp200_MChi150_hbb/160905_195327");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp200_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp200_MChi1_hbb/160905_195443");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp200_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp200_MChi50_hbb/160905_195559");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp20_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp20_MChi1_hbb/160905_195716");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp295_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp295_MChi150_hbb/160905_195832");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp300_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp300_MChi1_hbb/160905_195947");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp300_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp300_MChi50_hbb/160905_200103");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp500_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp500_MChi150_hbb/160905_200224");

  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp500_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp500_MChi1_hbb/160905_200432");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp500_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp500_MChi500_hbb/160905_200603");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp50_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp50_MChi10_hbb/160905_200747");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp50_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp50_MChi1_hbb/160905_200905");
  boosted_xAna_monoHiggsBase("ZpHS_ZpHS_MZp50_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp50_MChi50_hbb/160905_201022");
  //boosted_xAna_monoHiggsBase("");
}
