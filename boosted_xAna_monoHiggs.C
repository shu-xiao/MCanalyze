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
#include <TCanvas.h>

void efferr(float nsig,float ntotal,float factor=1)
{
  float eff = nsig/ntotal;
  float err = sqrt( (1-eff)*eff/ntotal);
  cout << "efficiency = " << eff*factor << " +- " << err*factor << endl;

}


using namespace std;
void boosted_xAna_monoHiggs(std::string inputFile){

  //get TTree from file ...
  TreeReader data(inputFile.data());

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  TCanvas* c1 = new TCanvas("c1","",889*1.5,768);
 // TCanvas* c2 = new TCanvas("c2","",889*1.5,768);
  TH1F* h_pfMet = new TH1F("h_pfMet", "phMet", 22,-50,2150);
  TH1F* h_higgsPt = new TH1F("h_higgsPt", "higgs Pt", 15,-50,1450);
  TH1F* h_higgsEta = new TH1F("h_higgsEta", "higgs Eta", 30,-3.1,2.9);

  TH1F* deltaR_0 = new TH1F("h_deltaR_0", "deltaR_0", 15,-0.05,1.45);
  TH1F* deltaR_1 = new TH1F("h_deltaR_1", "deltaR_1", 15,-0.05,1.45);
  TH1F* deltaR_subjet = new TH1F("h_deltaR_subjet", "deltaR_subjet", 20,-0.05,1.95);
  TH1F* extraEle = new TH1F("h_extraEle", "extra electrons", 6,-0.5,5.5);
  TH1F* extraMuo = new TH1F("h_extraMuo", "extra muons", 6,-0.5,5.5);
  TH1F* extraTau = new TH1F("h_extraTau", "extra tau", 6,-0.5,5.5);
  TH1F* extrabj = new TH1F("h_extraBJet", "extra b jet", 8,-0.5,7.5);
  TH1F* extraAK4j = new TH1F("h_extraAK4Jet", "extra AK4 jet", 8,-0.5,7.5);

  TH1F* h_higgsJetM = new TH1F("h_higgsJetMass", "higgs jet Mass", 20,100-2.5,200-2.5);
  float mindphi = 99;
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
	    //	    cout << thisTrig << endl;
	    passTrigger=true;
	    break;
	  }


      }

    if(!passTrigger)continue;
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
    if( isData && !passFilter )continue;
    nPass[2]++;
	



    float pfMet = data.GetFloat("pfMetCorrPt");
    float pfMetPhi = data.GetFloat("pfMetCorrPhi");
    
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

    if(myEles.size()>0)  extraEle->Fill(myEles.size());
    if(myMuos.size()>0)  extraMuo->Fill(myMuos.size());
    if(myTaus.size()>0)  extraTau->Fill(myTaus.size());
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
    vector<float>   *subjetSDPy  =  data.GetPtrVectorFloat("FATsubjetSDPy", nFATJets);
    vector<float>   *subjetSDPz  =  data.GetPtrVectorFloat("FATsubjetSDPz", nFATJets);
    vector<float>   *subjetSDE   =  data.GetPtrVectorFloat("FATsubjetSDE", nFATJets);  
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

    TLorentzVector* higgsJet = (TLorentzVector*)fatjetP4->At(HIndex);
    float higgsPt, higgsEta, higgsM;
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
    TLorentzVector* thinJet = (TLorentzVector*)thinjetP4->At(HIndex);
  /*  TLorentzVector  bjet[2];
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

    if(nGoodTHINJets>1) extraAK4j->Fill(nGoodTHINJets); // extra AK4 jet
    if(nGoodTHINBJets>0) extrabj->Fill(nGoodTHINBJets); // extra b jet
    if(nGoodTHINJets>1)continue;
    if(nGoodTHINBJets>0)continue;
    nPass[8]++;
    nPass[9]++;

    bool passDphi=true;
    TLorentzVector AK4Bjet;
    for(unsigned int i=0; i<indexForDPhi.size(); i++)
      {
	int jetIndex=indexForDPhi[i];
	TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(jetIndex);
       AK4Bjet +=  *thisJet;
	double dphi=TVector2::Phi_mpi_pi(pfMetPhi-thisJet->Phi());
       float tem_mindfhp = 99;
       if (tem_mindfhp > fabs(dphi)) tem_mindfhp = fabs(dphi);
	if(fabs(dphi)<0.4)
	  {
	    passDphi=false;
	    //break;
	  }	
       if (tem_mindfhp < mindphi && passDphi && i == indexForDPhi.size()-1) mindphi = fabs(dphi);
      }
    h_higgsJetM->Fill(AK4Bjet.M());
    if(!passDphi)continue;
    nPass[10]++;

    
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


    h_pfMet->Fill(pfMet);
    h_higgsPt->Fill(higgsPt);
    h_higgsEta->Fill(higgsEta);
  //  deltaR_0->Fill(higgsJet->DeltaR(bjet[0]));
    //deltaR_1->Fill(higgsJet->DeltaR(bjet[1]));
    deltaR_subjet->Fill(thinJet->DeltaR(*higgsJet));
    h_higgsJetM->Fill(AK4Bjet.M());
  } // end of loop over entries


  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;


  efferr(nPass[11],nTotal);
  std::cout << "Minimum deltaPhi between the AK4 jets and MET = " << mindphi << endl;
  h_pfMet->Draw("hist");
  c1->SaveAs("h_phMet_boosted.png");
  h_higgsPt->Draw("hist");
  c1->SaveAs("higgsPt_boosted.png");
  h_higgsEta->Draw("hist");
  c1->SaveAs("higgsEta_boosted.png");

/*  deltaR_0->Draw("hist");
  c1->SaveAs("deltaR_0_boosted.png");
  deltaR_1->Draw("hist");
  c1->SaveAs("deltaR_1_boosted.png");*/
  deltaR_subjet->Draw("hist");
  c1->SaveAs("deltaR_subjet_boosted.png");

  extraEle->Draw("hist");
  c1->SaveAs("extraElectron_boosted.png");
  extraMuo->Draw("hist");
  c1->SaveAs("extraMuon_boosted.png");
  extraTau->Draw("hist");
  c1->SaveAs("extraTau_boosted.png");

  extrabj->Draw("hist");
  c1->SaveAs("extraBJet_boosted.png");
  extraAK4j->Draw("hist");
  c1->SaveAs("extraAK4jet_boosted.png");

  h_higgsJetM->Draw("hist");
  c1->SaveAs("h_higgsJetM_boosted.png");

  TFile* outFile = new TFile("plot_boosted.root","recreate");
  h_pfMet->Write();
  h_higgsPt->Write();
  h_higgsEta->Write();
  deltaR_0->Write();
  deltaR_1->Write();
  deltaR_subjet->Write();
  extraEle->Write();
  extraMuo->Write();
  extraTau->Write();
  extrabj->Write();
  extraAK4j->Write();
  h_higgsJetM->Write();

  outFile->Close();
  // TFile* outFile = new TFile("test.root","recreate");

  // h_hh->Write();
  // h_hh2->Write();
  // h_SD1->Write();
  // h_SD2->Write();

  // outFile->Close();



}
