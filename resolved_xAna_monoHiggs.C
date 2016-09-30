// example code to run 2015 mono-Higgs resolved selections on signal (EXO-16-012)

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
#include <TGraph.h>
#include <fstream>

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
void resolved_xAna_monoHiggsBase(std::string inputFile){
  //get TTree from file ...
  TreeReader data(Form("MonoH-%s/0000/NCUGlobalTuples_1.root",inputFile.data()));

  TString endfix;
  endfix=gSystem->GetFromPipe(Form("file=%s; test=${file%%/crab*}; echo \"${test}\"",inputFile.data()));
  outputFile=Form("%s_resolved.txt",endfix.Data());
  

  Long64_t nTotal=0;
  Long64_t nPass[20]={0};
  TCanvas* c1 = new TCanvas("c1","",889*1.5,768);
 // TCanvas* c2 = new TCanvas("c2","",889*1.5,768);
  //TH1F* h_pfMet_fin = new TH1F("h_pfMet_fin", "phMet_fin", 25,0,2500);
  TH1F* h_pfMet = new TH1F("h_pfMet", "phMet", 25,0,2500);  
  TH1F* h_higgsPt = new TH1F("h_higgsPt", "higgs Pt", 20,0,2000);
  TH1F* h_higgsEta = new TH1F("h_higgsEta", "higgs Eta", 30,-3,3);

  //TH1F* h_deltaR_0 = new TH1F("h_deltaR_0", "deltaR_0", 15,-0.05,1.45);
  //TH1F* h_deltaR_1 = new TH1F("h_deltaR_1", "deltaR_1", 15,-0.05,1.45);
  TH1F* h_deltaR_subjet = new TH1F("h_deltaR_subjet", "deltaR_subjet", 20,0.,2);
  TH1F* h_extraEle = new TH1F("h_extraEle", "extra electrons", 6,-0.5,5.5);
  TH1F* h_extraMuo = new TH1F("h_extraMuo", "extra muons", 6,-0.5,5.5);
  TH1F* h_extraTau = new TH1F("h_extraTau", "extra tau", 6,-0.5,5.5);
  TH1F* h_extrabj = new TH1F("h_extraBJet", "extra b jet", 8,-0.5,7.5);
  TH1F* h_extraAK4j = new TH1F("h_extraAK4Jet", "extra AK4 jet", 8,-0.5,7.5);
  TH1F* h_mindphi = new TH1F("h_mindphi", "Minimum delta phi", 17,0,3.4);
  TH1F* h_higgsJetM = new TH1F("h_higgsJetMass", "higgs jet Mass", 24,50,170);

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
    if(pfMet<170.)continue;
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
    if(myEles.size()>0)  continue;
    if(myMuos.size()>0)  continue;
    if(myTaus.size()>0)  continue;
    nPass[4]++;
    nPass[5]++;
    nPass[6]++;
				      

    //find a pair of b-jets that could be a Higgs candidate
    const int nTHINJets     = data.GetInt("THINnJet");
    TClonesArray* thinjetP4 = (TClonesArray*) data.GetPtrTObject("THINjetP4");
    float* thinJetCSV =  data.GetPtrFloat("THINjetCISVV2");
    vector<bool>& passThinJetLooseID = *((vector<bool>*) data.GetPtr("THINjetPassIDLoose"));
    vector<bool>& passThinJetPUID = *((vector<bool>*) data.GetPtr("THINisPUJetID"));    
    
    float maxHpt=-999;
    int Hindex[2]={-1,-1};

    for(int ij=0; ij < nTHINJets; ij++){
      TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
      if(thisJet->Pt()<30)continue;
      if(fabs(thisJet->Eta())>2.4)continue;
      if(!passThinJetLooseID[ij])continue;
      if(!passThinJetPUID[ij])continue;
      
      // for b-jet (medium ID)
      if(thinJetCSV[ij]<0.80)continue;

      for(int jj=0; jj < ij ; jj++){
	TLorentzVector* thatJet = (TLorentzVector*)thinjetP4->At(jj);
	if(thatJet->Pt()<30)continue;
	if(fabs(thatJet->Eta())>2.4)continue;
	if(!passThinJetLooseID[jj])continue;
	if(!passThinJetPUID[jj])continue;
	
	// for b-jet (medium ID)
	if(thinJetCSV[jj]<0.80)continue;

	float thisHpt = (*thisJet + *thatJet).Pt();
	float thisHmass = (*thisJet + *thatJet).M();
	if(thisHpt < 150.)continue;
      h_higgsJetM->Fill(thisHmass);
	if(thisHmass < 100)continue;
	if(thisHmass > 150)continue;
	
	if(thisHpt>maxHpt)
	  {
	    Hindex[0]=jj;
	    Hindex[1]=ij;
	    maxHpt   =thisHpt;
	  }

      } // end of inner loop jet


    } // end of outer loop jet


    if(Hindex[0]<0 || Hindex[1]<0)continue;
    nPass[7]++;

    TLorentzVector  bjet[2];
    for(int ib=0; ib<2;ib++)bjet[ib] = 
			      *((TLorentzVector*)thinjetP4->At(Hindex[ib]));
    TLorentzVector  higgsJet = bjet[0]+bjet[1];
    float higgsPt, higgsEta, higgsM;
    higgsPt = higgsJet.Pt();
    higgsEta = higgsJet.Eta();
    higgsM = higgsJet.M();
    //h_higgsJetM->Fill(higgsM);
    h_deltaR_subjet->Fill(bjet[0].DeltaR(bjet[1]));
    // cout << "Hindex [0] = " << Hindex[0] << endl;
    // cout << "Hindex [1] = " << Hindex[1] << endl;

    //veto n>=1 AK4 jets and extra AK4 b jet
    
    unsigned int nGoodTHINBJets=0;
    unsigned int nGoodTHINJets=0;
    int jetIndex=-1; // extra AK4 jets if there is one
    std::vector<int> indexForDPhi;
    indexForDPhi.clear();

    for(int ij=0; ij < nTHINJets; ij++){
      TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(ij);
      if(thisJet->Pt()<30)continue;
      if(fabs(thisJet->Eta())>4.5)continue;
      if(!passThinJetLooseID[ij])continue;
      if(!passThinJetPUID[ij])continue;

      indexForDPhi.push_back(ij);

      if(thisJet->DeltaR(bjet[0])<0.4)continue;
      if(thisJet->DeltaR(bjet[1])<0.4)continue;
      

      nGoodTHINJets++;
      jetIndex=ij;

      // for b-jet
      if(fabs(thisJet->Eta())>2.4)continue;
      if(thinJetCSV[ij]<0.46)continue;
      nGoodTHINBJets++;

    } // end of loop

    h_extraAK4j->Fill(nGoodTHINJets); // extra AK4Jet
    h_extrabj->Fill(nGoodTHINBJets); // extra b Jet
    if(nGoodTHINBJets>0)   continue; 
    if(nGoodTHINJets>1) continue;
    nPass[8]++;    
    nPass[9]++;


    bool passDphi=true;
    TLorentzVector AK4Bjet(0,0,0,0);
    float mindphi = 99;
    for(unsigned int i=0; i<indexForDPhi.size(); i++)
      {
	int jetIndex=indexForDPhi[i];
	TLorentzVector* thisJet = (TLorentzVector*)thinjetP4->At(jetIndex); 
       //AK4Bjet +=  *thisJet;
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


   // h_pfMet_fin->Fill(pfMet);
    h_higgsPt->Fill(higgsPt);
    h_higgsEta->Fill(higgsEta);
    //h_deltaR_0->Fill(higgsJet.DeltaR(bjet[0]));
    //h_deltaR_1->Fill(higgsJet.DeltaR(bjet[1]));

    //h_higgsJetM->Fill(higgsM);
  } // end of loop over entries


  std::cout << "nTotal    = " << nTotal << std::endl;
  for(int i=0;i<20;i++)
    if(nPass[i]>0)
      std::cout << "nPass[" << i << "]= " << nPass[i] << std::endl;


  efferr(nPass[10],nTotal);
  //std::cout << "Minimum deltaPhi between the AK4 jets and MET = " << mindphi << endl;
  // draw and save figure
  /*h_pfMet->Draw("hist");
  c1->SaveAs("h_phMet_resolved.png");
  h_pfMet_fin->Draw("hist");
  c1->SaveAs("h_phMet_fin_resolved.png");
  h_higgsPt->Draw("hist");
  c1->SaveAs("higgsPt_resolved.png");
  h_higgsEta->Draw("hist");
  c1->SaveAs("higgsEta_resolved.png");
*/
  /*h_deltaR_0->Draw("hist");
  c1->SaveAs("deltaR_0_resolved.png");
  h_deltaR_1->Draw("hist");
  c1->SaveAs("deltaR_1_resolved.png");*/
  /*h_deltaR_subjet->Draw("hist");
  c1->SaveAs("deltaR_subjet_resolved.png");

  h_extraEle->Draw("hist");
  c1->SaveAs("extraElectron_resolved.png");
  h_extraMuo->Draw("hist");
  c1->SaveAs("extraMuon_resolved.png");
  h_extraTau->Draw("hist");
  c1->SaveAs("extraTau_resolved.png");

  h_extrabj->Draw("hist");
  c1->SaveAs("extraBJet_resolved.png");
  h_extraAK4j->Draw("hist");
  c1->SaveAs("extraAK4jet_resolved.png");

  h_higgsJetM->Draw("hist");
  c1->SaveAs("h_higgsJetM_resolved.png");
*/
  outputRootFile=Form("%s_resolved.root",endfix.Data());

  TFile* outFile_r = new TFile(outputRootFile.data(),"recreate");
  h_pfMet->Write();
  h_higgsPt->Write();
  h_higgsEta->Write();
  //deltaR_0->Write();
  //deltaR_1->Write();
  h_deltaR_subjet->Write();
  h_extraEle->Write();
  h_extraMuo->Write();
  h_extraTau->Write();
  h_extrabj->Write();
  h_extraAK4j->Write();
  h_higgsJetM->Write();

  outFile_r->Close();

}

void resolved_xAna_monoHiggs(){

  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi1000_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi1000_hbb/160904_134335");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi10_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi10_hbb/160904_134511");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi150_hbb/160904_134422");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi1_hbb/160904_134649");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi500_hbb/160904_134558");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10000_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp10000_MChi50_hbb/160904_134825");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp1000_MChi1000_hbb/crab_MonoH-Scalar_Scalar_MZp1000_MChi1000_hbb/160904_134737");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp1000_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp1000_MChi150_hbb/160904_135001");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp1000_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp1000_MChi1_hbb/160904_134913");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi1000_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi1000_hbb/160904_135049");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi150_hbb/160904_135316");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi1_hbb/160904_135226");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi500_hbb/160904_135457");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp10_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp10_MChi50_hbb/160904_135404");

  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp15_MChi10_hbb/crab_MonoH-Scalar_Scalar_MZp15_MChi10_hbb/160904_135633");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp200_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp200_MChi150_hbb/160904_135544");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp200_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp200_MChi1_hbb/160904_135813");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp200_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp200_MChi50_hbb/160904_135721");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp20_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp20_MChi1_hbb/160904_135947");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp295_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp295_MChi150_hbb/160904_135859");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp300_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp300_MChi1_hbb/160904_140124");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp500_MChi150_hbb/crab_MonoH-Scalar_Scalar_MZp500_MChi150_hbb/160904_140035");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp500_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp500_MChi500_hbb/160904_140302");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp50_MChi10_hbb/crab_MonoH-Scalar_Scalar_MZp50_MChi10_hbb/160904_140211");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp50_MChi1_hbb/crab_MonoH-Scalar_Scalar_MZp50_MChi1_hbb/160904_140438");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp50_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp50_MChi50_hbb/160904_140348");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp95_MChi50_hbb/crab_MonoH-Scalar_Scalar_MZp95_MChi50_hbb/160904_140624");
  resolved_xAna_monoHiggsBase("Scalar_Scalar_MZp995_MChi500_hbb/crab_MonoH-Scalar_Scalar_MZp995_MChi500_hbb/160904_140534");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi1000_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi1000_hbb/160904_140712");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi150_hbb/160904_140804");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi500_hbb/160904_140852");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10000_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10000_MChi50_hbb/160904_140940");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp1000_MChi1000_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp1000_MChi1000_hbb/160904_141029");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp1000_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp1000_MChi150_hbb/160904_141118");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp100_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp100_MChi10_hbb/160904_141206");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi1000_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi1000_hbb/160904_141255");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi10_hbb/160904_141341");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi1_hbb/160904_141431");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi500_hbb/160904_141524");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp10_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp10_MChi50_hbb/160904_141638");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp15_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp15_MChi10_hbb/160904_141729");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp2000_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp2000_MChi1_hbb/160904_141817");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp2000_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp2000_MChi500_hbb/160904_141908");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp200_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp200_MChi150_hbb/160904_141956");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp200_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp200_MChi1_hbb/160904_142058");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp200_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp200_MChi50_hbb/160904_142150");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp20_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp20_MChi1_hbb/160904_142237");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp295_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp295_MChi150_hbb/160904_142327");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp300_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp300_MChi1_hbb/160904_142415");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp300_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp300_MChi50_hbb/160904_142504");
  
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp500_MChi150_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp500_MChi150_hbb/160904_142555");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp500_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp500_MChi1_hbb/160904_142645");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp500_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp500_MChi500_hbb/160904_142732");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp50_MChi10_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp50_MChi10_hbb/160904_142822");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp50_MChi1_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp50_MChi1_hbb/160904_142909");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp50_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp50_MChi50_hbb/160904_143002");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp95_MChi50_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp95_MChi50_hbb/160904_143102");
  resolved_xAna_monoHiggsBase("ZpBaryonic_ZpBaryonic_MZp995_MChi500_hbb/crab_MonoH-ZpBaryonic_ZpBaryonic_MZp995_MChi500_hbb/160904_143203");
  
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi1000_hbb/160905_192611");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi10_hbb/160905_192728");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi150_hbb/160905_192843");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi1_hbb/160905_193000");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi500_hbb/160905_193119");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10000_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp10000_MChi50_hbb/160905_193235");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1000_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp1000_MChi1000_hbb/160905_193352");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1000_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp1000_MChi150_hbb/160905_193509");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1000_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp1000_MChi1_hbb/160905_193625");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp100_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp100_MChi10_hbb/160905_193741");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp100_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp100_MChi1_hbb/160905_193859");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi1000_hbb/160905_194016");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi10_hbb/160905_194139");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi150_hbb/160905_194254");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi1_hbb/160905_194410");
  
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi500_hbb/160905_194527");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp10_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp10_MChi50_hbb/160905_194644");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp15_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp15_MChi10_hbb/160905_194805");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp1995_MChi1000_hbb/crab_MonoH-ZpHS_ZpHS_MZp1995_MChi1000_hbb/160905_194927");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp2000_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp2000_MChi1_hbb/160905_195046");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp2000_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp2000_MChi500_hbb/160905_195208");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp200_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp200_MChi150_hbb/160905_195327");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp200_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp200_MChi1_hbb/160905_195443");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp200_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp200_MChi50_hbb/160905_195559");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp20_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp20_MChi1_hbb/160905_195716");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp295_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp295_MChi150_hbb/160905_195832");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp300_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp300_MChi1_hbb/160905_195947");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp300_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp300_MChi50_hbb/160905_200103");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp500_MChi150_hbb/crab_MonoH-ZpHS_ZpHS_MZp500_MChi150_hbb/160905_200224");

  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp500_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp500_MChi1_hbb/160905_200432");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp500_MChi500_hbb/crab_MonoH-ZpHS_ZpHS_MZp500_MChi500_hbb/160905_200603");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp50_MChi10_hbb/crab_MonoH-ZpHS_ZpHS_MZp50_MChi10_hbb/160905_200747");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp50_MChi1_hbb/crab_MonoH-ZpHS_ZpHS_MZp50_MChi1_hbb/160905_200905");
  resolved_xAna_monoHiggsBase("ZpHS_ZpHS_MZp50_MChi50_hbb/crab_MonoH-ZpHS_ZpHS_MZp50_MChi50_hbb/160905_201022");
  //resolved_xAna_monoHiggsBase("");
}
