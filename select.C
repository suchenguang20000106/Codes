#define loose_cxx
#include "loose.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <iostream>

using namespace std;
using std:vector;

void loose::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L loose.C
//      Root > loose t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only Lcpair branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();
   //TH1F *mPi0 = new TH1F("mPi0","mPi0",100,0.11,0.16);

   TH1F *m_bc = new TH1F("Mbc","Mbc",100, 2.25, 2.3);
   m_bc->GetXaxis()->SetTitle("Mbc(GeV)");
   m_bc->GetYaxis()->SetTitle("Events");

   TH1F *mbc_sig = new TH1F("Mbc_sig","Mbc_signal",100, 2.25, 2.3);
   mbc_sig->GetXaxis()->SetTitle("Mbc(GeV)");
   mbc_sig->GetYaxis()->SetTitle("Events");

   TH1F *mbc_side = new TH1F("Mbc_side","Mbc_sideband",100, 2.25, 2.3);
   mbc_side->GetXaxis()->SetTitle("Mbc(GeV)");
   mbc_side->GetYaxis()->SetTitle("Events");

   TH1F *delta_E = new TH1F("delta_E","delta_E", 100, -0.1, 0.1);
   delta_E->GetXaxis()->SetTitle("deltaE(GeV)");
   delta_E->GetYaxis()->SetTitle("Events");
   
//   TFile root(rootname,"RECREATE");

   TObjArray hist(0);
   hist.Add(m_bc);
   hist.Add(mbc_sig);
   hist.Add(mbc_side);
   hist.Add(delta_E);

   //tree
   double MBC;
   double MBC_side;
   double MBC_sig;
   double DE;
   TTree *T = new TTree("Lcpair","Lcpair");
   T->Branch("Mbc", &MBC, "Mbc/D");
   T->Branch("deltaE", &DE, "deltaE/D");
   T->Branch("Mbc_side", &MBC_side, "Mbc_side/D");
   T->Branch("Mbc_sig", &MBC_sig, "Mbc_sig/D");

   for (Long64_t jentry=0; jentry<nentries; jentry++) 
   {  
       GetEntry(jentry);

       //select protons;
       vector<int> ip;
       for(Long64_t i = 0; i < n_P; i++)
       {
          if(p_Vr[i] > 0.2) continue;//Vr <= 0.2cm for protons
          ip.push_back(i);
       }
       if (ip.size() == 0) continue;
      
      //cout << "ip" << endl;

       //select photons;
       vector<int> igam;
       for(Long64_t i = 0; i < n_Gam; i++)
       {
          if(mdang[i] < 8) continue; //angle displacement for all charged tracks
          if(mpdang[i] < 20) continue; //angle displacement for antiprotons
          igam.push_back(i);
       }
       if(igam.size() == 0) continue;
      //cout << "igam" << endl;

       //select pi0
       vector<int> ipi0;
       for(Long64_t i = 0; i < igam.size()-1; i++)
       {
          for(Long64_t j = i + 1; j < igam.size(); j++)
          {  
             int sequence = igam[i]*n_Gam-igam[i]*(igam[i]+1)/2+igam[j]-igam[i]-1;
             if(chisquare[sequence] > 50)  continue; //chisquare requirement for karman fit
             if(Mpi0[sequence] > 0.152967 || Mpi0[sequence] < 0.116967) continue;//invariant mass requirement
             ipi0.push_back(sequence);
          }
       }
       if(ipi0.size() == 0) continue;

       //cout << "ipi0" << endl;
/*       for(Long64_t i = 0; i < n_Pi0; i++)
       {
          if(chisquare[i] > 50) continue;
          if(Mpi0[i] > 0.152967 || Mpi0[i] < 0.116967) continue;
          ipi0.push_back(i);
       }
*/

       //combination of pi0 and p
       const double ECMS = 2.3;
       double mmbc = 0.;
       double mdeltae = 0.2;
       for(Long64_t i = 0; i < ip.size(); i++)
       {
          for(Long64_t j = 0; j < ipi0.size(); j++)
          {
             //calculate p4 for Lc
             double p4lc[4];
             for (int k = 0; k < 4; k++)
             {
               //p4lc[k] = p4_P[i][k] + p4_Pi0[ipi0[j]][k];
               p4lc[k] = p4_Lc[ip[i]*n_Pi0 + ipi0[j]][k];
             }
             double deltae = p4lc[3] - ECMS;
             double mbc2 = ECMS*ECMS - p4lc[1]*p4lc[1] - p4lc[2]*p4lc[2] - p4lc[0]*p4lc[0];
             if(mbc2 < 0) continue;
             double mbc = sqrt(mbc2);
             if(fabs(deltae) < fabs(mdeltae) && mbc < 2.3 && mbc > 2.25) //deltae minimum and mbc in the cut range
             {
                mdeltae = deltae;
                mmbc = mbc;
             }
          }
       }
       
       if(mmbc < 2.25 || fabs(mdeltae) > 0.1) continue;

       m_bc->Fill(mmbc);
       delta_E->Fill(mdeltae);
       
       if(mdeltae < 0.0345291 && mdeltae > -0.0622625) 
       {
          mbc_sig->Fill(mmbc);
          MBC_sig = mmbc;
          MBC_side = 0.;
       }
       else
       {
          mbc_side->Fill(mmbc);
          MBC_side = mmbc;
          MBC_sig = 0.;
       }
       
       MBC = mmbc;
       DE = mdeltae;
       T->Fill();
   }
    
    TFile _1M("Lc.root","RECREATE");
    T->Write();
   _1M.Close();
   /*
    TFile root(rootname,"RECREATE");
    hist.Write();
    root.Close();
*/
   TCanvas *c = new TCanvas("c","c",1400,800);
   c->Divide(2,2);
   c->cd(1); m_bc->Draw();
   c->cd(2); delta_E->Draw();
   c->cd(3); mbc_sig->Draw();
   c->cd(4); mbc_side->Draw();
   
}
