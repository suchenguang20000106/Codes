#include "TH1.h"
#include "TStyle.h"
#include "TCanvas.h"
#include "TFile.h"
#include <vector>
#include <iostream>

void background(char hisname[])
{   
    TCanvas *bkg = new TCanvas("bkg","background",1500,800);
    bkg->Divide(3,2);

    bkg->cd(1);
    TFile* DD = new TFile("DD.root");
    TH1F *dd = (TH1F*)DD->Get(hisname);
    dd->SetFillColor(3);
    dd->Scale(7.26);
    dd->Rebin(4);
    dd->Draw();
    
    bkg->cd(2);
    TFile* ditau = new TFile("ditau.root");
    TH1F *Ditau = (TH1F*)ditau->Get(hisname);
    Ditau->SetFillColor(5);
    Ditau->Scale(3.4);
    Ditau->Rebin(4);
    Ditau->Draw();

    bkg->cd(3);
    TFile* DsDss = new TFile("DsDss.root");
    TH1F *dsdss = (TH1F*)DsDss->Get(hisname);
    dsdss->SetFillColor(0);
    dsdss->Scale(10);
    dsdss->Rebin(4);
    dsdss->Draw();

    bkg->cd(4);
    TFile* ISR = new TFile("ISR.root");
    TH1F *isr = (TH1F*)ISR->Get(hisname);
    isr->SetFillColor(9);
    isr->Scale(0.67);
    isr->Rebin(4);
    isr->Draw();

    bkg->cd(5);
    TFile* qqbar1M = new TFile("qqbar1M.root");
    TH1F *Qqbar = (TH1F*)qqbar1M->Get(hisname);
    Qqbar->SetFillColor(4);
    Qqbar->Scale(1.24);
    Qqbar->Rebin(4);
    Qqbar->Draw();
    
    bkg->cd(6);
    TFile* Lcpair = new TFile("Lcpair.root");
    TH1F *lcpair = (TH1F*)Lcpair->Get(hisname);
    lcpair->SetFillColor(2);
    lcpair->Scale(0.237);
    lcpair->Rebin(4);
    lcpair->Draw();


/*    lcpair->GetXaxis()->SetTitle("Mbc(GeV/c^{2}");
    lcpair->GetYaxis()->SetTitle("Events");
*/
    TCanvas *c = new TCanvas("c","c");
    THStack *bcg = new THStack("bcg","background");
    bcg->Add(dd);
    bcg->Add(Ditau);
    bcg->Add(dsdss);
    bcg->Add(isr);
    bcg->Add(Qqbar);
    bcg->Add(lcpair);
    bcg->Draw();
}