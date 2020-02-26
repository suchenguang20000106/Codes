#ifndef __CINT__
#include "RooGlobalFunc.h"
#endif

#include "RooRealVar.h"
#include "RooDataSet.h"
#include "RooGaussian.h"
#include "RooFFTConvPdf.h"
#include "TCanvas.h"
#include "RooKeysPdf.h"
#include "RooNDKeysPdf.h"
#include "RooProdPdf.h"
#include "RooPlot.h"
#include "TAxis.h"
#include "RooConstVar.h
#include "RooArgusBG.h"
using namespace RooFit ;

void argus()
{
    TFile *file = new TFile("1M.root");
    TTree *qqbar = (TTree*)file -> Get("Lcpair");



    RooRealVar Mbc("Mbc", "Mbc", 2.25, 2.3);
    RooPlot *Mbcframe = Mbc.frame(Title("Argus fit for all"));

    RooRealVar Mbc_side("Mbc_side", "Mbc_side", 2.25, 2.3);
    RooPlot *Mbc_sideframe = Mbc_side.frame(Title("Argus fit for sideband region"));

    RooRealVar Mbc_sig("Mbc_sig", "Mbc_sig", 2.25, 2.3);
    RooPlot *Mbc_sigframe = Mbc_sig.frame(Title("Argus fit for sigal region"));



    RooDataSet total("total", "total", qqbar, Mbc);
    total.plotOn(Mbcframe);

    RooDataSet sideband("sideband", "sideband", qqbar, Mbc_side);
    sideband.plotOn(Mbc_sideframe, MarkerStyle(24));

    RooDataSet signal("signal", "signal", qqbar, Mbc_sig);
    signal.plotOn(Mbc_sigframe);
    RooConstVar E0("E0", "E0", 2.3);



    //fit total
    RooRealVar c("c", "c", -10., -100., -1.);
    RooRealVar p("p", "p", 0.5, 0., 10.);
    RooArgusBG bg_total("bg_total", "bg_total", Mbc, E0, c, p);
    bg_total.fitTo(total);
    bg_total.plotOn(Mbcframe);

    //fit sideband region
    RooRealVar c_side("c_side", "c_side", -10., -100., -1.);
    RooRealVar p_side("p_side", "p_side", 0.5, 0., 10.);
    RooArgusBG bg_side("bg_side", "bg_side", Mbc_side, E0, c_side, p_side);
    bg_side.fitTo(sideband);
    bg_side.plotOn(Mbc_sideframe, LineColor(kRed));

    //fit signal region
    RooRealVar c_sig("c_sig", "c_sig", -10., -100., -1.);
    RooRealVar p_sig("p_sig", "p_sig", 0.5, 0., 10.);
    RooArgusBG bg_sig("bg_sig", "bg_sig", Mbc_sig, E0, c_sig, p_sig);
    bg_sig.fitTo(signal);
    bg_sig.plotOn(Mbc_sigframe);

    TCanvas* c1 = new TCanvas("c1", "c1", 1200, 1000);
    c1 -> Divide(2,2);
    c1 -> cd(1); Mbcframe -> DrawClone();
    c1 -> cd(2); Mbc_sideframe -> DrawClone();
    c1 -> cd(3); Mbc_sigframe -> DrawClone();
    c1 -> cd(4); Mbc_sigframe -> DrawClone(); Mbc_sideframe -> DrawClone("same");
    
}