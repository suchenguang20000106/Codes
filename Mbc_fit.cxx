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
using namespace RooFit ;
//using namespace std;


void cbs_cut()
{
 TTree* tree = (TTree*) gDirectory->Get("ppi0");
 
 RooRealVar mLcppi0("mLcppi0", "mass of L_c(GeV)", 2.25, 2.3);
 mLcppi0.setBins(1000,"range");
 RooDataSet data("data", "data of mass", tree, mLcppi0);
 
 RooPlot* mframe = mLcppi0.frame(Title("cbs for mbc"));
 data.plotOn(mframe,Name("data"));

 RooDataHist *hist = data.binnedClone();
 RooHistPdf histpdf("histpdf", "histpdf", mLcppi0, *hist, 9);
 histpdf.plotOn(mframe,RooFit::Precision(0.00001));

 mframe.Draw();
}