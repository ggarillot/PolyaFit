#include <iostream>
#include <sstream>
#include <set>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TText.h>
#include <TLatex.h>
#include <TPaveText.h>
#include <TStyle.h>

#include "PolyaFitter.h"
#include "GraphManager.h"

int main(int argc , char** argv)
{
	gStyle->SetOptStat(0) ;

	if ( argc != 4 )
	{
		std::cerr << "ERROR : problem with arguments passed for the program" << std::endl ;
		return -1 ;
	}

	int layer = atoi(argv[1]) ;
	int dif = atoi(argv[2]) ;
	int asic = atoi(argv[3]) ;

	GraphManager::AsicID asicID{layer,dif,asic} ;

	GraphManager a ;
	GraphManager b ;

	a.ProcessData( std::string("/home/guillaume/PolyaFit/json/SPS_Oct2015.json") ) ;
    b.ProcessSimu( std::string("/home/guillaume/files/MultiplicityMap/Cosmic/Geant4.10.04/QGSP_BERT/Eff_mu-_50GeV_Analog.root") ) ;

	auto resData = a.fitEffGraph(asicID) ;
	auto resSim = b.fitEffGraph(asicID) ;

	auto graphData = a.getGraph(asicID) ;
	auto graphSim = b.getGraph(asicID) ;



	std::stringstream toto ; toto << asicID ;
	TCanvas* canvas = new TCanvas(toto.str().c_str() , toto.str().c_str() , 900 , 900) ;
	canvas->cd() ;
	canvas->SetTopMargin(0.05f) ;
	canvas->SetRightMargin(0.05f) ;
	canvas->SetTicks() ;
	canvas->SetLogx() ;


	std::stringstream sstext ;
	sstext << asicID << ";Threshold (pC);Efficiency";
	std::stringstream sstext2 ;
	sstext2 << ";Threshold (pC);Efficiency";
	TH2D* range = new TH2D(sstext.str().c_str() , sstext2.str().c_str() , 1 , 0.1 , 27 , 1 , 0 , 1.1) ;
	range->GetXaxis()->SetTitleFont(62) ;
	range->GetXaxis()->SetTitleOffset(1.1f) ;
	range->GetYaxis()->SetTitleFont(62) ;
	range->Draw() ;

	sstext.str("") ;
	sstext << asicID ;
	TLatex* t = new TLatex(.525 , .96 , sstext.str().c_str()) ;
	t->SetNDC() ;
	t->SetTextAlign(21) ;
	t->SetTextSize(0.03f) ;
	t->Draw() ;

	graphData->SetMarkerStyle(20) ;
	graphData->Draw("P same") ;

	TF1* fit = new TF1("fit", PolyaFitter::polyaVsThr , 0 , 30 , 3) ;
	fit->SetParameters(resData.qbar , resData.delta , resData.eff0) ;
	fit->SetNpx(2000) ;
	fit->SetLineColor(kBlack) ;
	fit->SetLineStyle(2) ;
	fit->Draw("same") ;


	graphSim->SetMarkerStyle(20) ;
	graphSim->SetMarkerColor(kRed-4) ;
	graphSim->SetLineColor(kRed-4) ;
	graphSim->Draw("P same") ;

	TF1* fitSim = new TF1("fitSim", PolyaFitter::polyaVsThr , 0 , 30 , 3) ;
	fitSim->SetParameters(resSim.qbar , resSim.delta , resSim.eff0) ;
	fitSim->SetNpx(2000) ;
	fitSim->SetLineColor(kRed-4) ;
	fitSim->SetLineStyle(2) ;
	fitSim->Draw("same") ;


	TPaveText* pt = new TPaveText(0.16 , 0.32 , 0.44 , 0.47 , "NDC") ;
	pt->SetBorderSize(0) ;
	pt->SetFillColor(kWhite) ;
	pt->SetTextAlign(11) ;

	pt->AddText("SPS_Oct2015") ;

	std::stringstream titi ;
	titi.precision(4) ;
	titi.width(3) ;
	titi << "#bar{q} = " << resData.qbar << " #pm " << std::round(resData.qbarError*1000)/1000 << " pC" ;
	pt->AddText( titi.str().c_str() ) ;

	titi.str("") ;
	titi << "#delta = " << resData.delta << " #pm " << std::round(resData.deltaError*1000)/1000 << " pC^{-1}" ;
	pt->AddText( titi.str().c_str() ) ;

	titi.str("") ;
	titi << "eff0 = " << resData.eff0 << " #pm " << std::round(resData.eff0Error*1000)/1000+0.0001 ;
	pt->AddText( titi.str().c_str() ) ;

	pt->Draw() ;


	TPaveText* pt2 = new TPaveText(0.16 , 0.15 , 0.44 , 0.30 , "NDC") ;
	pt2->SetTextColor(kRed-4) ;
	pt2->SetBorderSize(0) ;
	pt2->SetFillColor(kWhite) ;
	pt2->SetTextAlign(11) ;

	pt2->AddText("Simulation") ;

	titi.str("") ;
	titi.precision(4) ;
	titi.width(3) ;
	titi << "#bar{q} = " << resSim.qbar << " #pm " << std::round(resSim.qbarError*1000)/1000 << " pC" ;
	pt2->AddText( titi.str().c_str() ) ;

	titi.str("") ;
	titi << "#delta = " << resSim.delta << " #pm " << std::round(resSim.deltaError*1000)/1000 << " pC^{-1}" ;
	pt2->AddText( titi.str().c_str() ) ;

	titi.str("") ;
	titi << "eff0 = " << resSim.eff0 << " #pm " << std::round(resSim.eff0Error*1000)/1000+0.0001 ;
	pt2->AddText( titi.str().c_str() ) ;


	pt2->Draw() ;


	canvas->SaveAs("polyaDataSim.root") ;
}
