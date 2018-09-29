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
	if ( argc != 4 )
	{
		std::cerr << "ERROR : problem with arguments passed for the program" << std::endl ;
		return -1 ;
	}

	int layer = atoi(argv[1]) ;
	int dif = atoi(argv[2]) ;
	int asic = atoi(argv[3]) ;

	GraphManager a ;
	GraphManager b ;

//	a.ProcessData( std::string("/home/garillot/files/MultiplicityMap/DATA/thrScan") ) ;
	a.ProcessData( std::string("/home/garillot/Code/PolyaFit/json/SPS_Oct2015.json") ) ;
	b.ProcessFile( std::string("/home/garillot/SDHCALMarlinProcessor/script/totoLog.root") ) ;
//	b.ProcessFile( std::string("/home/garillot/files/TUNINGTEMP/Eff_cosmicMuons50GeV_Analog.root") ) ;
//	a.ProcessFile( std::string("/home/garillot/files/PolyaScan/MulResults/map_7_0.5_0.3.root") ) ;

//	a.openGraphs("Data_graphs.root") ;
//	a.openGraphs("/home/garillot/files/PolyaScan/Graphs/20_16_0.2graphs.root") ;

//	a.openGraphs("/home/garillot/Code/PolyaFit/analogTest.root") ;
//	a.openGraphs("/home/garillot/SDHCALMarlinProcessor/analogTest.root") ;
//	a.openGraphs("/home/garillot/Code/PolyaFit/Data_graphs.root") ;
//	a.createGraphs(5 , 2) ;
	TGraphAsymmErrors* graphA = a.getGraph(layer,dif,asic) ;
	TGraphAsymmErrors* graphB = b.getGraph(layer,dif,asic) ;

	PolyaFitter::PolyaFitResult resA = a.fitGraph(layer,dif,asic) ;
	resA.print() ;
	PolyaFitter::PolyaFitResult resB = b.fitGraph(layer,dif,asic) ;
	resB.print() ;

	TApplication* app = new TApplication("app" , 0 , 0) ;
	TCanvas* c1 = new TCanvas("c1","c1",900,900) ;
	c1->cd() ;
	gStyle->SetOptStat(0) ;

	TF1* fit = new TF1("fit", PolyaFitter::baseFunc , 0 , 30 , 3) ;
	fit->SetParameters(resA.qbar , resA.delta , resA.eff0) ;
	fit->SetNpx(2000) ;

//	TF1* fit = new TF1("fit", PolyaFitter::baseFunc , 0 , 30 , 4) ;
//	fit->SetParameters(res.qbar , res.delta , res.eff0 , 0.967282) ;
//	fit->SetNpx(2000) ;

	std::cout << "MinimizerA status : " << resA.minimStatus << std::endl ;
	std::cout << "MinimizerB status : " << resB.minimStatus << std::endl ;

	std::stringstream sstext ;
	sstext << "Layer " << layer << ", Dif " << dif << ", Asic " << asic << ";Threshold (pC);Efficiency";
	TH2D* range = new TH2D("range" , sstext.str().c_str() , 1 , 0.1 , 27 , 1 , 0 , 1.1) ;
	range->Draw() ;

	graphA->SetMarkerStyle(20) ;
	graphA->Draw("P same") ;
	graphB->SetMarkerColor(kRed-4) ;
	graphB->SetLineColor(kRed-4) ;
	graphB->Draw("P same") ;

	fit->Draw("same") ;



	TPaveText* pt = new TPaveText(0.6 , 0.75 , 0.88 , 0.88) ;
	pt->ConvertNDCtoPad() ;
	pt->SetBorderSize(0) ;
	pt->SetFillColor(kWhite) ;
	pt->SetTextAlign(11) ;

	std::stringstream titi ;
	titi.precision(4) ;
	titi.width(3) ;
	titi << "#bar{q} = " << resA.qbar << " #pm " << resA.qbarError << " pC" ;
	pt->AddText( titi.str().c_str() ) ;

	titi.str("") ;
	titi << "#delta = " << resA.delta << " #pm " << resA.deltaError << " pC^{-1}" ;
	pt->AddText( titi.str().c_str() ) ;

	titi.str("") ;
	titi << "eff0 = " << resA.eff0 << " #pm " << resA.eff0Error ;
	pt->AddText( titi.str().c_str() ) ;


	pt->Draw() ;

	TPaveText* pt2 = new TPaveText(0.6 , 0.61 , 0.88 , 0.74) ;
	pt2->ConvertNDCtoPad() ;
	pt2->SetTextColor(kRed-4) ;
	pt2->SetBorderSize(0) ;
	pt2->SetFillColor(kWhite) ;
	pt2->SetTextAlign(11) ;

	std::stringstream titi2 ;
	titi2.precision(4) ;
	titi2.width(3) ;
	titi2 << "#bar{q} = " << resB.qbar << " #pm " << resB.qbarError << " pC" ;
	pt2->AddText( titi2.str().c_str() ) ;

	titi2.str("") ;
	titi2 << "#delta = " << resB.delta << " #pm " << resB.deltaError << " pC^{-1}" ;
	pt2->AddText( titi2.str().c_str() ) ;

	titi2.str("") ;
	titi2 << "eff0 = " << resB.eff0 << " #pm " << resB.eff0Error ;
	pt2->AddText( titi2.str().c_str() ) ;


	pt2->Draw() ;


	c1->SetLogx() ;
	c1->SaveAs("draw.png") ;
	app->Run() ;


	return 0 ;
}
