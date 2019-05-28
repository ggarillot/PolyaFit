#include <iostream>

#include <TFile.h>
#include <TDirectory.h>
#include <TCanvas.h>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TPaveText.h>
#include <TLatex.h>

#include <sstream>

#include "GraphManager.h"
#include "PolyaFitter.h"

int main()
{
	gStyle->SetOptStat(0) ;

	GraphManager a ;

	//	a.ProcessData( std::string("/home/garillot/files/PolyaScan/DATA") ) ;
	a.ProcessData( std::string("/home/guillaume/PolyaFit/json/SPS_Sept2018.json") ) ;
	//	a.ProcessData( std::string("/home/guillaume/PolyaFit/json/SPS_Oct2015.json") ) ;
	//	b.ProcessFile( std::string("/home/garillot/files/TUNINGTEMP/Eff_cosmicMuons_Analog.root") ) ;
	//	b.ProcessFile( std::string("/home/garillot/SDHCALMarlinProcessor/Eff_mu-_50GeV_Analog.root") ) ;

	auto aRes = a.fitAllGraphs() ;
	auto aMulRes = a.fitAllMulGraphs() ;

	std::map<GraphManager::AsicID , TCanvas*> canvasMap ;
	std::map<GraphManager::AsicID , TCanvas*> canvasMulMap ;

	//data
	for ( const auto& it : aRes )
	{
		if ( canvasMap.find(it.first) == canvasMap.end() )
		{
			std::stringstream toto ; toto << it.first ;
			TCanvas* c1 = new TCanvas(toto.str().c_str() , toto.str().c_str() , 900 , 900) ;
			c1->cd() ;
			c1->SetTopMargin(0.05f) ;
			c1->SetRightMargin(0.05f) ;
			c1->SetTicks() ;
			c1->SetLogx() ;

			std::stringstream sstext ;
			sstext << it.first << ";Threshold (pC);Efficiency";
			std::stringstream sstext2 ;
			sstext2 << ";Threshold (pC);Efficiency";
			TH2D* range = new TH2D(sstext.str().c_str() , sstext2.str().c_str() , 1 , 0.1 , 27 , 1 , 0 , 1.1) ;
			range->GetXaxis()->SetTitleFont(62) ;
			range->GetXaxis()->SetTitleOffset(1.1f) ;
			range->GetYaxis()->SetTitleFont(62) ;
			range->Draw() ;

			sstext.str("") ;
			sstext << it.first ;
			TLatex* t = new TLatex(.525 , .96 , sstext.str().c_str()) ;
			t->SetNDC() ;
			t->SetTextAlign(21) ;
			t->SetTextSize(0.03f) ;
			t->Draw() ;

			canvasMap.insert( { it.first , c1 } ) ;
		}

		auto jt = canvasMap.find(it.first) ;
		jt->second->cd() ;

		auto graph = a.getGraph(it.first) ;
		graph->SetMarkerStyle(20) ;
		graph->Draw("P same") ;

		TF1* fit = new TF1("fit", PolyaFitter::polyaVsThr , 0 , 30 , 3) ;
		fit->SetParameters(it.second.qbar , it.second.delta , it.second.eff0) ;
		fit->SetNpx(2000) ;
		fit->SetLineColor(kBlack) ;
		fit->SetLineStyle(2) ;
		fit->Draw("same") ;

		TPaveText* pt = new TPaveText(0.16 , 0.20 , 0.44 , 0.33 , "NDC") ;
		pt->SetBorderSize(0) ;
		pt->SetFillColor(kWhite) ;
		pt->SetTextAlign(11) ;

		std::stringstream titi ;
		titi.precision(4) ;
		titi.width(3) ;
		titi << "#bar{q} = " << it.second.qbar << " #pm " << std::round(it.second.qbarError*1000)/1000 << " pC" ;
		pt->AddText( titi.str().c_str() ) ;

		titi.str("") ;
		titi << "#delta = " << it.second.delta << " #pm " << std::round(it.second.deltaError*1000)/1000 << " pC^{-1}" ;
		pt->AddText( titi.str().c_str() ) ;

		titi.str("") ;
		titi << "eff0 = " << it.second.eff0 << " #pm " << std::round(it.second.eff0Error*1000)/1000+0.0001 ;
		pt->AddText( titi.str().c_str() ) ;

		pt->Draw() ;
	}

	for ( const auto& it : aMulRes )
	{
		if ( canvasMulMap.find(it.first) == canvasMulMap.end() )
		{
			std::stringstream toto ; toto << it.first << "_Mul";
			TCanvas* c1 = new TCanvas(toto.str().c_str() , toto.str().c_str() , 900 , 900) ;
			c1->SetTopMargin(0.05f) ;
			c1->SetRightMargin(0.05f) ;
			c1->SetTicks() ;

			std::stringstream sstext ;
			sstext << it.first << ";Threshold (pC);Multiplicity";
			std::stringstream sstext2 ;
			sstext2 << ";Threshold (pC);Multiplicity";
			TH2D* range = new TH2D(sstext.str().c_str() , sstext2.str().c_str() , 1 , 0.1 , 0.549 , 1 , 1.01 , 3.1) ;
			range->GetXaxis()->SetTitleFont(62) ;
			range->GetXaxis()->SetTitleOffset(1.1f) ;
			range->GetYaxis()->SetTitleFont(62) ;
			range->Draw() ;

			sstext.str("") ;
			sstext << it.first ;
			TLatex* t = new TLatex(.525 , .96 , sstext.str().c_str()) ;
			t->SetNDC() ;
			t->SetTextAlign(21) ;
			t->SetTextSize(0.03f) ;
			t->Draw() ;

			canvasMulMap.insert( { it.first , c1 } ) ;
		}

		auto jt = canvasMulMap.find(it.first) ;
		jt->second->cd() ;

		auto graph = a.getMulGraph(it.first) ;
		graph->SetMarkerStyle(20) ;
		graph->Draw("P same") ;

		TF1* fit = new TF1("fit", MultiplicityFitter::multiplicityVsThr , 0 , 30.6 , 3) ;
		fit->SetParameters(it.second.f , it.second.p , it.second.c) ;
		fit->SetNpx(2000) ;
		fit->SetLineColor(kBlack) ;
		fit->SetLineStyle(2) ;
		fit->Draw("same") ;

		TPaveText* pt = new TPaveText(0.6 , 0.75 , 0.88 , 0.88 , "NDC") ;
		pt->SetBorderSize(0) ;
		pt->SetFillColor(kWhite) ;
		pt->SetTextAlign(11) ;

		std::stringstream titi ;
		titi.precision(4) ;
		titi.width(3) ;
		titi << "factor = " << it.second.f << " pC^{-1}" ;
		pt->AddText( titi.str().c_str() ) ;

		titi.str("") ;
		titi << "power = " << it.second.p ;
		pt->AddText( titi.str().c_str() ) ;

		titi.str("") ;
		titi << "constant = " << it.second.c ;
		pt->AddText( titi.str().c_str() ) ;

		pt->Draw() ;
	}



	TFile* file = new TFile("draw2018.root" , "RECREATE") ;
	file->cd() ;

	TDirectoryFile* dir = new TDirectoryFile("dir","dir") ;
	TDirectoryFile* dirMul = new TDirectoryFile("dirMul","dirMul") ;

	dir->cd() ;

	std::map<int,TDirectoryFile*> layerMap = {} ;

	for ( const auto& it : canvasMap )
	{
		dir->cd() ;
		if ( !layerMap.count(it.first.layerID) && it.first.layerID != -1 )
		{
			std::stringstream toto ; toto << "Layer" << it.first.layerID ;
			TDirectoryFile* dira = new TDirectoryFile(toto.str().c_str(),toto.str().c_str()) ;
			layerMap.insert({it.first.layerID , dira}) ;
		}

		if ( it.first.layerID != -1 )
			layerMap.at(it.first.layerID)->cd() ;

		std::cout << "Write canvas " << it.first.layerID << "," << it.first.difID << "," << it.first.asicID << std::endl ;
		it.second->Write() ;
	}

	dirMul->cd() ;

	layerMap.clear() ;

	for ( const auto& it : canvasMulMap )
	{
		dirMul->cd() ;
		if ( !layerMap.count(it.first.layerID) && it.first.layerID != -1 )
		{
			std::stringstream toto ; toto << "Layer" << it.first.layerID ;
			TDirectoryFile* dira = new TDirectoryFile(toto.str().c_str(),toto.str().c_str()) ;
			layerMap.insert({it.first.layerID , dira}) ;
		}

		if ( it.first.layerID != -1 )
			layerMap.at(it.first.layerID)->cd() ;

		std::cout << "Write canvas " << it.first.layerID << "," << it.first.difID << "," << it.first.asicID << std::endl ;
		it.second->Write() ;
	}

	file->Close() ;

	return 0 ;
}
