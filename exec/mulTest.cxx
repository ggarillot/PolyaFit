#include <iostream>
#include <sstream>
#include <set>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TGraphAsymmErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>

#include "PolyaFitter.h"
#include "MultiplicityFitter.h"
#include "GraphManager.h"


int main(int argc , char** argv)
{
	int runstmp[] = { 730630, 730627, 730626, 730625, 730619, 730618,
					  730617, 730616, 730615, 730611, 730609, 730607,
					  730569, 730568, 730567, 730566, 730631, 730633, 730545, 730677 } ;
	int dac0tmp[] = { 188,199,210,221,232,243,254,265,276,287,299,310,321,332,343,354,365,376,387,170 } ;
	int dac1tmp[] = { 130,147,164,181,197,214,231,248,265,282,298,315,332,349,366,383,399,416,433,498 } ;
	int dac2tmp[] = { 168,185,202,220,237,254,271,288,305,323,340,357,374,391,408,425,443,460,477,342 } ;

	std::vector<int> runs ;
	std::vector<double> thr1Vec ;
	std::vector<double> thr2Vec ;
	std::vector<double> thr3Vec ;

	for ( unsigned int i = 0 ; i < 20 ; i++ )
	{
		runs.push_back(runstmp[i]) ;
		thr1Vec.push_back( (dac0tmp[i]-90)/700.0 ) ;
		thr2Vec.push_back( (dac1tmp[i]-98)/80.0 ) ;
		thr3Vec.push_back( (dac2tmp[i]-98)/16.3 ) ;
	}

	TGraphErrors* graph = new TGraphErrors ;

	for ( unsigned int iRun = 0 ; iRun < runs.size() ; ++iRun )
	{
		std::vector<double> thresholds ;
		thresholds.push_back( thr1Vec.at(iRun) ) ;
		thresholds.push_back( thr2Vec.at(iRun) ) ;
		thresholds.push_back( thr3Vec.at(iRun) ) ;

		std::stringstream filePath ;
		filePath << std::string("/home/garillot/files/MultiplicityMap/DATA/thrScan") << "/Eff_" << runs.at(iRun) << ".root" ;

		std::cout << "Process " << filePath.str() << std::endl ;
		TFile* file = new TFile( filePath.str().c_str() , "READ") ;
		TTree* tree = dynamic_cast<TTree*>( file->Get("tree") ) ;
		if ( !tree )
		{
			std::cout << "Error in ProcessData : tree not present in " << filePath.str() << std::endl ;
			file->Close() ;
			//			return ;
			continue ;
		}

		std::vector<double>* efficiencies = nullptr ;
		std::vector<double>* efficienciesLowerBound = nullptr ;
		std::vector<double>* efficienciesUpperBound = nullptr ;


		int difID , asicID , layerID , padID ;
		std::vector<double>* multiplicities = nullptr ;
		std::vector<double>* multiplicitiesError = nullptr ;
		std::vector<double>* position = nullptr ;

		int nTrack ;
		int globalNTrack = 0 ;

		tree->SetBranchAddress("LayerID" , &layerID) ;
		tree->SetBranchAddress("DifID" , &difID) ;
		tree->SetBranchAddress("AsicID" , &asicID) ;
		tree->SetBranchAddress("PadID" , &padID) ;
		tree->SetBranchAddress("Multiplicities" , &multiplicities) ;
		tree->SetBranchAddress("MultiplicitiesError" , &multiplicitiesError) ;
		tree->SetBranchAddress("Position" , &position) ;
		tree->SetBranchAddress("Ntrack" , &nTrack) ;

		int iEntry = 0 ;
		while ( tree->GetEntry(iEntry++) )
		{
			if ( !(layerID == 47 && difID == 1 && asicID == 48 && padID == -1) ) // because padID > -1 means stats for an individual pad
				continue ;


			for ( unsigned int i = 0 ; i < 3 ; ++i )
			{
				if ( i > 0 )
					continue ;
				std::cout << "thr : " << thresholds.at(i) << "  " << multiplicities->at(i) << "   " << multiplicitiesError->at(i) << std::endl ;
				if ( multiplicities->at(i) < std::numeric_limits<double>::epsilon() )
					continue ;
				graph->SetPoint(graph->GetN() , thresholds.at(i) , multiplicities->at(i) ) ;
				graph->SetPointError(graph->GetN()-1 , 0 , multiplicitiesError->at(i) ) ;
			}

			//			std::cout << "thr : " << thresholds.at(0) << "  " << multiplicities->at(0) << "   " << multiplicitiesError->at(0) << std::endl ;

			//			graph->SetPoint(graph->GetN() , thresholds.at(0) , multiplicities->at(0) ) ;
			//			graph->SetPointError(graph->GetN()-1 , 0 , multiplicitiesError->at(0) ) ;
		}
	}

	MultiplicityFitter a ;
	a.getPoints(graph) ;

	a.setParams() ;
	a.minimize() ;

	auto res = a.getFitResult() ;
	res.print() ;


	TCanvas* c1 = new TCanvas("c1" , "c1" , 900 , 900) ;
	c1->cd() ;
	graph->SetMarkerStyle(20) ;
	graph->Draw("AP") ;

	TF1* fit = new TF1("fit", MultiplicityFitter::baseFunc , 0 , 30 , 3) ;
	fit->SetParameters(res.f , res.p , res.c) ;
	fit->SetNpx(2000) ;
	fit->SetLineColor(kBlack) ;
	fit->SetLineStyle(2) ;
	fit->Draw("same") ;

	c1->SaveAs("test.root") ;
	/*
	GraphManager a ;

//	a.ProcessData( std::string("/home/garillot/files/PolyaScan/DATA") ) ;
	a.ProcessData( std::string("/home/garillot/files/MultiplicityMap/DATA/thrScan") ) ;
//	a.ProcessFile( std::string("/home/garillot/files/PolyaScan/MulResults/20_16_0.2.root") ) ;

	a.fitAllGraphs() ;

//	a.writeResultTree("resSim.root") ;
	a.writeResultTree("resData.root") ;
*/

	return 0 ;
}
