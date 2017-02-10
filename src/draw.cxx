#include <iostream>
#include <sstream>
#include <set>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <TCanvas.h>
#include <TApplication.h>

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

//	a.ProcessData( std::string("/home/garillot/files/PolyaScan/DATA") ) ;

	a.ProcessFile( std::string("/home/garillot/files/PolyaScan/MulResults/map_7_0.5_0.3.root") ) ;

//	a.openGraphs("Data_graphs.root") ;
//	a.openGraphs("/home/garillot/files/PolyaScan/Graphs/20_16_0.2graphs.root") ;

//	a.openGraphs("/home/garillot/Code/PolyaFit/analogTest.root") ;
//	a.openGraphs("/home/garillot/SDHCALMarlinProcessor/analogTest.root") ;
//	a.openGraphs("/home/garillot/Code/PolyaFit/Data_graphs.root") ;
//	a.createGraphs(5 , 2) ;
	TGraphErrors* graph = a.getGraph(layer,dif,asic) ;

	PolyaFitter::PolyaFitResult res = a.fitGraph(layer,dif,asic) ;
	res.print() ;

	TApplication* app = new TApplication("app" , 0 , 0) ;
	TCanvas* c1 = new TCanvas("c1","c1",900,900) ;
	c1->cd() ;

	TF1* fit = new TF1("fit", PolyaFitter::baseFunc , 0 , 30 , 3) ;
	fit->SetParameters(res.qbar , res.delta , res.eff0) ;
	fit->SetNpx(2000) ;

//	TF1* fit = new TF1("fit", PolyaFitter::baseFunc , 0 , 30 , 4) ;
//	fit->SetParameters(res.qbar , res.delta , res.eff0 , 0.967282) ;
//	fit->SetNpx(2000) ;

	std::cout << "Minimizer status : " << res.minimStatus << std::endl ;

	graph->Draw("AP") ;
	fit->Draw("same") ;
	c1->SetLogx() ;
	c1->SaveAs("draw.root") ;
	app->Run() ;


	return 0 ;
}
