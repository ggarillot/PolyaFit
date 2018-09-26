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

#include "PolyaFitter.h"
#include "GraphManager.h"


int main(int argc , char** argv)
{
//	if ( argc != 3 )
//	{
//		std::cerr << "ERROR : problem with arguments passed for the program" << std::endl ;
//		return -1 ;
//	}

//	double qbar = atof(argv[1]) ;
//	double delta = atof(argv[2]) ;

//	GraphManager a ;
//	GraphManager b ;

//	a.createGraphsData() ;
//	a.writeGraphsInFile("Data_graphs.root") ;

//	a.openGraphs("/home/garillot/Code/PolyaFit/analogTest.root") ;
//	b.openGraphs("/home/garillot/Code/PolyaFit/Data_graphs.root") ;

//	a.openGraphs("/home/garillot/SDHCALMarlinProcessor/analogTest.root") ;

//	a.fitAllGraphs() ;
//	b.fitAllGraphs() ;
//	a.writeResultTree("resSim.root") ;
//	a.writeResultTree("resData.root") ;


	GraphManager a ;

//	a.ProcessData( std::string("/home/garillot/files/PolyaScan/DATA") ) ;
	a.ProcessData( std::string("/home/garillot/files/MultiplicityMap/DATA/thrScan") ) ;
//	a.ProcessFile( std::string("/home/garillot/files/PolyaScan/MulResults/20_16_0.2.root") ) ;

	a.fitAllGraphs() ;

//	a.fitMulGraph(0,94,12) ;
	a.fitAllMulGraphs() ;

//	a.writeResultTree("resSim.root") ;
	a.writeResultTree("resData.root") ;


	return 0 ;
}
