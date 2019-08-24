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
	if ( argc != 2 )
	{
		std::cerr << "ERROR : problem with arguments passed for the program" << std::endl ;
		return -1 ;
	}

	std::string inputFile = argv[1] ;
	std::string outputFile = std::string("Fit.root") ;

//	double qbar = atof(argv[1]) ;
//	double delta = atof(argv[2]) ;
//	double d = atof(argv[3]) ;

	GraphManager a ;

//	std::stringstream filePath ; filePath << "/home/garillot/files/PolyaScan/MulResults/map_" << qbar << "_" << delta << "_" << d << ".root" ;
//	a.ProcessFile( filePath.str() ) ;
    a.ProcessSimu( inputFile ) ;

	a.fitAllEffGraphs() ;
	a.fitAllMulGraphs() ;

	std::cout << "write tree" << std::endl ;
//	std::stringstream outFilePath ; outFilePath << "/home/garillot/files/PolyaScan/Fits/Fits_" << qbar << "_" << delta << "_" << d << ".root" ;
//	a.writeResultTree( outFilePath.str() ) ;
	a.writeResultTree( outputFile ) ;

	return 0 ;
}
