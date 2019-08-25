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
		return 1 ;
	}

	std::string inputFile = argv[1] ;
	std::string outputFile = std::string("Fit.root") ;

	GraphManager a ;

    a.ProcessSimu( inputFile ) ;

	a.fitAllEffGraphs() ;
	a.fitAllMulGraphs() ;

	a.writeResultTree( outputFile ) ;

	return 0 ;
}
