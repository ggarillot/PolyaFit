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

	std::string jsonFile = argv[1] ;

	GraphManager a ;

	a.ProcessData( std::string(jsonFile) ) ;

	a.fitAllEffGraphs() ;
	a.fitAllMulGraphs() ;

	a.writeResultTree("resData.root") ;

	return 0 ;
}
