#include <iostream>
#include <sstream>
#include <fstream>
#include <set>
#include <ctime>

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <TF1.h>
#include <TH1D.h>
#include <TH2.h>
#include <TCanvas.h>

#include "PolyaFitter.h"
#include "GraphManager.h"


int main(int argc , char** argv)
{
	if ( argc != 4 )
	{
		std::cerr << "ERROR : problem with arguments passed for the program" << std::endl ;
		return -1 ;
	}

	double mul_ = atof(argv[1]) ;
	double eff2_ = atof(argv[2]) ;
	double eff3_ = atof(argv[3]) ;


	TFile* file = new TFile("/home/garillot/Code/PolyaFit/resData.root") ;
	TTree* tree = dynamic_cast<TTree*>(file->Get("tree")) ;

	int layerID ;
	int difID ;
	int asicID ;

	double f = 0 ;
	double p = 0 ;
	double c = 0 ;
	double qbar = 0 ;
	double delta = 0 ;
	double eff0 = 0 ;

	int minimStatus = 0 ;

	tree->SetBranchAddress("factor" , &f) ;
	tree->SetBranchAddress("power" , &p) ;
	tree->SetBranchAddress("constant" , &c) ;
	tree->SetBranchAddress("qbar" , &qbar) ;
	tree->SetBranchAddress("delta" , &delta) ;
	tree->SetBranchAddress("eff0" , &eff0) ;

	tree->SetBranchAddress("LayerID" , &layerID) ;
	tree->SetBranchAddress("DifID" , &difID) ;
	tree->SetBranchAddress("AsicID" , &asicID) ;

	tree->SetBranchAddress("minimStatus" , &minimStatus) ;


	TFile* fileOut = new TFile("/home/garillot/Code/PolyaFit/targetThr.root" , "RECREATE") ;
	TTree* treeOut = new TTree("tree" , "tree") ;
	treeOut->SetDirectory(0) ;

	double thr1 = 0 ;
	double thr2 = 0 ;
	double thr3 = 0 ;
	int dac1 = 0 ;
	int dac2 = 0 ;
	int dac3 = 0 ;

	double eff1 = 0 ;


	treeOut->Branch("LayerID" , &layerID) ;
	treeOut->Branch("DifID" , &difID) ;
	treeOut->Branch("AsicID" , &asicID) ;

	treeOut->Branch("thr1" , &thr1) ;
	treeOut->Branch("thr2" , &thr2) ;
	treeOut->Branch("thr3" , &thr3) ;

	treeOut->Branch("dac1" , &dac1) ;
	treeOut->Branch("dac2" , &dac2) ;
	treeOut->Branch("dac3" , &dac3) ;

	treeOut->Branch("eff1" , &eff1) ;

	int iEntry = 0 ;
	while ( tree->GetEntry(iEntry++) )
	{
		double mul = mul_ ;
		double eff2 = eff2_ ;
		double eff3 = eff3_ ;

		if ( layerID == 9 )
		{
			mul = 0.8*mul_ ;
			eff2 = 0.8*eff2_ ;
			eff3 = 0.19*eff3_ ;
		}

		if ( asicID == -1 )
			continue ;
		if ( minimStatus > 0 )
			continue ;

		thr1 = std::exp( std::log( (mul-c)/f )/p) ;

		dac1 = static_cast<int>( thr1*700 + 90 ) ;
		dac1 = std::min(dac1,440) ;
		dac1 = std::max(dac1,200) ;

		double params[] = {qbar,delta,eff0} ;

		eff1 = PolyaFitter::polyaVsThr(thr1 , params) ;

		double a = 0 ;
		double b = 30 ;
		thr2 = 0.5*(b-a) ;

		while ( b-a > 0.01 )
		{
			double effCurrent = PolyaFitter::polyaVsThr(thr2 , params) ;

			if ( effCurrent > eff2 )
				a = thr2 ;
			else
				b = thr2 ;

			thr2 = 0.5*(b-a) + a ;
		}

		dac2 = static_cast<int>( thr2*80 + 98 ) ;
		dac2 = std::min(dac2 , 500) ;
		dac2 = std::max(dac2 , 130) ;

		a = 0 ;
		b = 40 ;
		thr3 = 0.5*(b-a) ;
		while ( b-a > 0.01 )
		{
			double effCurrent = PolyaFitter::polyaVsThr(thr3 , params) ;

			if ( effCurrent > eff3 )
				a = thr3 ;
			else
				b = thr3 ;

			thr3 = 0.5*(b-a) + a ;
		}

		dac3 = static_cast<int>( thr3*16.3 + 98 ) ;

		dac3 = std::min(dac3 , 700) ;
		dac3 = std::max(dac3 , 130) ;

		treeOut->Fill() ;
	}

	file->Close() ;
	fileOut->cd() ;
	treeOut->Write() ;
	fileOut->Close() ;


	TFile* filer = new TFile("/home/garillot/Code/PolyaFit/targetThr.root" , "READ") ;
	TTree* treer = dynamic_cast<TTree*>(filer->Get("tree")) ;

	treer->SetBranchAddress("LayerID" , &layerID) ;
	treer->SetBranchAddress("DifID" , &difID) ;
	treer->SetBranchAddress("AsicID" , &asicID) ;

	treer->SetBranchAddress("dac1" , &dac1) ;
	treer->SetBranchAddress("dac2" , &dac2) ;
	treer->SetBranchAddress("dac3" , &dac3) ;

	std::ofstream txtFile ;
	txtFile.open ("changeThr.py") ;
	txtFile << "import OracleAccess as oa" << std::endl ;
	txtFile << "s = oa.OracleAccess(\"Dome_42chambres_Reference_v4_244\")" << std::endl ;

	iEntry = 0 ;
	while ( treer->GetEntry(iEntry++) )
	{
		txtFile << "print 'change " << difID << " , " << asicID << "'" << std::endl ;
		txtFile << "s.ChangeThreshold(" << dac1 << "," << dac2 << "," << dac3 << "," << difID << "," << asicID <<")" << std::endl ;

	}

	txtFile << "s.uploadChanges()" << std::endl ;
	txtFile << "s.dumpStateNames()" << std::endl ;

	txtFile.close() ;

	return 0 ;
}
