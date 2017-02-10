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


struct ResultPair
{
		PolyaFitter::PolyaFitResult moi ;
		PolyaFitter::PolyaFitResult arnaud ;
		bool isMoi ;
		bool isArnaud ;
} ;

int main()
{
	GraphManager a ;

	a.openGraphs("Data_graphs.root") ;
	std::map<GraphManager::AsicID,PolyaFitter::PolyaFitResult> map = a.fitAllGraphs() ;
	a.writeResultTree("res.root") ;

	TFile* file = new TFile("res.root" , "READ") ;
	TTree* tree = dynamic_cast<TTree*> ( file->Get("tree") ) ;

	GraphManager::AsicID id ;
	PolyaFitter::PolyaFitResult res ;

	float qbar ;
	float qbarError ;
	float delta ;
	float deltaError ;
	float eff0 ;
	float chi2 ;

	tree->SetBranchAddress("layerID" , &id.layerID) ;
	tree->SetBranchAddress("difID" , &id.difID) ;
	tree->SetBranchAddress("asicID" , &id.asicID) ;
	tree->SetBranchAddress("qbar" , &res.qbar) ;
	tree->SetBranchAddress("qbarError" , &res.qbarError) ;
	tree->SetBranchAddress("delta" , &res.delta) ;
	tree->SetBranchAddress("deltaError" , &res.deltaError) ;
	tree->SetBranchAddress("eff0" , &res.eff0) ;
	tree->SetBranchAddress("eff0Error" , &res.eff0Error) ;
	tree->SetBranchAddress("chi2" , &res.chi2) ;

	std::map<GraphManager::AsicID , ResultPair> resMap ;

	int iEntry = 0 ;
	while( tree->GetEntry(iEntry++) )
	{
		if ( id == GraphManager::AsicID(-1,0,0) )
			continue ;

		std::map<GraphManager::AsicID , ResultPair>::iterator it = resMap.find(id) ;
		if ( it == resMap.end() )
		{
			ResultPair pair ;
			pair.isArnaud = false ;
			pair.isMoi = true ;
			pair.moi = res ;
			resMap.insert( std::make_pair(id , pair) ) ;
		}
		else
		{
			it->second.moi = res ;
			it->second.isMoi = true ;
		}
	}

	file->Close() ;


	for ( int i = 0 ; i < 48 ; i++ )
	{
		if ( i==1 || i ==34 )
			continue ;

		int layerID = i ;
		std::stringstream fileName ; fileName << "/home/garillot/duringTB/dataLayers/" << "layer" << i << ".root" ;
		TFile* file = new TFile(fileName.str().c_str() , "READ") ;
		TTree* tree = dynamic_cast<TTree*> ( file->Get("tree") ) ;

		tree->SetBranchAddress("dif" , &id.difID) ;
		tree->SetBranchAddress("asic" , &id.asicID) ;
		tree->SetBranchAddress("polya_qbar" , &qbar) ;
		tree->SetBranchAddress("polya_qbar_error" , &qbarError) ;
		tree->SetBranchAddress("polya_delta" , &delta) ;
		tree->SetBranchAddress("polya_delta_error" , &deltaError) ;
		tree->SetBranchAddress("polya_eff0" , &eff0) ;
		tree->SetBranchAddress("efficiency_chi2" , &chi2) ;

		int iEntry = 0 ;
		while( tree->GetEntry(iEntry++) )
		{
			id.layerID = layerID ;
			res.qbar = qbar ;
			res.delta = delta ;
			res.eff0 = eff0 ;
			res.chi2 = chi2 ;
			res.qbarError = qbarError ;
			res.deltaError = deltaError ;
			res.eff0Error = 0.0 ;

			std::map<GraphManager::AsicID , ResultPair>::iterator it = resMap.find(id) ;
			if ( it == resMap.end() )
			{
				ResultPair pair ;
				pair.arnaud = res ;
				pair.isMoi = false ;
				pair.isArnaud = true ;
				resMap.insert( std::make_pair(id , pair) ) ;
			}
			else
			{
				it->second.arnaud = res ;
				it->second.isArnaud = true ;
			}
		}

		file->Close() ;
	}

	TFile* outFile = new TFile("compar.root" , "RECREATE") ;
	tree = new TTree("tree","tree") ;

	GraphManager::AsicID iD ;
	PolyaFitter::PolyaFitResult resMoi ;
	PolyaFitter::PolyaFitResult resArnaud ;

	tree->Branch("layerID" , &iD.layerID) ;
	tree->Branch("difID" , &iD.difID) ;
	tree->Branch("asicID" , &iD.asicID) ;

	tree->Branch("qbarMoi" , &resMoi.qbar) ;
	tree->Branch("deltaMoi" , &resMoi.delta) ;
	tree->Branch("eff0Moi" , &resMoi.eff0) ;
	tree->Branch("chi2Moi" , &resMoi.chi2) ;

	tree->Branch("qbarArnaud" , &resArnaud.qbar) ;
	tree->Branch("deltarnaud" , &resArnaud.delta) ;
	tree->Branch("eff0rArnaud" , &resArnaud.eff0) ;
	tree->Branch("chi2Arnaud" , &resArnaud.chi2) ;


	for ( std::map<GraphManager::AsicID , ResultPair>::iterator it = resMap.begin() ; it != resMap.end() ; ++it )
	{
		if ( !(it->second.isMoi && it->second.isArnaud) )
		{
			if (it->second.isMoi && !(it->second.isArnaud) )
				resArnaud = PolyaFitter::PolyaFitResult(-1,-1,-1,-1,-1,-1,-1,-1) ;
			else
				continue ;
		}


		iD = it->first ;
		resMoi = it->second.moi ;
		if ( it->second.isArnaud )
			resArnaud = it->second.arnaud ;
		tree->Fill() ;

	}

	outFile->cd() ;
	tree->Write() ;
	return 0 ;
}
