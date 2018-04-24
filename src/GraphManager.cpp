#include "GraphManager.h"

#include "PolyaFitter.h"

#include <iostream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <set>
#include <algorithm>

#include <TDirectory.h>
#include <TList.h>

GraphManager::GraphManager()
{
}


GraphManager::~GraphManager()
{
	//	for ( std::map<AsicID,TGraphErrors*>::const_iterator it = graphMap.begin() ; it != graphMap.end() ; ++it )
	//	{
	//		if ( it->second )
	//			delete it->second ;
	//	}
}

void GraphManager::reset()
{
	//	for ( std::map<AsicID,TGraphErrors*>::const_iterator it = graphMap.begin() ; it != graphMap.end() ; ++it )
	//	{
	//		if ( it->second )
	//			delete it->second ;
	//	}
	graphMap.clear() ;
	resultMap.clear() ;
}

bool GraphManager::AsicID::operator<(const AsicID& b) const
{
	if ( this->layerID != b.layerID )
		return this->layerID < b.layerID ;
	else if ( this->difID != b.difID )
		return this->difID < b.difID ;
	else
		return this->asicID < b.asicID ;
}

bool GraphManager::AsicID::operator==(const AsicID& b) const
{
	return ( this->layerID == b.layerID ) && ( this->difID == b.difID ) && ( this->asicID == b.asicID ) ;
}


void GraphManager::openGraphs(std::string fileName)
{
	std::cout << "Open Graphs in " << fileName << std::endl ;
	graphMap.clear() ;

	TFile* file = new TFile(fileName.c_str() , "READ") ;
	file->cd() ;

	TDirectoryFile* graphsDir = dynamic_cast<TDirectoryFile*>( file->GetDirectory("Graphs") ) ;
	if (!graphsDir)
	{
		std::cerr << "ERROR : \"Graphs\" directory not present" << std::endl ;
		return ;
	}

	TList* layerList = graphsDir->GetListOfKeys() ;
	TIter iter(layerList) ;

	TObject* obj = NULL ;
	while ( (obj = iter()) )
	{
		std::string name( obj->GetName() ) ;
		if (name == std::string("Global") )
		{
			TGraphErrors* globalGraph = dynamic_cast<TGraphErrors*>( graphsDir->Get("Global") ) ;
			graphMap.insert( std::make_pair( AsicID(-1,-1,-1) , globalGraph) ) ;
			continue ;
		}

		std::string ref("Layer") ;
		std::size_t found = name.find(ref) ;
		if ( found == std::string::npos )
			continue ;

		TDirectoryFile* layerDir = dynamic_cast<TDirectoryFile*>( graphsDir->GetDirectory( name.c_str() ) ) ;
		openGraphsInLayer(layerDir) ;
	}

	std::cout << graphMap.size() << " graphs opened" << std::endl ;
}


void GraphManager::openGraphsInLayer(TDirectoryFile* layerDir)
{
	std::string layerIDstr = layerDir->GetName() ;
	std::string ref("Layer") ;
	std::size_t found = layerIDstr.find(ref) ;

	layerIDstr.replace( found , ref.length() , "" ) ;
	int layerID = atoi( layerIDstr.c_str() ) ;

	std::cout << "Open Layer " << layerID << "..." << std::endl ;

	TList* graphList = layerDir->GetListOfKeys() ;
	TIter iter(graphList) ;

	TObject* obj = NULL ;
	while ( (obj = iter()) )
	{
		std::string name( obj->GetName() ) ;

		std::size_t found = name.find(",") ;
		if ( found == std::string::npos )
			continue ;

		const long pos = static_cast<const long>(found) ;

		int difID = atoi( std::string( name.begin() , name.begin()+pos ).c_str() ) ;
		int asicID = atoi( std::string( name.begin()+pos+1 , name.end() ).c_str() ) ;

		TGraphErrors* graph = dynamic_cast<TGraphErrors*>( layerDir->Get( name.c_str() ) ) ;
		graphMap.insert( std::make_pair( AsicID(layerID , difID , asicID) , graph) ) ;
	}
}

void GraphManager::ProcessFile(std::string fileName)
{
	std::cout << "ProcessFile" << std::endl ;
	TFile* file = new TFile(fileName.c_str() , "READ") ;
	TTree* tree = dynamic_cast<TTree*>( file->Get("tree") ) ;
	if (!tree)
	{
		std::cout << "Tree not present in file" << std::endl ;
		return ;
	}

	std::vector<double>* efficiencies = NULL ;
	std::vector<double>* efficienciesError = NULL ;

	std::vector<double>* thresholds = reinterpret_cast< std::vector<double>* >( file->Get("Thresholds") ) ;

	int difID , asicID , layerID , padID ;
	std::vector<double>* multiplicities = NULL ;
	std::vector<double>* multiplicitiesError = NULL ;
	std::vector<double>* position = NULL ;

	tree->SetBranchAddress("LayerID" , &layerID) ;
	tree->SetBranchAddress("DifID" , &difID) ;
	tree->SetBranchAddress("AsicID" , &asicID) ;
	tree->SetBranchAddress("PadID" , &padID) ;
	tree->SetBranchAddress("Efficiencies" , &efficiencies ) ;
	tree->SetBranchAddress("EfficienciesError" , &efficienciesError ) ;
	tree->SetBranchAddress("Multiplicities" , &multiplicities) ;
	tree->SetBranchAddress("MultiplicitiesError" , &multiplicitiesError) ;
	tree->SetBranchAddress("Position" , &position) ;


	int iEntry = 0 ;
	while ( tree->GetEntry(iEntry++) )
	{
		if ( padID > -1 ) // because padID > -1 means stats for an individual pad
			continue ;

		AsicID asicKey(layerID , difID , asicID) ;

		//		if ( eff1 < 0.7 )
		if ( (*std::max_element(efficiencies->begin() , efficiencies->end() )) < 0.1 )
			continue ;

		TGraphErrors* graph = NULL ;
		graph = new TGraphErrors( static_cast<int>( thresholds->size() ) , &(*thresholds)[0] , &(*efficiencies)[0] , NULL , &(*efficienciesError)[0] ) ;
		graph->SetMarkerStyle(20) ;
		graphMap.insert( std::make_pair(asicKey , graph) ) ;

		mulMap.insert( std::make_pair(asicKey , multiplicities->at(0)) ) ;
		mulErrMap.insert( std::make_pair(asicKey , multiplicitiesError->at(0)) ) ;
		posMap.insert( std::make_pair( asicKey , std::vector<double>(*position) ) ) ;
	}

	file->Close() ;
}

void GraphManager::ProcessData(std::string dataPath)
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

	for ( unsigned int iRun = 0 ; iRun < runs.size() ; ++iRun )
	{
		std::vector<double> thresholds ;
		thresholds.push_back( thr1Vec.at(iRun) ) ;
		thresholds.push_back( thr2Vec.at(iRun) ) ;
		thresholds.push_back( thr3Vec.at(iRun) ) ;

		std::stringstream filePath ;
		filePath << dataPath << "/map_" << runs.at(iRun) << ".root" ;

		std::cout << "Process " << filePath.str() << std::endl ;
		TFile* file = new TFile( filePath.str().c_str() , "READ") ;
		TTree* tree = dynamic_cast<TTree*>( file->Get("tree") ) ;
		if ( !tree )
		{
			std::cout << "Error in ProcessData : tree not present in " << filePath.str() << std::endl ;
			file->Close() ;
			return ;
		}

		std::vector<double>* efficiencies = NULL ;
		std::vector<double>* efficienciesError = NULL ;

		int difID , asicID , layerID , padID ;
		double multiplicity , multiplicityError ;
		std::vector<double>* position = NULL ;
		//FIXME multiplicies is a vector now but still old data files need to reprocess data and update code
		tree->SetBranchAddress("LayerID" , &layerID) ;
		tree->SetBranchAddress("DifID" , &difID) ;
		tree->SetBranchAddress("AsicID" , &asicID) ;
		tree->SetBranchAddress("PadID" , &padID) ;
		tree->SetBranchAddress("Efficiencies" , &efficiencies) ;
		tree->SetBranchAddress("EfficienciesError" , &efficienciesError) ;
		tree->SetBranchAddress("Multiplicity" , &multiplicity) ;
		tree->SetBranchAddress("MultiplicityError" , &multiplicityError) ;
		tree->SetBranchAddress("Position" , &position) ;

		AsicID globalKey(-1,-1,-1) ;
		std::vector<double> globalEff(3 , 0.0) ;
		std::vector<double> globalEffErr(3 , 0.0) ;
		int nOkAsicsGlobal = 0 ;


		int iEntry = 0 ;
		while ( tree->GetEntry(iEntry++) )
		{
			if ( padID > -1 ) // because padID > -1 means stats for an individual pad
				continue ;

			AsicID asicKey(layerID , difID , asicID) ;

			if ( (*std::max_element(efficiencies->begin() , efficiencies->end() )) < 0.1 )
				continue ;

			if ( runs.at(iRun) == 730677 )
			{
				mulMap.insert( std::make_pair( asicKey , multiplicity) ) ;
				mulErrMap.insert( std::make_pair( asicKey , multiplicity) ) ;

				posMap.insert( std::make_pair( asicKey , std::vector<double>(*position) ) ) ;
				continue ;
			}

			if ( asicID == -1 )
				continue ;
			if ( layerID == 1 || layerID == 34 ) //Dead layers
				continue ;


			std::map<AsicID,TGraphErrors*>::const_iterator it = graphMap.find(asicKey) ;
			if ( it == graphMap.end() )
			{
				TGraphErrors* graph = new TGraphErrors ;
				graphMap.insert( std::make_pair(asicKey , graph) ) ;
				it = graphMap.find(asicKey) ;
			}

			for ( unsigned int i = 0 ; i < 3 ; ++i )
			{
				addPoint(it->second, thresholds.at(i), efficiencies->at(i), 0.0 , efficienciesError->at(i) ) ;

				globalEff.at(i) += efficiencies->at(i) ;
				globalEffErr.at(i) += 1.0/( efficienciesError->at(i)*efficienciesError->at(i) ) ;
			}
			nOkAsicsGlobal++ ;
		}

		file->Close() ;

		if ( runs.at(iRun) == 730677 )
			continue ;



		TGraphErrors* graph = new TGraphErrors ;
		graphMap.insert( std::make_pair(globalKey , graph) ) ;
		std::map<AsicID,TGraphErrors*>::const_iterator it = graphMap.find(globalKey) ;

		for ( unsigned int i = 0 ; i < 3 ; ++i )
		{
			globalEff.at(i) /= nOkAsicsGlobal ;
			globalEffErr.at(i) = std::sqrt( 1.0/globalEffErr.at(i) ) ;

			addPoint(it->second , thresholds.at(i) , globalEff.at(i) , 0.0 , globalEffErr.at(i) ) ;
		}

	}

}

void GraphManager::writeGraphsInFile(std::string fileName)
{
	TFile* file = new TFile(fileName.c_str() , "RECREATE") ;
	file->cd() ;
	TDirectory* dir = file->mkdir("Graphs") ;
	dir->cd() ;

	for ( std::map<AsicID,TGraphErrors*>::iterator it = graphMap.begin() ; it != graphMap.end() ; ++it )
	{
		std::stringstream graphName ;

		if ( it->first == AsicID(-1,-1,-1) )
		{
			dir->cd() ;
			graphName << "Global" ;
		}
		else
		{
			std::stringstream dirName ; dirName << "Layer" << it->first.layerID ;
			TDirectory* layerDir = dir->GetDirectory( dirName.str().c_str() ) ;
			if ( !layerDir )
				layerDir = dir->mkdir( dirName.str().c_str() ) ;

			layerDir->cd() ;

			graphName << it->first.difID << "," << it->first.asicID ;
		}

		it->second->SetMarkerStyle(20) ;
		it->second->Write(graphName.str().c_str() ) ;
	}

	file->Close() ;
}


PolyaFitter::PolyaFitResult GraphManager::fitGraph(int layer , int dif , int asic)
{
	AsicID id(layer,dif,asic) ;
	std::map<AsicID,TGraphErrors*>::iterator it = graphMap.find( id ) ;
	if ( it == graphMap.end() )
	{
		std::cerr << "ERROR : graph not present" << std::endl ;
		return PolyaFitter::PolyaFitResult() ;
	}

	std::cout << "Fit graph " ; id.print() ;
	PolyaFitter b ;
	b.getPoints( it->second ) ;
	b.setParams() ;
	b.minimize() ;
	return b.getFitResult() ;
}

std::map<GraphManager::AsicID,PolyaFitter::PolyaFitResult> GraphManager::fitAllGraphs()
{
	resultMap.clear() ;

	for ( std::map<AsicID,TGraphErrors*>::iterator it = graphMap.begin() ; it != graphMap.end() ; ++it )
		resultMap.insert( std::make_pair(it->first , fitGraph(it->first.layerID , it->first.difID , it->first.asicID) ) ) ;

	return resultMap ;
}

void GraphManager::writeResultTree(std::string fileName)
{
	TFile* file = new TFile(fileName.c_str() , "RECREATE") ;
	TTree* tree = new TTree("tree","tree") ;

	int layerID ;
	int difID ;
	int asicID ;
	double qbar , qbarError ;
	double delta , deltaError ;
	double eff0 , eff0Error ;
	double mul , mulError ;
	double chi2 ;
	int minimStatus ;
	std::vector<double> position ;

	tree->Branch("LayerID" , &layerID) ;
	tree->Branch("DifID" , &difID) ;
	tree->Branch("AsicID" , &asicID) ;
	tree->Branch("mul" , &mul) ;
	tree->Branch("mulErr" , &mulError) ;
	tree->Branch("qbar" , &qbar) ;
	tree->Branch("qbarError" , &qbarError) ;
	tree->Branch("delta" , &delta) ;
	tree->Branch("deltaError" , &deltaError) ;
	tree->Branch("eff0" , &eff0) ;
	tree->Branch("eff0Error" , &eff0Error) ;
	tree->Branch("chi2" , &chi2) ;
	tree->Branch("minimStatus" , &minimStatus) ;
	tree->Branch("Position" , &position) ;

	for ( std::map<AsicID,PolyaFitter::PolyaFitResult>::const_iterator it = resultMap.begin() ; it != resultMap.end() ; ++it )
	{
		AsicID id = it->first ;
		PolyaFitter::PolyaFitResult res = it->second ;

		if ( mulMap.find(id) == mulMap.end() )
		{
			continue ;
			//			layerID = id.layerID ;
			//			difID = id.difID ;
			//			asicID = id.asicID ;

			//			qbar = -1 ;
			//			qbarError = -1 ;
			//			delta = -1 ;
			//			deltaError = -1 ;
			//			eff0 = 0 ;
			//			eff0Error = 0 ;
			//			chi2 = -1 ;

			//			minimStatus = -1 ;

			//			mul = 0 ;
			//			mulError = 0 ;
			//			position = posMap[it->first] ;

			//			tree->Fill() ;

			//			continue ;
		}

		layerID = id.layerID ;
		difID = id.difID ;
		asicID = id.asicID ;

		qbar = res.qbar ;
		qbarError = res.qbarError ;
		delta = res.delta ;
		deltaError = res.deltaError ;
		eff0 = res.eff0 ;
		eff0Error = res.eff0Error ;
		chi2 = res.chi2 ;

		minimStatus = res.minimStatus ;

		mul = mulMap[it->first] ;
		mulError = mulErrMap[it->first] ;

		position = posMap[it->first] ;

		tree->Fill() ;
	}

	file->cd() ;
	tree->Write("tree") ;

	file->Close() ;
}

void GraphManager::writeResultTree(double qbar , double delta)
{
	std::stringstream fileName ; fileName << qbar << "_" << delta << "_results.root" ;
	writeResultTree( fileName.str() ) ;
}


void GraphManager::addPoint(TGraphErrors* graph, double x, double y, double ex , double ey)
{
	if (!graph)
	{
		std::cerr << "ERROR in GraphManager::addPoint : graph ptr = NULL" << std::endl ;
		return ;
	}
	int point = graph->GetN() ;
	graph->SetPoint(point , x , y) ;
	graph->SetPointError(point , ex , ey) ;
}

TGraphErrors* GraphManager::getGraph(AsicID id) const
{
	std::map<AsicID,TGraphErrors*>::const_iterator it = graphMap.find( id ) ;
	if ( it == graphMap.end() )
		return nullptr ;
	else
		return it->second ;
}

TGraphErrors* GraphManager::getGraph(int layer , int dif , int asic) const
{
	return getGraph( AsicID(layer,dif,asic) ) ;
}

TGraphErrors* GraphManager::getGlobalGraph() const
{
	return getGraph(-1,-1,-1) ;
}



