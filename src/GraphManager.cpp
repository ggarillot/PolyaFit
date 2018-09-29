#include "GraphManager.h"

#include "PolyaFitter.h"

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <sstream>
#include <set>
#include <algorithm>
#include <cassert>

#include <TDirectory.h>
#include <TList.h>
#include <TEfficiency.h>

void GraphManager::reset()
{
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

	TObject* obj = nullptr ;
	while ( (obj = iter()) )
	{
		std::string name( obj->GetName() ) ;
		if (name == std::string("Global") )
		{
			TGraphAsymmErrors* globalGraph = dynamic_cast<TGraphAsymmErrors*>( graphsDir->Get("Global") ) ;
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

	TObject* obj = nullptr ;
	while ( (obj = iter()) )
	{
		std::string name( obj->GetName() ) ;

		std::size_t found2 = name.find(",") ;
		if ( found2 == std::string::npos )
			continue ;

		const long pos = static_cast<const long>(found2) ;

		int difID = atoi( std::string( name.begin() , name.begin()+pos ).c_str() ) ;
		int asicID = atoi( std::string( name.begin()+pos+1 , name.end() ).c_str() ) ;

		TGraphAsymmErrors* graph = dynamic_cast<TGraphAsymmErrors*>( layerDir->Get( name.c_str() ) ) ;
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

	std::vector<double>* efficiencies = nullptr ;
	std::vector<double>* efficienciesLowerBound = nullptr ;
	std::vector<double>* efficienciesUpperBound = nullptr ;

	std::vector<double>* thresholds = reinterpret_cast< std::vector<double>* >( file->Get("Thresholds") ) ;

	int difID , asicID , layerID , padID ;
	std::vector<double>* multiplicities = nullptr ;
	std::vector<double>* multiplicitiesError = nullptr ;
	std::vector<double>* position = nullptr ;

	tree->SetBranchAddress("LayerID" , &layerID) ;
	tree->SetBranchAddress("DifID" , &difID) ;
	tree->SetBranchAddress("AsicID" , &asicID) ;
	tree->SetBranchAddress("PadID" , &padID) ;
	tree->SetBranchAddress("Efficiencies" , &efficiencies ) ;
	tree->SetBranchAddress("EfficienciesLowerBound" , &efficienciesLowerBound) ;
	tree->SetBranchAddress("EfficienciesUpperBound" , &efficienciesUpperBound) ;
	tree->SetBranchAddress("Multiplicities" , &multiplicities) ;
	tree->SetBranchAddress("MultiplicitiesError" , &multiplicitiesError) ;
	tree->SetBranchAddress("Position" , &position) ;


	int iEntry = 0 ;
	while ( tree->GetEntry(iEntry++) )
	{
		if ( padID > -1 ) // because padID > -1 means stats for an individual pad
			continue ;

		AsicID asicKey(layerID , difID , asicID) ;

		if ( (*std::max_element(efficiencies->begin() , efficiencies->end() )) < 0.1 )
			continue ;

		posMap.insert( std::make_pair( asicKey , std::vector<double>(*position) ) ) ;

		TGraphAsymmErrors* graph = nullptr ;

		std::vector<double> low( efficiencies->size() , 0 ) ;
		std::vector<double> high( efficiencies->size() , 0 ) ;

		for ( unsigned int i = 0 ; i < efficiencies->size() ; ++i )
		{
			low[i] = efficiencies->at(i) - efficienciesLowerBound->at(i) ;
			high[i] = efficienciesUpperBound->at(i) - efficiencies->at(i) ;
		}
		graph = new TGraphAsymmErrors( static_cast<int>( thresholds->size() ) , &(*thresholds)[0] , &(*efficiencies)[0] , nullptr , nullptr , &(low)[0] , &(high)[0] ) ;
		graph->SetMarkerStyle(20) ;
		graphMap.insert( std::make_pair(asicKey , graph) ) ;

		mulMap.insert( std::make_pair(asicKey , multiplicities->at(0)) ) ;
		mulErrMap.insert( std::make_pair(asicKey , multiplicitiesError->at(0)) ) ;
	}

	file->Close() ;
}

void GraphManager::ProcessData(std::string jsonFileName)
{
	std::ifstream jsonFile(jsonFileName) ;
	auto json = nlohmann::json::parse(jsonFile) ;

	std::string directory = json.at("directory") ;

	std::vector<std::string> fileList = {} ;
	std::vector<std::array<int,3>> thrList = {} ;
	unsigned int mulRef = std::numeric_limits<unsigned int>::max() ;

	auto list = json.at("runs") ;
	for ( const auto& i : list )
	{
		auto files = i.at("files") ;
		for ( const auto& file : files )
		{
			fileList.push_back( file ) ;
			thrList.push_back( i.at("thresholds") ) ;
		}
//		fileList.push_back( i.at("file") ) ;
//		thrList.push_back( i.at("thresholds") ) ;
		if ( i.count("mulRef") )
		{
			if ( mulRef != std::numeric_limits<unsigned int>::max() )
				std::cout << "WARNING : multiple mulRef defined" << std::endl ;
			mulRef = static_cast<unsigned int>( fileList.size()-1 ) ;
		}
	}

	assert( fileList.size() == thrList.size() ) ;

	for ( unsigned int i = 0 ; i < fileList.size() ; i++ )
	{
		std::vector<double> thresholds ;
		thresholds.push_back( (thrList.at(i)[0]-90)/700.0 ) ;
		thresholds.push_back( (thrList.at(i)[1]-98)/80.0 ) ;
		thresholds.push_back( (thrList.at(i)[2]-98)/16.3 ) ;

		std::stringstream fileName ;
		fileName << directory << "/" << fileList.at(i) ;
		std::cout << "Process " << fileName.str() << std::endl ;
		TFile* file = new TFile( fileName.str().c_str() , "READ") ;
		TTree* tree = dynamic_cast<TTree*>( file->Get("tree") ) ;
		if ( !tree )
		{
			std::cout << "Error in ProcessData : tree not present in " << fileName.str() << std::endl ;
			file->Close() ;
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
		tree->SetBranchAddress("Efficiencies" , &efficiencies) ;
		tree->SetBranchAddress("EfficienciesLowerBound" , &efficienciesLowerBound) ;
		tree->SetBranchAddress("EfficienciesUpperBound" , &efficienciesUpperBound) ;
		tree->SetBranchAddress("Multiplicities" , &multiplicities) ;
		tree->SetBranchAddress("MultiplicitiesError" , &multiplicitiesError) ;
		tree->SetBranchAddress("Position" , &position) ;
		tree->SetBranchAddress("Ntrack" , &nTrack) ;

		AsicID globalKey(-1,-1,-1) ;
		std::vector<double> globalEff(3 , 0.0) ;

		std::vector<double> globalMul(3 , 0.0) ;
		std::vector<double> globalMulSq(3 , 0.0) ;

		//		double globalMul = 0.0 ;
		//		double globalMulSq = 0.0 ;

		int nOkAsicsGlobal = 0 ;


		int iEntry = 0 ;
		while ( tree->GetEntry(iEntry++) )
		{
			if ( padID > -1 ) // because padID > -1 means stats for an individual pad
				continue ;

			AsicID asicKey(layerID , difID , asicID) ;

			if ( (*std::max_element(efficiencies->begin() , efficiencies->end() )) < 0.1 )
				continue ;

			posMap.insert( std::make_pair( asicKey , std::vector<double>(*position) ) ) ;

			if ( i == mulRef )
			{
				mulMap.insert( std::make_pair( asicKey , multiplicities->at(0) ) ) ;
				mulErrMap.insert( std::make_pair( asicKey , multiplicitiesError->at(0) ) ) ;
			}

			std::map<AsicID,TGraphAsymmErrors*>::const_iterator it = graphMap.find(asicKey) ;
			std::map<AsicID,TGraphErrors*>::const_iterator itMul = graphMulMap.find(asicKey) ;
			if ( it == graphMap.end() )
			{
				TGraphAsymmErrors* graph = new TGraphAsymmErrors ;
				graphMap.insert( std::make_pair(asicKey , graph) ) ;
				it = graphMap.find(asicKey) ;
			}

			if ( itMul == graphMulMap.end() )
			{
				TGraphErrors* graph = new TGraphErrors ;
				graphMulMap.insert( std::make_pair(asicKey , graph) ) ;
				itMul = graphMulMap.find(asicKey) ;
			}

			for ( unsigned int j = 0 ; j < 3 ; ++j )
			{
				auto errLow = efficiencies->at(j) - efficienciesLowerBound->at(j) ;
				auto errHigh = efficienciesUpperBound->at(j) - efficiencies->at(j) ;

				assert( errLow >= 0 && errHigh >= 0 ) ;
				addPoint(it->second, thresholds.at(j), efficiencies->at(j) , errLow , errHigh ) ;

				globalEff.at(j) += efficiencies->at(j) ;

				if ( j > 0 )
					continue ;
				if ( multiplicities->at(j) < std::numeric_limits<double>::epsilon() )
					continue ;
				addPoint(itMul->second, thresholds.at(j), multiplicities->at(j) , multiplicitiesError->at(j) ) ;
				globalMul.at(j) += multiplicities->at(j) ;
				globalMulSq.at(j) += multiplicities->at(j)*multiplicities->at(j) ;

			}

			//			if ( !(multiplicities->at(0) < std::numeric_limits<double>::epsilon()) )
			//			{
			//				//				continue ;
			//				addPoint(itMul->second, thresholds.at(0), multiplicities->at(0) , multiplicitiesError->at(0) ) ;
			//			}
			//			globalMul += multiplicities->at(0) ;
			//			globalMulSq += multiplicities->at(0)*multiplicities->at(0) ;

			globalNTrack += nTrack ;
			nOkAsicsGlobal++ ;
		}

		file->Close() ;

		//		if ( runs.at(iRun) == 730677 )
		//			continue ;



//		TGraphAsymmErrors* graph = new TGraphAsymmErrors ;
//		graphMap.insert( std::make_pair(globalKey , graph) ) ;
//		std::map<AsicID,TGraphAsymmErrors*>::const_iterator it = graphMap.find(globalKey) ;

//		TGraphErrors* graphMul = new TGraphErrors ;
//		graphMulMap.insert( std::make_pair(globalKey , graphMul) ) ;
//		std::map<AsicID,TGraphErrors*>::const_iterator itMul = graphMulMap.find(globalKey) ;

//		for ( unsigned int j = 0 ; j < 3 ; ++j )
//		{
//			globalEff.at(j) /= nOkAsicsGlobal ;


//			constexpr double level = 0.683 ;

//			double a = globalEff.at(j)*globalNTrack + 1 ;
//			double b = globalNTrack - globalEff.at(j)*globalNTrack + 1 ;

//			double lowerBound = 0 ;
//			double upperBound = 0 ;
//			TEfficiency::BetaShortestInterval( level , a , b , lowerBound , upperBound ) ;

//			auto errLow = globalEff.at(j) - lowerBound ;
//			auto errHigh = upperBound - globalEff.at(j) ;

//			assert( errLow > 0 && errHigh > 0 ) ;

//			addPoint(it->second , thresholds.at(j) , globalEff.at(j) , errLow , errHigh ) ;


//			if ( j > 0 )
//				continue ;
//			double mulErr = 0.0 ;
//			double var = globalMulSq.at(j)/globalNTrack - (globalMul.at(j)/globalNTrack)*(globalMul.at(j)/globalNTrack) ;

//			if ( var < std::numeric_limits<double>::epsilon() )
//				var = 1.0/( std::sqrt(12*globalNTrack) ) ;

//			mulErr = sqrt( var/(globalNTrack-1.0) ) ;

//			globalMul.at(j) /= nOkAsicsGlobal ;
//			addPoint(itMul->second , thresholds.at(j) , globalMul.at(j) , mulErr ) ;
//		}


		//		double mulErr = 0.0 ;
		//		double var = globalMulSq/globalNTrack - (globalMul/globalNTrack)*(globalMul/globalNTrack) ;

		//		if ( var < std::numeric_limits<double>::epsilon() )
		//			var = 1.0/( std::sqrt(12*globalNTrack) ) ;

		//		mulErr = sqrt( var/(globalNTrack-1.0) ) ;

		//		globalMul /= nOkAsicsGlobal ;
		//		addPoint(itMul->second , thresholds.at(0) , globalMul , mulErr ) ;
	}

}

void GraphManager::writeGraphsInFile(std::string fileName)
{
	TFile* file = new TFile(fileName.c_str() , "RECREATE") ;
	file->cd() ;
	TDirectory* dir = file->mkdir("Graphs") ;
	dir->cd() ;

	for ( std::map<AsicID,TGraphAsymmErrors*>::iterator it = graphMap.begin() ; it != graphMap.end() ; ++it )
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
	std::map<AsicID,TGraphAsymmErrors*>::iterator it = graphMap.find( id ) ;
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

	for ( std::map<AsicID,TGraphAsymmErrors*>::iterator it = graphMap.begin() ; it != graphMap.end() ; ++it )
		resultMap.insert( std::make_pair(it->first , fitGraph(it->first.layerID , it->first.difID , it->first.asicID) ) ) ;

	return resultMap ;
}


MultiplicityFitter::MulFitResult GraphManager::fitMulGraph(int layer , int dif , int asic)
{
	AsicID id(layer,dif,asic) ;
	std::map<AsicID,TGraphErrors*>::iterator it = graphMulMap.find( id ) ;
	if ( it == graphMulMap.end() )
	{
		std::cerr << "ERROR : graph not present" << std::endl ;
		return MultiplicityFitter::MulFitResult() ;
	}

	std::cout << "Fit graph " ; id.print() ;
	MultiplicityFitter b ;
	b.getPoints( it->second ) ;
	b.setParams() ;
	b.minimize() ;
	return b.getFitResult() ;
}

std::map<GraphManager::AsicID,MultiplicityFitter::MulFitResult> GraphManager::fitAllMulGraphs()
{
	resultMulMap.clear() ;

	for ( std::map<AsicID,TGraphErrors*>::iterator it = graphMulMap.begin() ; it != graphMulMap.end() ; ++it )
		resultMulMap.insert( std::make_pair(it->first , fitMulGraph(it->first.layerID , it->first.difID , it->first.asicID) ) ) ;

	return resultMulMap ;
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

	double factor , factorError ;
	double power , powerError ;
	double constant , constantError ;
	double chi2Mul ;
	int minimStatusMul ;

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

	tree->Branch("factor" , &factor) ;
	tree->Branch("factorError" , &factorError) ;
	tree->Branch("power" , &power) ;
	tree->Branch("powerError" , &powerError) ;
	tree->Branch("constant" , &constant) ;
	tree->Branch("constantError" , &constantError) ;
	tree->Branch("chi2Mul" , &chi2Mul) ;
	tree->Branch("minimStatusMul" , &minimStatusMul) ;

	tree->Branch("Position" , &position) ;

	for ( std::map<AsicID,PolyaFitter::PolyaFitResult>::const_iterator it = resultMap.begin() ; it != resultMap.end() ; ++it )
	{
		AsicID id = it->first ;
		PolyaFitter::PolyaFitResult res = it->second ;

		if ( posMap.find(id) == posMap.end() )
			continue ;

		MultiplicityFitter::MulFitResult resMul = resultMulMap.at(id) ;
		if ( resultMulMap.find(id) == resultMulMap.end() )
		{
			std::cout << "WTF" << std::endl ;
			continue ;
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

		factor = resMul.f ;
		factorError = resMul.fErr ;
		power = resMul.p ;
		powerError = resMul.pErr ;
		constant = resMul.c ;
		constantError = resMul.cErr ;

		chi2Mul = resMul.chi2 ;
		minimStatusMul = resMul.minimStatus ;


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

void GraphManager::addPoint(TGraphErrors* graph , double x , double y , double ey)
{
	if (!graph)
	{
		std::cerr << "ERROR in GraphManager::addPoint : graph ptr = NULL" << std::endl ;
		return ;
	}
	int point = graph->GetN() ;
	graph->SetPoint(point , x , y) ;
	graph->SetPointError(point , 0 , ey) ;
}

void GraphManager::addPoint(TGraphAsymmErrors* graph , double x , double y , double ey)
{
	if (!graph)
	{
		std::cerr << "ERROR in GraphManager::addPoint : graph ptr = NULL" << std::endl ;
		return ;
	}
	int point = graph->GetN() ;
	graph->SetPoint(point , x , y) ;
	graph->SetPointError(point , 0 , 0 , ey , ey) ;
}

void GraphManager::addPoint(TGraphAsymmErrors* graph , double x , double y , double eylow , double eyhigh)
{
	if (!graph)
	{
		std::cerr << "ERROR in GraphManager::addPoint : graph ptr = NULL" << std::endl ;
		return ;
	}
	int point = graph->GetN() ;
	graph->SetPoint(point , x , y) ;
	graph->SetPointError(point , 0 , 0 , eylow , eyhigh) ;
}

TGraphAsymmErrors* GraphManager::getGraph(AsicID id) const
{
	std::map<AsicID,TGraphAsymmErrors*>::const_iterator it = graphMap.find( id ) ;
	if ( it == graphMap.end() )
		return nullptr ;
	else
		return it->second ;
}

TGraphAsymmErrors* GraphManager::getGraph(int layer , int dif , int asic) const
{
	return getGraph( AsicID(layer,dif,asic) ) ;
}

TGraphAsymmErrors* GraphManager::getGlobalGraph() const
{
	return getGraph(-1,-1,-1) ;
}


TGraphErrors* GraphManager::getMulGraph(AsicID id) const
{
	std::map<AsicID,TGraphErrors*>::const_iterator it = graphMulMap.find( id ) ;
	if ( it == graphMulMap.end() )
		return nullptr ;
	else
		return it->second ;
}

TGraphErrors* GraphManager::getMulGraph(int layer , int dif , int asic) const
{
	return getMulGraph( AsicID(layer,dif,asic) ) ;
}

TGraphErrors* GraphManager::getMulGlobalGraph() const
{
	return getMulGraph(-1,-1,-1) ;
}



