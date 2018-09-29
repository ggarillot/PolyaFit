#ifndef GraphManager_h
#define GraphManager_h

#include <TFile.h>
#include <TTree.h>

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>

#include <iostream>
#include <map>

#include "json.hpp"

#include "PolyaFitter.h"
#include "MultiplicityFitter.h"

class GraphManager
{
	public :
		GraphManager() = default ;
		~GraphManager()  = default ;

		struct AsicID
		{
				AsicID() : layerID(-1) , difID(0) , asicID(0) {}
				AsicID(int lID , int dID , int aID) : layerID(lID) , difID(dID) , asicID(aID) {}
				int layerID ;
				int difID ;
				int asicID ;
				bool operator<(const AsicID& b) const ;
				bool operator==(const AsicID& b) const ;
				void print() const
				{
					std::cout << "Layer " << layerID << " , Dif " << difID << " , Asic " << asicID << std::endl ;
				}

		} ;
		friend std::ostream& operator<<(std::ostream& os , const AsicID& asicID)
		{
			os << "Layer " << asicID.layerID << " , Dif " << asicID.difID << " , Asic " << asicID.asicID ;
			return os ;
		}


		void openGraphs(std::string fileName) ;

		//		void createGraphsData() ;

		void writeGraphsInFile(std::string fileName) ;

		PolyaFitter::PolyaFitResult fitGraph(int layer , int dif , int asic) ;
		std::map<AsicID,PolyaFitter::PolyaFitResult> fitAllGraphs() ;

		MultiplicityFitter::MulFitResult fitMulGraph(int layer , int dif , int asic) ;
		std::map<AsicID,MultiplicityFitter::MulFitResult> fitAllMulGraphs() ;

		void writeResultTree(std::string fileName) ;
		void writeResultTree(double qbar , double delta) ;

		void ProcessFile(std::string fileName) ;
		void ProcessData(std::string jsonFileName) ;

		TGraphAsymmErrors* getGraph(int layer , int dif , int asic) const ;
		TGraphAsymmErrors* getGraph(AsicID id) const ;
		TGraphAsymmErrors* getGlobalGraph() const ;

		TGraphErrors* getMulGraph(int layer , int dif , int asic) const ;
		TGraphErrors* getMulGraph(AsicID id) const ;
		TGraphErrors* getMulGlobalGraph() const ;

		void reset() ;

	protected :

		static void addPoint(TGraphErrors* graph , double x , double y , double ey) ;
		static void addPoint(TGraphAsymmErrors* graph , double x , double y , double ey) ;
		static void addPoint(TGraphAsymmErrors* graph , double x , double y , double eylow , double eyhigh) ;

		void openGraphsInLayer(TDirectoryFile* layerDir) ;


		std::map<AsicID,TGraphAsymmErrors*> graphMap = {} ;
		std::map<AsicID,PolyaFitter::PolyaFitResult> resultMap = {} ;

		std::map<AsicID,TGraphErrors*> graphMulMap = {} ;
		std::map<AsicID,MultiplicityFitter::MulFitResult> resultMulMap = {} ;


		std::map<AsicID,double> mulMap = {} ;
		std::map<AsicID,double> mulErrMap = {} ;
		std::map<AsicID, std::vector<double> > posMap = {} ;

} ;

#endif //GraphManager_h
