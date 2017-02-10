#ifndef GraphManager_h
#define GraphManager_h

#include <TFile.h>
#include <TTree.h>
#include <TGraphErrors.h>

#include <iostream>
#include <map>

#include "PolyaFitter.h"

class GraphManager
{
	public :
		GraphManager() ;
		~GraphManager() ;

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

		void writeResultTree(std::string fileName) ;
		void writeResultTree(double qbar , double delta) ;

		void ProcessFile(std::string fileName) ;
		void ProcessData(std::string dataPath = "") ;

		TGraphErrors* getGraph(int layer , int dif , int asic) const ;
		TGraphErrors* getGlobalGraph() const ;

		void reset() ;

	protected :

		static void addPoint(TGraphErrors* graph , double x , double y , double ex , double ey) ;

		void openGraphsInLayer(TDirectoryFile* layerDir) ;


		std::map<AsicID,TGraphErrors*> graphMap ;
		std::map<AsicID,PolyaFitter::PolyaFitResult> resultMap ;
		std::map<AsicID,double> mulMap ;
		std::map<AsicID,double> mulErrMap ;

} ;

#endif //GraphManager_h
