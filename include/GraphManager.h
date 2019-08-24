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
		~GraphManager() = default ;

		struct AsicID
		{
                AsicID(int lID = -1 , int dID = 0 , int aID = 0)
                    : layerID(lID) , difID(dID) , asicID(aID) {}
				int layerID ;
				int difID ;
				int asicID ;
				bool operator<(const AsicID& b) const ;
				bool operator==(const AsicID& b) const ;
				void print() const
				{
					std::cout << "Layer " << layerID << " , DIF " << difID << " , ASIC " << asicID << std::endl ;
				}

		} ;
		friend std::ostream& operator<<(std::ostream& os , const AsicID& asicID)
		{
			os << "Layer " << asicID.layerID << " , DIF " << asicID.difID << " , ASIC " << asicID.asicID ;
			return os ;
		}


		void openGraphs(std::string fileName) ;

		void writeGraphsInFile(std::string fileName) ;

        PolyaFitResult fitEffGraph(int layer , int dif , int asic) const ;
        PolyaFitResult fitEffGraph(AsicID id) const ;
        std::map<AsicID,PolyaFitResult> fitAllEffGraphs() ;

		MultiplicityFitResult fitMulGraph(int layer , int dif , int asic) const ;
		MultiplicityFitResult fitMulGraph(AsicID id) const ;
		std::map<AsicID,MultiplicityFitResult> fitAllMulGraphs() ;

		void writeResultTree(std::string fileName) ;
		void writeResultTree(double qbar , double delta) ;

        void ProcessSimu(std::string fileName) ;
		void ProcessData(std::string jsonFileName) ;

		TGraphAsymmErrors* getGraph(int layer , int dif , int asic) const ;
		TGraphAsymmErrors* getGraph(AsicID id) const ;
		TGraphAsymmErrors* getGlobalGraph() const ;

		TGraphErrors* getMulGraph(int layer , int dif , int asic) const ;
		TGraphErrors* getMulGraph(AsicID id) const ;
		TGraphErrors* getMulGlobalGraph() const ;

		void reset() ;

        void setColor(Color_t col) { color = col ; }

	protected :

        Color_t color = kBlack ;

		void openGraphsInLayer(TDirectoryFile* layerDir) ;


		std::map<AsicID,TGraphAsymmErrors*> graphMap = {} ;
		std::map<AsicID,PolyaFitResult> resultMap = {} ;

		std::map<AsicID,TGraphErrors*> graphMulMap = {} ;
		std::map<AsicID,MultiplicityFitResult> resultMulMap = {} ;


		std::map<AsicID,double> mulMap = {} ;
		std::map<AsicID,double> mulErrMap = {} ;
		std::map<AsicID, std::vector<double> > posMap = {} ;

} ;

#endif //GraphManager_h
