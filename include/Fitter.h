#ifndef Fitter_h
#define Fitter_h

#include "Minimisation.h"

#include <vector>
#include <string>
#include <limits>

#include <TGraph.h>

class Fitter : public Minimisation
{
	public :
		Fitter(unsigned int nParam_) ;
		virtual ~Fitter() = default ;

		void getPointsFromGraph(const TGraph* graph) ;

		virtual double eval(double x , const double* params) = 0 ;

		double functionToMinimize(const double* param) ;
		double minimize() ;

		virtual void createResultStruct() = 0 ;

		const std::vector<double>& getValues() const { return values ; }

	protected :
		unsigned int nPoints = 0 ;

		std::vector<double> xVec = {} ;
		std::vector<double> values = {} ;
		std::vector<double> lowerErrors = {} ;
		std::vector<double> upperErrors = {} ;
} ;

#endif // Fitter_h
