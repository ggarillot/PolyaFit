#include "Fitter.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <utility>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/GSLMinimizer.h>
#include <Math/GSLSimAnMinimizer.h>
#include <Math/Functor.h>

using namespace std ;

Fitter::Fitter(unsigned int nParam_)
	: Minimisation(nParam_)
{
}

void Fitter::getPointsFromGraph(const TGraph* graph)
{
	nPoints = static_cast<unsigned int>( graph->GetN() ) ;

	xVec = std::vector<double>( nPoints , 0 ) ;
	values = std::vector<double>( nPoints , 0 ) ;

	lowerErrors = std::vector<double>( nPoints , 0 ) ;
	upperErrors = std::vector<double>( nPoints , 0 ) ;

	//if an error == 0, it will 'destroy' the fit

	for ( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		xVec[i] = graph->GetX()[i] ;
		values[i] = graph->GetY()[i] ;

		if ( graph->GetErrorYlow( static_cast<int>(i) ) < zeroLimit && values[i] > zeroLimit )
		{
			std::cout << "errlow" << i << " = 0 , set to 1e-12 instead" << std::endl ;
			lowerErrors[i] = 1e-12 ;
		}
		else
			lowerErrors[i] = graph->GetErrorYlow( static_cast<int>(i) ) ;


		if ( graph->GetErrorYhigh( static_cast<int>(i) ) < zeroLimit && (1.0-values[i]) > zeroLimit )
		{
			std::cout << "errhigh " << i << " = 0 , set to 1e-12 instead" << std::endl ;
			upperErrors[i] = 1e-12 ;
		}
		else
			upperErrors[i] = graph->GetErrorYhigh( static_cast<int>(i) ) ;
	}
}

double Fitter::functionToMinimize(const double* params)
{
	double chi2_ = 0.0 ;

	for( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		double x = xVec[i] ;
		double fitValue = eval(x , params) ;

		double error = 0 ;
		if ( fitValue < values[i] )
			error = lowerErrors[i] ;
		else
			error = upperErrors[i] ;

		if ( std::abs(fitValue - 1) < zeroLimit ) //if fit == 1, take lowError
			error = lowerErrors[i] ;
		if ( std::abs(fitValue) < zeroLimit ) //if fit == 0, take upError
			error = upperErrors[i] ;

		chi2_ += (values[i] - fitValue)*(values[i] - fitValue)/(error*error) ;
	}

	chi2_ /= (nPoints-1) ;
	return chi2_ ;
}

double Fitter::minimize()
{
	auto val = Minimisation::minimize() ;
	createResultStruct() ;
	return val ;
}
