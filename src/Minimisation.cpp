#include "Minimisation.h"

#include <iostream>
#include <sstream>
#include <cmath>
#include <cassert>
#include <limits>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/GSLMinimizer.h>
#include <Math/GSLSimAnMinimizer.h>
#include <Math/Functor.h>

constexpr double Minimisation::minValue ;
constexpr double Minimisation::maxValue ;

Minimisation::Minimisation(unsigned int nParam_)
	: nParameters(nParam_)
{
	parameters = std::vector<double>(nParameters , 0) ;
	parametersErrors = std::vector<double>(nParameters , 0) ;

	parametersLimits = std::vector< std::pair<double,double> >( nParameters , {minValue,maxValue} ) ;
	parametersNames = std::vector<std::string>(nParameters , "") ;

	for ( unsigned int i = 0 ; i < nParameters ; ++i )
	{
		std::stringstream toto ; toto << "param" << i ;
		parametersNames[i] = toto.str() ;
	}
}

double Minimisation::chi2Calc(const std::vector<double>& param)
{
	return functionToMinimize( &param[0] ) ;
}


double Minimisation::minimize()
{
	assert( parameters.size() == parametersLimits.size() ) ;
	assert( parameters.size() == nParameters ) ;
	assert( parameters.size() == parametersErrors.size() ) ;

	if ( printLevel > 0 )
	{
		std::cout << "Launching minimizing" << std::endl ;
		double prevMin = chi2Calc(parameters) ;
		std::cout << "prevMin : " << prevMin << std::endl ;
	}

	ROOT::Minuit2::Minuit2Minimizer min ;
	min.SetMaxFunctionCalls(400000) ;
	min.SetMaxIterations(2) ;
	min.SetTolerance(precision) ;
	min.SetPrintLevel(printLevel) ;

	ROOT::Math::Functor f(this , &Minimisation::functionToMinimize , nParameters) ;

	min.SetFunction(f) ;

	minimStatus = -1 ;
	int nTry = 0 ;

	if ( min.PrintLevel() > 0 )
		printParam() ;

	while ( minimStatus != 0 && nTry < nTries )
	{
		if ( nTry > 0 && printLevel > 0 )
			std::cout << "Minimisation did not converge : another try" << std::endl ;

		for ( unsigned int i = 0 ; i < nParameters ; i++ )
		{
			auto minLimit = parametersLimits[i].first ;
			auto maxLimit = parametersLimits[i].second ;

			parameters[i] = std::min(parameters[i] , maxLimit) ;
			parameters[i] = std::max(parameters[i] , minLimit) ;

			min.SetLimitedVariable(i , parametersNames[i].c_str() , parameters[i] , step , minLimit , maxLimit) ;
		}

		min.Minimize() ;

		const double* xs = min.X() ;
		const double* xsErr = min.Errors() ;

		minimStatus = min.Status() ;
		chi2 = functionToMinimize(xs) ;

		for ( unsigned int i = 0 ; i < nParameters ; i++ )
		{
			parameters[i] = xs[i] ;
			parametersErrors[i] = xsErr[i] ;
		}

		nTry++ ;
	}

	return chi2Calc(parameters) ;
}

void Minimisation::printParam() const
{
	std::cout << "\"parameters\": [" ;
	for  (unsigned int i = 0 ; i < nParameters ; i++ )
	{
		std::cout << parameters[i] ;
		if ( i < nParameters - 1 )
			std::cout << "," ;
	}
	std::cout << "]" << std::endl ;
}

void Minimisation::setParam(unsigned int i , double value)
{
	assert( i < nParameters ) ;
	parameters[i] = value ;
}
void Minimisation::setParams(std::vector<double> values)
{
	assert( values.size() == nParameters ) ;
	parameters = values ;
}

void Minimisation::resetParams()
{
	for ( auto& i : parameters )
		i = 0.0 ;
}
