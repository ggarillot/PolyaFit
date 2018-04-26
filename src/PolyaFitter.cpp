#include "PolyaFitter.h"

#include <iostream>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/Functor.h>

PolyaFitter::PolyaFitter()
	: fitResult()
{
	nParam = 3 ;
	param = std::vector<double>(nParam , 0) ;
	setParams() ;
}

PolyaFitter::~PolyaFitter()
{
	deletePoints() ;
}

void PolyaFitter::deletePoints()
{

}

double PolyaFitter::baseFunc(const double* x , const double* params)
{
	return polyaEff(*x , params) ;
}

double PolyaFitter::polyaEff(double x , const double* params)
{
	double alpha = params[0]/params[1] ;
	double delta = params[1] ;

	return params[2]*ROOT::Math::gamma_cdf_c(x , alpha , delta , 0.0) ;
}

void PolyaFitter::getPoints(TGraphAsymmErrors* graph)
{
	nPoints = static_cast<unsigned int>( graph->GetN() ) ;

	q = std::vector<double>( nPoints , 0 ) ;
	values = std::vector<double>( nPoints , 0 ) ;

	lowerBound = std::vector<double>( nPoints , 0 ) ;
	upperBound = std::vector<double>( nPoints , 0 ) ;

	//if an error == 0, it will 'destroy' the fit

	for ( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		q[i] = graph->GetX()[i] ;
		values[i] = graph->GetY()[i] ;

		if ( graph->GetErrorYlow( static_cast<int>(i) ) < zeroLimit && values[i] > zeroLimit )
		{
			std::cout << "errlow" << i << " = 0 , set to 1e-12 instead" << std::endl ;
			lowerBound[i] = 1e-12 ;
		}
		else
			lowerBound[i] = graph->GetErrorYlow( static_cast<int>(i) ) ;


		if ( graph->GetErrorYhigh( static_cast<int>(i) ) < zeroLimit && (1.0-values[i]) > zeroLimit )
		{
			std::cout << "errhigh" << i << " = 0 , set to 1e-12 instead" << std::endl ;
			upperBound[i] = 1e-12 ;
		}
		else
			upperBound[i] = graph->GetErrorYhigh( static_cast<int>(i) ) ;
	}
}

void PolyaFitter::setParams(const double* params)
{
	if (!params)
	{
		param[0] = 4.5 ;
		param[1] = 2 ;

		if ( ! values.empty() )
		{
			param[2] = values[0] ;

			if (param[2] < 0.5)
				param[1] = 0 ;
		}
		else
			param[2] = 0.96 ;

		//		param[3] = param[2] ;
		return ;
	}

	for( unsigned int i = 0 ; i < nParam ; i++ )
		param[i] = params[i] ;
}

double PolyaFitter::functionToMinimize(const double* params)
{
	double chi2 = 0.0 ;

	for( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		if ( lowerBound[i] < 1e-100 && values[i] > zeroLimit )
		{
			std::cout << "errlow" << i << " = 0 , set to 1e-12" << std::endl ;
			lowerBound[i] = 1e-12 ;
		}

		if ( upperBound[i] < 1e-100 && (1.0-values[i]) > zeroLimit )
		{
			std::cout << "errhigh" << i << " = 0 , set to 1e-12" << std::endl ;
			upperBound[i] = 1e-12 ;
		}

		double x = q[i] ;
		double fitValue = polyaEff(x , params) ;

		double error = 0 ;
		if ( fitValue < values[i] )
			error = lowerBound[i] ;
		else
			error = upperBound[i] ;

		if ( std::abs(fitValue - 1) < zeroLimit )
			error = lowerBound[i] ;

		chi2 += (values[i] - fitValue)*(values[i] - fitValue)/(error*error) ;
	}

	chi2 /= (nPoints-1) ;
	return chi2 ;
}

void PolyaFitter::minimize()
{
	ROOT::Minuit2::Minuit2Minimizer min ;
	min.SetMaxFunctionCalls(1000000) ;
	min.SetMaxIterations(100000) ;
	min.SetTolerance(1e-4) ;
	min.SetPrintLevel(0) ;
	double step = 0.1 ;

	ROOT::Math::Functor f(this , &PolyaFitter::functionToMinimize , nParam) ;

	min.SetFunction(f) ;

	int minimizerStatus = -1 ;
	int nTry = 0 ;

	if ( min.PrintLevel() > 0)
		std::cout << "Init Params : " << param[0] << " , " << param[1] << " , "	<< param[2] << std::endl ;

	while ( minimizerStatus != 0 && nTry < 10 )
	{
		if ( nTry>0 )
			std::cout << "Minimisation did not converge : another try" << std::endl ;

		min.SetLowerLimitedVariable(0 , "qbar" , param[0] , step , 0) ;
		min.SetLowerLimitedVariable(1 , "delta" , param[1] , step , 0) ;
		min.SetLimitedVariable(2 , "e0" , param[2] , step , 0 , 1) ;
		//		min.SetLimitedVariable(3 , "c" , param[3] , step , 0 , 1) ;

		min.Minimize() ;
		minimizerStatus = min.Status() ;

		const double* xs = min.X() ;
		const double* xsErr = min.Errors() ;
		double chi2 = functionToMinimize(xs) ;

		for(unsigned int i = 0 ; i < nParam ; i++)
			param[i] = xs[i] ;

		fitResult = PolyaFitResult(param[0] , param[1] , param[2] , chi2 , xsErr[0] , xsErr[1] , xsErr[2] , minimizerStatus) ;

		nTry++ ;
	}
}
