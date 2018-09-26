#include "MultiplicityFitter.h"

#include <iostream>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/Functor.h>

MultiplicityFitter::MultiplicityFitter()
	: fitResult()
{
	nParam = 3 ;
	param = std::vector<double>(nParam , 0) ;
	setParams() ;
}

MultiplicityFitter::~MultiplicityFitter()
{
	deletePoints() ;
}

void MultiplicityFitter::deletePoints()
{

}

double MultiplicityFitter::baseFunc(const double* x , const double* params)
{
	return mulVsThr(*x , params) ;
}

double MultiplicityFitter::mulVsThr(double x , const double* params)
{
	return params[0]*std::pow(x , params[1]) + params[2] ;
}

void MultiplicityFitter::getPoints(TGraphErrors* graph)
{
	nPoints = static_cast<unsigned int>( graph->GetN() ) ;

	q = std::vector<double>( nPoints , 0 ) ;
	values = std::vector<double>( nPoints , 0 ) ;

	errors = std::vector<double>( nPoints , 0 ) ;

	//if an error == 0, it will 'destroy' the fit

	for ( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		q[i] = graph->GetX()[i] ;
		values[i] = graph->GetY()[i] ;

		if ( graph->GetErrorY( static_cast<int>(i) ) < zeroLimit )
		{
			std::cout << "errlow" << i << " = 0 , set to 1e-12 instead" << std::endl ;
			errors[i] = 1e-12 ;
		}
		else
			errors[i] = graph->GetErrorYlow( static_cast<int>(i) ) ;
	}
}

void MultiplicityFitter::setParams(const double* params)
{
	if (!params)
	{
		param[0] = 0.4 ;
		param[1] = -0.3 ;
		param[2] = 1 ;

		return ;
	}

	for ( unsigned int i = 0 ; i < nParam ; i++ )
		param[i] = params[i] ;
}

double MultiplicityFitter::functionToMinimize(const double* params)
{
	double chi2 = 0.0 ;

	for( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		double x = q[i] ;
		double fitValue = mulVsThr(x , params) ;

		double error = errors[i] ;

		chi2 += (values[i] - fitValue)*(values[i] - fitValue)/(error*error) ;
	}

	chi2 /= (nPoints-1) ;
	return chi2 ;
}

void MultiplicityFitter::minimize()
{
//	for ( const auto i : q )
//		std::cout << i << " "  ;
//	std::cout << std::endl ;
//	for ( const auto i : values )
//		std::cout << i << " "  ;
//	std::cout << std::endl ;
//	for ( const auto i : errors )
//		std::cout << i << " "  ;
//	std::cout << std::endl ;

	ROOT::Minuit2::Minuit2Minimizer min ;
	min.SetMaxFunctionCalls(1000000) ;
	min.SetMaxIterations(100000) ;
	min.SetTolerance(1e-2) ;
	min.SetPrintLevel(0) ;
	double step = 0.1 ;

	ROOT::Math::Functor f(this , &MultiplicityFitter::functionToMinimize , nParam) ;

	min.SetFunction(f) ;

	int minimizerStatus = std::numeric_limits<int>::max() ;
	int nTry = 0 ;

	if ( min.PrintLevel() > 0)
		std::cout << "Init Params : " << param[0] << " , " << param[1] << " , "	<< param[2] << std::endl ;

	while ( minimizerStatus > 1 && nTry < 10 )
	{
		if ( nTry>0 )
			std::cout << "Minimisation did not converge : another try" << std::endl ;

				min.SetLowerLimitedVariable(0 , "factor" , param[0] , step , 0) ;
				min.SetUpperLimitedVariable(1 , "power" , param[1] , step , 0) ;
				min.SetLowerLimitedVariable(2 , "constant" , param[2] , step , 0) ;

//		min.SetLimitedVariable(0 , "factor" , param[0] , step , 0 , 2) ;
//		min.SetLimitedVariable(1 , "power" , param[1] , step , -1 , 0) ;
//		min.SetLimitedVariable(2 , "constant" , param[2] , step , 0 , 2) ;

		min.Minimize() ;
		minimizerStatus = min.Status() ;

		const double* xs = min.X() ;
		const double* xsErr = min.Errors() ;
		double chi2 = functionToMinimize(xs) ;

		for(unsigned int i = 0 ; i < nParam ; i++)
			param[i] = xs[i] ;

		fitResult = MulFitResult(param[0] , param[1] , param[2] , chi2 , xsErr[0] , xsErr[1] , xsErr[2] , minimizerStatus) ;

		nTry++ ;
	}
}
