#include "PolyaFitter.h"

#include <iostream>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/Functor.h>

PolyaFitter::PolyaFitter()
{
	nParam = 3 ;
	//	nParam = 4 ;
	q = NULL ;
	values = NULL ;
	errors = NULL ;

	param = new double[nParam] ;
	setParams() ;
}

PolyaFitter::~PolyaFitter()
{
	if (param)
		delete[] param ;
	deletePoints() ;
}

void PolyaFitter::deletePoints()
{
	if (q)
		delete[] q ;
	if (values)
		delete[] values ;
	if (errors)
		delete[] errors ;

	q = NULL ;
	values = NULL ;
	errors = NULL ;
}

double PolyaFitter::baseFunc(double* x , const double* params)
{
	return polyaEff(*x , params) ;
}

double PolyaFitter::polyaEff(double x , const double* params)
{
	//	double alpha = params[1] + 1 ;
	//	double theta = params[0]/(params[1] + 1) ;

	//	return params[2]*ROOT::Math::gamma_cdf_c(x , alpha , theta , 0.0) ;
	//	return params[2] - params[3]*ROOT::Math::gamma_cdf(x , alpha , theta , 0.0) ;


	double alpha = params[0]/params[1] ;
	double delta = params[1] ;

	return params[2]*ROOT::Math::gamma_cdf_c(x , alpha , delta , 0.0) ;
}

void PolyaFitter::getPoints(TGraphErrors* graph)
{
	nPoints = static_cast<unsigned int>( graph->GetN() ) ;

	q = new double[nPoints] ;
	values = new double[nPoints] ;
	errors = new double[nPoints] ;

	//if an error == 0, it will 'destroy' the fit, so if an error == 0, i put the minimum error of all other points instead
	double minError = -1 ;
	for ( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		if ( minError < 0 && graph->GetEY()[i] > 0 )
		{
			minError = graph->GetEY()[i] ;
			continue ;
		}

		if ( graph->GetEY()[i] < minError )
			minError = graph->GetEY()[i] ;
	}

	if (minError < 0)
		minError = 1.0 ;

	for ( unsigned int i = 0 ; i < nPoints ; i++ )
	{
		q[i] = graph->GetX()[i] ;
		values[i] = graph->GetY()[i] ;

		if ( graph->GetEY()[i] < 1e-100 )
		{
			std::cout << "err" << i << " = 0 , set to " << minError << " instead" << std::endl ;
			errors[i] = minError ;
		}
		else
			errors[i] = graph->GetEY()[i] ;
	}
}

void PolyaFitter::setParams(const double* params)
{
	if (!params)
	{
		param[0] = 4.5 ;
		param[1] = 2 ;

		if (values)
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
		if ( errors[i] < 1e-100 )
		{
			std::cout << "err" << i << " = 0 , set to 1e-12" << std::endl ;
			errors[i] = 1e-12 ;
		}

		double x = q[i] ;
		double fitValue = polyaEff(x , params) ;
		chi2 += (values[i] - fitValue)*(values[i] - fitValue)/(errors[i]*errors[i]) ;
	}

	chi2 /= (nPoints-1) ;
	return chi2 ;
}

void PolyaFitter::minimize()
{
	ROOT::Minuit2::Minuit2Minimizer min ;
	min.SetMaxFunctionCalls(1000000) ;
	min.SetMaxIterations(100000) ;
	min.SetTolerance(1e-8) ;
	min.SetPrintLevel(0) ;
	double step = 0.1 ;

	ROOT::Math::Functor f(this , &PolyaFitter::functionToMinimize , nParam) ;

	min.SetFunction(f) ;

	int minimizerStatus = -1 ;
	int nTry = 0 ;

	if ( min.PrintLevel() > 0)
		std::cout << "Init Params : " << param[0] << " , " << param[1] << " , "	<< param[2] << std::endl ;
	//			std::cout << "Init Params : " << param[0] << " , " << param[1] << " , "	<< param[2] << " , " << param[3] << std::endl ;

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
