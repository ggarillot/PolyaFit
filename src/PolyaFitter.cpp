#include "PolyaFitter.h"

#include <iostream>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/Functor.h>

PolyaFitter::PolyaFitter()
	: Fitter(3)
{
	setParams( {4.5,2,0.96} ) ;
	parametersNames = {"qbar" , "delta" , "eff0"} ;
	parametersLimits = { {zeroLimit,maxValue} , {zeroLimit,maxValue} , {0.0,1.0} } ;
}

double PolyaFitter::polyaVsThr(double x , const double* params)
{
	double a[1] = {x} ;
	return polyaVsThr(a,params) ;
}

double PolyaFitter::polyaVsThr(const double* x , const double* params)
{
	double alpha = params[0]/params[1] ;
	double delta = params[1] ;

	return params[2]*ROOT::Math::gamma_cdf_c(x[0] , alpha , delta , 0.0) ;
}

double PolyaFitter::eval(double x , const double* params)
{
	double xt[1] = {x} ;
	return polyaVsThr(xt,params) ;
}

void PolyaFitter::createResultStruct()
{
	fitResult = PolyaFitResult{parameters[0] , parameters[1] , parameters[2] , chi2 , parametersErrors[0] , parametersErrors[1] , parametersErrors[2] , minimStatus} ;
}
