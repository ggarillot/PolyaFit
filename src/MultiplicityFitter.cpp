#include "MultiplicityFitter.h"

#include <iostream>

#include <Minuit2/Minuit2Minimizer.h>
#include <Math/ProbFuncMathCore.h>
#include <Math/Functor.h>

MultiplicityFitter::MultiplicityFitter()
	: Fitter(3)
{
	setParams( {0.4,-0.3,1.0} ) ;
	parametersNames = {"f" , "p" , "c"} ;
	parametersLimits = { {0.0,maxValue} , {minValue,0.0} , {1.0,maxValue} } ;

	precision = 1e-2 ;
}

double MultiplicityFitter::multiplicityVsThr(const double* x , const double* params)
{
	return params[0]*std::pow(x[0] , params[1]) + params[2] ;
}

double MultiplicityFitter::eval(double x , const double* params)
{
	double xt[1] = {x} ;
	return multiplicityVsThr(xt,params) ;
}

void MultiplicityFitter::createResultStruct()
{
	fitResult = MultiplicityFitResult{parameters[0] , parameters[1] , parameters[2] , chi2 , parametersErrors[0] , parametersErrors[1] , parametersErrors[2] , minimStatus} ;
}
