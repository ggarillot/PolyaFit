#ifndef MultiplicityFitter_h
#define MultiplicityFitter_h

#include "Fitter.h"

#include <TGraphErrors.h>
#include <iostream>
#include <vector>

struct MultiplicityFitResult
{
		MultiplicityFitResult()
			: f(0) , p(0) , c(0) , chi2(0) , fErr(0) , pErr(0) , cErr(0) , minimStatus(0)
		{}
		MultiplicityFitResult(double factor , double power , double cons , double c2 , double fE , double pE , double cE , int m)
			: f(factor) , p(power) , c(cons) , chi2(c2) , fErr(fE) , pErr(pE) , cErr(cE) , minimStatus(m)
		{}

		double f ;
		double p ;
		double c ;
		double chi2 ;
		double fErr ;
		double pErr ;
		double cErr ;

		int minimStatus ;

		void print() const
		{
			std::cout << "Factor :  " << f << " += " << fErr << std::endl ;
			std::cout << "Power :  " << p << " += " << pErr << std::endl ;
			std::cout << "Constante :  " << c << " += " << cErr << std::endl ;
		}
} ;

class MultiplicityFitter : public Fitter
{
	public :

		MultiplicityFitter() ;
		~MultiplicityFitter() = default ;

		static double multiplicityVsThr(const double* x , const double* params) ;
		double eval(double x , const double* params) ;


		MultiplicityFitResult getFitResult() const { return fitResult ; }

		void createResultStruct() ;

	protected :
		MultiplicityFitResult fitResult = {} ;
} ;

#endif //MultiplicityFitter_h
