#ifndef MultiplicityFitter_h
#define MultiplicityFitter_h

#include <TGraphAsymmErrors.h>
#include <TGraphErrors.h>
#include <iostream>
#include <vector>

class MultiplicityFitter
{
	public :

		MultiplicityFitter() ;
		~MultiplicityFitter() ;
		static double baseFunc(const double* x , const double* params) ;
		static double mulVsThr(double x , const double* params) ;

		void getPoints(TGraphErrors* graph) ;

		double functionToMinimize(const double* params) ;

		void minimize() ;

		void setParams(const double* params = nullptr) ;

		std::vector<double> getParams() const { return param ; }

		struct MulFitResult
		{
				MulFitResult()
					: f(0) , p(0) , c(0) , chi2(0) , fErr(0) , pErr(0) , cErr(0) , minimStatus(0)
				{}
				MulFitResult(double factor , double power , double cons , double c2 , double fE , double pE , double cE , int m)
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

		MulFitResult getFitResult() const { return fitResult ; }

		static constexpr double zeroLimit = std::numeric_limits<double>::epsilon() ;

	protected :

		unsigned int nPoints = 0 ;

		std::vector<double> q = {} ;
		std::vector<double> values = {} ;
		std::vector<double> errors = {} ;

		unsigned int nParam = 0 ;
		std::vector<double> param = {} ;

		MulFitResult fitResult ;

		void deletePoints() ;

} ;

#endif //MultiplicityFitter_h
