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
					: p(0) , c(0) , chi2(0)  , pErr(0) , cErr(0) , minimStatus(0)
				{}
				MulFitResult(double power , double cons , double c2 , double pE , double cE , int m)
					: p(power) , c(cons) , chi2(c2) , pErr(pE) , cErr(cE) , minimStatus(m)
				{}

				double p ;
				double c ;
				double chi2 ;
				double pErr ;
				double cErr ;

				int minimStatus ;

				void print() const
				{
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
		std::vector<double> lowerBound = {} ;
		std::vector<double> upperBound = {} ;

		unsigned int nParam = 0 ;
		std::vector<double> param = {} ;

		MulFitResult fitResult ;

		void deletePoints() ;

} ;

#endif //MultiplicityFitter_h
