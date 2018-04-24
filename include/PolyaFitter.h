#ifndef PolyaFitter_h
#define PolyaFitter_h

#include <TGraphErrors.h>
#include <iostream>
#include <vector>

class PolyaFitter
{
	public :

		PolyaFitter() ;
		~PolyaFitter() ;
		static double baseFunc(const double* x , const double* params) ;
		static double polyaEff(double x , const double* params) ;

		void getPoints(TGraphErrors* graph) ;

		double functionToMinimize(const double* params) ;

		void minimize() ;

		void setParams(const double* params = NULL) ;

		const double* getParams() const { return param ; }

		struct PolyaFitResult
		{
				PolyaFitResult()
					: qbar(0) , delta(0) , eff0(0) , chi2(0) , qbarError(0) , deltaError(0) , eff0Error(0) , minimStatus(0)
				{}
				PolyaFitResult(double q , double d , double e , double c , double qe , double de , double ee , int st)
					: qbar(q) , delta(d) , eff0(e) , chi2(c) , qbarError(qe) , deltaError(de) , eff0Error(ee) , minimStatus(st)
				{}

				double qbar ;
				double delta ;
				double eff0 ;
				double chi2 ;
				double qbarError ;
				double deltaError ;
				double eff0Error ;

				int minimStatus ;

				void print() const
				{
					std::cout << "Qbar :  " << qbar << " += " << qbarError << std::endl ;
					std::cout << "Delta :  " << delta << " += " << deltaError << std::endl ;
					std::cout << "Eff0 :  " << eff0 << " += " << eff0Error << std::endl ;
					std::cout << "Chi2 : " << chi2 << std::endl ;
				}
		} ;

		PolyaFitResult getFitResult() const { return fitResult ; }


	protected :

		unsigned int nPoints ;
		double* q ;
		double* values ;
		double* errors ;

		std::vector<double> lowerBound = {} ;
		std::vector<double> upperBound = {} ;

		unsigned int nParam ;
		double* param ;

		PolyaFitResult fitResult ;

		void deletePoints() ;

} ;

#endif //PolyaFitter_h
