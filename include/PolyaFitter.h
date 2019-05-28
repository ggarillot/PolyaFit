#ifndef PolyaFitter_h
#define PolyaFitter_h

#include "Fitter.h"

#include <TGraphAsymmErrors.h>
#include <iostream>
#include <vector>

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

class PolyaFitter : public Fitter
{
	public :

		PolyaFitter() ;
		~PolyaFitter() = default ;

		static double polyaVsThr(const double* x , const double* params) ;
		double eval(double x , const double* params) ;


		PolyaFitResult getFitResult() const { return fitResult ; }

		void createResultStruct() ;

	protected :
		PolyaFitResult fitResult = {} ;
} ;



#endif //PolyaFitter_h
