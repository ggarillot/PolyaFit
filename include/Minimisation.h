#ifndef Minimisation_h
#define Minimisation_h

#include <vector>
#include <string>
#include <limits>

class Minimisation
{
	public :
		Minimisation(unsigned int nParam_) ;
		virtual ~Minimisation() = default ;

		virtual double functionToMinimize(const double* param) = 0 ;
		double chi2Calc(const std::vector<double>& param) ;
		virtual double minimize() ;

		virtual void printParam() const ;

		std::vector<double> getParameters() const { return parameters ; }

		void setParam(unsigned int i , double value) ;
		void setParams(std::vector<double> values) ;

		void resetParams() ;

		void setPrecision(double prec) { precision = prec ; }
		void setStep(double s) { step = s ; }
		void setPrintLevel(int p) { printLevel = p ; }
		void setNTries(int n) { nTries = n ; }

		static constexpr double zeroLimit = std::numeric_limits<double>::epsilon() ;

		static constexpr double minValue = -5e6 ; //for MINUIT
		static constexpr double maxValue = 5e6 ; //for MINUIT

	protected :
		unsigned int nParameters ;
		std::vector<double> parameters = {} ;
		std::vector<double> parametersErrors = {} ;

		int minimStatus = -1 ;
		double chi2 = maxValue ;

		std::vector< std::pair<double,double> > parametersLimits = {} ;
		std::vector<std::string> parametersNames = {} ;

		double step = 0.1 ;
		double precision = 1e-4 ;
		int nTries = 10 ;

		int printLevel = 0 ;
} ;

#endif // Minimisation_h
