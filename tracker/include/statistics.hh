#include "TMinuit.h"
#include "physics.hh"
#include "LinearAlgebra.hh"
#include <iostream>
#include <fstream>
#include <math.h>
#include <TMath.h>

#ifndef STAT_HH
#define STAT_HH

using par = std::pair<int, int>;

class TrackFitter
{
public:
	//pars is x0, y0, z0, t0, vx, vy, vz
	static void chi2_error(int &npar, double *gin, double &f, double *par, int iflag);
	static void timeless_chi2_error(int &npar, double *gin, double &f, double *par, int iflag);
	static std::vector<physics::digi_hit *> digi_list;
	static std::vector<double> parameters;
	static std::vector<double> parameter_errors;
	const static int npar = 7;
	static double cov_matrix[npar][npar];

	bool fit(std::vector<physics::digi_hit *> _digi_list, std::vector<double> arg_guess = {})
	{
		digi_list = _digi_list;
		parameters.resize(npar);
		parameter_errors.resize(npar);

		std::vector<double> guess = arg_guess;

		TMinuit minimizer(npar);
		int ierflg = 0;
		minimizer.SetFCN(chi2_error);

		double first_step_size = 0.1;
		auto maxcalls = 500000.0;
		auto tolerance = 0.1;
		double arglist[2];
		arglist[0] = maxcalls;
		arglist[1] = tolerance;

		int quiet_mode = -1;
		int normal = 0;

		minimizer.SetPrintLevel(quiet_mode);

		double vel_y_direction = arg_guess[4]; //use the direction of the seed to help fit the track

		//we find the lowest t, and take the associated y-value and fix

		double y0_fix;

		double min_t = 9999999.99;
		for (auto hit : digi_list)
		{
			if (hit->t < min_t)
			{
				min_t = hit->t;
				y0_fix = hit->y;
			}
		}

		minimizer.mnparm(0, "x0", guess[0], first_step_size, detector::x_min, detector::x_max, ierflg);
		minimizer.mnparm(1, "y0", y0_fix, first_step_size, 0, 0, ierflg);
		minimizer.mnparm(2, "z0", guess[2], first_step_size, detector::z_min_wall, detector::z_max, ierflg);
		minimizer.mnparm(3, "vx", guess[3], first_step_size, -1.0 * constants::c, constants::c, ierflg);
		if (vel_y_direction >= 0)
			minimizer.mnparm(4, "vy", guess[4], first_step_size, 0.0, constants::c, ierflg);
		if (vel_y_direction < 0)
			minimizer.mnparm(4, "vy", guess[4], first_step_size, -1.0 * constants::c, 0.0, ierflg);
		minimizer.mnparm(5, "vz", guess[5], first_step_size, -1.0 * constants::c, constants::c, ierflg);
		minimizer.mnparm(6, "t0", guess[6], first_step_size, 0.0, pow(10.0, 10.0), ierflg);

		minimizer.FixParameter(1);

		minimizer.mnexcm("MIGRAD", arglist, 2, ierflg);

		for (int k = 0; k < npar; k++)
		{
			minimizer.GetParameter(k, parameters[k], parameter_errors[k]);
		}

		minimizer.mnemat(&cov_matrix[0][0], npar);

		for (int k = 0; k < parameter_errors.size(); k++)
			parameter_errors[k] = TMath::Sqrt(cov_matrix[k][k]);

		//GETTTING STATUS:
		double fmin = 0.0;
		double fedm = 0.0;
		double errdef = 0.0;
		int npari = 0;
		int nparx = 0;
		/////////////////////////////////////////////////////
		int istat = 0; //this is the one we really care about
		// 0 - no covariance
		// 1 - not accurate, diagonal approximation
		// 2 - forced positive definite
		// 3 - full accurate matrix--succesful convergence

		minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);

		return (istat >= 2) ? true : false;
	}

	double rand_guess()
	{
		return 0.0;
	}

}; //class TrackFinder

class VertexFitter
{
public:
	static void nll(int &npar, double *gin, double &f, double *par, int iflag); //NEGATIVE LOG LIKLIHOOD FOR FITTING VERTICES
	
	static void cost_adaptive(int &npar, double *gin, double &f, double *par, int iflag); //NEGATIVE LOG LIKLIHOOD FOR FITTING VERTICES
	static double calc_adaptive_weight(double chi2){
		auto a = TMath::Exp(-chi2/(2*adaptive_temperature));
		auto b = TMath::Exp(-adaptive_chi2_cutoff/(2*adaptive_temperature));
		return (a/(a+b));
	}
	static int adaptive_iterations; //	int adaptive_iterations = 10;
	static double adaptive_annealing_factor; //	double adaptive_annealing_factor = 0.8;
	static double adaptive_temperature; //	double adaptive_temperature = 3;
	static double adaptive_chi2_cutoff; //	double adaptive_chi2_cutoff = 60;

	static std::vector<physics::track *> track_list;
	static std::vector<double> parameters;
	static std::vector<double> parameter_errors;
	const static int npar = 4;
	static bool bad_fit;
	static double cov_matrix[npar][npar];
	static double _merit;
	double ndof;
	std::vector<double> delta_chi2_list;


	double merit(){
		return _merit / ndof;
	}

	std::vector<double> delta_chi2()
	{
		delta_chi2_list.clear();

		for (auto track : track_list)
		{
			delta_chi2_list.push_back(track->vertex_residual(parameters));
		}

		// double ndof = static_cast<double>(npar * (track_list.size() - 1));
		// double ndof = static_cast<double>(3 * (track_list.size())- npar);
		// ndof = ndof > 1. ? ndof : 1.0;
		// std::cout<<"   NDOF"<<ndof<<std::endl;

		return delta_chi2_list;
	}


	// Fit, with adjustable tolerance. IMINUIT DOC: The default tolerance is 0.1, and the minimization will stop when the estimated vertical distance to the minimum (EDM) is less than 0.001*[tolerance]*UP
	bool fit(std::vector<physics::track *> _track_list, std::vector<double> arg_guess = {}, double tolerance=0.1)
	{

		track_list = _track_list;
		ndof = static_cast<double>(3 * (track_list.size())- npar);
		ndof = ndof > 1. ? ndof : 1.0;

		bad_fit = false;

		parameters.resize(npar);
		parameter_errors.resize(npar);

		std::vector<double> guess = arg_guess;

		TMinuit minimizer(npar);
		int ierflg = 0;
		minimizer.SetFCN(nll);

		double first_step_size = 0.010;
		auto maxcalls = 500000.0;
		// auto tolerance = 0.1;
		double arglist[2];
		arglist[0] = maxcalls;
		arglist[1] = tolerance;

		int quiet_mode = -1;
		int normal = 0;

		minimizer.SetPrintLevel(quiet_mode);

		double min_y = detector::y_min - 10.0 * units::cm;
		double max_y = detector::y_max + 10.0 * units::cm;

		minimizer.mnparm(0, "x", guess[0], 1.0, 0, 0, ierflg);
		minimizer.mnparm(1, "y", guess[1], 1.0, 0, 0, ierflg);
		minimizer.mnparm(2, "z", guess[2], 1.0, 0, 0, ierflg);
		minimizer.mnparm(3, "t", guess[3], first_step_size, 0, 0, ierflg);

		minimizer.mnexcm("MIGRAD", arglist, 2, ierflg);

		//GETTTING STATUS:
		double fmin = 0.0;
		double fedm = 0.0;
		double errdef = 1; //0: NLL, 1:Least-square
		int npari = 0;
		int nparx = 0;
		int istat = 0; //this is the one we really care about

		minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);
		_merit = fmin;

		//while (ierflg) minimizer.mnexcm("MIGRAD", arglist ,2,ierflg);
		for (int k = 0; k < npar; k++)
		{
			minimizer.GetParameter(k, parameters[k], parameter_errors[k]);
		}

		minimizer.mnemat(&cov_matrix[0][0], npar);

		return (istat >= 2) ? true : false;
	}

	// Fit that calls the adaptive cost function
	bool fit_adaptive(std::vector<physics::track *> _track_list, std::vector<double> arg_guess = {}, double tolerance=0.1)
	{
		std::cout <<"  -Adaptive fit: initial parameters: " << arg_guess[0] << ", " << arg_guess[1] << ", " << arg_guess[2] << ", " << arg_guess[3] << std::endl;
		bool fit_status = fit_adaptive_oneiteration(_track_list, arg_guess, tolerance);
		std::vector<double> parameter_pre = parameters;
		double total_interations = adaptive_iterations;
		while (adaptive_iterations>0)
		{
			adaptive_iterations-=1;
			adaptive_temperature = adaptive_temperature*adaptive_annealing_factor;
			fit_status = fit_adaptive_oneiteration(_track_list, parameter_pre, tolerance);
			if (fit_status==false){
				break;
				std::cout << " -Adaptive fit failed" << std::endl;
			}
			if (TMath::Abs(parameters[0]-parameter_pre[0])<0.1 &&TMath::Abs(parameters[1]-parameter_pre[1])<0.1 &&TMath::Abs(parameters[2]-parameter_pre[2])<0.1 && TMath::Abs(parameters[3]-parameter_pre[3])<0.1){
				break;
			}
			parameter_pre = parameters;
			// if (par_handler->par_map["debug_vertex"]==1)
				std::cout <<"  -Adaptive fit iteration " << total_interations-adaptive_iterations << ", parameters: " << parameters[0] << ", " << parameters[1] << ", " << parameters[2] << ", " << parameters[3] << std::endl;
		}

		return fit_status;
	}

	bool fit_adaptive_oneiteration(std::vector<physics::track *> _track_list, std::vector<double> arg_guess = {}, double tolerance=0.1)
	{

		track_list = _track_list;
		// Don't update ndof. It is already taken care of in the cost function
		// ndof = static_cast<double>(3 * (track_list.size())- npar);
		// ndof = ndof > 1. ? ndof : 1.0;

		bad_fit = false;

		parameters.resize(npar);
		parameter_errors.resize(npar);

		std::vector<double> guess = arg_guess;

		TMinuit minimizer(npar);
		int ierflg = 0;
		minimizer.SetFCN(cost_adaptive);

		double first_step_size = 0.010;
		auto maxcalls = 500000.0;
		// auto tolerance = 0.1;
		double arglist[2];
		arglist[0] = maxcalls;
		arglist[1] = tolerance;

		int quiet_mode = -1;
		int normal = 0;

		minimizer.SetPrintLevel(quiet_mode);

		double min_y = detector::y_min - 10.0 * units::cm;
		double max_y = detector::y_max + 10.0 * units::cm;

		minimizer.mnparm(0, "x", guess[0], 1.0, 0, 0, ierflg);
		minimizer.mnparm(1, "y", guess[1], 1.0, 0, 0, ierflg);
		minimizer.mnparm(2, "z", guess[2], 1.0, 0, 0, ierflg);
		minimizer.mnparm(3, "t", guess[3], first_step_size, 0, 0, ierflg);

		minimizer.mnexcm("MIGRAD", arglist, 2, ierflg);

		//GETTTING STATUS:
		double fmin = 0.0;
		double fedm = 0.0;
		double errdef = 1; //0: NLL, 1:Least-square
		int npari = 0;
		int nparx = 0;
		int istat = 0; //this is the one we really care about

		minimizer.mnstat(fmin, fedm, errdef, npari, nparx, istat);
		_merit = fmin;

		//while (ierflg) minimizer.mnexcm("MIGRAD", arglist ,2,ierflg);
		for (int k = 0; k < npar; k++)
		{
			minimizer.GetParameter(k, parameters[k], parameter_errors[k]);
		}

		minimizer.mnemat(&cov_matrix[0][0], npar);

		return (istat >= 2) ? true : false;
	}	

}; //class VertexFitter

class Stat_Funcs
{
public:

	double chi_prob(double chi2, double ndof) {
		/* P-value for a chi squared approximation that is computationally efficient
		formula is eqn (8) of the following paper

		Beh, E. (2018). Exploring How to Simply Approximate the P-value of a
		Chi-squared Statistic. Austrian Journal of Statistics, 47(3), 63–75.
		https://doi.org/10.17713/ajs.v47i3.757 */

		std::vector<double> c = {-1.37266,
					1.06807,
					2.13161,
					-0.04859}; // fitting constants

		double chi = std::sqrt(chi2);
		double nu = std::sqrt(ndof);

		if ( chi < (c[0] + c[1] * nu) ) return 1;

		else {
			double expon = ( chi - (c[0] + c[1] * nu) )
					/ (c[2] + c[3] * nu);

			return std::pow( 0.1 , expon*expon );
		};


	}

	double chi_prob_eld(double chi2, double ndof) {
		/* P-value for a chi squared approximation, from the following paper

		Elderton, W. P. (1902). Tables for Testing the Goodness of Fit of Theory to Observation.
		Biometrika, 1(2), 155–163. https://doi.org/10.2307/2331485

		assumes all ndof <= 30, since ndof = 4 * hits - 6 and hits <= 9, we're good!*/

		double sum = 0;
		double chi = std::sqrt(chi2);

		if (ndof == 1) return 1 - erf(std::sqrt(2) * chi);

		bool odd = ((int)(ndof) % 2 == 1);

		for (int i = 1 ? odd : 0; i < ndof / 2; i++) {

			sum += std::pow(chi, 2 * i) * std::pow(2, i) * tgamma(i + 1) / tgamma(2.0 * i + 1); // chi^i / (2 i - 1)!!

		};

		sum *= exp(- chi2 / 2);

		if (odd) {
			sum *= std::sqrt(2 / M_PI);
			sum += 1 - erf(std::sqrt(2) * chi);
		};
		return sum;

	}

};

#endif
