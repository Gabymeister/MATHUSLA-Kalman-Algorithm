#include<iostream>
#include <TLorentzVector.h>
#include "Geometry.hh"
#include "globals.hh"
#include "TFitter.h"
#include "TMatrix.h"
#include "TMatrixD.h"
#include "LinearAlgebra.hh"
#include "physics.hh"
#include <Eigen/Dense>


using Vector = vector::Vector;

namespace physics {


	track::track(std::vector<double> params, std::vector<double> par_errors){
			x0 = params[0]; ex0 = par_errors[0];
			y0 = params[1]; ey0 = par_errors[1];
			z0 = params[2]; ez0 = par_errors[2];
			vx = params[3]; evx = par_errors[3];
			vy = params[4]; evy = par_errors[4];
			vz = params[5]; evz = par_errors[5];
			if (params.size() == 7) {
				t0 = params[6];
				et0 = par_errors[6];
				with_t = true;
			}

	}

	std::vector<double> track::parameters(){
			std::vector<double> p = { x0, y0, z0, vx, vy, vz };
			return p;
	}

	void track::parameters(std::vector<double> pars){
			x0 = pars[0];
			y0 = pars[1];
			z0 = pars[2];
			vx = pars[3];
			vy = pars[4];
			vz = pars[5];
			if (pars.size() == 7) {
				t0 = pars[6];
				with_t = true;
			}
		}

	void track::par_errors(std::vector<double> epars){
			ex0 = epars[0];
			ey0 = epars[1];
			ez0 = epars[2];
			evx = epars[3];
			evy = epars[4];
			evz = epars[5];
			if (epars.size() == 7) {
				et0 = epars[6];
				with_t = true;
			}
		}

	double track::chi2(){
			double chi_2 = 0.0;
			double t;
			double res_t, res_x, res_z;
			for (auto hit : hits){
				t = (hit->y - y0)/vy;

				res_x = ((x0 + vx*t) - hit->x)/hit->ex;
				res_z = ((z0 + vz*t) - hit->z)/hit->ez;

				if  (!with_t) chi_2 +=  res_x*res_x + res_z*res_z;
				else {
					res_t = ((t0 + t) - hit->t)/hit->et;
					chi_2 +=  res_x*res_x + res_z*res_z + res_t*res_t;

					//std::cout << "zero pos are " << x0 << ", " << t0 << ", " << z0 << ", " << y0 << std::endl;
					//std::cout << "zero res are " << res_x << ", " << res_t << ", " << res_z << std::endl;

				}
			}

			return chi_2;
		}

    double track::chi2_per_dof(){
      	int n_track_params = 6;
      	int ndof = (4.0*hits.size() - n_track_params );
      	if (ndof < 1) ndof = 1;
      	return chi2()/ndof; //FIX ME
      }

    double track::beta(){
      	return TMath::Sqrt( vx*vx + vy*vy + vz*vz  )/constants::c;
      }

    double track::beta_err(){
        double norm = beta()*(constants::c * constants::c);
        std::vector<double> derivatives = { 0., 0., vx/norm, vy/norm, vz/norm, 0.};

      	double error = 0.0;

        for (int i = 0; i < derivatives.size(); i++){
            for (int j = 0; j < derivatives.size(); j++){

                error += derivatives[i]*derivatives[j]*_CovMatrix[i][j];

            }
        }

        return TMath::Sqrt(error);
      }


    std::vector<int> track::layers(){
    	std::vector<int> layer_indices;
    	for (auto hit : hits){
    		bool new_layer = true;
    		int layer_index = hit->det_id.layerIndex;
    		for (int layer_n : layer_indices){
    			if (layer_n == layer_index){
    				new_layer = false;
    				break;
    			}
    		}
	    	if (new_layer) layer_indices.push_back(layer_index);
    	}
        std::sort(layer_indices.begin(), layer_indices.end());
    	return layer_indices;
    } //nlayers


    int track::nlayers(){
    	//returns the number of layers that a track has hits in
    	return layers().size();
    } //nlayers


    Vector track::Position_at_Y(double y){
    	double delta_t = (y-y0)/vy;
    	return Vector(x0, y0, z0) + Vector(vx, vy, vz).Scale(delta_t);
    }

    Vector track::VelVector(){
    	return Vector(vx, vy, vz);
    }

    Vector track::direction(){
    	double velocity = beta()*constants::c;
    	return VelVector().Scale(1.0/velocity);
    }

    Vector track::P0Vector(){
    	return Vector(x0, y0, z0);
    }

    Vector track::position(double t){ //global time t
    	double dt = t-t0;
    	return P0Vector() + VelVector().Scale(dt);
    }

    double track::distance_to(Vector point, double t){
    	auto pos = position(t);
    	return (pos - point).Magnitude();
    }

    double track::err_distance_to(Vector point, double t){

    	std::vector<double> derivatives = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //partial derivates of position_at
    	//note, we don't include the derivative with respect to y. y0 is fixed, and so it isnt included in the covariance matrix
    	auto x = point.x;
    	auto y = point.y;
    	auto z = point.z;

    	double dist = distance_to(point, t);

    	// derivatives[0] = ((t-t0)*vx - x + x0)/dist;
    	// //derivatives[1] = detector::scintillator_height;
    	// derivatives[1] = ((t-t0)*vz - z + z0)/dist;
    	// derivatives[2] = (t-t0)*((t-t0)*vx - x + x0)/dist;
    	// derivatives[3] = (t-t0)*((t-t0)*vy - y + y0)/dist;
    	// derivatives[4] = (t-t0)*((t-t0)*vz - z + z0)/dist;
    	// derivatives[5] = -1.0*(vx*((t-t0)*vx - x + x0) + vy*((t-t0)*vy - y + y0) + vz*((t-t0)*vz - z + z0))/dist;

    	derivatives[0] = ((t-t0)*vx - x + x0)/dist;
    	derivatives[1] = ((t-t0)*vy - y + y0)/dist;
    	derivatives[2] = ((t-t0)*vz - z + z0)/dist;
    	derivatives[3] = (t-t0)*((t-t0)*vx - x + x0)/dist;
    	derivatives[4] = (t-t0)*((t-t0)*vy - y + y0)/dist;
    	derivatives[5] = (t-t0)*((t-t0)*vz - z + z0)/dist;
    	derivatives[6] = -1.0*(vx*((t-t0)*vx - x + x0) + vy*((t-t0)*vy - y + y0) + vz*((t-t0)*vz - z + z0))/dist;		

    	//now we calculate the actual error
    	double error = 0.0;

    	for (int i = 0; i < derivatives.size(); i++){
    		for (int j = 0; j < derivatives.size(); j++){
				// std::cout<<"i: "<<i<<" ,j: "<<j<<std::endl;
    			error += derivatives[i]*derivatives[j]*_CovMatrix[i][j];

    		}
    	}

    	return TMath::Sqrt(error);
    }

    double track::err_distance_to_mod(Vector point, double t){

    	std::vector<double> derivatives = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //partial derivates of position_at
    	std::vector<double> derivatives_time = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; //partial derivates of position_at
    	//note, we don't include the derivative with respect to y. y0 is fixed, and so it isnt included in the covariance matrix
    	auto x = point.x;
    	auto y = point.y;
    	auto z = point.z;

    	double dist = distance_to(point, t);

		auto dy = y-y0;
		auto dt = t0+dy/vy-t;
		auto err_x = x0 + vx*dy/vy -x;
		auto err_z = z0 + vz*dy/vy -z;

    	derivatives[0] = err_x/dist;
    	derivatives[1] = -(err_x*vx/vy + err_z*vz/vy)/dist;
    	derivatives[2] = err_z/dist;
    	derivatives[3] = err_x*dy/vy/dist;
    	derivatives[4] = -(err_x*vx*dy/(vy*vy) + err_z*vz*dy/(vy*vy))/dist;
    	derivatives[5] = err_z*dy/vy/dist;
    	derivatives[6] = 0;

		derivatives_time[0] = 0;
		derivatives_time[1] = -1/vy;
		derivatives_time[2] = 0;
		derivatives_time[3] = 0;
		derivatives_time[4] = -dy/vy/vy;
		derivatives_time[5] = 0;
		derivatives_time[6] = 1;

    	//now we calculate the actual error
    	double error = 0.0;
    	double error_time = 0.0;

    	for (int i = 0; i < derivatives.size(); i++){
    		for (int j = 0; j < derivatives.size(); j++){
				// std::cout<<"i: "<<i<<" ,j: "<<j<<std::endl;
    			error += derivatives[i]*derivatives[j]*_CovMatrix[i][j];
    			error_time += derivatives_time[i]*derivatives_time[j]*_CovMatrix[i][j];
    		}
    	}

		auto nll = 0.5 * (dist*dist / error) + 0.5 * TMath::Log(error)
				 + 0.5 * (dt*dt/error_time) + 0.5 * TMath::Log(error_time) 
				 + 2*0.5* TMath::Log(2*3.141592653589);

    	return nll;
    }	

	double track::chi2_distance_to(Vector point, double t){

    	auto _x = point.x;
    	auto _y = point.y;
    	auto _z = point.z;
		double dy = _y - y0;

		// Jacobian from (_x,_z,t) to (x0,y0,z0, vx,vy,vz, t)
		Eigen::MatrixXd jac;
		jac = Eigen::MatrixXd::Zero(3, 7);
		jac << 	1, -vx/vy,	0,	 	dy/vy,	-vx*dy/(vy*vy),		0,		0,
			   	0, -vz/vy,  1,			0,	-vz*dy/(vy*vy),	dy/vy,		0,
				0,  -1/vy, 	0,			0,	   -dy/(vy*vy),		0, 		1;
		// Copy the covariance matrix of (x0,y0,z0, vx,vy,vz, t) into an Eigen matrix
		Eigen::MatrixXd CovMatrix_track;
		CovMatrix_track = Eigen::MatrixXd::Zero(7, 7);
		for (int ii=0; ii<7; ii++){
			for (int jj=0; jj<7; jj++){
				CovMatrix_track(ii,jj) = _CovMatrix[ii][jj];
			}
		}		
		// Calculate the final uncertainty matrix of (_x,_z,t)
		auto CovMatrix_vertex = jac * CovMatrix_track * jac.transpose();

    	//now we calculate the Chi2
		Eigen::VectorXd residual_vector(3);
		residual_vector<<		_x - (x0+vx*dy/vy), 	_z- (z0+vz*dy/vy),		t- (t0+dy/vy);
    	auto error = residual_vector.transpose()*CovMatrix_vertex.inverse()*residual_vector;

    	return error;
    }	


    double track::vertex_residual(std::vector<double> params){

    	double x = params[0];
    	double y = params[1];
    	double z = params[2];
    	double t = params[3];

    	return chi2_distance_to(Vector(x, y, z), t);
    }


    //SEE docs/Closest_Approach_Of_Parametric_Vectors.pdf for the math
    double track::closest_approach(track* tr2){
    	using namespace vector;
        if (tr2->index == index) return 0.00;

        Vector rel_v = tr2->VelVector() - this->VelVector();
    	double rel_v2 = rel_v ^ rel_v ; //that's the dot product :/

    	Vector displacement = this->P0Vector() - tr2->P0Vector();

    	double t_ca = (  (displacement ^ rel_v) - ( (this->VelVector().Scale(this->t0) - tr2->VelVector().Scale(tr2->t0)) ^ rel_v)  )/rel_v2;

    	Vector pos1 = this->P0Vector() + this->VelVector().Scale(t_ca - this->t0);
    	Vector pos2 =  tr2->P0Vector() +  tr2->VelVector().Scale(t_ca -  tr2->t0);

    	return (pos1 - pos2).Magnitude();
    }


    //SEE docs/Closest_Approach_Of_Parametric_Vectors.pdf for the math
    Vector track::closest_approach_midpoint(track* tr2){

    	using namespace vector;
        // if (tr2->index == index) return P0Vector();
        if (tr2->index == index) {
			// std::cout<< "seed using same track"<<std::endl;
			return P0Vector();
		}

        Vector rel_v = tr2->VelVector() - this->VelVector();
    	double rel_v2 = rel_v ^ rel_v ; //that's the dot product :/

    	Vector displacement = this->P0Vector() - tr2->P0Vector();

    	double t_ca = (  (displacement ^ rel_v) - ( (this->VelVector().Scale(this->t0) - tr2->VelVector().Scale(tr2->t0)) ^ rel_v)  )/rel_v2;

    	Vector pos1 = this->P0Vector() + this->VelVector().Scale(t_ca - this->t0);
    	Vector pos2 =  tr2->P0Vector() +  tr2->VelVector().Scale(t_ca -  tr2->t0);
    
    	return (pos1 + pos2).Scale(0.5);
    }

    std::vector<double> track::closest_approach_midpoint_4d(track* tr2){

    	using namespace vector;
        // if (tr2->index == index) return P0Vector();
        if (tr2->index == index) {
			// std::cout<< "seed using same track"<<std::endl;
			return {tr2->x0,tr2->y0,tr2->z0,tr2->t0};
		}

        Vector rel_v = tr2->VelVector() - this->VelVector();
    	double rel_v2 = rel_v ^ rel_v ; //that's the dot product :/

    	Vector displacement = this->P0Vector() - tr2->P0Vector();

    	double t_ca = (  (displacement ^ rel_v) - ( (this->VelVector().Scale(this->t0) - tr2->VelVector().Scale(tr2->t0)) ^ rel_v)  )/rel_v2;

    	Vector pos1 = this->P0Vector() + this->VelVector().Scale(t_ca - this->t0);
    	Vector pos2 =  tr2->P0Vector() +  tr2->VelVector().Scale(t_ca -  tr2->t0);
    
    	Vector mid =  (pos1 + pos2).Scale(0.5);
		std::vector<double> mid_t = {mid.x,mid.y,mid.z, t_ca};
		return mid_t;
    }	

    Eigen::VectorXd track::ca_midpoint_kalman(track* tr2){

        using namespace vector;
//        if (tr2->index == index) return P0Vector();

        Vector rel_v = tr2->VelVector() - this->VelVector();
        double rel_v2 = rel_v ^ rel_v ; //that's the dot product :/

        Vector displacement = this->P0Vector() - tr2->P0Vector();

        double t_ca = (  (displacement ^ rel_v) - ( (this->VelVector().Scale(this->t0) - tr2->VelVector().Scale(tr2->t0)) ^ rel_v)  )/rel_v2;

        Vector pos1 = this->P0Vector() + this->VelVector().Scale(t_ca - this->t0);
        Vector pos2 =  tr2->P0Vector() +  tr2->VelVector().Scale(t_ca -  tr2->t0);

        Vector mid = (pos1 + pos2).Scale(0.5);

	Eigen::VectorXd mid_t(4);

	mid_t[0] = mid.x;
	mid_t[1] = mid.y;
	mid_t[2] = mid.z;
	mid_t[3] = t_ca;

	return mid_t;
    }


    double track::cos_angle_from_ip(){

        using namespace detector;

        Vector ip_direction = Vector(x0 - ip_x, y0 - ip_y, z0 - ip_z );
        double track_vel = beta()*constants::c;

        return (ip_direction ^ VelVector())/(track_vel*ip_direction.Magnitude());
    }

    // calculate shortest dist. from point to line
    double track::shortDistance()
    {
        //      std::cout << "222: " << line_point2 << " 111: " << line_point1 << "\n";
        Vector point = Vector(detector::ip_x, detector::ip_y, detector::ip_z );
        vector::Vector AB = vector::Vector(vx, vy, vz);
        vector::Vector AC = point - vector::Vector(x0, y0, z0);
        double area = vector::Vector(AB * AC).Magnitude();
        double CD = area / AB.Magnitude();
        return CD;
    }
    





}; //namespcae physics
