#include "globals.hh"
#include "VertexFinder.hh"
#include "statistics.hh"
#include "physics.hh"
#include "LinearAlgebra.hh"
#include "kalman.hh"
#include "kalman-test.hh"
#include <iostream>
#include <fstream>
#include <TMath.h>
#include "Math/ProbFunc.h"

#include <Eigen/Dense>


void VertexFinder::Seed_k_m()
{
	for (int n1 = 0; n1 < tracks_k_m.size(); n1++)
	{
		for (int n2 = n1 + 1; n2 < tracks_k_m.size(); n2++)
		{

			auto tr1 = tracks_k_m[n1];
			auto tr2 = tracks_k_m[n2];
			auto closest_dist = tr1->closest_approach(tr2);

			if (closest_dist < par_handler->par_map["vertex_seed_dist"])
			{

				// Find chi2 and midpoint of the seed.
				VertexFitter fitter_seed;
				auto current_seed = vertex_seed(tr1, tr2);
				current_seed.closest_dist = closest_dist;
				current_seed.seed_midpoint = current_seed.guess();
				std::vector<physics::track *> seed_tracks = {current_seed.tracks.first, current_seed.tracks.second};
				auto status0 = fitter_seed.fit(seed_tracks, current_seed.seed_midpoint);
				if (status0 == true) {
					current_seed.seed_midpoint = fitter_seed.parameters;
					current_seed.seed_midpoint_err = fitter_seed.parameter_errors;
					current_seed.chi2 = fitter_seed._merit;
					current_seed.closest_dist = current_seed.tracks.first->distance_to(Vector(current_seed.seed_midpoint[0],current_seed.seed_midpoint[1],current_seed.seed_midpoint[2]),current_seed.seed_midpoint[3]) + current_seed.tracks.second->distance_to(Vector(current_seed.seed_midpoint[0],current_seed.seed_midpoint[1],current_seed.seed_midpoint[2]),current_seed.seed_midpoint[3]);
				}	


				if(current_seed.closest_dist < par_handler->par_map["vertex_seed_dist"] && current_seed.chi2< par_handler->par_map["vertex_seed_chi2"]){
					seeds_k_m.push_back(current_seed);
					// if(par_handler->par_map["debug_vertex"]==1)
						// std::cout <<"  seed chi2: "<<fitter_seed._merit<<"\t. Seed dist: "<<closest_dist << "\t,new dist "<< current_seed.closest_dist<< ".\t Layers: " << current_seed.tracks.first->chi_s.size()+current_seed.tracks.second->chi_s.size()<<std::endl;
				}


			}
		} //n2
	}	  //n1


	for (int i1 = 0; i1 < seeds_k_m.size(); i1++){
		auto sd = seeds_k_m[i1];
		for (int n1 = 0; n1 < tracks_k_m.size(); n1++){
			auto tr1 = tracks_k_m[n1];
			auto dist = tr1->distance_to(Vector(sd.seed_midpoint[0],sd.seed_midpoint[1],sd.seed_midpoint[2]),sd.seed_midpoint[3]);
			if (dist< par_handler->par_map["vertex_add_max_distance"])
				sd.compatible_tracks+=1;
		}
		seeds_k_m[i1] = sd;
	}

	std::sort(seeds_k_m.begin(), seeds_k_m.end(), [](vertex_seed a, vertex_seed b) -> bool
			  { return a.score() < b.score(); });

	// if(par_handler->par_map["debug_vertex"]==1){
	// 	for (auto sd:seeds_k_m)
	// 		std::cout <<"  seed chi2: "<<sd.chi2<<".\t Seed dist: "<<sd.closest_dist << ",\t N tracks "<< sd.compatible_tracks<< ".\t Layers: " << sd.tracks.first->chi_s.size()+sd.tracks.second->chi_s.size()<<std::endl;
	// }			  

} //VF:Seed

void VertexFinder::FindVertices_k_m_hybrid()
{
	if (seeds_k_m.size() < 1){
		noSeeds+=1;
		return; // no seeds
	}


	while (seeds_k_m.size() > 0 and tracks_k_m.size() > 0)
	{
		auto current_seed = seeds_k_m[0];
		seeds_k_m.erase(seeds_k_m.begin());
		// Remove seeds with tracks that are already used in vertex
		// First, get a list of ID of remaining tracks
		std::vector <int> track_ids;
		for (auto tk : tracks_k_m){
			track_ids.push_back(tk->index);
		}
		// Then, see if the seed can be found in the remaining tracks
		std::vector<int>::iterator find1, find2;
		find1 = std::find(track_ids.begin(), track_ids.end(),current_seed.tracks.first->index);
		find2 = std::find(track_ids.begin(), track_ids.end(),current_seed.tracks.second->index);
		// Skip this seed, if either hit in the seeds is not found in the remaining hits:
		if ((find1 == track_ids.end()) || (find2 == track_ids.end())){
			continue;
		}

		std::vector<physics::track *> used_tracks = {};
		std::vector<physics::track *> unused_tracks = {};

		// Find the best midpoint of the seed
		if(par_handler->par_map["debug_vertex"]==1){
			std::cout <<"  New seed--------------------"<<std::endl;
			std::cout <<"   seed chi2: "<<current_seed.chi2<<std::endl;
		}
		auto seed_midpoint = current_seed.seed_midpoint;
		auto parameter_errors = current_seed.seed_midpoint_err;
		auto chi2_current = current_seed.chi2;
		VertexFitter fitter_find;		

		// Add tracks to the track list for fit
		// First, add tracks from seed
		used_tracks.push_back(current_seed.tracks.first);
		used_tracks.push_back(current_seed.tracks.second);
		// Second, calculate the distance from each track to the seed
		// , and sort by the distance (small->large)
		// for (auto tr : tracks_k_m)
		// 	tr->distance = tr->distance_to(Vector(seed_midpoint[0],seed_midpoint[1],seed_midpoint[2]),seed_midpoint[3]);
		// std::sort(tracks_k_m.begin(), tracks_k_m.end(), [](physics::track* a, physics::track* b) -> bool
		// 	  { return a->distance < b->distance; });
		
		// Then deal with the rest 
		// for (auto tr : tracks_k_m)
		auto ntracks = tracks_k_m.size();
		for (auto itrack=0; itrack<ntracks; itrack++)
		{	
			// sort tracks by distance to seed_midpoint. This is inside the for loop to keep it updated with changing seed_midpoint
			for (auto tr_temp : tracks_k_m){
				// tr_temp->distance = tr_temp->distance_to(Vector(seed_midpoint[0],seed_midpoint[1],seed_midpoint[2]),seed_midpoint[3]);
				tr_temp->distance_weighted = tr_temp->chi2_distance_to_pointerror(Vector(seed_midpoint[0],seed_midpoint[1],seed_midpoint[2]),seed_midpoint[3], parameter_errors);
			}
			std::sort(tracks_k_m.begin(), tracks_k_m.end(), [](physics::track* a, physics::track* b) -> bool
				{ return a->distance_weighted < b->distance_weighted; });
			auto tr = tracks_k_m.front();
			tr->distance = tr->distance_to(Vector(seed_midpoint[0],seed_midpoint[1],seed_midpoint[2]),seed_midpoint[3]);

			if(par_handler->par_map["debug_vertex"]==1)
				// std::cout << "                   -parameters now: " << seed_midpoint[0] << ", " << seed_midpoint[1] << ", " << seed_midpoint[2] << ", " << seed_midpoint[3] << std::endl;
				std::cout << "    track #"<<tr->index<<",\t" << tr->distance_weighted << " [sigma] " <<  tr->distance << " [cm] to estimated vertex.";

			// Special treatment for tracks in seed: add them directly even if they do not pass the distance cut.
			if (tr->index==current_seed.tracks.first->index || tr->index==current_seed.tracks.second->index){
				if(par_handler->par_map["debug_vertex"]==1)
					std::cout <<"\t- added (seed)"<<std::endl;
			}
			// Add the track to pool if distance is below threshold
			// if (current_seed.closest_approach(tr) < par_handler->par_map["vertex_add_max_distance"]) // Tom: I think this is wrong. We should calculate the distance to the best guess of the seed, not the individual tracks in the seed. 
			else if (tr->distance < par_handler->par_map["vertex_add_max_distance"])
			{

				// auto chi2_delta = tr->chi2_distance_to_pointerror(Vector(seed_midpoint[0],seed_midpoint[1],seed_midpoint[2]),seed_midpoint[3], current_seed.seed_midpoint_err);
				// if (chi2_delta<par_handler->par_map["vertex_chi2_add"]){
				// 	if(par_handler->par_map["debug_vertex"]==1)
				// 		std::cout <<"\t- added, delta-Chi2: " <<chi2_delta<<std::endl;
				// 	used_tracks.push_back(tr);
				// }
				// else {
				// 	unused_tracks.push_back(tr);
				// 	if(par_handler->par_map["debug_vertex"]==1)
				// 		std::cout <<"\t  - dropped, delta-Chi2: " <<chi2_delta<<std::endl;
				// }


				// Do a fit 
				used_tracks.push_back(tr);
				auto status = fitter_find.fit(used_tracks, seed_midpoint, 1.0);
				auto chi2_delta = fitter_find._merit - chi2_current;
				// Erase the hit from used list if delta chi2 is too large.
				if(status == false || chi2_delta>par_handler->par_map["vertex_chi2_add"]){
					unused_tracks.push_back(tr);
					used_tracks.erase(used_tracks.end()-1);
					if (par_handler->par_map["debug_vertex"]==1){
						if (status==true)
							std::cout <<"\t  - dropped, failed vertex_chi2_add cut. delta-Chi2: " <<chi2_delta<<std::endl;
						else
							std::cout <<"\t  - dropped, fit failed." <<std::endl;
					}
				}
				// Else, update the new midpoint and chi2
				else{
					seed_midpoint = fitter_find.parameters;
					parameter_errors = fitter_find.parameter_errors;
					chi2_current = fitter_find._merit;
					if (par_handler->par_map["debug_vertex"]==1)
						std::cout <<"\t- added, delta-Chi2: " <<chi2_delta<<std::endl;
				}
			}			
			else
			{
				if(par_handler->par_map["debug_vertex"]==1) std::cout<< std::endl;
				unused_tracks.push_back(tr);
			}

			tracks_k_m.erase(tracks_k_m.begin());
		}

		// Do the fit on all tracks added to the seed
		VertexFitter::adaptive_iterations = par_handler->par_map["vertex_adaptive_niters"];
		VertexFitter::adaptive_annealing_factor = par_handler->par_map["vertex_adaptive_r"];
		VertexFitter::adaptive_temperature = par_handler->par_map["vertex_adaptive_T0"];
		VertexFitter::adaptive_chi2_cutoff = par_handler->par_map["vertex_adaptive_chi2"];		
		VertexFitter fitter;

		auto status = fitter.fit(used_tracks, seed_midpoint, 1.0);
		auto delta_chi2_list = fitter.delta_chi2();
		if(par_handler->par_map["debug_vertex"]==1)
			std::cout <<"    -parameters before dropping: " << seed_midpoint[0] << ", " << seed_midpoint[1] << ", " << seed_midpoint[2] << ", " << seed_midpoint[3] << std::endl;

		// Drop bad tracks iteratively (if using the normal chi2 fit)
		if (par_handler->par_map["vertex_use_adaptive_fit"]==0)
		{

			// for(int i=delta_chi2_list.size()-1; i>=0; i--){
			// 	if (delta_chi2_list[i]> par_handler->par_map["vertex_chi2_drop"]){
			// 		unused_tracks.push_back(used_tracks.at(i));
			// 		used_tracks.erase(used_tracks.begin()+i);
			// 		if(par_handler->par_map["debug_vertex"]==1)
			// 			std::cout<<"     -track #"<<unused_tracks.back()->index <<" dropped from fit with delta-Chi2 "<<delta_chi2_list[i] <<std::endl;
			// 	}
			// }

			// Drop one bad tracks each time
			auto max_chi2 = max_element(delta_chi2_list.begin(), delta_chi2_list.end());
			while (*max_chi2> par_handler->par_map["vertex_chi2_drop"]){
				auto i_track_max_chi2 = max_chi2 - delta_chi2_list.begin();
				unused_tracks.push_back(used_tracks.at(i_track_max_chi2));
				used_tracks.erase(used_tracks.begin()+i_track_max_chi2);	
				if(par_handler->par_map["debug_vertex"]==1)
					std::cout<<"     -track #"<<unused_tracks.back()->index <<" dropped from fit with delta-Chi2 "<<*max_chi2<<std::endl;

				if (used_tracks.size()<2){
					break; // Stop dropping, if there are not enough tracks
				}					

				// Fit again
				auto status = fitter.fit(used_tracks, seed_midpoint, 1.0);	
				delta_chi2_list = fitter.delta_chi2();
				max_chi2 = max_element(delta_chi2_list.begin(), delta_chi2_list.end());	
				seed_midpoint = fitter.parameters;
				if(par_handler->par_map["debug_vertex"]==1)
					std::cout <<"    -parameters updated: " << seed_midpoint[0] << ", " << seed_midpoint[1] << ", " << seed_midpoint[2] << ", " << seed_midpoint[3] << std::endl;

			}

			if (used_tracks.size()<2){
				for (auto track : used_tracks)
					unused_tracks.push_back(track);
				tracks_k_m = unused_tracks;
				if(par_handler->par_map["debug_vertex"]==1) std::cout <<"   FAILED (no enough tracks) "<<std::endl;
				continue; // Skip, if there are not enough tracks
			}			
			// Fit again with higher accuracy
			status = fitter.fit(used_tracks, seed_midpoint, 0.1);

		}
		else{
			status = fitter.fit_adaptive(used_tracks, seed_midpoint);
		}
		delta_chi2_list = fitter.delta_chi2();


		if (used_tracks.size()<2){
			for (auto track : used_tracks)
				unused_tracks.push_back(track);
			tracks_k_m = unused_tracks;
			if(par_handler->par_map["debug_vertex"]==1) std::cout <<"   FAILED (no enough tracks) "<<std::endl;
			continue; // Skip, if there are not enough tracks
		}			


		auto chi2_sum = fitter._merit;
		auto chi2_dof = fitter.ndof;
		// if (status == false or (ROOT::Math::chisquared_cdf(chi2_sum, chi2_dof) >= par_handler->par_map["vertex_chi2_final_pval"]))
		if (status == false or chi2_sum/chi2_dof >= par_handler->par_map["vertex_chi2_final"])
		{
			if (status == false){
				if(par_handler->par_map["debug_vertex"]==1)
					std::cout <<"   FAILED (fit did not converge) "<<std::endl;
				noConverge += 1;
			}
			else{
				if(par_handler->par_map["debug_vertex"]==1)
					std::cout <<"   FAILED (rejected by chi2 cut). Track chi2/dof: "<<chi2_sum << "/" << chi2_dof << " = " << chi2_sum/chi2_dof <<std::endl;
				missedChi2 += 1;
			}
			for (auto track : used_tracks)
				unused_tracks.push_back(track);
			tracks_k_m = unused_tracks;			
			continue;
		}

		if(par_handler->par_map["debug_vertex"]==1)
			std::cout <<"   ** SUCCEED. Vertex #"<<vertices_k_m.size() <<" found with "<< used_tracks.size()<<" tracks"<<std::endl;
		

		double cos_opening_angle = -1.0;
		if (used_tracks.size() == 2)
		{

			auto tr1 = used_tracks[0];
			auto tr2 = used_tracks[1];

			cos_opening_angle = tr1->vx * tr2->vx + tr1->vy * tr2->vy + tr1->vz * tr2->vz;
			cos_opening_angle = cos_opening_angle / (tr1->beta() * tr2->beta() * constants::c * constants::c);
		}

		auto good_vertex = new physics::vertex(fitter.parameters, cos_opening_angle);

		for (auto track : used_tracks)
		{
			good_vertex->track_indices.push_back(track->index);
		}

		good_vertex->CovMatrix(fitter.cov_matrix, fitter.npar);
		good_vertex->merit(fitter.merit());
		good_vertex->delta_chi2_list = delta_chi2_list;
		vertices_k_m.push_back(good_vertex);
		tracks_k_m = unused_tracks;
	}

}

void VertexFinder::FindVertices_k()
{
	if (seeds_k_m.size() < 1)
		return;

	while (seeds_k_m.size() > 0 and tracks_k_m.size() > 0)
	{
		std::vector<physics::track *> used_tracks = {};
		std::vector<physics::track *> unused_tracks = {};

		auto current_seed = seeds_k_m[0];

		seeds_k_m.erase(seeds_k_m.begin());

		used_tracks = tracks_k_m;

		double drops = -1;

		int i = 0;

		while (drops != 0)
		{

		kalman_vertex kfv;
		kfv.dropping = true;
		kfv.vertexer(used_tracks, &current_seed);

		if (kfv.status != 2) {
			break;
		}

		used_tracks = {}; // clear used tracks, only push back good ones

		drops = 0;

		// drop tracks that don't make the beta cut
		for (int k=0; k < kfv.added_tracks.size(); k++) {
//			double v = 0;
//			for (int i=0; i < 3; i++) v += std::pow(kfv.added_tracks[k]->q_s[i],2);
//			v = std::sqrt(v);

//			if (std::abs(kfv.pulls_v[k]) < kalman::pull_cut_drop) {

			double beta = (kfv.added_tracks[k]->q_s).norm() / constants::c;

//			if (!(kalman::v_cut_drop[0] < beta && beta < kalman::v_cut_drop[1])) {
			if (!(par_handler->par_map["v_cut_drop[0]"] < beta && beta < par_handler->par_map["v_cut_drop[1]"])) {
//			if (!(kalman::v_cut_drop[0] < kfv.pulls_v[k] && kfv.pulls_v[k] < kalman::v_cut_drop[1])) {
				unused_tracks.push_back(kfv.added_tracks[k]);
				drops++;
			}
			else {
				used_tracks.push_back(kfv.added_tracks[k]);
			}
		}

		if (used_tracks.size() < 2) {
			break;
			//continue;
		}

		//if (failed) break;

		if (i == 25) {std::cout << "i break" << std::endl; break;}
		i++;

		if (drops != 0) continue; // only evaluate vertex if no tracks were dropped, otherwise keep going

		//kalman_vertex kfv_2;
		//kfv_2.dropping = false;
		//kfv_2.vertexer(used_tracks, &current_seed);

		int n_track_params = 4;
        	int ndof = ((4.0 + 3.0) * kfv.added_tracks.size() - n_track_params ); // 4 x_k and 3 q_k parameters
	        if (ndof < 1) ndof = 1;

//		if (kfv.chi_v / ndof > cuts::kalman_vertex_chi) {
		if (kfv.chi_v / ndof > par_handler->par_map["kalman_vertex_chi"]) {
			break;
//			continue;
		}

		auto good_vertex = new physics::vertex(kfv.x_s, -2.0);

                for (auto track : used_tracks)
                {
                        good_vertex->track_indices.push_back(track->index);

			good_vertex->q_f.push_back(track->q_f);
			good_vertex->D_f.push_back(track->D_f);

			good_vertex->q_s.push_back(track->q_s);
			good_vertex->D_s.push_back(track->D_s);
                }

		vertices_k_m.push_back(good_vertex);

                tracks_k_m = unused_tracks;

		} // drop while loop

	}
}

std::vector<physics::track *> VertexFitter::track_list = {};
std::vector<double> VertexFitter::parameters = {};
std::vector<double> VertexFitter::parameter_errors = {};
double VertexFitter::_merit = 0.;
double VertexFitter::cov_matrix[VertexFitter::npar][VertexFitter::npar];

int VertexFitter::adaptive_iterations = 10;
double VertexFitter::adaptive_annealing_factor = 0.8;
double VertexFitter::adaptive_temperature = 3;
double VertexFitter::adaptive_chi2_cutoff = 60;


bool VertexFitter::bad_fit = false;
void VertexFitter::nll(int &npar, double *gin, double &f, double *pars, int iflag)
{
	// std::ofstream file;
	// file.open("print.txt", std::ios_base::app);

	using Vector = vector::Vector;
	double _x = pars[0];
	double _y = pars[1];
	double _z = pars[2];
	double _t = pars[3];

	double error = 0.0;

	int n = 0;

	for (auto track : VertexFitter::track_list)
	{

		double dist = track->distance_to(Vector(_x, _y, _z), _t);


		// -- Use Negative Log Likelihood
		// double err = track->err_distance_to(Vector(_x, _y, _z), _t);
		// error += 0.5 * (dist / err) * (dist / err) + TMath::Log(err);
		// error += track->err_distance_to_mod(Vector(_x, _y, _z), _t);

		// -- Use chi2 instead
		error += track->chi2_distance_to(Vector(_x, _y, _z), _t);


		

		if (isnan(error))
		{
			bad_fit = true;
			//std::cout << " Bad Vertex fit! " << std::endl;
			//std::cout << dist << " " << err << std::endl;
//			track->CovMatrix().Print();
			return;
		}
	}

	f = error;
}

void VertexFitter::cost_adaptive(int &npar, double *gin, double &f, double *pars, int iflag)
{
// Adaptive cost funtion that deweight outliers.

	// std::ofstream file;
	// file.open("print.txt", std::ios_base::app);

	using Vector = vector::Vector;
	double _x = pars[0];
	double _y = pars[1];
	double _z = pars[2];
	double _t = pars[3];

	double error = 0.0;
	// ndof = 0;

	for (auto track : VertexFitter::track_list)
	{
		// -- Use weighted chi2
		auto delta_chi2 = track->chi2_distance_to(Vector(_x, _y, _z), _t);
		auto weight_i = calc_adaptive_weight(delta_chi2);

		// add to total chi2
		error += delta_chi2*weight_i;

		// add to DOF
		// ndof += 3*weight_i;

		if (isnan(error))
		{
			bad_fit = true;
			//std::cout << " Bad Vertex fit! " << std::endl;
			//std::cout << dist << " " << err << std::endl;
			// track->CovMatrix().Print();
			return;
		}
	}
	// ndof = ndof - 4;
	f = error;
}
