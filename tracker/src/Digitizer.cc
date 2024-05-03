#include "Digitizer.hh"
//#include "NoiseMaker.hh"
#include "physics.hh"
#include "globals.hh"
#include <TRandom.h>
#include <TRandom3.h>
#include <time.h>
#include <stdlib.h>
#include <algorithm>
#include <functional>
#include <array>
#include <iostream>
#include <vector>
#include <sys/time.h>
//#include "units.hh"
//DIGITIZATION ALGORITHM
namespace physics{
	bool time_sort(physics::sim_hit *hit1, physics::sim_hit *hit2) {return (hit1->t < hit2->t);}
};

std::vector<physics::digi_hit*> Digitizer::Digitize(){

	//this is the vector of digi_hits we will return at the end of the function
	std::vector<physics::digi_hit*> digis;


	//looping through each detector ID
	std::vector<physics::sim_hit*> current_hits;
	std::vector<physics::sim_hit*> current_remaining_hits = hits;
	std::vector<physics::sim_hit*> next_remaining_hits;

	while (current_remaining_hits.size() > 0){

		//current detector id which we are working in
		auto current_id = (current_remaining_hits[0])->det_id;


		//taking out all hits with the same detector id to be digitized, leaving the remaing for the next iteration
		for (auto hit : current_remaining_hits){
			if (hit->det_id.IsNull()) continue;
			if (hit->det_id == current_id){ current_hits.push_back(hit);}
			else {next_remaining_hits.push_back(hit);}

		}

		//time sorting current hits

		std::sort(current_hits.begin(), current_hits.end(), &physics::time_sort);


		// going through all hits until they are either all added to digis, or dropped

		while (current_hits.size() > 0){

			std::vector<physics::sim_hit*> used_hits;
			std::vector<physics::sim_hit*> unused_hits;

			double t0 = (current_hits[0])->t;

			double e_sum = 0;

			for (auto hit : current_hits){
				if ( hit->t < t0 + cuts::digi_spacing ){
					e_sum += hit->e;
					used_hits.push_back(hit);
				} else { unused_hits.push_back(hit);}
			}
            
            // ignoring all hits in ignored floors/walls or above wall y cut
            
            auto current_center = current_id.GetCenter();
            
			if (e_sum > cuts::SiPM_energy_threshold){
				physics::digi_hit* current_digi = new physics::digi_hit();
				current_digi->det_id = current_id;
				for (auto hit : used_hits){current_digi->AddHit(hit);}
				current_digi->index = ( digis.size() );
				digis.push_back(current_digi);
				current_hits = unused_hits;
			}  else { current_hits.erase(current_hits.begin());}

		} // while (current_hits.size() > 0)

		//resetting all the sorting vectors, and assigning the next remianing hits to the next iteration for current remaining
		current_remaining_hits.clear();
		current_remaining_hits = next_remaining_hits;
		next_remaining_hits.clear();
		current_hits.clear();

	} // while (current_remaining_hits.size() > 0)

	//at this point, all of the digi_hits in the digi_vector have the hits which will make them up. However, they don't have any of their energy, position, or timing information added.
	//Below, we compute the energy, time, and position of all of the digi hits
	//We incoorporate the time and position smearing into this calculation as well

	int counter = 0;

	struct timeval curTime;
	gettimeofday(&curTime, NULL);
	long int micro_sec = curTime.tv_usec;
	srand( micro_sec );

	// Now we throw out hits in the floor and wall to simulate reduced detector efficiency
	std::vector<physics::digi_hit*> digis_not_dropped;

	// now manage hits in the floor and wall
	for (auto digi : digis){
		auto current_id = digi->det_id;

		auto center = current_id.GetCenter();
		auto long_direction_index = current_id.GetLongIndex();
		auto uncertainty = current_id.uncertainty();

	   //first IDs are the veto layers	

		double e_sum = 0;
		double long_direction_sum = 0.0;
		double t_sum = 0;
		double x_sum = 0;
		double y_sum = 0;

		for (auto hit : digi->hits){
			e_sum += hit->e;
			t_sum += hit->t * hit->e;
			y_sum += hit->y * hit->e;
			x_sum += hit->x * hit->e;
			if (long_direction_index == 0){
				long_direction_sum += hit->x * hit->e;
			} else if (long_direction_index == 1) {
				long_direction_sum += hit->y * hit->e;
			} else if (long_direction_index == 2) {
				long_direction_sum += hit->z * hit->e;
			}
		}

		digi->e = e_sum;
		digi->t = t_sum/e_sum;
		//digi->y = y_sum/e_sum;
		digi->ex = uncertainty[0];
		digi->ey = uncertainty[1];
		digi->ez = uncertainty[2];

		// --------------------Digitize---------------------------------------------
		// Wall: to the center of each square
		// Floor and tracking layers: to the center of each bar
		//note: et is the same for all of them and is set in the digi class defintion
		// if (current_id.isFloorElement || current_id.isWallElement){
		if (current_id.allignment == 1){//vertical
			digi->z = center[1];
			if (long_direction_index == 0){
				digi->x = long_direction_sum/e_sum;
				digi->y = center[0];
			} else {// y is long direction
				digi->x = center[0];
				digi->y = long_direction_sum/e_sum;
			}			

	    } else {//horizontal
		    if (long_direction_index == 0){
				digi->x = long_direction_sum/e_sum;
				digi->y = center[0];
				digi->z = center[1];
			} else {// z is long direction
				digi->x = center[0];
				digi->y = center[1];
				digi->z= long_direction_sum/e_sum;
			}
		}
		digi->long_direction_index = long_direction_index;

		//TIME AND POSITION SMEARING!!!!!!!!!!!!!!!
		//we see the random number generator with a number that should be completly random:
		//the clock time times the layer index times the number of digis

		// Time smearing
		digi->t += generator.Gaus(0.0, digi->et);



		//--Position smearing for Wall hits 		
		//--2023-06-30 Tom: turn off
		if (current_id.allignment == 1) {
			if (long_direction_index == 0) {
				double smeared_x = digi->x + generator.Gaus(0.0, digi->ex);
				digi->x = smeared_x;
			} else if (long_direction_index == 1){
				double smeared_y = digi->y + generator.Gaus(0.0, digi->ey);
				digi->y = smeared_y;
			}
			digis_not_dropped.push_back(digi);

		// Position smearing for Tracker && Floor hits
		} else if (long_direction_index == 0) {
			double smeared_x = digi->x + generator.Gaus(0.0, digi->ex);
			digi->x = smeared_x;
		} else if (long_direction_index == 2) {
			double smeared_z = digi->z + generator.Gaus(0.0, digi->ez);
			digi->z = smeared_z;
		}


		digis_not_dropped.push_back(digi); // it's a tracking / trigger layer hit, so it will never be dropped
	}

	digis = digis_not_dropped; // only keep hits not dropped by inefficiency in the floor or wall
	//setting digi indices
	int k = 0;
	for (auto digi : digis) {
		digi->index = k++;		
	}

	return digis;



}
