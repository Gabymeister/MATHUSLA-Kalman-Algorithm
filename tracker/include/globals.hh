#include<cstdlib>
#include<vector>
#include <fstream>
#include <string>
#include <iostream>


#ifndef UNITS_HH
#define UNITS_HH

// A python-like print function
template<typename T, typename ...TAIL>
void print(const T &t, TAIL... tail);



namespace units{
	const double mm = 10.0;
	const double cm = 1.00;
	const double ns = 1.00;
	const double MeV = 1.00;
};

namespace detector{ 
	using namespace units;
	const double ip_x = 0.0;
	const double ip_y = 0.0;
	const double ip_z = 0.0;


const std::vector<std::vector<double>> MODULE_X = {{-1950.0*cm,-1050.0*cm},
{-950.0*cm,-50.0*cm},
{50.0*cm,950.0*cm},
{1050.0*cm,1950.0*cm},
};

const std::vector<std::vector<double>> MODULE_Z = {{7000.0*cm,7900*cm},
{8000.0*cm,8900.0*cm},
{9000.0*cm,9900.0*cm},
{10000.0*cm,10900.0*cm},
};

	const std::vector<double> COSMIC_SHIFT = {0.0, 547*cm, 11950.0*cm}; // shift of sim cosmic -> main coordinates
	const double surface_height = 8547*units::cm;

const int n_modules = 16;
const int n_floors = 2;
	const double scintillator_length = 450.0*units::cm;
	const double scintillator_width = 4.50*units::cm;
const double scintillator_height = 1.6*units::cm;
const double scintillator_thickness = 1.0*units::cm; //Just the sensitive part
	const double time_resolution = 1.0*units::ns;

	// 2023-04-25: changing to 6 TOP layers. Numbers updated from simulation by print. 
const std::vector<std::vector<double>> LAYERS_Y={{8550.0*cm,8551.6*cm},
{8631.6*cm,8633.2*cm},
{9893.2*cm,9894.8*cm},
{9974.8*cm,9976.4*cm},
{10056.4*cm,10058.0*cm},
{10138.0*cm,10139.6*cm},
};

	const int n_layers = LAYERS_Y.size();

	const double x_min = MODULE_X[0][0];
	const double y_min = LAYERS_Y[0][0];
	const double z_min = MODULE_Z[0][0];

	const double x_max = MODULE_X[MODULE_X.size()-1][1];
	const double y_max = LAYERS_Y[LAYERS_Y.size()-1][1];
	const double z_max = MODULE_Z[MODULE_Z.size()-1][1];

	//FLOOR TILE WIDTHS

	const double floor_x_width = scintillator_length; // 50.0*units::cm;
	const double floor_z_width = scintillator_width; // 50.0*units::cm;

    //WALL TILE WIDTHS

    const double wall_x_width = scintillator_width; // 50.0*units::cm;
    const double wall_y_width = scintillator_length; // 50.0*units::cm;

    //FRONT WALL PARAMETERS

const double wall_gap = 1.0*units::cm; //gap on each side of wall
const double wall_gap2 = 100.0*units::cm; //gap between two walls
const double wall_height = 1343.2*units::cm;
const double wall_start_y = y_min - 3*units::cm; //min y value of wall (casing included)
const std::vector<std::vector<double>> FRONT_WALL_Z{{9895.8*cm,9897.4*cm},
{9997.4*cm, 6999.0*cm}
};
const int n_walls = FRONT_WALL_Z.size();
    //FOR statistics.hh ONLY - NEW MIN Z WITH WALL
    const double z_min_wall = z_min - wall_gap - scintillator_height;

    //BACK WALL PARAMETERS
    const double back_wall_top = LAYERS_Y[0][0];
    const double back_wall_bottom = LAYERS_Y[0][0] - (MODULE_X[0][1] - MODULE_X[0][0]); //1 module thickness down from bottom of tracking layers
    const std::vector<std::vector<double>> BACK_WALL_Z{{13900.0*cm,13901.6*cm},
    {13981.6*cm,13983.2*cm},
    {14063.2*cm,14064.8*cm},
    {14144.8*cm,14146.4*cm},
    };
    const double n_back_walls = BACK_WALL_Z.size();
};

namespace constants{

	const double c = 29.97*units::cm/units::ns;
	const double optic_fiber_n = 1.580; //2023.5.3 Tom: change from 1.5->1.58 according to sim. estimate for the optical fiber index of refraction

};

namespace cuts{

	//digi cuts and constants
	const double digi_spacing = 20.0*units::ns;
	const double SiPM_energy_threshold = 0.65*units::MeV;

	//seeding
	const double seed_ds2 = 5.0; //sigma
	const double seed_residual = 10.0; //sigma

	//tracking
	const double residual_drop = 12.0; //sigma
	const double residual_add = 12.0; //sigma
	const double track_chi2 = 5.0;
	const int track_nlayers = 4; // used to be 3
	const int nseed_hits = 4;
	const double time_difference_drop = 12.0; //sigma
	const double seed_time_difference = 10.0; //ns
	const int ntrack_hits = 4;
	const double distance_to_hit = 75.0*units::cm;

	//kalman tracker
        const double kalman_chi_s = 150.0;
        const double kalman_chi_add = 200.0;
        const double kalman_track_chi = 15.0;
        const std::vector<double> kalman_v_add = {0.800000,1.200000};
        const std::vector<double> kalman_v_drop = {0.900000,1.100000};

	//merging
        const double merge_cos_theta = 0.998;
        const double merge_distance = 25.0*units::cm;

	//cleaning step
	const double chi2_add = 10.0;
	const double chi2_drop = 16.0;
	const int cleaning_nhits = 6;

	//vertexing
        const double vertex_seed_dist = 100.0*units::cm;
        const double vertex_chi2 = 15.0;
        const double vertex_add_max_distance = 100.0*units::cm;
        const double kalman_vertex_chi_add = 100000.0;
        const double kalman_vertex_chi = 100.0;

	//run options
        const int start_ev = 0;
        const int end_ev = 50;

    //digi hit cuts for floors and wall
    const std::vector<bool> include_floor = { true, true }; //ith index for ith floor from bottom

    const bool include_wall = true;
    const double wall_y_cut = detector::wall_start_y + detector::wall_height; //all digi hits above this will be thrown out
};


namespace kalman{

	const double sigma_ms_p = 5.01; // [rad MeV]
        const double p = 500.0; // [MeV] representative momentum

//        const std::vector<double> v_cut_add = {-2.0,2.0};
//        const std::vector<double> v_cut_drop = {-2.0,2.0};
        const std::vector<double> v_cut_add = {0.0,2.0};
        const std::vector<double> v_cut_drop = {0.5,1.5};

	const double pull_cut_add = 100.0;
	const double pull_cut_drop = 100.0;
}



#endif
