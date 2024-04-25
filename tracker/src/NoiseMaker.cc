
#include "NoiseMaker.hh"

double NoiseMaker::window = 1000*units::ns;

//------------------------------------------------------------------//

NoiseMaker::NoiseMaker(std::vector<physics::digi_hit*>& digis ){

get_real_hits(digis);
int index = 0;

hit_generator;
hit_generator.SetSeed( rand()*rand()*rand() % rand());
double average = window*rate_of_hits;
int total_hits = 0;

for(detID id:detID_list){ 		
	if(ts) std::cout<<"-----------------------for detID number: "<<index<<std::endl;
	int hit_q = hit_generator.Poisson(average);
	total_hits+= hit_q;
	if(hit_q !=0){
		std::vector<double> noise_times = event_timing(id,hit_q);
		
		make_digis(noise_times, id);
	}
	index++;
}


}

//------------------------------------------------------------------//

void NoiseMaker::preDigitizer(){
ParHandler hndlr;
hndlr.Handle();
double noise_hz = hndlr.par_map["noise_hz"];
run = false;
if(noise_hz>0){
	bool run = true;
	hits_per_second = noise_hz;
	rate_of_hits = hits_per_second/(1000000000*units::ns);	

	_geometry = new Geometry;

	if(ts) std::cout<<"Start of layer, wall and floor detID_list "<<std::endl;
	layer_detIDs(detID_list);
	wall_detIDs(detID_list);
	floor_detIDs(detID_list);
	back_detIDs(detID_list);

	detID_q = detID_list.size();
	double average_hits = window*rate_of_hits*detID_q;
	std::cout<<"-------------------------------------------------------------------------"<<std::endl;
	std::cout<<"NoiseMaker is turned on with "<<hits_per_second<<" hz noise rate over a window of "<<window<<" nanoseconds."<<std::endl;
	std::cout<<"The average number of noise digis per event is "<<average_hits<<std::endl;
	std::cout<<"-------------------------------------------------------------------------"<<std::endl;
}else if(noise_hz<0){
	std::cout<<"--------------------------------------ERROR-----------------------------------"<<std::endl;
	std::cout<<"The noisemaker will not run because the noise rate of "<<noise_hz<<" is invalid."<<std::endl;
}


}

std::vector<physics::digi_hit*> NoiseMaker::return_digis(){

return digi_hits;

}

////////////////////////////////////////////////////////////////////// Private Functions Below

void NoiseMaker::get_real_hits(std::vector<physics::digi_hit*> digis){
if(ts) if(digis.size() != 0){ std::cout<<"digis.size() grabbed "<<digis.size()<<std::endl;}
for(auto digi:digis){
        auto current_id = digi->det_id;
        double time = digi->t;
        bool id_already_exists = false;

        std::vector<real_hit_times>::iterator row;
        for(row = real_hits.begin(); row !=real_hits.end(); row++){
                if( current_id == row->id){ (row->hit_times).push_back(time); id_already_exists=true;}
        }
        if(!id_already_exists){
                std::vector<double> newtime{time};
                for(auto moment : newtime){
                        if(ts) std::cout<<"Time being inserted "<<moment<<std::endl;
                }
                real_hits.push_back(real_hit_times{current_id,newtime});
        }
}

if(ts)std::cout<<"---------------------------------Size of real_hits: "<<real_hits.size()<<std::endl;

/*
 std::vector<real_hit_times>::iterator row;
 int iteration = 0;
 for(row = real_hits.begin(); row !=real_hits.end(); row++){
 	detID id = row->id;
 	std::vector<double> times = row->hit_times;
 	std::cout<<"Iteration: "<<iteration<<" number of times: "<<times.size()<<" and detID"<<std::endl;
 	for(auto time: times){
 		std::cout<<"Time: "<<time<<std::endl;
        }
        id.Print();
        iteration++;
 
 
 }*/
}

//------------------------------------------------------------------//

std::vector<double> NoiseMaker::event_timing(detID id, int hit_q){
int hit_r = hit_q; //remaining hits
std::vector<double> hit_times;

while(hit_r>0){
        bool add_hit = false;
        double current_time = hit_generator.Uniform(window);
        if(hit_times.size() !=0){
                for(double previous_time:hit_times){
                        double difference = previous_time - current_time;
                        if(abs(difference)>= 25*units::ns){
                                add_hit = true;
                        }
                }
        }else{ add_hit = true;}

        std::vector<double> real_hit_times;
        bool get_times= get_detID_specific_hit_times( &real_hit_times, id);
        if(get_times){
                for(auto real_time : real_hit_times){
                        double difference = real_time - current_time;
                        if(abs(difference)< 25*units::ns){ add_hit = false;}
        }
        }

        if(add_hit){
        hit_times.push_back(current_time);
        hit_r--;
        //std::cout<<"hit_r: "<<hit_r<<" time: "<<current_time<<std::endl;
        }
}

return hit_times;
}

//------------------------------------------------------------------//

bool NoiseMaker::get_detID_specific_hit_times(std::vector<double>* times, detID id){
bool times_exist =false;
std::vector<real_hit_times>::iterator row;
int realhits = 0;
for(row = real_hits.begin(); row !=real_hits.end(); row++){
                if( id == row->id){ times = &(row->hit_times); times_exist = true;
                        if(ts) std::cout<<"Real hit(s) alread exist in this detID"<<std::endl;
                        for(auto time: *(times)){
                        if(ts) std::cout<<"Time: "<<time<<std::endl;
                        }                }
        }

return times_exist;

}

//------------------------------------------------------------------//

void NoiseMaker::make_digis(std::vector<double> times, detID id){
	for(double time : times){
	physics::digi_hit* digi = new physics::digi_hit();
	
	digi->det_id = id;
	digi->t = time;
	std::vector<double> location = set_hit_location(id);
	digi->x = location[0];
	digi->y = location[1];
	digi->z = location[2];
	int layer_index = id.yModule;;
	auto uncertainty = id.uncertainty();
	digi->ex = uncertainty[0];
	digi->ey = uncertainty[1];
	digi->ez = uncertainty[2];
	//dummy variables
	digi->e = 999;
	digi->px = 999;
	digi->py = 999;
	digi->pz = 999;
	digi->particle_mass = 999;
	digi->particle_energy = 999;
	digi->pdg = 999;
	

	digi_hits.push_back(digi);
	}
}

//------------------------------------------------------------------//

std::vector<double> NoiseMaker::set_hit_location(detID id){
    std::vector<double> location = {0,0,0};
    location = _geometry->GetCenter(id);

    int layer_index =id.yModule;//only modify hit location for regular layers, all others leave as center
    if(layer_index >= detector::n_floors){ //is regular layer detID
        std::vector<double> layer_widths = _geometry->GetDimensions(id);
        int long_direction_index = id.GetLongIndex(); //if 0 along x if 1 along z
	double length;
	if (long_direction_index == 0) length = layer_widths[0];
	else length = layer_widths[2];
            double randomlength = hit_generator.Uniform(-length/2, length/2);
            if(long_direction_index ==0){//length along x
                location[0]+=randomlength;
            }
            else{//length along z
                location[2]+=randomlength;
            }
    }
if(location == std::vector<double> {0,0,0}){std::cout<<"location is unmodified for a detID"<<std::endl;}
return location;
}

//------------------------------------------------------------------//

void NoiseMaker::layer_detIDs(std::vector<detID>& _detID_list){
    if(ts) std::cout<<"NoiseMaker::layer_detIDs"<<std::endl;
    int layerids = 0;
    for (int yModule = 0; yModule < detector::LAYERS_Y.size(); yModule++){
        if(yModule < 2 ){continue;} //ensure layers aren't floor layers
        for (int xModule = 0; detector::MODULE_X.size(); xModule++){
	    for (int zModule = 0; detector::MODULE_Z.size(); zModule++){
	        double length = detector::MODULE_X[0][1] - detector::MODULE_X[0][0];
		int num_x, num_z;
		if (yModule % 2 == 0) {
		   num_x = (int)(length / detector::scintillator_length);
		   num_z = (int)(length / detector::scintillator_width);
		} else {
		   num_x = (int)(length / detector::scintillator_width);
		   num_z = (int)(length / detector::scintillator_length);
		}
		for (int xIndex = 0; xIndex < num_x; xIndex++) {
		    for (int zIndex = 0; zIndex < num_z; zIndex++) {
		        detID _id = detID(xModule, xIndex, yModule, 0, zModule + detector::n_walls, zIndex);
			_detID_list.push_back(_id);
		    }
		}
	    }
        }
    }
}

//------------------------------------------------------------------//

void NoiseMaker::back_detIDs(std::vector<detID>& _detID_list){
    if(ts) std::cout<<"NoiseMaker::layer_detIDs"<<std::endl;
    for (int zModule = 0; zModule < detector::n_back_walls; zModule++){
        for (int xModule = 0; detector::MODULE_X.size(); xModule++){
	    double length = detector::MODULE_X[0][1] - detector::MODULE_X[0][0];
	    int num_x, num_y;
	    if (zModule % 2 == 0) {
		num_x = (int)(length / detector::scintillator_length);
		num_y = (int)(length / detector::scintillator_width);
	    } else {
		num_x = (int)(length / detector::scintillator_width);
		num_y = (int)(length / detector::scintillator_length);
	    }
	    for (int xIndex = 0; xIndex < num_x; xIndex++) {
		for (int yIndex = 0; yIndex < num_y; yIndex++) {
		    detID _id = detID(xModule, xIndex, 0, yIndex, zModule + detector::MODULE_Z.size() + detector::n_walls, 0);
		    _detID_list.push_back(_id);
		}
	    }
        }
    }
}

//------------------------------------------------------------------//

void NoiseMaker::wall_detIDs(std::vector<detID>& _detID_list){
    if(ts) std::cout<<"NoiseMaker::wall_detIDs"<<std::endl;
    for (int zModule = 0; zModule < detector::n_walls; zModule++){
        double length = detector::MODULE_X[0][1] - detector::MODULE_X[0][0];
	int num_x, num_y;
	if (zModule % 2 == 0) {
	    num_x = (int)(length / detector::scintillator_length);
	    num_y = (int)(length / detector::scintillator_width);
	} else {
	    num_x = (int)(length / detector::scintillator_width);
	    num_y = (int)(length / detector::scintillator_length);
	}
	for (int xIndex = 0; xIndex < num_x; xIndex++) {
	    for (int yIndex = 0; yIndex < num_y; yIndex++) {
	        detID _id = detID(0, xIndex, 0, yIndex, zModule, 0);
	        _detID_list.push_back(_id);
	    }
        }
    }
}

//------------------------------------------------------------------//

void NoiseMaker::floor_detIDs(std::vector<detID>& _detID_list){
    if(ts) std::cout<<"NoiseMaker::floor_detIDs"<<std::endl;
    for (int yModule = 0; yModule < detector::n_floors; yModule++){
        double length = detector::MODULE_X[0][1] - detector::MODULE_X[0][0];
	int num_x, num_z;
	if (yModule % 2 == 0) {
	    num_x = (int)(length / detector::scintillator_length);
	    num_z = (int)(length / detector::scintillator_width);
	} else {
	    num_x = (int)(length / detector::scintillator_width);
	    num_z = (int)(length / detector::scintillator_length);
	}
	for (int xIndex = 0; xIndex < num_x; xIndex++) {
	    for (int zIndex = 0; zIndex < num_z; zIndex++) {
	        detID _id = detID(0, xIndex, yModule, 0, detector::n_walls, zIndex);
	        _detID_list.push_back(_id);
	    }
        }
    }
}

//------------------------------------------------------------------//
