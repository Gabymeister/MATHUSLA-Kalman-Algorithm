
#include "NoiseMaker.hh"

double NoiseMaker::window = 1000*units::ns;

//------------------------------------------------------------------//

NoiseMaker::NoiseMaker(std::vector<physics::digi_hit*>& digis){
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
			make_digis(noise_times, id, longCenters[index]);
		}
		index++;
	}
}

//------------------------------------------------------------------//

void NoiseMaker::preDigitizer(GeometryHandler* _geometry){
	Geometry = _geometry;
	ParHandler hndlr;
	hndlr.Handle();
	double noise_hz = hndlr.par_map["noise_hz"];
	if(noise_hz>0){
		run = true;
		hits_per_second = noise_hz;
		rate_of_hits = hits_per_second/(1000000000*units::ns);	

		if(ts) std::cout<<"Start of layer, wall and floor detID_list "<<std::endl;
		layer_detIDs(detID_list, longCenters);
		wall_detIDs(detID_list, longCenters);
		floor_detIDs(detID_list, longCenters);
		back_detIDs(detID_list, longCenters);

		detID_q = detID_list.size();
		double average_hits = window*rate_of_hits*detID_q;
		std::cout<<"-------------------------------------------------------------------------"<<std::endl;
		std::cout<<"NoiseMaker is turned on with "<<hits_per_second<<" hz noise rate over a window of "<<window<<" nanoseconds."<<std::endl;
		std::cout<<"The average number of noise digis per event is "<<average_hits<<std::endl;
		std::cout<<"-------------------------------------------------------------------------"<<std::endl;
	} else if(noise_hz<0) {
		std::cout<<"--------------------------------------ERROR-----------------------------------"<<std::endl;
		std::cout<<"The noisemaker will not run because the noise rate of "<<noise_hz<<" is invalid."<<std::endl;
	}
}

//---------------------------------------------------------------//

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
                        if(ts) std::cout<<"Real hit(s) already exist in this detID"<<std::endl;
                        for(auto time: *(times)){
                        if(ts) std::cout<<"Time: "<<time<<std::endl;
                        }                }
        }

return times_exist;

}

//------------------------------------------------------------------//

void NoiseMaker::make_digis(std::vector<double> times, detID id, double _longCenter){
	for(double time : times){
	physics::digi_hit* digi = new physics::digi_hit();
	
	digi->det_id = id;
	digi->t = time;
	std::vector<double> location = set_hit_location(id, _longCenter);
	digi->x = location[0];
	digi->y = location[1];
	digi->z = location[2];
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

std::vector<double> NoiseMaker::set_hit_location(detID id, double _longCenter){
    std::vector<double> location;
	double length = Geometry->GetLength();
	double randomlength = hit_generator.Uniform(-length/2, length/2);
	std::vector<double> centers = id.GetCenter();
	if (id.GetLongIndex() == 0) {
		location = {_longCenter, centers[0], centers[1]};
		location[0]+=randomlength;
	} else if (id.GetLongIndex() == 1) {
		location = {centers[0], _longCenter, centers[1]};
		location[1]+=randomlength;
	} else if (id.GetLongIndex() == 2) {
		location = {centers[0], centers[1], _longCenter};
		location[2]+=randomlength;
	}

	if(location == std::vector<double> {0,0,0}){std::cout<<"location is unmodified for a detID"<<std::endl;}
	return location;
}

//------------------------------------------------------------------//

void NoiseMaker::layer_detIDs(std::vector<detID>& _detID_list, std::vector<double>& _longCenters){
    if(ts) std::cout<<"NoiseMaker::layer_detIDs"<<std::endl;
	double length = Geometry->GetXModule()[1] - Geometry->GetXModule()[0];//Length of module
	double distance = Geometry->GetXModule()[2] - Geometry->GetXModule()[1];//distance between modules
	int num_xMods = (int)(Geometry->GetXModule().size()/2); //Number of xModules
	int num_zMods = (int)(Geometry->GetZModule().size()/2);// number of zModules
	int num_x, num_z, bar_direction;
	double shortCenter, longCenter;
    for (int yModule = 0; yModule < (int)(Geometry->GetYLayers().size()/2); yModule++){
    	int layerid = (int)(Geometry->GetYFloor().size()/2) + 
					  (int)(Geometry->GetZWalls().size()/2) + yModule;
		double yCenter = (Geometry->GetYLayers()[yModule*2] + Geometry->GetYLayers()[yModule*2+1])/2;
		if (yModule % 2 == 0) {
		   	num_x = (int)(length / Geometry->GetLength());//num scintillators per module
		   	num_z = (int)(length / Geometry->GetWidth());
			bar_direction = 0;
		} else {
			num_x = (int)(length / Geometry->GetWidth());
			num_z = (int)(length / Geometry->GetLength());
			bar_direction = 2;
		}
        for (int xModule = 0; xModule < (int)(Geometry->GetXModule().size()/2); xModule++){
	    	for (int zModule = 0; zModule < (int)(Geometry->GetZModule().size()/2); zModule++){
		for (int xIndex = 0; xIndex < num_x; xIndex++) {
		    for (int zIndex = 0; zIndex < num_z; zIndex++) {
				if (yModule % 2 == 0) {
					shortCenter = zModule*(length+distance) + Geometry->GetWidth()*(zIndex+0.5);
					shortCenter += Geometry->GetZStart();
					longCenter = xModule*(length+distance) + Geometry->GetLength()*(xIndex+0.5);
					longCenter += Geometry->GetXStart();
				} else {
					shortCenter = xModule*(length+distance) + Geometry->GetWidth()*(xIndex+0.5);
					shortCenter += Geometry->GetXStart();
					longCenter = zModule*(length+distance) + Geometry->GetLength()*(zIndex+0.5);
					longCenter += Geometry->GetZStart();
		    	}
				int idNumber = xIndex + xModule*num_x + (num_xMods*num_x)*(zModule*num_z + zIndex);
				idNumber = idNumber*100 + layerid;
				detID _id;
		        if (yModule%2==0) _id = detID(yCenter, shortCenter, bar_direction, 0, idNumber);
				else _id = detID(shortCenter, yCenter, bar_direction, 0, idNumber);
				_detID_list.push_back(_id);
				_longCenters.push_back(longCenter);
			}
	    }
        }
    	}
	}
}

//------------------------------------------------------------------//

void NoiseMaker::back_detIDs(std::vector<detID>& _detID_list, std::vector<double>& _longCenters){
    if(ts) std::cout<<"NoiseMaker::back_detIDs"<<std::endl;
	double length = Geometry->GetXModule()[1] - Geometry->GetXModule()[0];//Length of module
	double distance = Geometry->GetXModule()[2] - Geometry->GetXModule()[1];//Distance between
	int num_xMods = (int)(Geometry->GetXModule().size()/2);
	int num_x, num_y, bar_direction;
	double shortCenter, longCenter;
    for (int zModule = 0; zModule < (int)(Geometry->GetZBack().size()/2); zModule++){
    	int layerid = (int)(Geometry->GetYFloor().size()/2) + 
					  (int)(Geometry->GetZWalls().size()/2) + 
					  (int)(Geometry->GetYLayers().size()/2) + zModule;
		double zCenter = (Geometry->GetZBack()[zModule*2] + Geometry->GetZBack()[zModule*2+1])/2;
		if (zModule % 2 == 0) {
		   	num_x = (int)(length / Geometry->GetLength());//num scintillators per module
			num_y = (int)(length / Geometry->GetWidth());
			bar_direction = 0;
		} else {
			num_x = (int)(length / Geometry->GetWidth());
			num_y = (int)(length / Geometry->GetLength());
			bar_direction = 1;
		}
        for (int xModule = 0; xModule < (int)(Geometry->GetXModule().size()/2); xModule++){
		for (int xIndex = 0; xIndex < num_x; xIndex++) {
		    for (int yIndex = 0; yIndex < num_y; yIndex++) {
				if (zModule % 2 == 0) {
					shortCenter = Geometry->GetWidth()*(yIndex+0.5) + Geometry->GetBackStart();
					longCenter = xModule*(length+distance) + Geometry->GetLength()*(xIndex+0.5);
					longCenter += Geometry->GetXStart();
				} else {
					shortCenter = xModule*(length+distance) + Geometry->GetWidth()*(xIndex+0.5);
					shortCenter += Geometry->GetXStart();
					longCenter = Geometry->GetLength()*(yIndex+0.5) + Geometry->GetBackStart();
		    	}
				int idNumber = xIndex + xModule*num_x + (num_xMods*num_x)*yIndex;
				idNumber = idNumber*100 + layerid;
		        detID _id = detID(shortCenter, zCenter, bar_direction, 1, idNumber);
				_detID_list.push_back(_id);
				_longCenters.push_back(longCenter);
			}
	    }
    	}
	}
}

//------------------------------------------------------------------//

void NoiseMaker::wall_detIDs(std::vector<detID>& _detID_list, std::vector<double>& _longCenters){
    if(ts) std::cout<<"NoiseMaker::wall_detIDs"<<std::endl;
	double length = Geometry->GetXEnd() - Geometry->GetXStart();//Length of wall
	double height = Geometry->GetWallEnd() - Geometry->GetWallStart();//height of wall
	int num_x, num_y, bar_direction;
	double shortCenter, longCenter;
    for (int zModule = 0; zModule < (int)(Geometry->GetZWalls().size()/2); zModule++){
    	int layerid = zModule;
		double zCenter = (Geometry->GetZWalls()[zModule*2] + Geometry->GetZWalls()[zModule*2+1])/2;
		if (zModule % 2 == 0) {
		   	num_x = (int)(length / Geometry->GetLength());//num scintillators in xdirecion
			num_y = (int)(height / Geometry->GetWidth());//num scintillators in ydirection
			bar_direction = 0;
		} else {
			num_x = (int)(length / Geometry->GetWidth());
			num_y = (int)(height / Geometry->GetLength()) + 1;
			bar_direction = 1;
			std::cout << num_y << std::endl;
		}
		for (int xIndex = 0; xIndex < num_x; xIndex++) {
			for (int yIndex = 0; yIndex < num_y; yIndex++) {
				if (zModule % 2 == 0) {
					shortCenter = Geometry->GetWidth()*(yIndex+0.5) + Geometry->GetWallStart();
					longCenter = Geometry->GetLength()*(xIndex+0.5) + Geometry->GetXStart();
				} else {
					shortCenter = Geometry->GetWidth()*(xIndex+0.5) + Geometry->GetXStart();
					longCenter = Geometry->GetLength()*(yIndex+0.5) + Geometry->GetWallStart();
		    	}
				int idNumber = xIndex + num_x*yIndex;
				idNumber = idNumber*100 + layerid;
				detID _id = detID(shortCenter, zCenter, bar_direction, 1, idNumber);
				_detID_list.push_back(_id);
				_longCenters.push_back(longCenter);
			}
		}
	}
}

//------------------------------------------------------------------//

void NoiseMaker::floor_detIDs(std::vector<detID>& _detID_list, std::vector<double>& _longCenters){
    if(ts) std::cout<<"NoiseMaker::floor_detIDs"<<std::endl;
	double length = Geometry->GetXEnd() - Geometry->GetXStart();//Length of floor
	int num_x, num_z, bar_direction;
	double shortCenter, longCenter;
    for (int yModule = 0; yModule < (int)(Geometry->GetYFloor().size()/2); yModule++){
    	int layerid = (int)(Geometry->GetZWalls().size()/2) + yModule;
		double yCenter = (Geometry->GetYFloor()[yModule*2] + Geometry->GetYFloor()[yModule*2+1])/2;
		if (yModule % 2 == 0) {
		   	num_x = (int)(length / Geometry->GetLength());//num scintillators in xdirecion
			num_z = (int)(length / Geometry->GetWidth());//num scintillators in zdirection
			bar_direction = 0;
		} else {
			num_x = (int)(length / Geometry->GetWidth());
			num_z = (int)(length / Geometry->GetLength());
			bar_direction = 2;
		}
		for (int xIndex = 0; xIndex < num_x; xIndex++) {
			for (int zIndex = 0; zIndex < num_z; zIndex++) {
				if (yModule % 2 == 0) {
					shortCenter = Geometry->GetWidth()*(zIndex+0.5) + Geometry->GetZStart();
					longCenter = Geometry->GetLength()*(xIndex+0.5) + Geometry->GetXStart();
				} else {
					shortCenter = Geometry->GetWidth()*(xIndex+0.5) + Geometry->GetXStart();
					longCenter = Geometry->GetLength()*(zIndex+0.5) + Geometry->GetZStart();
		    	}
				int idNumber = xIndex + num_x*zIndex;
				idNumber = idNumber*100 + layerid;
				detID _id;
				if (yModule%2==0) _id = detID(yCenter, shortCenter, bar_direction, 0, idNumber);
				else _id = detID(shortCenter, yCenter, bar_direction, 0, idNumber);
				_detID_list.push_back(_id);
				_longCenters.push_back(longCenter);
			}
		}
	}
}

//------------------------------------------------------------------//
