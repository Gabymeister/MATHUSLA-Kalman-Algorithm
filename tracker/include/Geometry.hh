#include "globals.hh"
#include "LinearAlgebra.hh"
#include <cmath>
#include <iostream>

#ifndef GEOMETRY_HH
#define GEOMETRY_HH

class detID{
public:
	int xModule;
	int xIndex;
	int yModule;
	int yIndex;
	int zModule;
	int zIndex;
	bool _null = false;
	bool isFloorElement = false;
    	bool isWallElement = false;
	bool isBackElement = false;

	detID(){_null = true;}

	void Print(){

		std::cout << "***************Printing DetID*************" << std::endl;
		if (_null){
			std::cout << "NOT IN KNOWN DETECTOR ELEMENT" << std::endl;
			return;
		}
		std::cout << "xModule: " << xModule << std::endl;
		std::cout << "yModule: " << yModule << std::endl;
		std::cout << "zModule: " << zModule << std::endl;
		std::cout << "x Index: " << xIndex << std::endl;
		std::cout << "y Index: " << yIndex << std::endl;
		std::cout << "z Index: " << zIndex << std::endl;
	}

	detID(int _xModule, int _x_index, int _yModule, int _y_index, int _zModule, int _z_index){
	    	xModule = _xModule;
		yModule = _yModule;
		zModule = _zModule;
		xIndex = _x_index;
		yIndex = _y_index;
		zIndex = _z_index;
		_null = false;
		//Setting floor/wall/back
		//NOTE: Only checks range!
		if (_zModule < 2) isWallElement = true;//2 walls
		else if (_yModule < 2) isFloorElement = true; //2 floors
		else if (_zModule > 2 + detector::MODULE_Z.size()) isBackElement = true;
		else {//main detector layers
			isFloorElement = false;
        		isWallElement = false;
			isBackElement = false;
		}
	}
	


	bool operator==(const detID &detID2){
		if (_null) return false;
		if (isFloorElement != detID2.isFloorElement) return false;
        if (isWallElement != detID2.isWallElement) return false;
        
	if (isFloorElement){
		if (xIndex != detID2.xIndex) return false;
		if (zIndex != detID2.zIndex) return false;
		if (yModule != detID2.yModule) return false;
		return true;
	}
        if (isWallElement) {
            if (xIndex != detID2.xIndex) return false;
            if (yIndex != detID2.yIndex) return false;
            if (zModule != detID2.zModule) return false;
            return true;
        }
	if (isBackElement) {
	    if (zModule != detID2.zModule) return false;
	    if (xModule != detID2.xModule) return false;
	    if (xIndex != detID2.xIndex) return false;
	    if (yIndex != detID2.yIndex) return false;
	}

		if (xModule != detID2.xModule) return false;
		if (zModule != detID2.zModule) return false;
		if (yModule != detID2.yModule) return false;
		if (xIndex != detID2.xIndex) return false;
		if (zIndex != detID2.zIndex) return false;
		return true;
	}
	
	bool IsNull() { return _null; }

	std::vector<double> uncertainty(){
	    if (isWallElement || isBackElement) {//vertical layers
	        if (zModule % 2 == 0) {
	            return {detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2),
		    detector::scintillator_width/sqrt(12.), detector::scintillator_thickness/sqrt(12.)};
		} else {
		    return {detector::scintillator_width/sqrt(12.),
		    detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2),
		    detector::scintillator_thickness/sqrt(12.)};
		}
	    } //horizontal layers
	    if (yModule % 2 == 0){
	        return {detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2),
		detector::scintillator_thickness/sqrt(12.), detector::scintillator_width/sqrt(12.)};
	    } else {
	        return {detector::scintillator_width/sqrt(12.), detector::scintillator_thickness/sqrt(12.),
		detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2)};
	    }
	    return {-1., -1. -1.};
	}

	//0 implies x is long direction for both walls and layers
	int GetLongIndex(){
	    if (isWallElement || isBackElement) {
	        return (zModule % 2);
	    } else {
	        return (yModule % 2);
	    }
	}
};

class Geometry{
public:

	~Geometry(){}

	//Assumes you know it is in the Floor
	detID GetDetIDFloor(double x, double y, double z){
	    int yModule = -1;
	    for (int i = 0; i < detector::n_floors; i++) {
	        if (y < detector::LAYERS_Y[i][1] && y > detector::LAYERS_Y[i][0])
		    yModule = i; 
	    }

	    if (yModule < 0) return detID(); //not in y_ranges
	    if (x < detector::x_min or x > detector::x_max) return detID();
	    if (z < detector::z_min or z > detector::z_max) return detID();

	    //indexing from far corner
	    int x_index, z_index;
	    if (yModule % 2 == 0) {
	        x_index = (int)((x - detector::x_min) / detector::scintillator_length); 
		z_index = (int)((z - detector::z_min)/ detector::scintillator_width);
	    } else {
	        x_index = (int)((x - detector::x_min) / detector::scintillator_width); 
		z_index = (int)((z - detector::z_min) / detector::scintillator_length);
	    }
	    return detID(0, x_index, yModule, 0, detector::n_walls, z_index);
	}

    	detID GetDetIDWall(double x, double y, double z) {
	    int zModule = -1;

	    for (int i = 0; i < detector::n_walls; i++) {
	        if (z < detector::FRONT_WALL_Z[i][1] && z > detector::FRONT_WALL_Z[i][0])
		    zModule = i; 
	    }

	    if (zModule < 0) return detID(); //not in wall ranges
	    if (x < detector::x_min or x > detector::x_max) return detID();
	    if (y < detector::wall_start_y or y > detector::wall_height + detector::wall_start_y) return detID();

	    //indexing from far corner
	    int x_index, y_index;
	    if (zModule % 2 == 0) {
	        x_index = (int)((x - detector::x_min) / detector::scintillator_length); 
		y_index = (int)((y - detector::wall_start_y) / detector::scintillator_width);
	    } else {
	        x_index = (int)((x - detector::x_min) / detector::scintillator_width); 
		y_index = (int)((y - detector::wall_start_y) / detector::scintillator_length);
	    }
            return detID(0, x_index, 0, y_index, zModule, 0);   
           
    	}

	detID GetDetIDBackWall(double x, double y, double z) {
	    int zModule = -1;
	    int xModule = -1;
	     
	    for (int i = 0; i < detector::n_back_walls; i++) {
	        if (z < detector::BACK_WALL_Z[i][1] && z > detector::BACK_WALL_Z[i][0])
		    zModule = i + detector::MODULE_Z.size() + detector::n_walls; //z_module starts front walls 
	    }
	    for (int i = 0; i < detector::MODULE_X.size(); i++) {
	        if (x < detector::MODULE_X[i][1] && x > detector::MODULE_X[i][0])
		    xModule = i; 
	    }

	    if (zModule < 0 || xModule < 0) return detID(); //not in wall ranges
	    if (y < detector::wall_start_y or y > detector::wall_height + detector::wall_start_y) return detID();

	    //indexing from far corner
	    int x_index, y_index;
	    if (zModule % 2 == 0) {
	        x_index = (int)((x - detector::MODULE_X[xModule][0]) / detector::scintillator_length); 
		y_index = (int) ((y - detector::back_wall_bottom)/ detector::scintillator_width);
	    } else {
	        x_index = (int)((x - detector::MODULE_X[xModule][0]) / detector::scintillator_width); 
		y_index = (int)((y - detector::back_wall_bottom) / detector::scintillator_length);
	    }
            return detID(xModule, x_index, 0, y_index, zModule, 0);   
           
    	}


	detID GetDetIDTrackingLayer(double x, double y, double z) {
	    int zModule = -1;
	    int xModule = -1;
	    int yModule = -1;
	     
	    for (int i = 0; i < detector::MODULE_Z.size(); i++) {
	        if (z < detector::MODULE_Z[i][1] && z > detector::BACK_WALL_Z[i][0])
		    zModule = i + detector::n_walls; //z_module starts front walls 
	    }
	    for (int i = 0; i < detector::MODULE_X.size(); i++) {
	        if (x < detector::MODULE_X[i][1] && x > detector::MODULE_X[i][0])
		    xModule = i; 
	    }
	    for (int i = detector::n_floors; i < detector::n_layers; i++) {
	    	if (y < detector::LAYERS_Y[i][1] && y > detector::LAYERS_Y[i][0])
		    yModule = i;
	    }

	    if (zModule < 0 || xModule < 0 || yModule < 0) return detID(); //not in wall ranges

	    //indexing from far corner
	    int x_index, z_index;
	    if (yModule % 2 == 0) {
	        x_index = (int)((x - detector::MODULE_X[xModule][0]) / detector::scintillator_length); 
		z_index = (int)((z - detector::MODULE_Z[zModule - detector::n_walls][0])/detector::scintillator_width);
	    } else {
	        x_index = (int)((x - detector::MODULE_X[xModule][0]) / detector::scintillator_width); 
		z_index = (int)((z - detector::MODULE_Z[zModule - detector::n_walls][0])/detector::scintillator_length);
	    }
            return detID(xModule, x_index, yModule, 0, zModule, z_index);   
           
    	}

	detID GetDetID(double x, double y, double z){

            if (z < detector::z_min) return GetDetIDWall(x, y, z);	
	    else if (z > detector::z_max) return GetDetIDBackWall(x,y,z);
	    else if (y < detector::LAYERS_Y[detector::n_floors][0]) return GetDetIDFloor(x,y,z);
	    else return GetDetIDTrackingLayer(x,y,z);

	}

	//pointer Hit
	template<class Hit>
	detID GetDetID(Hit _hit){
		return GetDetID(_hit->x, _hit->y, _hit->z);
	}

	//for vector::Vector of position
	detID GetDetID(vector::Vector _hit){
		return GetDetID(_hit.x, _hit.y, _hit.z);
	}

	detID GetDetID(std::vector<double> _hit){
		return GetDetID(_hit[0], _hit[1], _hit[2]);
	}

	std::vector<double> GetCenterFloor(detID _id){
	    double x,y,z;
	    y = detector::LAYERS_Y[_id.yModule][0] + detector::scintillator_height / 2;
	    if (_id.yModule % 2 == 0) {
	        x = (0.5 + _id.xIndex)*detector::scintillator_length + detector::x_min;
		z = (0.5 + _id.zIndex)*detector::scintillator_width + detector::z_min;
	    } else {
	        x = (0.5 + _id.xIndex)*detector::scintillator_width + detector::x_min;
		z = (0.5 + _id.zIndex)*detector::scintillator_length + detector::z_min;
	    } 
	    return {x,y,z};
	}

    	std::vector<double> GetCenterWall(detID _id){
	    double x,y,z;
	    z = detector::FRONT_WALL_Z[_id.zModule][0] + detector::scintillator_height / 2;
	    if (_id.zModule % 2 == 0) {
	        x = (0.5 + _id.xIndex)*detector::scintillator_length + detector::x_min;
		y = (0.5 + _id.yIndex)*detector::scintillator_width + detector::wall_start_y;
	    } else {
	        x = (0.5 + _id.xIndex)*detector::scintillator_width + detector::x_min;
		y = (0.5 + _id.yIndex)*detector::scintillator_length + detector::wall_start_y;
	    } 
	    return {x,y,z};
    	}

    	std::vector<double> GetCenterBackWall(detID _id){
	    double x,y,z;
	    int BackzIndex = _id.zModule - detector::MODULE_Z.size() - detector::n_walls;
	    z = detector::BACK_WALL_Z[BackzIndex][0] + detector::scintillator_height / 2;
	    if (_id.zModule % 2 == 0) {
	        x = (0.5 + _id.xIndex)*detector::scintillator_length + detector::MODULE_X[_id.xModule][0];
		y = (0.5 + _id.yIndex)*detector::scintillator_width + detector::back_wall_bottom;
	    } else {
	        x = (0.5 + _id.xIndex)*detector::scintillator_width + detector::MODULE_X[_id.xModule][0];
		y = (0.5 + _id.yIndex)*detector::scintillator_length + detector::back_wall_bottom;
	    } 
	    return {x,y,z};
    	}

    	std::vector<double> GetCenterTrackingLayer(detID _id){
	    double x,y,z;
	    int TrackingzIndex = _id.zModule - detector::n_walls;
	    y = detector::LAYERS_Y[_id.yModule][0] + detector::scintillator_height / 2;
	    if (_id.yModule % 2 == 0) {
	        x = (0.5 + _id.xIndex)*detector::scintillator_length + detector::MODULE_X[_id.xModule][0];
		z = (0.5 + _id.zIndex)*detector::scintillator_width + detector::MODULE_Z[TrackingzIndex][0];
	    } else {
	        x = (0.5 + _id.xIndex)*detector::scintillator_width + detector::MODULE_X[_id.xModule][0];
		z = (0.5 + _id.zIndex)*detector::scintillator_length + detector::MODULE_Z[TrackingzIndex][0];
	    } 
	    return {x,y,z};
    	}

	std::vector<double> GetCenter(detID _id){
	    if (_id.IsNull()) {
		std::cout << "detID is null" << std::endl;
		return {};
	    }
	    if (_id.isWallElement) return GetCenterWall(_id);
	    if (_id.isFloorElement) return GetCenterFloor(_id);
	    if (_id.isBackElement) return GetCenterBackWall(_id);
	    return GetCenterTrackingLayer(_id);
	}


	bool inBox(double x, double y, double z) {

		if (detector::x_min < x < detector::x_max) {
			if (detector::y_min < y < detector::y_max) {
				if (detector::z_min < z < detector::z_max) {
					return true;
				};
			};
		};

		return false;
	}

	std::vector<double> GetDimensions(detID _id){
	    double x,y,z;
	    if (_id.xModule % 2 == 0) {//x is long
	        if (_id.isWallElement || _id.isBackElement) {//vertical layers
		    return {detector::scintillator_length, detector::scintillator_width, detector::scintillator_thickness};
		} else {//horizontal layers
		    return {detector::scintillator_length, detector::scintillator_thickness, detector::scintillator_width};
		}
	    } else {
	        if (_id.isWallElement || _id.isBackElement) {//vertical layers
		    return {detector::scintillator_width, detector::scintillator_length, detector::scintillator_thickness};
		} else {//horizontal layers
		    return {detector::scintillator_width, detector::scintillator_thickness, detector::scintillator_length};
		}
	    }
	}





}; //Geometry Class




#endif
