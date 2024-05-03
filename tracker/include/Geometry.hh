#include "globals.hh"
#include <cmath>
#include <iostream>

#ifndef GEOMETRY_HH
#define GEOMETRY_HH

class detID{
public:
	double center1;
	double center2;
	//0-x, 1-y, 2-z
	int bar_direction;
	// 0 is layer, 1 is wall
	int allignment;
	int detectorID;
	//front_walls, floor, layers, back wall
	int layerID;
	bool _null;

	void Print(){

		std::cout << "***************Printing DetID*************" << std::endl;
		if (_null){
			std::cout << "NOT IN KNOWN DETECTOR ELEMENT" << std::endl;
			return;
		}
		if (bar_direction == 0) {
			std::cout << "yCenter: " << center1  << std::endl;
			std::cout << "zCenter: " << center2 << std::endl;
			std::cout << "bar direction: x" << std::endl;
		} else if (bar_direction == 1) {
			std::cout << "xCenter: " << center1  << std::endl;
			std::cout << "zCenter: " << center2 << std::endl;
			std::cout << "bar direction: y" << std::endl;
		} else {
			std::cout << "xCenter: " << center1  << std::endl;
			std::cout << "yCenter: " << center2 << std::endl;
			std::cout << "bar direction: z" << std::endl;
		}
		if (allignment == 0) {
			std::cout << "Wall Hit" << std::endl;
		} else { 
			std::cout << "Layer Hit" << std::endl;
		}
	}

	detID() {_null = true;}

	detID(double _center1, double _center2, int _bar_direction, int _allignment, int _detectorID){
		center1 = _center1;
		center2 = _center2;
		bar_direction = _bar_direction;
		allignment = _allignment;
		detectorID = _detectorID;
		layerID = _detectorID % 100;
	}

	bool IsNull() {return _null;}

	bool operator==(const detID &detID2){
		return (detectorID == detID2.detectorID);
	}
	
	std::vector<double> uncertainty(){
	    if (allignment) {//vertical layers
	        if (bar_direction == 0) { // x long
	            return {detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2),
		    detector::scintillator_width/sqrt(12.), detector::scintillator_thickness/sqrt(12.)};
			} else {
		    	return {detector::scintillator_width/sqrt(12.),
		    	detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2),
		    	detector::scintillator_thickness/sqrt(12.)};
			}
	    } //horizontal layers
	    if (bar_direction == 0){
	        return {detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2),
		detector::scintillator_thickness/sqrt(12.), detector::scintillator_width/sqrt(12.)};
	    } else {
	        return {detector::scintillator_width/sqrt(12.), detector::scintillator_thickness/sqrt(12.),
		detector::time_resolution*(constants::c/constants::optic_fiber_n)/sqrt(2)};
	    }
	}

	//0 implies x is long direction for both walls and layers
	bool isVertical(){
		return allignment == 1;
	}

	std::vector<double> GetCenter(){
		return {center1, center2};
	}

	int GetLongIndex() {
		return bar_direction;
	}
};




#endif
