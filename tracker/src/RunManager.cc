#include <iostream>
#include "RunManager.hh"
#include "TreeHandler.hh"
#include "NoiseMaker.hh"
#include "Digitizer.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include "par_handler.hh"

int RunManager::StartTracking()
{
	TreeHandler _handler(_InputTree_Name, _InputFile_Name, _OutputTree_Name, OutFileName());
	if (_handler.IsNull()) {
		std::cout << "Sorry, I couldn't open that file" << std::endl;
		return 0;
	} 

	GeometryHandler _geometry(_geomTree_Name, _InputFile_Name);
	_geometry.LoadGeometry();

	TH = &_handler;
	int events_handled = 0;

	int made_its_k = 0;

	int events_w_k_tracks = 0;

	int verts_k_m = 0;

	int dropped_hits = 0;
	int floor_wall_hits = 0;

	std::vector<int> zeros(9, 0);
	std::vector<int> failures_k = zeros;

	ParHandler hndlr;
	hndlr.Handle();

	if (hndlr.par_map["branch"] == 1.0) std::cout << "Running in Cosmic Mode" << std::endl;

	_digitizer->par_handler = &hndlr;

	NoiseMaker::preDigitizer(&_geometry);
	std::cout << "Finished preDigitizer" <<std::endl;

	_digitizer->InitGenerators();

	while (TH->Next() >= 0)
	{
		if (events_handled >= hndlr.par_map["end_ev"]) //cuts::end_ev)
		{
			break;
		}
		if (events_handled >= hndlr.par_map["start_ev"]) //cuts::start_ev)
		{

			if (hndlr.par_map["debug"] == 1) 
				std::cout << "\n"<< std::endl;
			if ((events_handled) % 1000 == 0 || hndlr.par_map["debug"] == 1 || hndlr.par_map["debug_vertex"] == 1 )
				std::cout << "=== Event is " << events_handled <<" ==="<< std::endl;

			TotalEventsProcessed++;
			_digitizer->clear();

			// copying the data to the new tree, and loading all the variables, incrementing index
			TH->LoadEvent();

			//adding all hits of the tree into the digitizer
			for (int n_hit = 0; n_hit < TH->sim_numhits; n_hit++)
			{
				physics::sim_hit *current = new physics::sim_hit(TH, n_hit);
				if (hndlr.par_map["branch"] == 1.0) {
					current->x += detector::COSMIC_SHIFT[0];
					current->y += detector::COSMIC_SHIFT[1];
					current->z += detector::COSMIC_SHIFT[2];
				}
				_digitizer->AddHit(current);
			}

			_digitizer->ev_num = events_handled;
			std::vector<physics::digi_hit *> digi_list = _digitizer->Digitize();
			
			std::vector<physics::digi_hit*> noise_digis;
			if(NoiseMaker::run){
                NoiseMaker* noise = new NoiseMaker(digi_list);
                noise_digis = noise->return_digis();
				if (noise_digis.size() > 0) {
				}
                for(auto digi:noise_digis){
                        digi_list.push_back(digi);
                }
        	}

			TH->ExportDigis(noise_digis, _digitizer->seed);

			TH->Fill();

			dropped_hits += _digitizer->dropped_hits;
			floor_wall_hits += _digitizer->floor_wall_hits;

		}

		events_handled++;
	}

	


	if (hndlr.file_opened) {

		TH->Write();

	}

	return 0;
}
