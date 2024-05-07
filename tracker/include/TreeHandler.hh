#include <TTree.h>
#include <TFile.h>
#include <TROOT.h>
#include <Eigen/Dense>
#include "globals.hh"

#ifndef TH_DEFINE
#define TH_DEFINE

class TreeHandler{ 
public:
	//PUT ALL INPUT AND OUTPUT BRANCHES HERE
	TTree* OutputTree;
	TTree* InputTree;
	TFile* OutputFile;
	TFile* InputFile;
	int index = -1;
	int NumEntries;
	bool _Null = false;
	bool IsNull(){return _Null;}

	int Next(){
    index++;
		if (index >= NumEntries) return -1;
		return index;
	}

	int LoadEvent(){
		if (InputTree == nullptr) return -1;
		InputTree->GetEvent(index);
		return 0;
	}

	void Fill(){
		OutputFile->cd();
		gROOT->cd();
		InputTree->GetEvent(index);
		OutputTree->Fill();
	}

	void Write(){
    InputFile->Close();
		OutputFile->cd();
		OutputFile->Write();
		OutputFile->Close();
	}

	template<class digi_hit>
	void ExportDigis(std::vector<digi_hit*>, long long int digi_seed);

	TreeHandler(TString input_tree_name, TString input_file_name, TString output_tree_name, TString outfile_name) {
		InputFile = TFile::Open(input_file_name);
		if (! InputFile) {
      		_Null = true;
      		return;
    	}

		InputTree = (TTree*) InputFile->Get(input_tree_name);

    	if (! InputTree){
      		_Null = true;
      		return;
    	}

		InputTree->SetBranchAddress("NumHits", &sim_numhits);
 		InputTree->SetBranchAddress("Hit_energy", &sim_hit_e);
 		InputTree->SetBranchAddress("Hit_time", &sim_hit_t);
 		InputTree->SetBranchAddress("Hit_Center1", &sim_hit_center1);
 		InputTree->SetBranchAddress("Hit_Center2", &sim_hit_center2);
 		InputTree->SetBranchAddress("Hit_bar_direction", &sim_hit_bar_direction);
 		InputTree->SetBranchAddress("Hit_allignment", &sim_hit_allignment);
 		InputTree->SetBranchAddress("Hit_LayerID", &sim_hit_layerID);

 		InputTree->SetBranchAddress("Hit_particlePdgId", &sim_hit_particlePdgId);
 		InputTree->SetBranchAddress("Hit_G4TrackId", &sim_hit_G4TrackId);
 		InputTree->SetBranchAddress("Hit_G4ParentTrackId", &sim_hit_G4ParentTrackId);
 		InputTree->SetBranchAddress("Hit_x", &sim_hit_x);
 		InputTree->SetBranchAddress("Hit_y", &sim_hit_y);
 		InputTree->SetBranchAddress("Hit_z", &sim_hit_z);
 		InputTree->SetBranchAddress("Hit_particleEnergy", &sim_hit_particleEnergy);
 		InputTree->SetBranchAddress("Hit_particlePx", &sim_hit_px);
 		InputTree->SetBranchAddress("Hit_particlePy", &sim_hit_py);
		InputTree->SetBranchAddress("Hit_particlePz", &sim_hit_pz);
		InputTree->SetBranchAddress("Hit_weight", &sim_hit_weight);
		InputTree->SetBranchAddress("GenParticle_energy", &sim_GenParticle_energy);
		InputTree->SetBranchAddress("GenParticle_pdgid", &sim_GenParticle_pdgid);
		InputTree->SetBranchAddress("GenParticle_index", &sim_GenParticle_index);
		InputTree->SetBranchAddress("GenParticle_G4index", &sim_GenParticle_G4index);
	 	InputTree->SetBranchAddress("GenParticle_time", &sim_GenParticle_time);
	 	InputTree->SetBranchAddress("GenParticle_x", &sim_GenParticle_x);
	 	InputTree->SetBranchAddress("GenParticle_y", &sim_GenParticle_y);
	 	InputTree->SetBranchAddress("GenParticle_z", &sim_GenParticle_z);
	 	InputTree->SetBranchAddress("GenParticle_px", &sim_GenParticle_px);
	 	InputTree->SetBranchAddress("GenParticle_py", &sim_GenParticle_py);
	 	InputTree->SetBranchAddress("GenParticle_pz", &sim_GenParticle_pz);
	 	InputTree->SetBranchAddress("GenParticle_mass", &sim_GenParticle_mass);


   		InputTree->SetBranchStatus("NumGenParticles", 0);
   		InputTree->SetBranchStatus("GenParticle_status", 0);
   		InputTree->SetBranchStatus("GenParticle_pt", 0);
   		InputTree->SetBranchStatus("GenParticle_eta", 0);
   		InputTree->SetBranchStatus("GenParticle_phi", 0);

   		InputTree->SetBranchStatus("GenParticle_mo1", 0);
   		InputTree->SetBranchStatus("GenParticle_mo2", 0);
   		InputTree->SetBranchStatus("GenParticle_dau1", 0);
   		InputTree->SetBranchStatus("GenParticle_dau2", 0);
		InputTree->SetBranchStatus("COSMIC_EVENT_ID", 0);
 		InputTree->SetBranchStatus("COSMIC_CORE_X", 0);
 		InputTree->SetBranchStatus("COSMIC_CORE_Y", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_PRIMARY_ENERGY", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_THETA", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_PHI", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_FIRST_HEIGHT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_ELECTRON_COUNT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_MUON_COUNT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_HADRON_COUNT", 0);
 		InputTree->SetBranchStatus("COSMIC_GEN_PRIMARY_ID", 0);
 		InputTree->SetBranchAddress("EXTRA_12", &sim_EXTRA_12);
 		InputTree->SetBranchAddress("EXTRA_13", &sim_EXTRA_13);
 		InputTree->SetBranchAddress("EXTRA_14", &sim_EXTRA_14);
 		InputTree->SetBranchAddress("EXTRA_15", &sim_EXTRA_15);

 		NumEntries = InputTree->GetEntries();

 		OutputFile = new TFile(outfile_name, "RECREATE");
		OutputTree = new TTree(output_tree_name, "MATHUSLA Tree");

		OutputTree->Branch("NumHits", &sim_numhits);
		OutputTree->Branch("Hit_energy", "std::vector<double>", sim_hit_e);
 		OutputTree->Branch("Hit_time", "std::vector<double>", sim_hit_t);
 		OutputTree->Branch("Hit_Center1", sim_hit_center1);
 		OutputTree->Branch("Hit_Center2", sim_hit_center2);
 		OutputTree->Branch("Hit_bar_direction", sim_hit_bar_direction);
 		OutputTree->Branch("Hit_allignment", sim_hit_allignment);
 		OutputTree->Branch("Hit_LayerID", sim_hit_layerID);
 		OutputTree->Branch("Hit_particlePdgId", "std::vector<double>", sim_hit_particlePdgId);
 		OutputTree->Branch("Hit_G4TrackId", "std::vector<double>", sim_hit_G4TrackId);
 		OutputTree->Branch("Hit_G4ParentTrackId", "std::vector<double>", sim_hit_G4ParentTrackId);
 		OutputTree->Branch("Hit_x", "std::vector<double>", sim_hit_x);
 		OutputTree->Branch("Hit_y", "std::vector<double>", sim_hit_y);
 		OutputTree->Branch("Hit_z", "std::vector<double>", sim_hit_z);
 		OutputTree->Branch("Hit_particleEnergy", "std::vector<double>", sim_hit_particleEnergy);
 		OutputTree->Branch("Hit_particlePx", "std::vector<double>", sim_hit_px);
 		OutputTree->Branch("Hit_particlePy", "std::vector<double>", sim_hit_py);
		OutputTree->Branch("Hit_particlePz", "std::vector<double>", sim_hit_pz);
 
		OutputTree->Branch("Digi_numHits", &Digi_numHits);
		OutputTree->Branch("Digi_time", &digi_hit_t);
        OutputTree->Branch("Digi_x", &digi_hit_x);
        OutputTree->Branch("Digi_y", &digi_hit_y);
        OutputTree->Branch("Digi_z", &digi_hit_z);
        OutputTree->Branch("Digi_energy", &digi_hit_e);
        OutputTree->Branch("Digi_px", &digi_hit_px);
        OutputTree->Branch("Digi_py", &digi_hit_py);
        OutputTree->Branch("Digi_pz", &digi_hit_pz);
        OutputTree->Branch("Digi_particle_energy", &digi_particle_energy);
        OutputTree->Branch("Digi_pdg_id", &digi_pdg);
        OutputTree->Branch("Digi_track_id", &digi_track_id);

		OutputTree->Branch("Digi_center1", &digi_center1);
		OutputTree->Branch("Digi_center2", &digi_center2);
		OutputTree->Branch("Digi_bar_direction", &digi_bar_direction);
		OutputTree->Branch("Digi_allignment", &digi_allignment);
		OutputTree->Branch("Digi_layer_id", &digi_LayerID);
		OutputTree->Branch("Digi_det_id", &digi_detID);
        OutputTree->Branch("Digi_seed", &digi_seed, "Digi_seed/L");
        OutputTree->Branch("Digi_hitIndices", &digi_hit_indices);

 		OutputTree->Branch("GenParticle_index", "std::vector<double>", sim_GenParticle_index);
 		OutputTree->Branch("GenParticle_G4index", "std::vector<double>", sim_GenParticle_G4index);
 		OutputTree->Branch("GenParticle_pdgid", "std::vector<double>", sim_GenParticle_pdgid);
 		OutputTree->Branch("GenParticle_time", "std::vector<double>", sim_GenParticle_time);
 		OutputTree->Branch("GenParticle_x", "std::vector<double>", sim_GenParticle_x);
 		OutputTree->Branch("GenParticle_y", "std::vector<double>", sim_GenParticle_y);
 		OutputTree->Branch("GenParticle_z", "std::vector<double>", sim_GenParticle_z);
 		OutputTree->Branch("GenParticle_energy", "std::vector<double>", sim_GenParticle_energy);
 		OutputTree->Branch("GenParticle_px", "std::vector<double>", sim_GenParticle_px);
 		OutputTree->Branch("GenParticle_py", "std::vector<double>", sim_GenParticle_py);
 		OutputTree->Branch("GenParticle_pz", "std::vector<double>", sim_GenParticle_pz);
 		OutputTree->Branch("GenParticle_mass", "std::vector<double>", sim_GenParticle_mass);
 		OutputTree->Branch("EXTRA_12", "std::vector<double>", sim_EXTRA_12);
 		OutputTree->Branch("EXTRA_13", "std::vector<double>", sim_EXTRA_13);
 		OutputTree->Branch("EXTRA_14", "std::vector<double>", sim_EXTRA_14);
 		OutputTree->Branch("EXTRA_15", "std::vector<double>", sim_EXTRA_15);
	}
//____________________________________________________________________________________________

//___Make Sim Branches_________________________________________________________________________
 	Double_t sim_numhits;
 	std::vector<double> *sim_hit_e = nullptr;
 	std::vector<double> *sim_hit_t = nullptr;
 	std::vector<double> *sim_hit_center1 = nullptr; 
 	std::vector<double> *sim_hit_center2 = nullptr; 
 	std::vector<double> *sim_hit_bar_direction = nullptr; 
 	std::vector<double> *sim_hit_allignment = nullptr;
 	std::vector<double> *sim_hit_layerID = nullptr;
 	std::vector<double> *sim_hit_particlePdgId = nullptr;
 	std::vector<double> *sim_hit_G4TrackId = nullptr;
 	std::vector<double> *sim_hit_G4ParentTrackId = nullptr;
 	std::vector<double> *sim_hit_x = nullptr;
 	std::vector<double> *sim_hit_y = nullptr;
 	std::vector<double> *sim_hit_z = nullptr;
 	std::vector<double> *sim_hit_particleEnergy = nullptr;
 	std::vector<double> *sim_hit_px = nullptr;
 	std::vector<double> *sim_hit_py = nullptr;
 	std::vector<double> *sim_hit_pz = nullptr;
 	std::vector<double> *sim_hit_weight = nullptr;
 	Double_t sim_NumGenParticles;
 	std::vector<double> *sim_GenParticle_index = nullptr;
 	std::vector<double> *sim_GenParticle_G4index = nullptr;
 	std::vector<double> *sim_GenParticle_pdgid = nullptr;
 	std::vector<double> *sim_GenParticle_status = nullptr;
 	std::vector<double> *sim_GenParticle_time = nullptr;
 	std::vector<double> *sim_GenParticle_x = nullptr;
 	std::vector<double> *sim_GenParticle_y = nullptr;
 	std::vector<double> *sim_GenParticle_z = nullptr;
 	std::vector<double> *sim_GenParticle_energy = nullptr;
 	std::vector<double> *sim_GenParticle_px = nullptr;
 	std::vector<double> *sim_GenParticle_py = nullptr;
 	std::vector<double> *sim_GenParticle_pz = nullptr;
 	std::vector<double> *sim_GenParticle_mo1 = nullptr;
 	std::vector<double> *sim_GenParticle_mo2 = nullptr;
 	std::vector<double> *sim_GenParticle_dau1 = nullptr;
 	std::vector<double> *sim_GenParticle_dau2 = nullptr;
 	std::vector<double> *sim_GenParticle_mass = nullptr;
 	std::vector<double> *sim_GenParticle_pt = nullptr;
 	std::vector<double> *sim_GenParticle_eta = nullptr;
 	std::vector<double> *sim_GenParticle_phi = nullptr;
	std::vector<double> *sim_COSMIC_EVENT_ID = nullptr;
 	std::vector<double> *sim_COSMIC_CORE_X = nullptr;
 	std::vector<double> *sim_COSMIC_CORE_Y = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_PRIMARY_ENERGY = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_THETA = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_PHI = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_FIRST_HEIGHT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_ELECTRON_COUNT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_MUON_COUNT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_HADRON_COUNT = nullptr;
 	std::vector<double> *sim_COSMIC_GEN_PRIMARY_ID = nullptr;
 	double_t G4KLong_code;
 	std::vector<double> *sim_EXTRA_12 = nullptr;
 	std::vector<double> *sim_EXTRA_13 = nullptr;
 	std::vector<double> *sim_EXTRA_14 = nullptr;
 	std::vector<double> *sim_EXTRA_15 = nullptr;

  //___Make Digi Branches_____________________________________________________________________
  	std::vector<double> digi_hit_t;
  	std::vector<double> digi_hit_x;
  	std::vector<double> digi_hit_y;
  	std::vector<double> digi_hit_z;
  	std::vector<double> digi_hit_e;
  	std::vector<double> digi_hit_px;
  	std::vector<double> digi_hit_py;
  	std::vector<double> digi_hit_pz;
    std::vector<double> digi_particle_energy;
    std::vector<int> digi_pdg;
    std::vector<int> digi_track_id;
	std::vector<double> digi_center1;
	std::vector<double> digi_center2;
	std::vector<int> digi_bar_direction;
	std::vector<int> digi_allignment;
	std::vector<int> digi_LayerID;
	std::vector<int> digi_detID;
  	std::vector<int> digi_hit_indices;
  	std::vector<int> Digi_numHits;
    long long int digi_seed;

}; //class TreeHandler

template<class digi_hit>
void TreeHandler::ExportDigis(std::vector<digi_hit*> digi_list, long long int seed){
	digi_hit_indices.clear();
	digi_hit_t.clear();
    digi_hit_x.clear();
    digi_hit_y.clear();
    digi_hit_z.clear();
    digi_hit_e.clear();
      digi_hit_px.clear();
      digi_hit_py.clear();
      digi_hit_pz.clear();
      Digi_numHits.clear();
      digi_particle_energy.clear();
      digi_pdg.clear();
	  digi_track_id.clear();
	  digi_center1.clear();
	  digi_center2.clear();
	  digi_bar_direction.clear();
	  digi_allignment.clear();
	  digi_LayerID.clear();
	  digi_detID.clear();
      
      digi_seed = seed;

      for (auto digi : digi_list){
        Digi_numHits.push_back(digi->hits.size());
        digi_hit_t.push_back(digi->t);
        digi_hit_x.push_back(digi->x);
        digi_hit_y.push_back(digi->y);
        digi_hit_z.push_back(digi->z);
        digi_hit_e.push_back(digi->e);
        digi_hit_px.push_back(digi->px);
        digi_hit_py.push_back(digi->py);
        digi_hit_pz.push_back(digi->pz);
        digi_particle_energy.push_back(digi->particle_energy);
        digi_pdg.push_back(digi->pdg);
		digi_track_id.push_back(digi->min_track_id);
		digi_center1.push_back(digi->det_id.center1);
		digi_center2.push_back(digi->det_id.center2);
		digi_bar_direction.push_back(digi->det_id.bar_direction);
		digi_allignment.push_back(digi->det_id.allignment);
		digi_LayerID.push_back(digi->det_id.layerID);
		digi_detID.push_back(digi->det_id.detectorID);
        for (auto hit : digi->hits){
          digi_hit_indices.push_back(hit->index);
        }
      }

 }

class GeometryHandler { 
public:
	//PUT ALL INPUT AND OUTPUT BRANCHES HERE
	TTree* InputTree;
	TFile* InputFile;
	int index = -1;
	int NumEntries;
	bool _Null = false;

//__BRANCHES______________________________________________________________________
	//All in CMS coordinates
	//starts and ends of the layer dimensions
	Double_t X_START;
	Double_t X_END;
	Double_t Z_START;
	Double_t Z_END;
	Double_t Wall_START;
	Double_t Wall_END;
	Double_t BACK_START;
	Double_t BACK_END;
	Double_t LENGTH;
	Double_t WIDTH;
	//starts/ends of each module in said dimension
	std::vector<double> *MODULE_X = nullptr;
	std::vector<double> *LAYERS_Y = nullptr;
	std::vector<double> *MODULE_Z = nullptr;
	std::vector<double> *FLOOR_Y = nullptr;
	std::vector<double> *WALL_Z = nullptr;
	std::vector<double> *BACK_Z = nullptr;


	bool IsNull(){return _Null;}

	int LoadGeometry(){
		if (InputTree == nullptr) return -1;
		InputTree->GetEntry(0);
		return 0;
	}

	double GetXStart(){return X_START;}
	double GetXEnd(){return X_END;}
	double GetZStart(){return Z_START;}
	double GetZEND(){return Z_END;}
	double GetWallStart(){return Wall_START;}
	double GetWallEnd(){return Wall_END;}
	double GetBackStart(){return BACK_START;}
	double GetBackEnd(){return BACK_END;}
	double GetLength(){return LENGTH;}
	double GetWidth(){return WIDTH;}
	std::vector<double>& GetXModule(){return *MODULE_X;}
	std::vector<double>& GetYLayers(){return *LAYERS_Y;}
	std::vector<double>& GetZModule(){return *MODULE_Z;}
	std::vector<double>& GetYFloor(){return *FLOOR_Y;}
	std::vector<double>& GetZWalls(){return *WALL_Z;}
	std::vector<double>& GetZBack(){return *BACK_Z;}

	GeometryHandler(TString input_tree_name, TString input_file_name) {
		InputFile = TFile::Open(input_file_name);
		if (! InputFile) {
      		_Null = true;
      		return;
    	}

		InputTree = (TTree*) InputFile->Get(input_tree_name);

    	if (! InputTree){
      		_Null = true;
      		return;
    	} 

		InputTree->SetBranchAddress("X_START", &X_START);
 		InputTree->SetBranchAddress("X_END", &X_END);
 		InputTree->SetBranchAddress("Z_START", &Z_START);
 		InputTree->SetBranchAddress("Z_END", &Z_END);
 		InputTree->SetBranchAddress("Wall_START", &Wall_START);
 		InputTree->SetBranchAddress("Wall_END", &Wall_END);
 		InputTree->SetBranchAddress("BACK_START", &BACK_START);
 		InputTree->SetBranchAddress("BACK_END", &BACK_END);
		InputTree->SetBranchAddress("LENGTH", &LENGTH);
		InputTree->SetBranchAddress("WIDTH", &WIDTH);

		InputTree->SetBranchAddress("MODULE_X", &MODULE_X);
		InputTree->SetBranchAddress("LAYERS_Y", &LAYERS_Y);
		InputTree->SetBranchAddress("MODULE_Z", &MODULE_Z);
		InputTree->SetBranchAddress("FLOOR_Y", &FLOOR_Y);
		InputTree->SetBranchAddress("WALL_Z", &WALL_Z);
		InputTree->SetBranchAddress("BACK_Z", &BACK_Z);
	}

};

#endif
