#include "TestModule.h"


// Preprocessor Options
#define COUT_TIMESTAMP 0
#define COUT_TIMESTAMP_WAIT 0  


// Standard library:
#include <sstream>
#include <stdexcept>
#include <cmath>


// Third party:
// - Bayeux/datatools:
//#include <datatools/service_manager.h>
// - Bayeux/geomtools:
//#include <geomtools/geometry_service.h>
//#include <geomtools/manager.h>
// - Bayeux/mctools:
#include <mctools/simulated_data.h>
//#include <mctools/utils.h>

// This project :
#include <falaise/snemo/datamodels/data_model.h>
//#include <falaise/snemo/processing/services.h>


// my own lib
#define UID_ENABLE 1
#if UID_ENABLE
    #include "/home/ecb/uid-assembler/uid_assembler.hpp"
#endif

// gen bb
#if CAFFE_ENABLE
    #include <bayeux/genbb_help/primary_particle.h>
#endif

// caffe
///home/ecb/caffe/include
#if CAFFE_ENABLE
    #include "caffe/util/io.hpp"
    #include "caffe/proto/caffe.pb.h"
    #define CPU_ONLY
#endif

DPP_MODULE_REGISTRATION_IMPLEMENT(TestModule, "TestModule")

bool TestModule::has_external_random() const {
    return _external_random_ != 0;
}

void TestModule::reset_external_random() {
    DT_THROW_IF(is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is already initialized ! ");
    _external_random_ = 0;
    return;
}

void TestModule::set_external_random(mygsl::rng& rng_) {
    DT_THROW_IF(is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is already initialized ! ");
    _external_random_ = &rng_;
    return;
}

mygsl::rng& TestModule::_get_random() {
    if (has_external_random()) return *_external_random_;
    return _random_;
}

void TestModule::initialize(const datatools::properties& myConfig,
                            datatools::service_manager& flServices,
                            dpp::module_handle_dict_type& moduleDict)
{
    std::cout << "initialize" << std::endl;
    
    DT_THROW_IF(is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is already initialized ! ");

    this->base_module::_common_initialize(/*setup_*/myConfig);
    
    if (_CD_label_.empty()) {
        if (/*setup_*/myConfig.has_key("CD_label")) {
            _CD_label_ = /*setup_*/myConfig.fetch_string("CD_label");
        }
    }
    if (_CD_label_.empty()) {
        _CD_label_ = snemo::datamodel::data_info::default_calibrated_data_label();
    }
    
    if (!has_external_random()) {
        int random_seed = 12345;
        if (/*setup_*/myConfig.has_key("random.seed")) {
            random_seed = /*setup_*/myConfig.fetch_integer("random.seed");
        }
        std::string random_id = "mt19937";
        if (/*setup_*/myConfig.has_key("random.id")) {
            random_id = /*setup_*/myConfig.fetch_string("random.id");
        }
        // Initialize the embedded random number generator:
        _random_.init(random_id, random_seed);
    }
    
    // Configuable data member
    // Extract the filename_out key from the supplied config, if
    // the key exists. datatools::properties throws an exception if
    // the key isn't in the config, so catch this if thrown and don't do
    // anything
    #if PYTHON_ANALYSIS
        try
        {
            myConfig.fetch("filename_out", this->filename_output_);
        }
        catch(std::logic_error& e)
        {
        }
    #endif

    // Other configurable data member
    #if CPLUSPLUS_ANALYSIS
        try
        {
            myConfig.fetch("filename_out_cpp", this->filename_output_cpp_);
        }
        catch(std::logic_error& e)
        {
        }
    #endif
    
    // ROOT
    std::cout << "In INIT: create TFile " << std::endl;
    
    #if PYTHON_ANALYSIS
        // Variables for Python analysis
        hfile_ = new TFile(filename_output_.c_str(), "RECREATE", "Output file of Simulation data");
        hfile_->cd();
    
        tree_ = new TTree("TSD", "TSD"); // timestamp data
        tree_->SetDirectory(hfile_);
    #endif
    
    // header data
    //tree_->Branch("header.runnumber",&header_.runnumber_);
    //tree_->Branch("header.eventnumber",&header_.eventnumber_);
    //tree_->Branch("header.date",&header_.date_);
    //tree_->Branch("header.runtype",&header_.runtype_);
    //tree_->Branch("header.simulated",&header_.simulated_);
    
    // generator data
    //tree_->Branch("truth.vertex_x", &gen_.vertex_x_);
    //tree_->Branch("truth.vertex_y", &gen_.vertex_y_);
    //tree_->Branch("truth.vertex_z", &gen_.vertex_z_);
    
    // timestamp data    
    #if PYTHON_ANALYSIS
        tree_->Branch("timestamp.anodic_t0", &timestamp_.anodic_t0);
        tree_->Branch("timestamp.anodic_t1", &timestamp_.anodic_t1);
        tree_->Branch("timestamp.anodic_t2", &timestamp_.anodic_t2);
        tree_->Branch("timestamp.anodic_t3", &timestamp_.anodic_t3);
        tree_->Branch("timestamp.anodic_t4", &timestamp_.anodic_t4);
        tree_->Branch("timestamp.cathodic_t5", &timestamp_.cathodic_t5);
        tree_->Branch("timestamp.cathodic_t6", &timestamp_.cathodic_t6);
        
        tree_->Branch("timestamp.cell_x", &timestamp_.cell_x);
        tree_->Branch("timestamp.cell_y", &timestamp_.cell_y);
        
        tree_->Branch("truth.n_gamma", &gen_.n_gamma_);
        tree_->Branch("truth.n_positron", &gen_.n_positron_);
        tree_->Branch("truth.n_electron", &gen_.n_electron_);
        tree_->Branch("truth.n_alpha", &gen_.n_alpha_);
        tree_->Branch("truth.caffe_category", &gen_.caffe_category_);
    #endif
    
    // Variables for C++ analysis
    #if CPLUSPLUS_ANALYSIS
        // Open output file and create tree for C++ analysis code
        hfile_cpp_ = new TFile(filename_output_cpp_.c_str(), "RECREATE", "Output file of Simulation data"); // TODO: change name
        hfile_cpp_->cd();

        tree_cpp_ = new TTree("histo", "histo"); // TODO: check name
        tree_cpp_->SetDirectory(hfile_cpp_);

        // Additional variables required for the C++ analysis code
        tree_cpp_->Branch("time", &store_.time);
        tree_cpp_->Branch("delay", &store_.delay);
        tree_cpp_->Branch("delay_since_good_trigger", &store_.delay_since_good_trigger);
        tree_cpp_->Branch("duration", &store_.duration);
        tree_cpp_->Branch("plasma_propagation_time", &store_.plasma_propagation_time);
        tree_cpp_->Branch("good_trigger", &store_.good_trigger);
        tree_cpp_->Branch("prev_good_trigger", &store_.prev_good_trigger);
        tree_cpp_->Branch("with_cathode", &store_.with_cathode);
        tree_cpp_->Branch("anode_peak", &store_.anode_peak);
        tree_cpp_->Branch("anode_time", &store_.anode_time);
        tree_cpp_->Branch("cathode_peak", &store_.cathode_peak);
        tree_cpp_->Branch("cathode_time", &store_.cathode_time);
        tree_cpp_->Branch("position", &store_.position);
        tree_cpp_->Branch("half_position", &store_.half_position);
        tree_cpp_->Branch("stop1", &store_.stop1);
        tree_cpp_->Branch("stop1_peak", &store_.stop1_peak);
        tree_cpp_->Branch("stop1_type", &store_.stop1_type);
        tree_cpp_->Branch("stop2", &store_.stop2);
        tree_cpp_->Branch("stop2_peak", &store_.stop2_peak);
        tree_cpp_->Branch("stop2_type", &store_.stop2_type);
        tree_cpp_->Branch("stopA", &store_.stopA);
        tree_cpp_->Branch("deriv_rms", &store_.deriv_rms);
        tree_cpp_->Branch("feast_t0", &store_.feast_t0);
        tree_cpp_->Branch("feast_t1", &store_.feast_t1);
        tree_cpp_->Branch("feast_t2", &store_.feast_t2);
        tree_cpp_->Branch("feast_t3", &store_.feast_t3);
        tree_cpp_->Branch("feast_t4", &store_.feast_t4);
        tree_cpp_->Branch("anode_histo", &store_.anode_histo);
        tree_cpp_->Branch("deriv_histo", &store_.deriv_histo);
        tree_cpp_->Branch("cathode_histo", &store_.cathode_histo);
    #endif
   
    // C++ analysis: Set variables to sensible default values
    // (most are unused)
    #if CPLUSPLUS_ANALYSIS
        store_.time = 0.0;
        store_.time = 0.0;
        store_.delay = 0.0;
        store_.delay_since_good_trigger = 0.0;
        store_.duration = 0;
        store_.plasma_propagation_time = 0;
        store_.good_trigger = 0;
        store_.prev_good_trigger = 0;
        store_.with_cathode = 0;
        store_.anode_peak = 0.0;
        store_.anode_time = 0.0;
        store_.cathode_peak = 0.0;
        store_.cathode_time = 0.0;
        store_.position = 0.0;
        store_.half_position = 0.0;
        store_.stop1 = 0.0;
        store_.stop1_peak = 0.0;
        store_.stop1_type = 0.0;
        store_.stop2 = 0.0;
        store_.stop2_peak = 0.0;
        store_.stop2_type = 0.0;
        store_.stopA = 0;
        store_.deriv_rms = 0.0;
        store_.feast_t0 = 0.0;
        store_.feast_t1 = 0.0;
        store_.feast_t2 = 0.0;
        store_.feast_t3 = 0.0;
        store_.feast_t4 = 0.0;
        store_.anode_histo = (TH1F*)0;
        store_.deriv_histo = (TH1F*)0;
        store_.cathode_histo = (TH1F*)0;
    #endif

    this->base_module::_set_initialized(true);
    return;
}

void TestModule::reset()
{
    DT_THROW_IF(!is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    if (!has_external_random()) {
        // Reset the random number generator:
        _random_.reset();
    }
    _external_random_ = 0;
    _set_defaults();
    
    #if PYTHON_ANALYSIS
        // write the output, finished streaming
        hfile_->cd();
        tree_->Write();
        hfile_->Close();
        std::cout << "In reset: finished conversion, file closed " << std::endl;
    
        // clean up
        //delete tree_;
        delete hfile_;
    #endif

    // Clean up (C++ analysis code)
    #if CPLUSPLUS_ANALYSIS
        hfile_cpp_->cd();
        tree_cpp_->Write();
        hfile_cpp_->Close();
        std::cout << "In reset: finished conversion, cpp file closed" << std::endl;
        delete hfile_cpp_;
    #endif

    // Reset default filenames
    #if PYTHON_ANALYSIS
        filename_output_ = "default.root";
    #endif

    #if CPLUSPLUS_ANALYSIS
        filename_output_cpp_ = "default_cpp.root";
    #endif

    this->base_module::_set_initialized(false);
    
    return;
}

void TestModule::_set_defaults() {
    _external_random_ = 0;
    _CD_label_.clear();
    return;
}

TestModule::TestModule() : dpp::base_module()
{
    std::cout << "TestModule" << std::endl;
    
    _external_random_ = 0;
    _set_defaults();
    
    #if PYTHON_ANALYSIS
        filename_output_ = "default.root";
    #endif

    #if CPLUSPLUS_ANALYSIS
        filename_output_cpp_ = "default_cpp.root";
    #endif

    return;
}

TestModule::~TestModule()
{
    std::cout << "~TestModule" << std::endl;
    
    if (is_initialized()) TestModule::reset();
    // MUST reset module at destruction
    //if(is_initialized()) // TODO
    //this->reset();
    return;
}

dpp::base_module::process_status
TestModule::process(datatools::things& workItem)
{

    //std::cout << "process" << std::endl;
    
    
    
    #if CPLUSPLUS_ANALYSIS
        // TODO: some values are junk/ignored
    
        // Tempoary local storage 
        // For C++ analysis program    
        Float_t time;
        Float_t delay;
        Float_t delay_since_good_trigger;
        Int_t duration;
        Float_t plasma_propagation_time;
        Bool_t good_trigger;
        Bool_t prev_good_trigger;
        Bool_t with_cathode;
        Float_t anode_peak;
        Float_t anode_time;
        Float_t cathode_peak;
        Float_t cathode_time;
        Float_t position;
        Float_t half_position;
        Float_t stop1;
        Float_t stop1_peak;
        Float_t stop1_type;
        Float_t stop2;
        Float_t stop2_peak;
        Float_t stop2_type;
        Int_t stopA;
        Float_t deriv_rms;
        Float_t feast_t0;
        Float_t feast_t1;
        Float_t feast_t2;
        Float_t feast_t3;
        Float_t feast_t4;
        TH1F *anode_histo;
        TH1F *deriv_histo;
        TH1F *cathode_hist;

        // Set class local storage to default values
        // TODO: not needed anymore?
        store_.time = 0;
        store_.delay = 0;
        store_.delay_since_good_trigger = 0;
        store_.duration = 0;
        store_.plasma_propagation_time = 0;
        store_.good_trigger = 0;
        store_.prev_good_trigger = 0;
        store_.with_cathode = 0;
        store_.anode_peak = 0;
        store_.anode_time = 0;
        store_.cathode_peak = 0;
        store_.cathode_time = 0;
        store_.position = 0;
        store_.half_position = 0;
        store_.stop1 = 0;
        store_.stop1_peak = 0;
        store_.stop1_type = 0;
        store_.stop2 = 0;
        store_.stop2_peak = 0;
        store_.stop2_type = 0;
        store_.stopA = 0;
        store_.deriv_rms = 0;
        store_.feast_t0 = 0;
        store_.feast_t1 = 0;
        store_.feast_t2 = 0;
        store_.feast_t3 = 0;
        store_.feast_t4 = 0;
        store_.anode_histo = (TH1F*)0;
        store_.deriv_histo = (TH1F*)0;
        store_.cathode_histo = (TH1F*)0;
        
    #endif

    // Tempoary local storage
    // timestamp event data
    // NOTE: required for both C++ and Python analysis modes
    std::vector<double> anodic_t0;
    std::vector<double> anodic_t1;
    std::vector<double> anodic_t2;
    std::vector<double> anodic_t3;
    std::vector<double> anodic_t4;
    std::vector<double> cathodic_t5;
    std::vector<double> cathodic_t6;
    std::vector<double> cell_x;
    std::vector<double> cell_y;
    std::vector<double> anodic_plasma_propagation_time;

    DT_THROW_IF(!is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");
    
    // Process CD bank to obtains timestamp data
    if(workItem.has("CD"))
    {
        //std::cout << "CD" << std::endl;
    
        const snemo::datamodel::calibrated_data & CD = workItem.get<snemo::datamodel::calibrated_data>("CD");
        //hits_.nofhits_ = CD.calibrated_tracker_hits().size();
        
        if(CD.has_calibrated_tracker_hits())
        {
            
            //timestamp.count = CD.get_number_of_calibrated_tracker_hits();
            //timestamp_.count = 0;
            
            const snemo::datamodel::calibrated_data::tracker_hit_collection_type & THCT = CD.calibrated_tracker_hits();
            
            timestamp_.count = THCT.size();
            std::cout << "timestamp_.count = " << timestamp_.count << std::endl;
            //int timestamp_count = THCT.size();
            //std::cout << "timestamp_count = " << timestamp_count << std::endl;
            
            //for(int ix = 0; ix < timestamp_count; ++ ix)
            for(int ix = 0; ix < THCT.size() /*hits_.nofhits_*/; ++ ix)
            {
                const snemo::datamodel::calibrated_data::tracker_hit_handle_type & THHT = THCT.at(ix);
                if(THHT.has_data())
                {
                    const snemo::datamodel::calibrated_tracker_hit & CTH = THHT.get();
                    if(CTH.has_anode_time())
                    {
                        double plasma_propagation_time_ = CTH.get_plasma_propagation_time();
                        
                        //for(size_t ix{0}; ix < 
                        const double anode_time = CTH.get_anode_time(); // microsecond?
                        
                        double anode_t1, anode_t2, anode_t3, anode_t4;
                        
                        /// NOTE: C++ analysis code (with TestTank data) used unit of
                        // microsecond throughout code
                        // Falaise uses unknown unit (presumed to be SI) and uses
                        // CLHEP::microsecond to convert to correct unit

                        // parameters from tracker-signals fit
                        // _cor -> taken from the correlation plot in the C++ analysis code
                        // TODO: check codes are the same (there are 2 lots in the C++ code?)
                        const double x0_cor = CLHEP::microsecond * -0.183547; // probably micro-second
                        const double y0_cor = CLHEP::microsecond * 0.179927;
                        const double sigma_x_cor = CLHEP::microsecond * 0.0458669; // microsecond ?
                        const double sigma_y_cor = CLHEP::microsecond * 0.0190218;
                        const double theta_cor = -0.363398; // rad?
                        
                        const double u1 = _get_random().gaussian(0.0, sigma_x_cor);
                        const double v1 = _get_random().gaussian(0.0, sigma_y_cor);
                        const double x1 = x0_cor + u1 * std::cos(-theta_cor) - v1 * std::sin(-theta_cor);
                        const double y1 = y0_cor + u1 * std::sin(-theta_cor) + v1 * std::cos(-theta_cor);
                        
                        const double u2 = _get_random().gaussian(0.0, sigma_x_cor);
                        const double v2 = _get_random().gaussian(0.0, sigma_y_cor);
                        const double x2 = x0_cor + u2 * std::cos(-theta_cor) - v2 * std::sin(-theta_cor);
                        const double y2 = y0_cor + u2 * std::sin(-theta_cor) + v2 * std::cos(-theta_cor);
                        
                        //anode_t1 = CTH.get_top_cathode_time() + x1;
                        //anode_t2 = CTH.get_top_cathode_time() + y1;
                        
                        //anode_t3 = CTH.get_bottom_cathode_time() + x2;
                        //anode_t4 = CTH.get_bottom_cathode_time() + y2;
                        

                        // timestamps T1, T3 go together
                        anode_t1 = CTH.get_top_cathode_time() + x1;
                        anode_t3 = CTH.get_top_cathode_time() + y1;
                        
                        // T2, T4 go together
                        anode_t2 = CTH.get_bottom_cathode_time() + x2;
                        anode_t4 = CTH.get_bottom_cathode_time() + y2;

                        //std::cout << "anode_t1=" << anode_t1 << std::endl;
                        
                        // NOTE TO SELF: To make data from this simulation match up with the test tank data,
                        // we replace the anode T0 offset with the mean value as obtained from the TestTank
                        // data
                        // DO NOT PUSH BACK THE anode_time !!!
                        // Push back the mean anode offset
                        // TODO: GET THE VARIANCE BACK HERE
                        const double testtank_mean_anode_time = CLHEP::microsecond * 4.808655335; // std::cout.precision(10)
                        const double falaise_mean_anode_time = CLHEP::microsecond * 1.438549573;
                        const double testtank_anode_time = testtank_mean_anode_time + anode_time - falaise_mean_anode_time;
                        anodic_t0.push_back(testtank_anode_time);
                        //anodic_t0.push_back(testtank_mean_anode_time);
                        /*timestamp_.*///anodic_t0.push_back(anode_time); // TODO: is this an absolute time, or relative to something?
                        /*timestamp_.*/anodic_t1.push_back(anode_t1); // If cathode top/bottom times are absolute this is correct
                        /*timestamp_.*/anodic_t2.push_back(anode_t2);
                        /*timestamp_.*/anodic_t3.push_back(anode_t3);
                        /*timestamp_.*/anodic_t4.push_back(anode_t4);
                        /*timestamp_.*/cathodic_t5.push_back(CTH.get_top_cathode_time()); // TODO: check if correct - might need to subtract anode time
                        /*timestamp_.*/cathodic_t6.push_back(CTH.get_bottom_cathode_time());
                        // NOTE: absolute time meaning relative to some EPOCH
                        
                        anodic_plasma_propagation_time.push_back(plasma_propagation_time_);
                        
                        //std::cout << anode_time << std::endl;
                        
                        // TIMESTAMP DATA IS OBTAINED HERE
                    }
                    
                    
                    int32_t module{CTH.get_module()};
                    int32_t side{CTH.get_side()};
                    int32_t layer{CTH.get_layer()};
                    int32_t row{CTH.get_row()};
                    double z_pos{CTH.get_z()}; // z position of hit use this to start with
                    // TODO: set other channels to zero
                    
                    if(side == 0)
                    {
                    
                    }
                    else if(side == 1)
                    {
                    
                    }
                    else
                    {
                        std::cout << "side=" << side << std::endl;
                    }
                    
                    if(layer < 0 || layer > 8)
                    {
                        std::cout << "layer=" << layer << std::endl;
                        std::cin.get();
                    }
                    if(row < 0 || row > 112)
                    {
                        std::cout << "row=" << row << std::endl;
                        std::cin.get();
                    }
                    
                    int32_t caffe_x{row};
                    int32_t caffe_y{side * 9 + layer};
                    int32_t ix_max{9 * 113};
                    int32_t caffe_ix{9 * caffe_y + caffe_x};
                    
                    // TODO: PROBLEM: this only sets the z_pos for a single cell, and yet we itterate over all
                    // cells - so need to create the data FIRST with BLANK (-1.0f) data and then SET the cells
                    // which have data/hits HERE without the loop (?)
                    #if CAFFE_ENABLE
                        google::protobuf::RepeatedField<float>* datumFloatData = datum.mutable_float_data();
                        for(int32_t ix{0}; ix < caffe_ix - 1; ++ ix)
                        {
                            datumFloatData->Add(-1.0f); // TODO! THIS LOOKS LIKE HITS ON z=0.0 !!! FIXME
                        }
                        // add at ix_ix
                        datumFloatData->Add(0.5 * (z_pos + 1.0)); // TODO: check
                        for(int32_t ix{caffe_ix + 1}; ix < ix_max; ++ ix)
                        {
                            datumFloatData->Add(-1.0f);
                        }
                        
                        // This is done in section below
                        //datum.set_label();
                        //datum.set_data(, );
                    #endif
                    
                    if(CTH.has_xy())
                    {
                        cell_x.push_back(CTH.get_x());
                        cell_y.push_back(CTH.get_y());
                        
                        // CELL X Y DATA IS OBTAINED HERE
                        
                        // DO NOT YET KNOW WHERE TRUTH INFORMATION IS OBTAINED FROM TODO
                    }
                    else
                    {
                        double invalid;
                        datatools::invalidate(invalid);
                        cell_x.push_back(invalid);
                        cell_y.push_back(invalid);
                    }
                    
                }
            
            }
        }
        else
        {
            timestamp_.count = 0;
        }
    }
    else
    {
        //hits_.nofhits_ = 0;
        std::cout << "DOES NOT HAVE CD" << std::endl;
    }
    
    timestamp_.anodic_t0 = &anodic_t0;
    timestamp_.anodic_t1 = &anodic_t1;
    timestamp_.anodic_t2 = &anodic_t2;
    timestamp_.anodic_t3 = &anodic_t3;
    timestamp_.anodic_t4 = &anodic_t4;
    timestamp_.cathodic_t5 = &cathodic_t5;
    timestamp_.cathodic_t6 = &cathodic_t6;
    timestamp_.cell_x = &cell_x;
    timestamp_.cell_y = &cell_y;
    timestamp_.plasma_propagation_time = &anodic_plasma_propagation_time;

    if(workItem.has("SD"))
    {
        const mctools::simulated_data& SD = workItem.get<mctools::simulated_data>("SD");
        //gen_.vertex_x_ = SD.get_vertex().x();
        //gen_.vertex_y_ = SD.get_vertex().y();
        //gen_.vertex_z_ = SD.get_vertex().z();
        
        // UID assembly
        #if UID_ENABLE
            uid_assembler<unsigned long long> uid;
            uid.add(genbb::primary_particle::particle_type::GAMMA);
            uid.add(genbb::primary_particle::particle_type::POSITRON);
            uid.add(genbb::primary_particle::particle_type::ELECTRON);
            uid.add(genbb::primary_particle::particle_type::ALPHA);
            uid.finalize();
            
            //TODO: all of these should have a CHECK before a GET
            const mctools::simulated_data::primary_event_type primary_evt{SD.get_primary_event()};
            for(unsigned int pix{0}; pix < primary_evt.get_number_of_particles(); ++ pix)
            {
                //if(primary_evt.has_particle())
                //{
                const genbb::primary_particle & pp{primary_evt.get_particle(pix)};
                
                if(pp.has_type())
                {
                    int pp_t{pp.get_type()};
                    
                    uid.increment_field(pp_t);
                    /*if(pp_t == particle_type::GAMMA)
                    {
                        uid.increment_field(particle_type::GAMMA);
                    }
                    else if(pp_t == particle_type::POSITRON)
                    {
                    
                    }
                    else if(pp_t == particle_type::ELECTRON)
                    {
                    
                    }
                    else if(pp_t == particle_type::ALPHA)
                    {
                    
                    }
                    */
                }
                //}
            }
            
            // set the count for each particle type
            gen_.n_gamma_ = uid.get(genbb::primary_particle::particle_type::GAMMA);
            gen_.n_positron_ = uid.get(genbb::primary_particle::particle_type::POSITRON);
            gen_.n_electron_ = uid.get(genbb::primary_particle::particle_type::ELECTRON);
            gen_.n_alpha_ = uid.get(genbb::primary_particle::particle_type::ALPHA);
            
            //std::cout << "**** UID ****" << std::endl;
            //uid.print_fields(std::cout); std::cout << std::endl;
            
            // set the Caffe uid
            gen_.caffe_category_ = uid.generate();
            //std::cout << gen_.caffe_category_ << std::endl;
        #else
            gen_.caffe_category_ = 0;
        #endif
        
        #if CAFFE_ENABLE
            datum.set_label(gen_.caffe_category_);
        #endif
        
    }
    else
    {
        std::cout << "Does not have SD!" << std::endl;
    }
    
    
    #if PYTHON_ANALYSIS
        tree_->Fill();
    #endif
    
    // At this point data is processed and stored in tree for Python analysis code
    // Re-format data for C++ analysis code
    
    // This loop converts to correct units

    for(size_t ix{0}; ix < anodic_t0.size(); ++ ix)
    {
        // Check ix is valid
        if((ix >= anodic_t1.size()) ||
           (ix >= anodic_t2.size()) ||
           (ix >= anodic_t3.size()) || 
           (ix >= anodic_t4.size()) ||
           (ix >= cathodic_t5.size()) || 
           (ix >= cathodic_t6.size())  )
        {
            std::cout << "WARNING: ix is out of range for other vectors" << std::endl;
            break;
        }

        // TODO: need to check same variables match up here
        
        // Get values and convert to correct unit
        anode_time = anodic_t0.at(ix) / CLHEP::microsecond;

        feast_t0 = anodic_t0.at(ix) / CLHEP::microsecond;
        feast_t1 = anodic_t1.at(ix) / CLHEP::microsecond;
        feast_t2 = anodic_t2.at(ix) / CLHEP::microsecond;
        feast_t3 = anodic_t3.at(ix) / CLHEP::microsecond;
        feast_t4 = anodic_t4.at(ix) / CLHEP::microsecond;

        cathode_time = cathodic_t5.at(ix) / CLHEP::microsecond;
        // only single side cathode time is stored

        plasma_propagation_time = anodic_plasma_propagation_time.at(ix) / CLHEP::microsecond;

        // Save values
        store_.anode_time = anode_time;
        
        store_.feast_t0 = feast_t0;
        store_.feast_t1 = feast_t1;
        store_.feast_t2 = feast_t2;
        store_.feast_t3 = feast_t3;
        store_.feast_t4 = feast_t4;

        store_.cathode_time = cathode_time;

        store_.plasma_propagation_time = plasma_propagation_time;


        // Print values
        #if COUT_TIMESTAMP
            std::cout << "cathode               : " << store_.cathode_time << "\n"\
                      << "t0                    : " << store_.feast_t0 << "\n"\
                      << "t1                    : " << store_.feast_t1 << "\n"\
                      << "t3                    : " << store_.feast_t3 << "\n"\
                      << "t2                    : " << store_.feast_t2 << "\n"\
                      << "t4                    : " << store_.feast_t4 << "\n"\
                      << "cathode + t0          : " << store_.cathode_time + store_.feast_t0 << "\n"\
                      << "t1 - (cathode + t0)   : " << store_.feast_t1 - (store_.cathode_time + store_.feast_t0) << "\n"\
                      << "t3 - (cathode + t0)   : " << store_.feast_t3 - (store_.cathode_time + store_.feast_t0) << "\n"\
                      << "t2 - (cathode + t0)   : " << store_.feast_t2 - (store_.cathode_time + store_.feast_t0) << "\n"\
                      << "t4 - (cathode + t0)   : " << store_.feast_t4 - (store_.cathode_time + store_.feast_t0) << "\n";
        
            std::cout.flush();

            #if COUT_TIMESTAMP_WAIT
                std::cin.get();
            #endif
        #endif


        // Fill tree
        tree_cpp_->Fill();
    }

    return PROCESS_OK;
}
