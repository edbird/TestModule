#include "TestModule.h"

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
//#include <mctools/simulated_data.h>
//#include <mctools/utils.h>

// This project :
#include <falaise/snemo/datamodels/data_model.h>
//#include <falaise/snemo/processing/services.h>

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
    
    this->base_module::_set_initialized(true);
    return;
}

void TestModule::reset()
{
    DT_THROW_IF(!is_initialized(), std::logic_error,
                "Module '" << get_name() << "' is not initialized !");

    this->base_module::_set_initialized(false);

    if (!has_external_random()) {
        // Reset the random number generator:
        _random_.reset();
    }
    _external_random_ = 0;
    _set_defaults();
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
    return;
}

TestModule::~TestModule()
{
    std::cout << "~TestModule" << std::endl;
    
    if (is_initialized()) TestModule::reset();
    return;
}

dpp::base_module::process_status
TestModule::process(datatools::things& workItem)
{

    //std::cout << "process" << std::endl;
    
    // Tempoary local storage
    // timestamp event data
    std::vector<double> anodic_t0;
    std::vector<double> anodic_t1;
    std::vector<double> anodic_t2;
    std::vector<double> anodic_t3;
    std::vector<double> anodic_t4;
    std::vector<double> cathodic_t5;
    std::vector<double> cathodic_t6;
    std::vector<double> cell_x;
    std::vector<double> cell_y;

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
                        //for(size_t ix{0}; ix < 
                        double anode_time = CTH.get_anode_time();
                        
                        double anode_t1, anode_t2, anode_t3, anode_t4;
                        
                        // parameters from tracker-signals fit
                        double x0_cor = -0.183547;
                        double y0_cor = 0.179927;
                        double sigma_x_cor = 0.0458669;
                        double sigma_y_cor = 0.0190218;
                        double theta_cor = -0.363398;
                        
                        double u1 = _get_random().gaussian(0.0, sigma_x_cor);
                        double v1 = _get_random().gaussian(0.0, sigma_y_cor);
                        double x1 = x0_cor + u1 * std::cos(-theta_cor) - v1 * std::sin(-theta_cor);
                        double y1 = y0_cor + u1 * std::sin(-theta_cor) + v1 * std::cos(-theta_cor);
                        
                        double u2 = _get_random().gaussian(0.0, sigma_x_cor);
                        double v2 = _get_random().gaussian(0.0, sigma_y_cor);
                        double x2 = x0_cor + u2 * std::cos(-theta_cor) - v2 * std::sin(-theta_cor);
                        double y2 = y0_cor + u2 * std::sin(-theta_cor) + v2 * std::cos(-theta_cor);
                        
                        anode_t1 = CTH.get_top_cathode_time() + x1;
                        anode_t2 = CTH.get_top_cathode_time() + y1;
                        
                        anode_t3 = CTH.get_bottom_cathode_time() + x2;
                        anode_t4 = CTH.get_bottom_cathode_time() + y2;
                        
                        //std::cout << "anode_t1=" << anode_t1 << std::endl;
                        
                        /*timestamp_.*/anodic_t0.push_back(anode_time);
                        /*timestamp_.*/anodic_t1.push_back(anode_t1);
                        /*timestamp_.*/anodic_t2.push_back(anode_t2);
                        /*timestamp_.*/anodic_t3.push_back(anode_t3);
                        /*timestamp_.*/anodic_t4.push_back(anode_t4);
                        /*timestamp_.*/cathodic_t5.push_back(CTH.get_top_cathode_time());
                        /*timestamp_.*/cathodic_t6.push_back(CTH.get_bottom_cathode_time());
                        
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
                    
                    if(layer < 0 || layer > 8) std::cout << "layer=" << layer << std::endl;
                    if(row < 0 || row > 112) std::cout << "row=" << row << std::endl;
                    std::cin.get();
                    
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
                    */
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
    
    
    
    if(workItem.has("SD"))
    {
        const mctools::simulated_data& SD = workItem.get<mctools::simulated_data>("SD");
        gen_.vertex_x_ = SD.get_vertex().x();
        gen_.vertex_y_ = SD.get_vertex().y();
        gen_.vertex_z_ = SD.get_vertex().z();
        
        // UID assembly
        #if UID_ENABLE
        uid_assembler<Long64_t> uid;
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
        
        // set the Caffe uid
        gen_.caffe_category_ = uid.generate();
        #else
        gen_.caffe_category_ = 0;
        #endif
        
        #if CAFFE_ENABLE
        datum.set_label(gen_.caffe_category_);
        #endif
        
    }
    
    
    
    
    return PROCESS_OK;
}