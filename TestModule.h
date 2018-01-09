#ifndef TESTMODULE_HH
#define TESTMODULE_HH

// Standard library:
#include <map>
#include <string>
#include <vector>

// Third party:
// - Bayeux/mygsl:
#include <mygsl/rng.h>
// - Bayeux/dpp:
#include <dpp/base_module.h>

// - Bayeux
//#include "bayeux/dpp/base_module.h"
//#include "bayeux/mygsl/rng.h"
//#include "bayeux/datatools/service_manager.h"
//#include "bayeux/geomtools/manager.h"
//#include "bayeux/geomtools/geometry_service.h"
//#include "bayeux/geomtools/line_3d.h"
//#include "bayeux/geomtools/helix_3d.h"

// - Falaise
//#include "falaise/snemo/datamodels/particle_track_data.h"
//#include "falaise/snemo/datamodels/particle_track.h"
//#include "falaise/snemo/datamodels/calibrated_calorimeter_hit.h"
#include "falaise/snemo/datamodels/calibrated_data.h"
//#include "falaise/snemo/datamodels/tracker_clustering_data.h"
//#include "falaise/snemo/datamodels/base_trajectory_pattern.h"
//#include "falaise/snemo/datamodels/event_header.h"



// TODO: clean
typedef struct GeneratorEventStorage
{
    double vertex_x_;
    double vertex_y_;
    double vertex_z_;
    int nofparticles_;
    std::vector<int>* pdgs_;
    std::vector<double>* px_;
    std::vector<double>* py_;
    std::vector<double>* pz_;
    
    // new
    int n_gamma_;
    int n_positron_;
    int n_electron_;
    int n_alpha_;
    
    // TODO: some events might have other particles in them - at the moment
    // i do not handle those cases
    
    // new - category flag for caffe
    unsigned long long caffe_category_;
} generatoreventstorage;



typedef struct TimestampStorage
{
    // 7 geiger timing
    int count;
    std::vector<double> * anodic_t0;
    std::vector<double> * anodic_t1;
    std::vector<double> * anodic_t2;
    std::vector<double> * anodic_t3;
    std::vector<double> * anodic_t4;
    std::vector<double> * cathodic_t5;
    std::vector<double> * cathodic_t6;
    
    // xy location of hit
    std::vector<double> * cell_x;
    std::vector<double> * cell_y;
    
} timestampstorage;


class TestModule : public dpp::base_module
{
public:
    /// Set the external PRNG
    void set_external_random(mygsl::rng& rng_);

    /// Reset the external PRNG
    void reset_external_random();

    /// Check if the module use an external PRNG
    bool has_external_random() const;

    
    TestModule();
    
    virtual ~TestModule();
    
    virtual void initialize(const datatools::properties& myConfig,
                            datatools::service_manager& flServices,
                            dpp::module_handle_dict_type& moduleDict);
    
    virtual dpp::base_module::process_status process(datatools::things& workItem);
    
    virtual void reset();
    
    
protected:
    /// Set default attributes values
    void _set_defaults();
    
    /// Getting random number generator
    mygsl::rng& _get_random();
    
    
private:
    std::string _module_category_;             //!< The geometry category of the SuperNEMO module
    mygsl::rng _random_;                       //!< internal PRN generator
    mygsl::rng* _external_random_;             //!< external PRN generator
    std::string _CD_label_;                    //!< The label of the calibrated data bank
    
    // Local storage
    TimestampStorage timestamp_;
    
    GeneratorEventStorage gen_;
    
    DPP_MODULE_REGISTRATION_INTERFACE(TestModule)
    
};

#endif