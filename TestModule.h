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
    
    DPP_MODULE_REGISTRATION_INTERFACE(TestModule)
    
};

#endif