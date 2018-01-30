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


// Third Party
#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>

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


// CAFFE CPU SWITCH
#define CAFFE_ENABLE 1
#if CAFFE_ENABLE
    #define CPU_ONLY
#endif

// caffe
///home/ecb/caffe/include
#if CAFFE_ENABLE
    #include "boost/scoped_ptr.hpp"
    #include "caffe/proto/caffe.pb.h"
    #include "caffe/util/db.hpp"
    #include "caffe/util/format.hpp"
    #include "caffe/util/io.hpp"
    
    // Note to self: line copied from
    // https://github.com/BVLC/caffe/blob/master/tools/convert_imageset.cpp
    //#include "caffe/util/db.hpp"
#endif


/******************************************************************************/
/* ANALYSIS SCRIPT SWITCH (OUTPUT MODE / TYPE SWITCH)                         */
/******************************************************************************/

// This switch changes between the origional C++ analysis code
// and a new Python analysis script (which is not debugged yet)

#define CPLUSPLUS_ANALYSIS 1
#define PYTHON_ANALYSIS 0

/******************************************************************************/
/*                                                                            */
/******************************************************************************/


// TODO: clean
// TODO: should remove this if CPLUSPLUS_ANALYSIS
//#if PYTHON_ANALYSIS
typedef struct GeneratorEventStorage
{
    //double vertex_x_;
    //double vertex_y_;
    //double vertex_z_;
    //int nofparticles_;
    //std::vector<int>* pdgs_;
    //std::vector<double>* px_;
    //std::vector<double>* py_;
    //std::vector<double>* pz_;
    
    // new
    int n_gamma_;
    int n_positron_;
    int n_electron_;
    int n_alpha_;
    
    // TODO: some events might have other particles in them - at the moment
    // i do not handle those cases
    
    // new - category flag for caffe
    unsigned long long caffe_category_;
}; //generatoreventstorage;
//#endif

// Notes:
// TestTankStorage is different from the other storage formats.
// Other storage formats use vectors to keep all data for each
// event separated.
// The C++ analysis code loops over a single loop and the data
// is formatted such that all events are "rolled together".
// In other words: Each event is a set of timestamps rather than
// a vector of timestamps.
// The regular output file is always created. The C++ analysis
// file is only created if the CPLUSPLUS_ANALSIS preprocessor
// flag is set.

#if CPLUSPLUS_ANALYSIS
    typedef struct TestTankStorage
    {
        Float_t time = 0;
        Float_t delay = 0;
        Float_t delay_since_good_trigger = 0;
        Int_t duration = 0;
        Float_t plasma_propagation_time = 0;
        Bool_t good_trigger = 0;
        Bool_t prev_good_trigger = 0;
        Bool_t with_cathode = 0;
        Float_t anode_peak = 0;
        Float_t anode_time = 0;
        Float_t cathode_peak = 0;
        Float_t cathode_time = 0;
        Float_t position = 0;
        Float_t half_position = 0;
        Float_t stop1 = 0;
        Float_t stop1_peak = 0;
        Float_t stop1_type = 0;
        Float_t stop2 = 0;
        Float_t stop2_peak = 0;
        Float_t stop2_type = 0;
        Int_t stopA = 0;
        Float_t deriv_rms = 0;
        Float_t feast_t0 = 0;
        Float_t feast_t1 = 0;
        Float_t feast_t2 = 0;
        Float_t feast_t3 = 0;
        Float_t feast_t4 = 0;
        TH1F *anode_histo = (TH1F*)0;
        TH1F *deriv_histo = (TH1F*)0;
        TH1F *cathode_histo = (TH1F*)0;
    };
#endif

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

    // plasma propagation time
    // vector format, as obtained from calibrated_tracker_hit class
    std::vector<double> * plasma_propagation_time;
    std::vector<double> * position;
    std::vector<double> * half_position;

}; //timestampstorage;


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
    
    #if CPLUSPLUS_ANALYSIS
        // Local storage
        // Matches data obtained from testtank
        // Additional variables for C++ analysis now in structure
        TestTankStorage store_;
    #endif
    
    // Local storage
    TimestampStorage timestamp_;
    GeneratorEventStorage gen_;
    
    // configurable data member
    #if PYTHON_ANALYSIS
        std::string filename_output_;
    #endif
    #if CPLUSPLUS_ANALYSIS
        std::string filename_output_cpp_;
    #endif

    // ROOT variables
    #if PYTHON_ANALYSIS
        File* hfile_;
        TTree* tree_;
    #endif
    #if CPLUSPLUS_ANALYSIS
        TFile* hfile_cpp_;
        TTree* tree_cpp_;
    #endif

    #if CAFFE_ENABLE
        // Note to self: not sure about these lines, copied from
        // https://github.com/BVLC/caffe/blob/master/tools/convert_imageset.cpp
        // Create new DB
        //DataParameter_DB_LMDB - don't know how to use this yet
        //boost::scoped_ptr<caffe::db::DB> db(caffe::db::GetDB(FLAGS_backend));
        caffe::db::DB *db_p{caffe::db::GetDB("lmdb")};
        //boost::scoped_ptr<caffe::db::DB> the_db(db_p);
        
        //boost::scoped_ptr<caffe::db::Transaction> the_txn;
        caffe::db::Transaction *txn_p;
        
        // TODO this goes in the init function
        //db->Open("temp_caffe_db.lmdb", db::NEW);
        //scoped_ptr<db::Transaction> txn(db->NewTransaction());
    #endif
    
    DPP_MODULE_REGISTRATION_INTERFACE(TestModule)
    
};

#endif
