// Copyright 2026, Diogo Costa (diogo.pinhodacosta@canada.ca)
// This file is part of OpenWQ model.

// This program, openWQ, is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

// =============================================================================
// OpenWQ Hydrolink Header for SUMMA
// =============================================================================
// This header defines the CLASSWQ_openwq class that provides the interface
// between SUMMA (Structure for Unifying Multiple Modeling Alternatives) and
// the OpenWQ water quality model.
//
// SUMMA Compartments:
//   - Canopy (1 layer)
//   - Snow (up to 5 layers, variable)
//   - Runoff (1 layer, tracked by OpenWQ)
//   - Soil (variable layers)
//   - Aquifer (1 layer)
//
// The interface provides five main methods:
//   1. decl()                 - Initialize OpenWQ and define compartments
//   2. openwq_run_time_start() - Called at the start of each timestep
//   3. openwq_run_space()      - Handle internal water/mass fluxes
//   4. openwq_run_space_in()   - Handle external water/mass fluxes (EWF)
//   5. openwq_run_time_end()   - Called at the end of each timestep
// =============================================================================

#ifndef OPENWQ_HYDROLINK_INCLUDED
#define OPENWQ_HYDROLINK_INCLUDED

// OpenWQ global module headers
#include "global/OpenWQ_hostModelConfig.hpp"
#include "global/OpenWQ_json.hpp"
#include "global/OpenWQ_wqconfig.hpp"
#include "global/OpenWQ_vars.hpp"

// OpenWQ functional module headers
#include "couplercalls/headerfile_CC.hpp"
#include "readjson/headerfile_nlohmann.hpp"
#include "initiate/headerfile_INIT.hpp"
#include "extwatflux_ss/headerfile_EWF_SS.hpp"
#include "units/headerfile_units.hpp"
#include "utils/headerfile_UTILS.hpp"
#include "compute/headerfile_compute.hpp"
#include "output/headerfile_OUT.hpp"

// OpenWQ model headers
#include "models_CH/headerfile_CH.hpp"
#include "models_TD/headerfile_TD.hpp"
#include "models_LE/headerfile_LE.hpp"
#include "models_SI/headerfile_SI.hpp"
#include "models_TS/headerfile_TS.hpp"

// Standard library headers
#include <iostream>
#include <time.h>
#include <vector>
#include <filesystem>
#include <memory>
#include <string>

// =============================================================================
// Global Compartment Indices
// =============================================================================
// These indices map SUMMA compartments to OpenWQ internal indices.

inline int canopy_index_openwq  = 0;  // Canopy water storage
inline int snow_index_openwq   = 1;  // Snow layers (up to max_snow_layers)
inline int runoff_index_openwq = 2;  // Surface runoff pool
inline int soil_index_openwq   = 3;  // Soil layers
inline int aquifer_index_openwq = 4;  // Aquifer/groundwater storage
inline int max_snow_layers     = 5;  // Maximum number of snow layers


// =============================================================================
// CLASSWQ_openwq Class Definition
// =============================================================================

class CLASSWQ_openwq {

private:
    // -------------------------------------------------------------------------
    // OpenWQ Module References (using smart pointers)
    // -------------------------------------------------------------------------

    // General modules
    std::unique_ptr<OpenWQ_hostModelconfig> OpenWQ_hostModelconfig_ref =
        std::make_unique<OpenWQ_hostModelconfig>();
    std::unique_ptr<OpenWQ_couplercalls> OpenWQ_couplercalls_ref =
        std::make_unique<OpenWQ_couplercalls>();
    std::unique_ptr<OpenWQ_json> OpenWQ_json_ref =
        std::make_unique<OpenWQ_json>();
    std::unique_ptr<OpenWQ_wqconfig> OpenWQ_wqconfig_ref =
        std::make_unique<OpenWQ_wqconfig>();
    std::unique_ptr<OpenWQ_units> OpenWQ_units_ref =
        std::make_unique<OpenWQ_units>();
    std::unique_ptr<OpenWQ_utils> OpenWQ_utils_ref =
        std::make_unique<OpenWQ_utils>();
    std::unique_ptr<OpenWQ_readjson> OpenWQ_readjson_ref =
        std::make_unique<OpenWQ_readjson>();
    std::unique_ptr<OpenWQ_initiate> OpenWQ_initiate_ref =
        std::make_unique<OpenWQ_initiate>();
    std::unique_ptr<OpenWQ_extwatflux_ss> OpenWQ_extwatflux_ss_ref =
        std::make_unique<OpenWQ_extwatflux_ss>();
    std::unique_ptr<OpenWQ_compute> OpenWQ_solver_ref =
        std::make_unique<OpenWQ_compute>();
    std::unique_ptr<OpenWQ_output> OpenWQ_output_ref =
        std::make_unique<OpenWQ_output>();

    // State variables (requires input from SUMMA)
    std::unique_ptr<OpenWQ_vars> OpenWQ_vars_ref;

    // Process models
    std::unique_ptr<OpenWQ_TD_model> OpenWQ_TD_model_ref =
        std::make_unique<OpenWQ_TD_model>();
    std::unique_ptr<OpenWQ_LE_model> OpenWQ_LE_model_ref =
        std::make_unique<OpenWQ_LE_model>();
    std::unique_ptr<OpenWQ_CH_model> OpenWQ_CH_model_ref =
        std::make_unique<OpenWQ_CH_model>();
    std::unique_ptr<OpenWQ_SI_model> OpenWQ_SI_model_ref =
        std::make_unique<OpenWQ_SI_model>();
    std::unique_ptr<OpenWQ_TS_model> OpenWQ_TS_model_ref =
        std::make_unique<OpenWQ_TS_model>();

    // -------------------------------------------------------------------------
    // Host Model Domain Information
    // -------------------------------------------------------------------------
    int num_HRU;              // Number of HRUs
    const float* hru_area;    // HRU areas
    long long* hruId;         // Array of HRU IDs

public:
    // Constructor and destructor
    CLASSWQ_openwq();
    ~CLASSWQ_openwq();

    // -------------------------------------------------------------------------
    // decl: Initialize OpenWQ
    // -------------------------------------------------------------------------
    // Sets up compartments, external fluxes, dependencies, and reads configuration.
    //
    // Parameters:
    //   num_HRU         - Number of HRUs
    //   nCanopy_2openwq - Number of canopy layers (fixed to 1)
    //   nSnow_2openwq   - Max number of snow layers (fixed to 5, varies at runtime)
    //   nSoil_2openwq   - Number of soil layers (variable)
    //   nRunoff_2openwq - Number of runoff layers (fixed to 1)
    //   nAquifer_2openwq - Number of aquifer layers (fixed to 1)
    //   nYdirec_2openwq - Number of y-direction layers (set to 1, unused in SUMMA)
    //   hruId           - Array of HRU identifiers
    //
    // Returns: 0 on success
    // -------------------------------------------------------------------------
    int decl(
        int num_HRU,
        int nCanopy_2openwq,
        int nSnow_2openwq,
        int nSoil_2openwq,
        int nRunoff_2openwq,
        int nAquifer_2openwq,
        int nYdirec_2openwq,
        long long hruId[]);

    // -------------------------------------------------------------------------
    // openwq_run_time_start: Begin timestep processing
    // -------------------------------------------------------------------------
    // Called at the start of each timestep to update water volumes and
    // dependency variables for each HRU.
    //
    // Parameters:
    //   last_hru_flag                    - True if this is the last HRU
    //   hru_index                        - Current HRU index (0-based)
    //   nSnow_2openwq                    - Current number of snow layers
    //   nSoil_2openwq                    - Number of soil layers
    //   simtime_summa                    - Time array [year, month, day, hour, minute]
    //   soilMoist_depVar_summa_frac      - Soil moisture fraction per layer
    //   soilTemp_depVar_summa_K          - Soil temperature [K] per layer
    //   airTemp_depVar_summa_K           - Air temperature [K]
    //   SWrad_depVar_summa_Wm2           - Shortwave radiation [W/m2]
    //   sweWatVol_stateVar_summa_m3      - Snow water equivalent volume [m3] per layer
    //   canopyWatVol_stateVar_summa_m3   - Canopy water volume [m3]
    //   soilWatVol_stateVar_summa_m3     - Soil water volume [m3] per layer
    //   aquiferWatVol_stateVar_summa_m3  - Aquifer water volume [m3]
    //
    // Returns: 0 on success
    // -------------------------------------------------------------------------
    int openwq_run_time_start(
        bool last_hru_flag,
        int hru_index,
        int nSnow_2openwq,
        int nSoil_2openwq,
        int simtime_summa[],
        double soilMoist_depVar_summa_frac[],
        double soilTemp_depVar_summa_K[],
        double airTemp_depVar_summa_K,
        double SWrad_depVar_summa_Wm2,
        double sweWatVol_stateVar_summa_m3[],
        double canopyWatVol_stateVar_summa_m3,
        double soilWatVol_stateVar_summa_m3[],
        double aquiferWatVol_stateVar_summa_m3);

    // -------------------------------------------------------------------------
    // openwq_run_space: Handle internal water/mass fluxes
    // -------------------------------------------------------------------------
    // Processes water flux between two compartment cells, transporting
    // dissolved chemicals proportionally to the water flux.
    //
    // Parameters:
    //   simtime_summa - Time array [year, month, day, hour, minute]
    //   source        - Source compartment index
    //   ix_s, iy_s, iz_s - Source cell coordinates (1-indexed from Fortran)
    //   recipient     - Recipient compartment index (-1 for outflow)
    //   ix_r, iy_r, iz_r - Recipient cell coordinates (1-indexed from Fortran)
    //   wflux_s2r     - Water flux from source to recipient [m3]
    //   wmass_source  - Water mass/volume in source cell [m3]
    //
    // Returns: 0 on success
    // -------------------------------------------------------------------------
    int openwq_run_space(
        int simtime_summa[],
        int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r, double wmass_source);

    // -------------------------------------------------------------------------
    // openwq_run_space_in: Handle external water/mass fluxes
    // -------------------------------------------------------------------------
    // Processes external water flux entering the domain (e.g., precipitation).
    //
    // Parameters:
    //   simtime_summa   - Time array [year, month, day, hour, minute]
    //   source_EWF_name - Name of external water flux source (e.g., "PRECIP")
    //   recipient       - Recipient compartment index
    //   ix_r, iy_r, iz_r - Recipient cell coordinates (1-indexed from Fortran)
    //   wflux_s2r       - Water flux entering the recipient [m3]
    //
    // Returns: 0 on success
    // -------------------------------------------------------------------------
    int openwq_run_space_in(
        int simtime_summa[],
        std::string source_EWF_name,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r);

    // -------------------------------------------------------------------------
    // openwq_run_time_end: End timestep processing
    // -------------------------------------------------------------------------
    // Called at the end of each timestep to solve equations and write outputs.
    //
    // Parameters:
    //   simtime_summa - Time array [year, month, day, hour, minute]
    //
    // Returns: 0 on success
    // -------------------------------------------------------------------------
    int openwq_run_time_end(
        int simtime_summa[]);

    // -------------------------------------------------------------------------
    // get_numHRU: Get number of HRUs
    // -------------------------------------------------------------------------
    int get_numHRU();
};

#endif // OPENWQ_HYDROLINK_INCLUDED
