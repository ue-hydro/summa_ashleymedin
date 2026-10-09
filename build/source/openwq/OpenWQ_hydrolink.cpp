// Copyright 2020, Diogo Costa (diogo.pinhodacosta@canada.ca)
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
// OpenWQ Hydrolink Implementation for SUMMA
// =============================================================================
// This file implements the CLASSWQ_openwq class methods that provide the
// interface between SUMMA and the OpenWQ water quality model.
//
// Key functionality:
//   - decl(): Initialize compartments, external fluxes, and dependencies
//   - openwq_run_time_start(): Update water volumes and dependency variables
//   - openwq_run_space(): Transport dissolved mass with internal water fluxes
//   - openwq_run_space_in(): Transport dissolved mass with external water fluxes
//   - openwq_run_time_end(): Solve equations and write outputs
// =============================================================================

#include "OpenWQ_hydrolink.h"
#include "OpenWQ_interface.h"

// =============================================================================
// Constructor and Destructor
// =============================================================================

CLASSWQ_openwq::CLASSWQ_openwq() {}

CLASSWQ_openwq::~CLASSWQ_openwq() {}

// =============================================================================
// decl: Initialize OpenWQ
// =============================================================================
// Sets up compartments, external fluxes, dependencies, and reads configuration.
// This method is called once during model initialization.
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
// =============================================================================
int CLASSWQ_openwq::decl(
    int num_HRU,
    int nCanopy_2openwq,
    int nSnow_2openwq,
    int nSoil_2openwq,
    int nRunoff_2openwq,
    int nAquifer_2openwq,
    int nYdirec_2openwq,
    long long hruId[]) {

    this->num_HRU = num_HRU;
    this->hruId = hruId;
    std::string msg_string;

    if (OpenWQ_hostModelconfig_ref->get_num_HydroComp() == 0) {

        // ---------------------------------------------------------------------
        // Define SUMMA compartments
        // Use capital letters for compartment names
        // ---------------------------------------------------------------------
        OpenWQ_hostModelconfig_ref->add_HydroComp(
            canopy_index_openwq, "SCALARCANOPYWAT",
            num_HRU, nYdirec_2openwq, nCanopy_2openwq);

        OpenWQ_hostModelconfig_ref->add_HydroComp(
            snow_index_openwq, "ILAYERVOLFRACWAT_SNOW",
            num_HRU, nYdirec_2openwq, max_snow_layers);

        OpenWQ_hostModelconfig_ref->add_HydroComp(
            runoff_index_openwq, "RUNOFF",
            num_HRU, nYdirec_2openwq, nRunoff_2openwq);

        OpenWQ_hostModelconfig_ref->add_HydroComp(
            soil_index_openwq, "ILAYERVOLFRACWAT_SOIL",
            num_HRU, nYdirec_2openwq, nSoil_2openwq);

        OpenWQ_hostModelconfig_ref->add_HydroComp(
            aquifer_index_openwq, "SCALARAQUIFER",
            num_HRU, nYdirec_2openwq, nAquifer_2openwq);

        // Transient pool with the water (and solute) delivered to the stream
        // in each step: surface runoff + exfiltration + soil baseflow +
        // aquifer baseflow. Filled in openwq_run_space, emptied in
        // openwq_run_time_end. One cell per HRU.
        OpenWQ_hostModelconfig_ref->add_HydroComp(
            stream_index_openwq, "RUNOFF_TO_STREAM",
            num_HRU, nYdirec_2openwq, 1);
        stream_vol_m3.assign(num_HRU, 0.0);
        stream_hru_area_m2.assign(num_HRU, 0.0);

        // ---------------------------------------------------------------------
        // Define external water fluxes (EWF)
        // Use capital letters for external flux names
        // ---------------------------------------------------------------------
        OpenWQ_hostModelconfig_ref->add_HydroExtFlux(
            0, "PRECIP", num_HRU, nYdirec_2openwq, 1);

        // ---------------------------------------------------------------------
        // Flux-concentration exports (openWQ prints these when the master-file
        // OUTPUT block selects them via FLUXES_CONC_TO_PRINT). The coupler fills
        // the flux through-volume at runtime.
        // Each export is named EXACTLY after the SUMMA variable it reports.
        // scalarRunoffVol_m3 draws its concentration from the RUNOFF
        // compartment (surface runoff only). averageRoutedRunoff and
        // scalarTotalRunoff draw it from the RUNOFF_TO_STREAM pool, so the
        // exported concentration is that of ALL the water delivered to the
        // stream, including the soil and aquifer baseflow. These two are the
        // ones mizuRoute ingests via EWF; they differ only in the flux
        // through-volume (used for mass output).
        // ---------------------------------------------------------------------
        OpenWQ_hostModelconfig_ref->add_FluxConcExport(
            scalarRunoffVol_fluxexp_openwq, "scalarRunoffVol_m3",
            runoff_index_openwq, num_HRU, nYdirec_2openwq, nRunoff_2openwq);
        OpenWQ_hostModelconfig_ref->add_FluxConcExport(
            averageRoutedRunoff_fluxexp_openwq, "averageRoutedRunoff",
            stream_index_openwq, num_HRU, nYdirec_2openwq, 1);
        OpenWQ_hostModelconfig_ref->add_FluxConcExport(
            scalarTotalRunoff_fluxexp_openwq, "scalarTotalRunoff",
            stream_index_openwq, num_HRU, nYdirec_2openwq, 1);

        // Initialize state variables container
        OpenWQ_vars_ref = std::make_unique<OpenWQ_vars>(
            OpenWQ_hostModelconfig_ref->get_num_HydroComp(),
            OpenWQ_hostModelconfig_ref->get_num_HydroExtFlux());

        // ---------------------------------------------------------------------
        // Define dependency variables for BGC kinetic expressions
        // These allow temperature/radiation-dependent reaction rates
        // ---------------------------------------------------------------------
        OpenWQ_hostModelconfig_ref->add_HydroDepend(
            0, "SM", num_HRU, nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq);

        OpenWQ_hostModelconfig_ref->add_HydroDepend(
            1, "Tair_K", num_HRU, nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq);

        OpenWQ_hostModelconfig_ref->add_HydroDepend(
            2, "Tsoil_K", num_HRU, nYdirec_2openwq, nSnow_2openwq + nSoil_2openwq);

        OpenWQ_hostModelconfig_ref->add_HydroDepend(
            3, "SWrad_Wm2", num_HRU, nYdirec_2openwq, 1);

        // Cell (HRU) area [m2] — used by the HBVSED sediment module to convert
        // water volumes (m3) to precipitation depth (mm) and areal sediment
        // densities (g/m2) to absolute mass (kg).
        OpenWQ_hostModelconfig_ref->add_HydroDepend(
            4, "cellArea_m2", num_HRU, nYdirec_2openwq, 1);

        // ---------------------------------------------------------------------
        // Map SUMMA element IDs (HRUs) to OpenWQ elements
        // This enables using hruId in SS/EWF JSON configuration files
        // ---------------------------------------------------------------------
        OpenWQ_hostModelconfig_ref->set_cellid_to_wqlabel("hruId");

        // Set cellid_to_wq values for referencing hostmodel element ids in outputs
        // IMPORTANT: This must be done BEFORE InitialConfig() so that SS/EWF JSON
        // files can use cell_id strings (e.g., "123456_z1") instead of indices
        std::vector<int> zdimension_cmp = {
            nCanopy_2openwq,
            max_snow_layers,
            nRunoff_2openwq,
            nSoil_2openwq,
            nAquifer_2openwq,
            1                    // RUNOFF_TO_STREAM pool
        };

        for (int cmp = 0; cmp < (int)OpenWQ_hostModelconfig_ref->get_num_HydroComp(); cmp++) {
            // Allocate the 3D structure for this compartment
            arma::Cube<double> domain_xyz(num_HRU, nYdirec_2openwq, zdimension_cmp[cmp]);
            OpenWQ_hostModelconfig_ref->set_cellid_to_wq_size(domain_xyz);

            // Set the cell_id values (format: "hruId_z<layer>")
            for (int x = 0; x < num_HRU; x++) {
                for (int z = 0; z < zdimension_cmp[cmp]; z++) {
                    OpenWQ_hostModelconfig_ref->set_cellid_to_wq_at(
                        cmp, x, 0, z,
                        std::to_string(static_cast<long long>(hruId[x])) + "_z" + std::to_string(z + 1));
                }
            }
        }

        // Set master JSON configuration file path
        OpenWQ_wqconfig_ref->set_OpenWQ_masterjson("openWQ_master.json");

        // Initialize OpenWQ with all configuration
        OpenWQ_couplercalls_ref->InitialConfig(
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_json_ref,
            *OpenWQ_wqconfig_ref,
            *OpenWQ_units_ref,
            *OpenWQ_utils_ref,
            *OpenWQ_readjson_ref,
            *OpenWQ_vars_ref,
            *OpenWQ_initiate_ref,
            *OpenWQ_TD_model_ref,
            *OpenWQ_LE_model_ref,
            *OpenWQ_CH_model_ref,
            *OpenWQ_SI_model_ref,
            *OpenWQ_TS_model_ref,
            *OpenWQ_extwatflux_ss_ref,
            *OpenWQ_output_ref);

        // Parse Sink/Source and External Water Flux data
        // Must be called AFTER InitialConfig (memory allocated)
        // and AFTER cell_id registration (so cell_id lookups succeed)
        OpenWQ_couplercalls_ref->ParseEWFandSS(
            *OpenWQ_json_ref,
            *OpenWQ_vars_ref,
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_wqconfig_ref,
            *OpenWQ_units_ref,
            *OpenWQ_utils_ref,
            *OpenWQ_output_ref,
            *OpenWQ_extwatflux_ss_ref);
    }

    return 0;
}

// =============================================================================
// openwq_run_time_start: Begin timestep processing
// =============================================================================
// Called at the start of each timestep to update water volumes and dependency
// variables for each HRU. The actual timestep processing is triggered when
// the last HRU is processed.
//
// Parameters:
//   last_hru_flag     - True if this is the last HRU
//   index_hru         - Current HRU index (0-based)
//   nSnow_2openwq     - Current number of snow layers
//   nSoil_2openwq     - Number of soil layers
//   simtime_summa     - Time array [year, month, day, hour, minute]
//   soilMoist_*       - Soil moisture fraction per layer
//   soilTemp_*        - Soil temperature [K] per layer
//   airTemp_*         - Air temperature [K]
//   SWrad_*           - Shortwave radiation [W/m2]
//   sweWatVol_*       - Snow water equivalent volume [m3] per layer
//   canopyWatVol_*    - Canopy water volume [m3]
//   soilWatVol_*      - Soil water volume [m3] per layer
//   aquiferWatVol_*   - Aquifer water volume [m3]
//
// Returns: 0 on success
// =============================================================================
int CLASSWQ_openwq::openwq_run_time_start(
    bool last_hru_flag,
    int index_hru,
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
    double aquiferWatVol_stateVar_summa_m3,
    double hru_area_m2) {

    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0],
        simtime_summa[1],
        simtime_summa[2],
        simtime_summa[3],
        simtime_summa[4],
        0);

    int runoff_vol = 0;

    // -------------------------------------------------------------------------
    // Update dependency variables (unlayered)
    // -------------------------------------------------------------------------
    OpenWQ_hostModelconfig_ref->set_dependVar_at(
        1, index_hru, 0, 0, airTemp_depVar_summa_K);

    OpenWQ_hostModelconfig_ref->set_dependVar_at(
        3, index_hru, 0, 0, SWrad_depVar_summa_Wm2);

    OpenWQ_hostModelconfig_ref->set_dependVar_at(
        4, index_hru, 0, 0, hru_area_m2);

    // -------------------------------------------------------------------------
    // Update water volumes (unlayered compartments)
    // -------------------------------------------------------------------------
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        canopy_index_openwq, index_hru, 0, 0, canopyWatVol_stateVar_summa_m3);

    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        runoff_index_openwq, index_hru, 0, 0, runoff_vol);

    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        aquifer_index_openwq, index_hru, 0, 0, aquiferWatVol_stateVar_summa_m3);

    // RUNOFF_TO_STREAM is a transient pool like RUNOFF: no water at the start
    // of the step, the delivered volume is set in openwq_set_fluxvol.
    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        stream_index_openwq, index_hru, 0, 0, 0.0);
    stream_vol_m3[index_hru] = 0.0;
    stream_hru_area_m2[index_hru] = hru_area_m2;

    // -------------------------------------------------------------------------
    // Update snow layer volumes. Layers beyond the CURRENT snowpack are zeroed
    // explicitly: when the snowpack melts (nSnow -> 0) the stored volumes would
    // otherwise keep their last nonzero values forever - which, among other
    // things, permanently trips HBVSED's snow erosion-inhibit gate.
    // -------------------------------------------------------------------------
    {
        const int nz_snow = (int) OpenWQ_hostModelconfig_ref->get_HydroComp_num_cells_z_at(snow_index_openwq);
        for (int z = 0; z < nz_snow; z++) {
            const double swe_vol = (z < nSnow_2openwq) ? sweWatVol_stateVar_summa_m3[z] : 0.0;
            OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
                snow_index_openwq, index_hru, 0, z, swe_vol);
        }
    }

    // -------------------------------------------------------------------------
    // Update soil layer volumes and dependency variables
    // -------------------------------------------------------------------------
    for (int z = 0; z < nSoil_2openwq; z++) {
        OpenWQ_hostModelconfig_ref->set_dependVar_at(
            0, index_hru, 0, z, soilMoist_depVar_summa_frac[z]);

        OpenWQ_hostModelconfig_ref->set_dependVar_at(
            2, index_hru, 0, z, soilTemp_depVar_summa_K[z]);

        OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
            soil_index_openwq, index_hru, 0, z, soilWatVol_stateVar_summa_m3[z]);
    }

    // -------------------------------------------------------------------------
    // Trigger timestep processing after all HRUs have been updated
    // -------------------------------------------------------------------------
    if (get_numHRU() - 1 == index_hru) {
        OpenWQ_couplercalls_ref->RunTimeLoopStart(
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_json_ref,
            *OpenWQ_wqconfig_ref,
            *OpenWQ_units_ref,
            *OpenWQ_utils_ref,
            *OpenWQ_readjson_ref,
            *OpenWQ_vars_ref,
            *OpenWQ_initiate_ref,
            *OpenWQ_TD_model_ref,
            *OpenWQ_LE_model_ref,
            *OpenWQ_CH_model_ref,
            *OpenWQ_SI_model_ref,
            *OpenWQ_TS_model_ref,
            *OpenWQ_extwatflux_ss_ref,
            *OpenWQ_solver_ref,
            *OpenWQ_output_ref,
            simtime);
    }

    return 0;
}

// =============================================================================
// openwq_run_space: Handle internal water/mass fluxes
// =============================================================================
// Processes water flux between two compartment cells, transporting dissolved
// chemicals proportionally to the water flux.
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
// =============================================================================
int CLASSWQ_openwq::openwq_run_space(
    int simtime_summa[],
    int source, int ix_s, int iy_s, int iz_s,
    int recipient, int ix_r, int iy_r, int iz_r,
    double wflux_s2r, double wmass_source) {

    // Convert Fortran 1-based indices to C++ 0-based indices
    // Preserve -1 for loss/outflow cases
    ix_s = std::max(-1, ix_s - 1);
    iy_s = std::max(-1, iy_s - 1);
    iz_s = std::max(-1, iz_s - 1);
    ix_r = std::max(-1, ix_r - 1);
    iy_r = std::max(-1, iy_r - 1);
    iz_r = std::max(-1, iz_r - 1);

    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0],
        simtime_summa[1],
        simtime_summa[2],
        simtime_summa[3],
        simtime_summa[4],
        0);

    // Every flux that SUMMA sends out of the domain (recipient == -1) is water
    // delivered to the stream: surface runoff, exfiltration, soil-layer
    // baseflow and aquifer baseflow. The transport below is left untouched, so
    // the source compartments (and the sediment and sorbed species that leave
    // with the surface runoff) behave exactly as before. The dissolved mass
    // that the call removes from the source is then copied into the
    // RUNOFF_TO_STREAM pool. The water volume of the pool is SUMMA's own total
    // runoff (see openwq_set_fluxvol), not the sum of the fluxes seen here, so
    // that concentration x runoff equals the mass delivered whatever set of
    // fluxes the chosen groundwater parameterization sends out of the domain.
    const bool to_stream = (recipient == -1)
        && (source != stream_index_openwq)
        && (ix_s >= 0) && (iy_s >= 0) && (iz_s >= 0);

    std::vector<double> d_source_before;
    if (to_stream) {
        const auto& d_src = (*OpenWQ_vars_ref->d_chemass_dt_transp_diss)(source);
        d_source_before.resize(d_src.n_elem);
        for (unsigned int chemi = 0; chemi < d_src.n_elem; chemi++)
            d_source_before[chemi] = d_src(chemi)(ix_s, iy_s, iz_s);
    }

    OpenWQ_couplercalls_ref->RunSpaceStep(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_TD_model_ref,
        *OpenWQ_TS_model_ref,
        *OpenWQ_LE_model_ref,
        *OpenWQ_CH_model_ref,
        *OpenWQ_SI_model_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime,
        source, ix_s, iy_s, iz_s,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r, wmass_source);

    if (to_stream) {
        const auto& d_src = (*OpenWQ_vars_ref->d_chemass_dt_transp_diss)(source);
        auto& d_pool = (*OpenWQ_vars_ref->d_chemass_dt_transp_diss)(stream_index_openwq);
        auto& mb = OpenWQ_vars_ref->mass_balance;
        for (unsigned int chemi = 0; chemi < d_src.n_elem; chemi++) {
            // dissolved mass that this flux took out of the source cell
            const double delivered = d_source_before[chemi] - d_src(chemi)(ix_s, iy_s, iz_s);
            if (delivered <= 0.0) continue;
            d_pool(chemi)(ix_s, 0, 0) += delivered;
            // It has not left the domain yet: it does when the pool is emptied
            // (openwq_run_time_end), where the tracker counts it again.
            if (mb.initialized && chemi < mb.num_species)
                mb.cumulative_out_flux[chemi] -= delivered;
        }
    }

    return 0;
}

// =============================================================================
// openwq_run_space_in: Handle external water/mass fluxes (EWF)
// =============================================================================
// Processes external water flux entering the domain (e.g., precipitation),
// adding dissolved chemicals according to the EWF configuration.
//
// Parameters:
//   simtime_summa   - Time array [year, month, day, hour, minute]
//   source_EWF_name - Name of external water flux source (e.g., "PRECIP")
//   recipient       - Recipient compartment index
//   ix_r, iy_r, iz_r - Recipient cell coordinates (1-indexed from Fortran)
//   wflux_s2r       - Water flux entering the recipient [m3]
//
// Returns: 0 on success
// =============================================================================
int CLASSWQ_openwq::openwq_run_space_in(
    int simtime_summa[],
    std::string source_EWF_name,
    int recipient, int ix_r, int iy_r, int iz_r,
    double wflux_s2r) {

    // Convert Fortran 1-based indices to C++ 0-based indices
    ix_r -= 1;
    iy_r -= 1;
    iz_r -= 1;

    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0],
        simtime_summa[1],
        simtime_summa[2],
        simtime_summa[3],
        simtime_summa[4],
        0);

    OpenWQ_couplercalls_ref->RunSpaceStep_IN(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_TD_model_ref,
        *OpenWQ_CH_model_ref,
        *OpenWQ_TS_model_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime,
        source_EWF_name,
        recipient, ix_r, iy_r, iz_r,
        wflux_s2r);

    return 0;
}

// =============================================================================
// openwq_update_runoff_vol: report the runoff through-volume of this step
// =============================================================================
// SUMMA's RUNOFF compartment is a transient routing pool: its volume is reset
// to zero at run_time_start. The coupling reports here the runoff volume
// accumulated during the space-step so that model_SI (sorption), the output
// concentration conversion, and other volume-dependent processes see the
// water that was actually routed through the compartment.
int CLASSWQ_openwq::openwq_update_runoff_vol(
    int index_hru,
    double runoff_vol_m3) {

    // Convert Fortran 1-based index to C++ 0-based
    index_hru -= 1;

    OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
        runoff_index_openwq, index_hru, 0, 0, runoff_vol_m3);

    return 0;
}

// =============================================================================
// openwq_set_fluxvol: fill the through-volume of a flux-concentration export
// =============================================================================
// Called by the coupler where the flux is computed. iflux is the 0-based export
// index (as registered via add_FluxConcExport); ix,iy,iz are Fortran 1-based
// cell indices. The output writer forms conc = chemass[src]/waterVol[src] and
// mass = conc * fluxVol, so the coupler simply reports the flux water volume.
int CLASSWQ_openwq::openwq_set_fluxvol(
    int iflux, int ix, int iy, int iz, double flux_vol_m3) {

    OpenWQ_hostModelconfig_ref->set_fluxVol_hydromodel_at(
        iflux, ix - 1, iy - 1, iz - 1, flux_vol_m3);

    // SUMMA's total runoff is the water delivered to the stream in this step:
    // it is the volume of the RUNOFF_TO_STREAM pool (reported to openWQ in
    // openwq_run_time_end).
    if (iflux == scalarTotalRunoff_fluxexp_openwq
            && ix >= 1 && ix <= (int) stream_vol_m3.size())
        stream_vol_m3[ix - 1] = (flux_vol_m3 > 0.0) ? flux_vol_m3 : 0.0;

    return 0;
}

// =============================================================================
// openwq_run_time_end: End timestep processing
// =============================================================================
// Called at the end of each timestep to solve equations and write outputs.
//
// Parameters:
//   simtime_summa - Time array [year, month, day, hour, minute]
//
// Returns: 0 on success
// =============================================================================
int CLASSWQ_openwq::openwq_run_time_end(
    int simtime_summa[]) {

    time_t simtime = OpenWQ_units_ref->convertTime_ints2time_t(
        *OpenWQ_wqconfig_ref,
        simtime_summa[0],
        simtime_summa[1],
        simtime_summa[2],
        simtime_summa[3],
        simtime_summa[4],
        0);

    // RUNOFF_TO_STREAM pool. Report the water delivered to the stream in this
    // step (SUMMA's total runoff), so that the concentration printed for the
    // pool and for the stream exports is delivered mass / delivered water, and
    // concentration x total runoff returns the delivered mass. Empty the pool: the
    // mass it holds at the start of the step is the delivery of the previous
    // step and now leaves the domain (flux = volume, so the advected fraction
    // is one).
    // As for RUNOFF, numerical dribbles (< 0.001 mm over the HRU) are reported
    // as no water instead of producing absurd concentrations.
    for (int ihru = 0; ihru < num_HRU; ihru++) {
        const double vol = stream_vol_m3[ihru];
        const double vol_reported =
            (vol >= 1.0e-6 * stream_hru_area_m2[ihru]) ? vol : 0.0;
        OpenWQ_hostModelconfig_ref->set_waterVol_hydromodel_at(
            stream_index_openwq, ihru, 0, 0, vol_reported);

        OpenWQ_couplercalls_ref->RunSpaceStep(
            *OpenWQ_hostModelconfig_ref,
            *OpenWQ_json_ref,
            *OpenWQ_wqconfig_ref,
            *OpenWQ_units_ref,
            *OpenWQ_utils_ref,
            *OpenWQ_readjson_ref,
            *OpenWQ_vars_ref,
            *OpenWQ_initiate_ref,
            *OpenWQ_TD_model_ref,
            *OpenWQ_TS_model_ref,
            *OpenWQ_LE_model_ref,
            *OpenWQ_CH_model_ref,
            *OpenWQ_SI_model_ref,
            *OpenWQ_extwatflux_ss_ref,
            *OpenWQ_solver_ref,
            *OpenWQ_output_ref,
            simtime,
            stream_index_openwq, ihru, 0, 0,
            -1, -1, -1, -1,
            1.0, 1.0);
    }

    OpenWQ_couplercalls_ref->RunTimeLoopEnd(
        *OpenWQ_hostModelconfig_ref,
        *OpenWQ_json_ref,
        *OpenWQ_wqconfig_ref,
        *OpenWQ_units_ref,
        *OpenWQ_utils_ref,
        *OpenWQ_readjson_ref,
        *OpenWQ_vars_ref,
        *OpenWQ_initiate_ref,
        *OpenWQ_TD_model_ref,
        *OpenWQ_LE_model_ref,
        *OpenWQ_CH_model_ref,
        *OpenWQ_SI_model_ref,
        *OpenWQ_TS_model_ref,
        *OpenWQ_extwatflux_ss_ref,
        *OpenWQ_solver_ref,
        *OpenWQ_output_ref,
        simtime);

    return 0;
}

// =============================================================================
// get_numHRU: Get number of HRUs
// =============================================================================
int CLASSWQ_openwq::get_numHRU() {
    return this->num_HRU;
}
