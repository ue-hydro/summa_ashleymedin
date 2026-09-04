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
// OpenWQ C Interface Header for SUMMA
// =============================================================================
// This header declares the C interface functions that allow SUMMA (Fortran)
// to call OpenWQ (C++) through ISO C bindings.
//
// The interface provides:
//   - Object creation and destruction
//   - Initialization (decl)
//   - Timestep start/end processing
//   - Internal and external water flux handling
//
// Implementation is in OpenWQ_interface.cpp
// =============================================================================

#ifndef OPENWQ_INTERFACE_H
#define OPENWQ_INTERFACE_H

#ifdef __cplusplus
extern "C" {
    class CLASSWQ_openwq;
    typedef CLASSWQ_openwq CLASSWQ_openwq;
#else
    typedef struct CLASSWQ_openwq CLASSWQ_openwq;
#endif

    // -------------------------------------------------------------------------
    // Object Management
    // -------------------------------------------------------------------------

    // Create OpenWQ object
    CLASSWQ_openwq* create_openwq();

    // Delete OpenWQ object
    void delete_openwq(CLASSWQ_openwq* openWQ);

    // -------------------------------------------------------------------------
    // Initialization
    // -------------------------------------------------------------------------

    // Initialize OpenWQ with SUMMA domain configuration
    int openwq_decl(
        CLASSWQ_openwq *openWQ,
        int hruCount,
        int nCanopy_2openwq,
        int nSnow_2openwq,
        int nSoil_2openwq,
        int nRunoff_2openwq,
        int nAquifer_2openwq,
        int nYdirec_2openwq,
        long long hruId[]);

    // -------------------------------------------------------------------------
    // Timestep Processing
    // -------------------------------------------------------------------------

    // Called at the start of each timestep to update water volumes
    int openwq_run_time_start(
        CLASSWQ_openwq *openWQ,
        bool last_hru_flag,
        int index_hru,
        int nSnow_2openwq,
        int nSoil_2openwq,
        int simtime_summa[],
        double soilMoist_depVar[],
        double soilTemp_K_depVar[],
        double airTemp_K_depVar,
        double SWrad_Wm2_depVar,
        double sweWatVol_stateVar[],
        double canopyWat,
        double soilWatVol_stateVar[],
        double aquiferStorage,
        double hru_area_m2);

    // Called at the end of each timestep to solve and write outputs
    int openwq_run_time_end(
        CLASSWQ_openwq *openWQ,
        int simtime_summa[]);

    // -------------------------------------------------------------------------
    // Water/Mass Flux Processing
    // -------------------------------------------------------------------------

    // Handle internal water fluxes between compartments
    int openwq_run_space(
        CLASSWQ_openwq *openWQ,
        int simtime_summa[],
        int source, int ix_s, int iy_s, int iz_s,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r, double wmass_source);

    // Handle external water fluxes (e.g., precipitation)
    int openwq_run_space_in(
        CLASSWQ_openwq *openWQ,
        int simtime_summa[],
        char* source_EWF_name,
        int recipient, int ix_r, int iy_r, int iz_r,
        double wflux_s2r);

    // Report the runoff through-volume of this step (see hydrolink)
    int openwq_update_runoff_vol(
        CLASSWQ_openwq *openWQ,
        int index_hru,
        double runoff_vol_m3);

    int openwq_set_fluxvol(
        CLASSWQ_openwq *openWQ,
        int iflux, int ix, int iy, int iz,
        double flux_vol_m3);

#ifdef __cplusplus
}
#endif

#endif // OPENWQ_INTERFACE_H
