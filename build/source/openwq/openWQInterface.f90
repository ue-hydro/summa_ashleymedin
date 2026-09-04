! Copyright 2020, Diogo Costa (diogo.pinhodacosta@canada.ca)
! This file is part of OpenWQ model.

! This program, openWQ, is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.

! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

! ==============================================================================
! OpenWQ Fortran-C Interface Bindings
! ==============================================================================
! This file contains the ISO C binding interface declarations that enable
! Fortran to call the C functions defined in OpenWQ_interface.h and
! implemented in OpenWQ_interface.cpp.
!
! These bindings provide the bridge between the Fortran CLASSWQ_openwq type
! and the C++ CLASSWQ_openwq class.
! ==============================================================================

interface

    ! ==========================================================================
    ! create_openwq_c: Create OpenWQ object
    ! ==========================================================================
    function create_openwq_c() bind(C, name="create_openwq")
        use iso_c_binding
        implicit none
        type(c_ptr) :: create_openwq_c
    end function

    ! ==========================================================================
    ! openwq_decl_c: Initialize OpenWQ
    ! ==========================================================================
    ! Sets up compartments, external fluxes, and dependencies.
    ! Returns: 0 (success) or -1 (failure)
    ! ==========================================================================
    function openwq_decl_c(     &
        openWQ,                 &
        num_hru,                &
        nCanopy_2openwq,        &
        nSnow_2openwq,          &
        nSoil_2openwq,          &
        nRunoff_2openwq,        &
        nAquifer_2openwq,       &
        y_direction,            &
        hruId) bind(C, name="openwq_decl")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_decl_c
        type(c_ptr), intent(in), value :: openWQ
        integer(c_int), intent(in), value  :: num_hru
        integer(c_int), intent(in), value  :: nCanopy_2openwq
        integer(c_int), intent(in), value  :: nSnow_2openwq
        integer(c_int), intent(in), value  :: nSoil_2openwq
        integer(c_int), intent(in), value  :: nAquifer_2openwq
        integer(c_int), intent(in), value  :: nRunoff_2openwq
        integer(c_int), intent(in), value  :: y_direction
        integer(c_long_long), intent(in) :: hruId(num_hru)

    end function

    ! ==========================================================================
    ! openwq_run_time_start_c: Begin timestep processing
    ! ==========================================================================
    ! Updates water volumes and dependency variables for each HRU.
    ! Returns: 0 (success) or -1 (failure)
    ! ==========================================================================
    function openwq_run_time_start_c(&
        openWQ,                             &
        last_hru_flag,                      &
        hru_index,                          &
        nSnow_2openwq,                      &
        nSoil_2openwq,                      &
        simtime_summa,                      &
        soilMoist_depVar_summa_frac,        &
        soilTemp_depVar_summa_K,            &
        airTemp_depVar_summa_K,             &
        SWrad_depVar_summa_Wm2,             &
        sweWatVol_stateVar_summa_m3,        &
        canopyWatVol_stateVar_summa_m3,     &
        soilWatVol_stateVar_summa_m3,       &
        aquiferWatVol_stateVar_summa_m3,    &
        hru_area_m2) bind(C, name="openwq_run_time_start")

        use iso_c_binding
        implicit none
        integer(c_int)                       :: openwq_run_time_start_c
        type(c_ptr),    intent(in), value    :: openWQ
        logical(c_bool),   intent(in)        :: last_hru_flag
        integer(c_int), intent(in), value    :: hru_index
        integer(c_int), intent(in), value    :: nSnow_2openwq
        integer(c_int), intent(in), value    :: nSoil_2openwq
        integer(c_int), intent(in)           :: simtime_summa(5)
        real(c_double), intent(in)           :: soilMoist_depVar_summa_frac(nSoil_2openwq)
        real(c_double), intent(in)           :: soilTemp_depVar_summa_K(nSoil_2openwq)
        real(c_double), intent(in), value    :: airTemp_depVar_summa_K
        real(c_double), intent(in), value    :: SWrad_depVar_summa_Wm2
        real(c_double), intent(in)           :: sweWatVol_stateVar_summa_m3(nSnow_2openwq)
        real(c_double), intent(in), value    :: canopyWatVol_stateVar_summa_m3
        real(c_double), intent(in)           :: soilWatVol_stateVar_summa_m3(nSoil_2openwq)
        real(c_double), intent(in), value    :: aquiferWatVol_stateVar_summa_m3
        real(c_double), intent(in), value    :: hru_area_m2

    end function

    ! ==========================================================================
    ! openwq_run_space_c: Handle internal water fluxes
    ! ==========================================================================
    ! Transports dissolved chemicals proportionally to water flux.
    ! Returns: 0 (success) or -1 (failure)
    ! ==========================================================================
    function openwq_run_space_c(&
        openWQ, &
        simtime, &
        source, ix_s, iy_s, iz_s, &
        recipient, ix_r, iy_r, iz_r, &
        wflux_s2r, &
        wmass_source) bind(C, name="openwq_run_space")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_space_c
        type(c_ptr),    intent(in), value      :: openWQ
        integer(c_int), intent(in)             :: simtime(5)
        integer(c_int), intent(in), value      :: source
        integer(c_int), intent(in), value      :: ix_s
        integer(c_int), intent(in), value      :: iy_s
        integer(c_int), intent(in), value      :: iz_s
        integer(c_int), intent(in), value      :: recipient
        integer(c_int), intent(in), value      :: ix_r
        integer(c_int), intent(in), value      :: iy_r
        integer(c_int), intent(in), value      :: iz_r
        real(c_double), intent(in), value      :: wflux_s2r
        real(c_double), intent(in), value      :: wmass_source

    end function

    ! ==========================================================================
    ! openwq_run_space_in_c: Handle external water fluxes (EWF)
    ! ==========================================================================
    ! Adds dissolved chemicals from external sources (e.g., precipitation).
    ! Returns: 0 (success) or -1 (failure)
    ! ==========================================================================
    function openwq_run_space_in_c( &
        openWQ, &
        simtime, &
        source_EWF_name, &
        recipient, ix_r, iy_r, iz_r, &
        wflux_s2r) bind(C, name="openwq_run_space_in")

        USE iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_space_in_c
        type(c_ptr), intent(in), value         :: openWQ
        integer(c_int), intent(in)             :: simtime(5)
        integer(c_int), intent(in), value      :: recipient
        integer(c_int), intent(in), value      :: ix_r
        integer(c_int), intent(in), value      :: iy_r
        integer(c_int), intent(in), value      :: iz_r
        real(c_double), intent(in), value      :: wflux_s2r
        character(c_char), intent(in)          :: source_EWF_name

    end function

    ! ==========================================================================
    ! openwq_update_runoff_vol_c: report runoff through-volume of this step
    ! ==========================================================================
    ! SUMMA's RUNOFF is a transient pool (start-of-step volume is zero);
    ! this reports the routed volume so sorption and concentration output work.
    ! Returns: 0 (success)
    ! ==========================================================================
    function openwq_update_runoff_vol_c(&
        openWQ, &
        index_hru, &
        runoff_vol_m3) bind(C, name="openwq_update_runoff_vol")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_update_runoff_vol_c
        type(c_ptr),    intent(in), value      :: openWQ
        integer(c_int), intent(in), value      :: index_hru
        real(c_double), intent(in), value      :: runoff_vol_m3

    end function

    ! ==========================================================================
    ! openwq_set_fluxvol_c: fill through-volume of a flux-concentration export
    ! ==========================================================================
    function openwq_set_fluxvol_c( &
        openWQ, iflux, ix, iy, iz, flux_vol_m3) bind(C, name="openwq_set_fluxvol")

        use iso_c_binding
        implicit none
        integer(c_int) :: openwq_set_fluxvol_c
        type(c_ptr),    intent(in), value      :: openWQ
        integer(c_int), intent(in), value      :: iflux
        integer(c_int), intent(in), value      :: ix
        integer(c_int), intent(in), value      :: iy
        integer(c_int), intent(in), value      :: iz
        real(c_double), intent(in), value      :: flux_vol_m3

    end function

    ! ==========================================================================
    ! openwq_run_time_end_c: End timestep processing
    ! ==========================================================================
    ! Solves equations and writes outputs.
    ! Returns: 0 (success) or -1 (failure)
    ! ==========================================================================
    function openwq_run_time_end_c( &
        openWQ, &
        simtime) bind(C, name="openwq_run_time_end")

        USE iso_c_binding
        implicit none
        integer(c_int) :: openwq_run_time_end_c
        type(c_ptr),    intent(in), value   :: openWQ
        integer(c_int), intent(in)          :: simtime(5)

    end function

end interface
