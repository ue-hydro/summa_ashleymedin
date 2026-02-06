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
! OpenWQ Fortran Module for SUMMA
! ==============================================================================
! This module provides the Fortran wrapper for the OpenWQ water quality model.
! It defines the CLASSWQ_openwq derived type that encapsulates a C pointer
! to the C++ OpenWQ object, enabling Fortran-C++ interoperability.
!
! The module provides type-bound procedures that map to C interface functions:
!   - decl: Initialize OpenWQ with SUMMA domain configuration
!   - openwq_run_time_start: Called at the start of each timestep
!   - openwq_run_space: Handle internal water fluxes
!   - openwq_run_space_in: Handle external water fluxes (e.g., precipitation)
!   - openwq_run_time_end: Called at the end of each timestep
! ==============================================================================

module openwq

   USE, intrinsic :: iso_c_binding
   USE nrtype
   private
   public :: CLASSWQ_openwq

   include "openWQInterface.f90"

   ! ===========================================================================
   ! CLASSWQ_openwq Derived Type
   ! ===========================================================================
   ! Encapsulates the C pointer to the OpenWQ C++ object.
   ! Type-bound procedures delegate to C interface functions.
   ! ===========================================================================
   type CLASSWQ_openwq
      private
      type(c_ptr) :: ptr  ! Pointer to C++ CLASSWQ_openwq object

   contains
      procedure :: decl => openWQ_init
      procedure :: openwq_run_time_start => openwq_run_time_start
      procedure :: openwq_run_space => openwq_run_space
      procedure :: openwq_run_space_in => openwq_run_space_in
      procedure :: openwq_run_time_end => openwq_run_time_end

   end type

   ! Constructor interface
   interface CLASSWQ_openwq
      procedure create_openwq
   end interface

contains

   ! ===========================================================================
   ! create_openwq: Constructor
   ! ===========================================================================
   ! Creates a new OpenWQ object by calling the C interface constructor.
   ! ===========================================================================
   function create_openwq()
      implicit none
      type(CLASSWQ_openwq) :: create_openwq
      create_openwq%ptr = create_openwq_c()
   end function

   ! ===========================================================================
   ! openWQ_init: Initialize OpenWQ (decl)
   ! ===========================================================================
   ! Initializes OpenWQ with SUMMA domain configuration including compartments,
   ! external fluxes, and dependencies.
   !
   ! Returns: 0 on success, -1 on failure
   ! ===========================================================================
   integer function openWQ_init( &
      this,                      &
      num_hru,                   &
      nCanopy_2openwq,           &
      nSnow_2openwq,             &
      nSoil_2openwq,             &
      nRunoff_2openwq,           &
      nAquifer_2openwq,          &
      nYdirec_2openwq,           &
      hruId)

      implicit none
      class(CLASSWQ_openwq) :: this
      integer(i4b), intent(in) :: num_hru
      integer(i4b), intent(in) :: nCanopy_2openwq
      integer(i4b), intent(in) :: nSnow_2openwq
      integer(i4b), intent(in) :: nSoil_2openwq
      integer(i4b), intent(in) :: nRunoff_2openwq
      integer(i4b), intent(in) :: nAquifer_2openwq
      integer(i4b), intent(in) :: nYdirec_2openwq
      integer(c_long_long), intent(in) :: hruId(num_hru)

      openWQ_init = openwq_decl_c(  &
         this%ptr,                  &
         num_hru,                   &
         nCanopy_2openwq,           &
         nSnow_2openwq,             &
         nSoil_2openwq,             &
         nRunoff_2openwq,           &
         nAquifer_2openwq,          &
         nYdirec_2openwq,           &
         hruId)

   end function

   ! ===========================================================================
   ! openwq_run_time_start: Begin timestep processing
   ! ===========================================================================
   ! Called at the start of each timestep to update water volumes and
   ! dependency variables for each HRU.
   !
   ! Returns: 0 on success, -1 on failure
   ! ===========================================================================
   integer function openwq_run_time_start(   &
      this,                                  &
      last_hru_flag,                         &
      hru_index,                             &
      nSnow_2openwq,                         &
      nSoil_2openwq,                         &
      simtime,                               &
      soilMoist_depVar_summa_frac,           &
      soilTemp_depVar_summa_K,               &
      airTemp_depVar_summa_K,                &
      SWrad_depVar_summa_Wm2,                &
      sweWatVol_stateVar_summa_m3,           &
      canopyWatVol_stateVar_summa_m3,        &
      soilWatVol_stateVar_summa_m3,          &
      aquiferWatVol_stateVar_summa_m3)

      implicit none
      class(CLASSWQ_openwq)      :: this
      logical(1), intent(in)     :: last_hru_flag
      integer(i4b), intent(in)   :: hru_index
      integer(i4b), intent(in)   :: nSnow_2openwq
      integer(i4b), intent(in)   :: nSoil_2openwq
      integer(i4b), intent(in)   :: simtime(5)
      real(rkind),  intent(in)   :: airTemp_depVar_summa_K
      real(rkind),  intent(in)   :: SWrad_depVar_summa_Wm2
      real(rkind),  intent(in)   :: soilTemp_depVar_summa_K(nSoil_2openwq)
      real(rkind),  intent(in)   :: soilMoist_depVar_summa_frac(nSoil_2openwq)
      real(rkind),  intent(in)   :: canopyWatVol_stateVar_summa_m3
      real(rkind),  intent(in)   :: sweWatVol_stateVar_summa_m3(nSnow_2openwq)
      real(rkind),  intent(in)   :: soilWatVol_stateVar_summa_m3(nSoil_2openwq)
      real(rkind),  intent(in)   :: aquiferWatVol_stateVar_summa_m3

      openwq_run_time_start = openwq_run_time_start_c( &
         this%ptr,                              &
         last_hru_flag,                         &
         hru_index,                             &
         nSnow_2openwq,                         &
         nSoil_2openwq,                         &
         simtime,                               &
         soilMoist_depVar_summa_frac,           &
         soilTemp_depVar_summa_K,               &
         airTemp_depVar_summa_K,                &
         SWrad_depVar_summa_Wm2,                &
         sweWatVol_stateVar_summa_m3,           &
         canopyWatVol_stateVar_summa_m3,        &
         soilWatVol_stateVar_summa_m3,          &
         aquiferWatVol_stateVar_summa_m3)

   end function

   ! ===========================================================================
   ! openwq_run_space: Handle internal water fluxes
   ! ===========================================================================
   ! Processes water flux between two compartment cells, transporting
   ! dissolved chemicals proportionally to the water flux.
   !
   ! Returns: 0 on success, -1 on failure
   ! ===========================================================================
   integer function openwq_run_space(  &
      this,                            &
      simtime,                         &
      source, ix_s, iy_s, iz_s,        &
      recipient, ix_r, iy_r, iz_r,     &
      wflux_s2r, wmass_source)

      implicit none
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5)
      integer(i4b), intent(in)   :: source
      integer(i4b), intent(in)   :: ix_s
      integer(i4b), intent(in)   :: iy_s
      integer(i4b), intent(in)   :: iz_s
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r
      real(rkind),  intent(in)   :: wmass_source

      openwq_run_space = openwq_run_space_c( &
         this%ptr,                           &
         simtime,                            &
         source, ix_s, iy_s, iz_s,           &
         recipient, ix_r, iy_r, iz_r,        &
         wflux_s2r, wmass_source)

   end function

   ! ===========================================================================
   ! openwq_run_space_in: Handle external water fluxes (EWF)
   ! ===========================================================================
   ! Processes external water flux entering the domain (e.g., precipitation),
   ! adding dissolved chemicals according to the EWF configuration.
   !
   ! Returns: 0 on success, -1 on failure
   ! ===========================================================================
   integer function openwq_run_space_in(  &
      this,                               &
      simtime,                            &
      source_EWF_name,                    &
      recipient, ix_r, iy_r, iz_r,        &
      wflux_s2r)

      implicit none
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5)
      integer(i4b), intent(in)   :: recipient
      integer(i4b), intent(in)   :: ix_r
      integer(i4b), intent(in)   :: iy_r
      integer(i4b), intent(in)   :: iz_r
      real(rkind),  intent(in)   :: wflux_s2r
      character(*), intent(in)   :: source_EWF_name

      openwq_run_space_in = openwq_run_space_in_c( &
         this%ptr,                                 &
         simtime,                                  &
         source_EWF_name,                          &
         recipient, ix_r, iy_r, iz_r,              &
         wflux_s2r)

   end function

   ! ===========================================================================
   ! openwq_run_time_end: End timestep processing
   ! ===========================================================================
   ! Called at the end of each timestep to solve equations and write outputs.
   !
   ! Returns: 0 on success, -1 on failure
   ! ===========================================================================
   integer function openwq_run_time_end(  &
      this,                               &
      simtime)

      implicit none
      class(CLASSWQ_openwq)      :: this
      integer(i4b), intent(in)   :: simtime(5)

      openwq_run_time_end = openwq_run_time_end_c( &
         this%ptr,                                 &
         simtime)

   end function

end module openwq
