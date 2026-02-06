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
! SUMMA-OpenWQ Integration Module
! ==============================================================================
! This module provides the high-level interface between SUMMA and OpenWQ.
! It handles the translation of SUMMA data structures to OpenWQ inputs and
! manages water quality calculations alongside hydrological simulations.
!
! Key Subroutines:
!   - openwq_init: Initialize OpenWQ with SUMMA domain configuration
!   - openwq_run_time_start: Update water volumes at timestep start
!   - openwq_run_space_step: Process all water fluxes between compartments
!   - openwq_run_time_end: Finalize timestep and write outputs
!
! SUMMA Compartments mapped to OpenWQ:
!   - Canopy (1 layer): Interception storage
!   - Snow (up to 5 layers): Snowpack water
!   - Runoff (1 layer): Surface runoff pool (tracked by OpenWQ)
!   - Soil (variable layers): Soil water storage
!   - Aquifer (1 layer): Groundwater storage
!
! Water Flux Pathways:
!   1. Precipitation -> Canopy/Snow/Runoff
!   2. Canopy -> Snow/Runoff (drainage, unloading)
!   3. Snow internal fluxes (between layers)
!   4. Snow -> Runoff (melt)
!   5. Runoff -> Soil (infiltration)
!   6. Runoff -> OUT (surface runoff)
!   7. Soil internal fluxes (between layers)
!   8. Soil -> OUT (baseflow, exfiltration)
!   9. Soil -> Aquifer (drainage)
!   10. Aquifer -> OUT (baseflow)
!
! NOTE: Evaporation/transpiration/sublimation do NOT transport dissolved
! chemicals. Mass stays in the remaining water, causing concentration
! enrichment. This happens automatically in OpenWQ.
! ==============================================================================

module summa_openwq
   USE nrtype
   USE openWQ, only: CLASSWQ_openwq
   USE data_types, only: gru_hru_doubleVec
   implicit none
   private

   ! Public subroutines
   public :: openwq_init
   public :: openwq_run_time_start
   public :: openwq_run_space_step
   public :: openwq_run_time_end

   ! Private subroutines
   private :: openWQ_run_time_start_inner

   ! ---------------------------------------------------------------------------
   ! Module-level data
   ! ---------------------------------------------------------------------------
   ! Copy of progStruct at timestep start for flux calculations
   type(gru_hru_doubleVec), save, public :: progStruct_timestep_start

   ! OpenWQ object instance
   type(CLASSWQ_openwq), save, public :: openwq_obj

contains

! ==============================================================================
! openwq_init: Initialize OpenWQ
! ==============================================================================
! Initializes the OpenWQ object and allocates space for timestep state tracking.
!
! Parameters:
!   err - Error code (output): 0 on success
! ==============================================================================
subroutine openwq_init(err)
   USE globalData, only: gru_struc
   USE globalData, only: prog_meta
   USE globalData, only: maxLayers, maxSnowLayers
   USE allocspace_progStuct_module, only: allocGlobal_porgStruct
   implicit none

   ! Output
   integer(i4b), intent(out) :: err

   ! Local variables
   integer(i4b) :: hruCount
   integer(i4b) :: hru_i
   integer(i8b), dimension(:), allocatable :: hruId
   integer(i4b) :: nSoil
   character(len=256) :: message

   ! OpenWQ dimensions (fixed for SUMMA)
   integer(i4b) :: nCanopy_2openwq  = 1   ! Canopy has only 1 layer
   integer(i4b) :: nRunoff_2openwq  = 1   ! Runoff has only 1 layer
   integer(i4b) :: nAquifer_2openwq = 1   ! Aquifer has only 1 layer
   integer(i4b) :: nYdirec_2openwq  = 1   ! Y-direction (not used in SUMMA)

   ! Create OpenWQ object
   openwq_obj = CLASSWQ_openwq()

   ! Calculate total HRU count and soil layers
   hruCount = sum(gru_struc(:)%hruCount)
   nSoil = maxLayers - maxSnowLayers

   ! Build HRU ID array for OpenWQ
   allocate(hruId(hruCount))
   do hru_i = 1, size(gru_struc)
      hruId(hru_i) = gru_struc(hru_i)%gru_id
   end do

   ! Initialize OpenWQ with SUMMA domain configuration
   err = openwq_obj%decl(    &
      hruCount,              &
      nCanopy_2openwq,       &
      maxSnowLayers,         &
      nSoil,                 &
      nRunoff_2openwq,       &
      nAquifer_2openwq,      &
      nYdirec_2openwq,       &
      hruId)

   ! Allocate storage for previous timestep state (needed for flux calculations)
   call allocGlobal_porgStruct(prog_meta, progStruct_timestep_start, maxSnowLayers, err, message)

end subroutine openwq_init


! ==============================================================================
! openwq_run_time_start: Begin timestep processing
! ==============================================================================
! Updates OpenWQ with current water volumes and dependency variables for all
! HRUs. Called at the start of each SUMMA timestep.
!
! Parameters:
!   summa1_struc - SUMMA master data structure (input)
! ==============================================================================
subroutine openwq_run_time_start(summa1_struc)
   USE summa_type, only: summa1_type_dec
   USE var_lookup, only: iLookINDEX
   implicit none

   ! Input
   type(summa1_type_dec), intent(in) :: summa1_struc

   ! Local variables
   integer(i4b) :: openWQArrayIndex
   integer(i4b) :: iGRU
   integer(i4b) :: iHRU
   integer(i4b) :: nHRU
   integer(i4b) :: nSoil
   integer(i4b) :: nSnow
   logical(1)   :: lastHRUFlag

   summaVars: associate(&
      progStruct => summa1_struc%progStruct, &
      timeStruct => summa1_struc%timeStruct, &
      attrStruct => summa1_struc%attrStruct, &
      indxStruct => summa1_struc%indxStruct, &
      nGRU       => summa1_struc%nGRU        &
   )

   openWQArrayIndex = 0
   lastHRUFlag = .false.

   ! Loop through all GRUs and HRUs
   do iGRU = 1, nGRU
      nHRU = size(progStruct%gru(iGRU)%hru(:))
      do iHRU = 1, nHRU
         ! Flag last HRU to trigger OpenWQ timestep processing
         if (iGRU == nGRU .and. iHRU == nHRU) then
            lastHRUFlag = .true.
         end if

         ! Get current layer counts
         nSnow = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1)
         nSoil = indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1)

         ! Update OpenWQ for this HRU
         call openwq_run_time_start_inner(openWQArrayIndex, iGRU, iHRU, &
                                          summa1_struc, nSnow, nSoil, lastHRUFlag)

         openWQArrayIndex = openWQArrayIndex + 1
      end do
   end do

   end associate summaVars
end subroutine openwq_run_time_start


! ==============================================================================
! openWQ_run_time_start_inner: HRU-level timestep start processing
! ==============================================================================
! Updates OpenWQ with water volumes and dependency variables for a single HRU.
! Also copies the current progStruct to progStruct_timestep_start for flux
! calculations that need the previous timestep's water volumes.
! ==============================================================================
subroutine openWQ_run_time_start_inner(openWQArrayIndex, iGRU, iHRU, &
                                       summa1_struc, nSnow, nSoil, last_hru_flag)
   USE summa_type, only: summa1_type_dec
   USE var_lookup, only: iLookPROG
   USE var_lookup, only: iLookATTR
   USE var_lookup, only: iLookINDEX
   USE var_lookup, only: iLookVarType
   USE var_lookup, only: iLookTIME
   USE var_lookup, only: iLookFORCE
   USE globalData, only: prog_meta
   USE globalData, only: realMissing
   USE multiconst, only: iden_water
   implicit none

   ! Input
   integer(i4b), intent(in) :: openWQArrayIndex
   integer(i4b), intent(in) :: iGRU
   integer(i4b), intent(in) :: iHRU
   type(summa1_type_dec), intent(in) :: summa1_struc
   integer(i4b), intent(in) :: nSnow
   integer(i4b), intent(in) :: nSoil
   logical(1), intent(in) :: last_hru_flag

   ! Local variables
   integer(i4b) :: simtime(5)
   real(rkind)  :: canopyWatVol_stateVar_summa_m3
   real(rkind)  :: sweWatVol_stateVar_summa_m3(nSnow)
   real(rkind)  :: soilTemp_depVar_summa_K(nSoil)
   real(rkind)  :: soilWatVol_stateVar_summa_m3(nSoil)
   real(rkind)  :: soilMoist_depVar_summa_frac(nSoil)
   real(rkind)  :: aquiferWatVol_stateVar_summa_m3
   real(rkind)  :: SWrad_summa_Wm2
   integer(i4b) :: ilay
   integer(i4b) :: iVar
   integer(i4b) :: iDat
   integer(i4b) :: index
   integer(i4b) :: offset
   integer(i4b) :: err

   summaVars: associate(&
      progStruct                  => summa1_struc%progStruct, &
      timeStruct                  => summa1_struc%timeStruct, &
      hru_area_m2                 => summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea), &
      Tair_summa_K                => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanairTemp)%dat(1), &
      SWRadAtm_summa_Wm2          => summa1_struc%forcStruct%gru(iGRU)%hru(iHRU)%var(iLookFORCE%SWRadAtm), &
      scalarCanopyWat_summa_kg_m2 => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1), &
      mLayerDepth_summa_m         => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(:), &
      mLayerVolFracWat_summa_frac => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(:), &
      Tsoil_summa_K               => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerTemp)%dat(:), &
      AquiferStorWat_summa_m      => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1) &
   )

   ! -------------------------------------------------------------------------
   ! Validate required inputs
   ! -------------------------------------------------------------------------
   if (Tair_summa_K == realMissing) then
      stop 'Error: OpenWQ requires air temperature (K)'
   endif

   ! -------------------------------------------------------------------------
   ! Process unlayered variables
   ! -------------------------------------------------------------------------

   ! Shortwave radiation [W/m2]
   if (SWRadAtm_summa_Wm2 == realMissing) then
      SWrad_summa_Wm2 = 0._rkind
   else
      SWrad_summa_Wm2 = SWRadAtm_summa_Wm2
   endif

   ! Canopy water volume [m3]
   ! Convert from kg/m2 to m3: multiply by area, divide by water density
   if (scalarCanopyWat_summa_kg_m2 == realMissing) then
      canopyWatVol_stateVar_summa_m3 = 0._rkind
   else
      canopyWatVol_stateVar_summa_m3 = scalarCanopyWat_summa_kg_m2 * hru_area_m2 / iden_water
   endif

   ! Aquifer water volume [m3]
   ! Convert from m to m3: multiply by area
   if (AquiferStorWat_summa_m == realMissing) then
      stop 'Error: OpenWQ requires aquifer storage (m3)'
   endif
   aquiferWatVol_stateVar_summa_m3 = AquiferStorWat_summa_m * hru_area_m2

   ! -------------------------------------------------------------------------
   ! Process snow layers
   ! -------------------------------------------------------------------------
   if (nSnow > 0) then
      do ilay = 1, nSnow
         ! Convert volumetric fraction to volume [m3]
         if (mLayerVolFracWat_summa_frac(ilay) /= realMissing) then
            sweWatVol_stateVar_summa_m3(ilay) = &
               mLayerVolFracWat_summa_frac(ilay) * mLayerDepth_summa_m(ilay) * hru_area_m2
         else
            sweWatVol_stateVar_summa_m3(ilay) = 0._rkind
         endif
      enddo
   endif

   ! -------------------------------------------------------------------------
   ! Process soil layers
   ! -------------------------------------------------------------------------
   do ilay = 1, nSoil
      ! Soil temperature [K]
      if (Tsoil_summa_K(nSnow + ilay) == realMissing) then
         stop 'Error: OpenWQ requires soil temperature (K)'
      endif
      soilTemp_depVar_summa_K(ilay) = Tsoil_summa_K(nSnow + ilay)

      ! Soil moisture fraction (placeholder - needs implementation)
      soilMoist_depVar_summa_frac(ilay) = 0

      ! Soil water volume [m3]
      if (mLayerVolFracWat_summa_frac(nSnow + ilay) == realMissing) then
         stop 'Error: OpenWQ requires soil water (m3)'
      endif
      soilWatVol_stateVar_summa_m3(ilay) = &
         mLayerVolFracWat_summa_frac(nSnow + ilay) * hru_area_m2 * mLayerDepth_summa_m(nSnow + ilay)
   enddo

   ! -------------------------------------------------------------------------
   ! Copy progStruct to timestep start storage (for flux calculations)
   ! -------------------------------------------------------------------------
   do iVar = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var)
      do iDat = 1, size(progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat)
         select case(prog_meta(iVar)%vartype)
            case(iLookVarType%ifcSoil)
               offset = 0
            case(iLookVarType%ifcToto)
               offset = 0
            case default
               offset = 1
         end select
         do index = offset, size(progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat) - 1 + offset
            progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iVar)%dat(index) = &
               progStruct%gru(iGRU)%hru(iHRU)%var(iVar)%dat(index)
         enddo
      end do
   end do

   ! -------------------------------------------------------------------------
   ! Build simulation time array and call OpenWQ
   ! -------------------------------------------------------------------------
   simtime(1) = timeStruct%var(iLookTIME%iyyy)   ! Year
   simtime(2) = timeStruct%var(iLookTIME%im)     ! Month
   simtime(3) = timeStruct%var(iLookTIME%id)     ! Day
   simtime(4) = timeStruct%var(iLookTIME%ih)     ! Hour
   simtime(5) = timeStruct%var(iLookTIME%imin)   ! Minute

   err = openwq_obj%openwq_run_time_start(&
      last_hru_flag, &
      openWQArrayIndex, &
      nSnow, &
      nSoil, &
      simtime, &
      soilMoist_depVar_summa_frac, &
      soilTemp_depVar_summa_K, &
      Tair_summa_K, &
      SWrad_summa_Wm2, &
      sweWatVol_stateVar_summa_m3, &
      canopyWatVol_stateVar_summa_m3, &
      soilWatVol_stateVar_summa_m3, &
      aquiferWatVol_stateVar_summa_m3)

   end associate summaVars

end subroutine openWQ_run_time_start_inner


! ==============================================================================
! openwq_run_space_step: Process all water fluxes
! ==============================================================================
! Handles all water flux pathways between SUMMA compartments, transporting
! dissolved chemicals proportionally with each flux.
!
! Key flux pathways:
!   1. Precipitation -> Canopy/Snow/Runoff
!   2. Canopy drainage -> Snow/Runoff
!   3. Snow internal fluxes and melt
!   4. Runoff -> Soil (infiltration) and OUT (surface runoff)
!   5. Soil internal fluxes, baseflow, and drainage
!   6. Aquifer baseflow
!
! NOTE: Evaporation/transpiration/sublimation do NOT call openwq_run_space
! because they don't transport dissolved chemicals. The concentration
! enrichment happens automatically as SUMMA updates water volumes.
! ==============================================================================
subroutine openwq_run_space_step(summa1_struc)
   USE var_lookup, only: iLookPROG
   USE var_lookup, only: iLookTIME
   USE var_lookup, only: iLookFLUX
   USE var_lookup, only: iLookATTR
   USE var_lookup, only: iLookINDEX
   USE var_lookup, only: iLookTYPE
   USE summa_type, only: summa1_type_dec
   USE data_types, only: var_dlength, var_i
   USE globalData, only: gru_struc
   USE globalData, only: data_step
   USE globalData, only: realMissing
   USE multiconst, only: iden_ice, iden_water
   USE module_sf_noahmplsm, only: isWater
   implicit none

   ! Input
   type(summa1_type_dec), intent(in) :: summa1_struc

   ! Local variables
   integer(i4b) :: hru_index
   integer(i4b) :: iHRU
   integer(i4b) :: iGRU
   integer(i4b) :: iLayer
   integer(i4b) :: simtime(5)
   integer(i4b) :: err

   ! OpenWQ compartment indices
   integer(i4b) :: canopy_index_openwq  = 0
   integer(i4b) :: snow_index_openwq    = 1
   integer(i4b) :: runoff_index_openwq  = 2
   integer(i4b) :: soil_index_openwq    = 3
   integer(i4b) :: aquifer_index_openwq = 4

   ! Flux routing variables
   integer(i4b) :: OpenWQindex_s, OpenWQindex_r
   integer(i4b) :: iy_r, iz_r, iy_s, iz_s
   real(rkind)  :: wflux_s2r, wmass_source

   ! Unit-converted SUMMA variables
   real(rkind)  :: scalarRainfall_summa_m3
   real(rkind)  :: scalarSnowfall_summa_m3
   real(rkind)  :: scalarThroughfallRain_summa_m3
   real(rkind)  :: scalarThroughfallSnow_summa_m3
   real(rkind)  :: canopyStorWat_kg_m3
   real(rkind)  :: scalarCanopySnowUnloading_summa_m3
   real(rkind)  :: scalarCanopyLiqDrainage_summa_m3
   real(rkind)  :: scalarCanopyTranspiration_summa_m3
   real(rkind)  :: scalarCanopyEvaporation_summa_m3
   real(rkind)  :: scalarCanopySublimation_summa_m3
   real(rkind)  :: scalarRunoffVol_m3
   real(rkind)  :: scalarSurfaceRunoff_summa_m3
   real(rkind)  :: scalarInfiltration_summa_m3
   real(rkind)  :: mLayerLiqFluxSnow_summa_m3
   real(rkind)  :: iLayerLiqFluxSoil_summa_m3
   real(rkind)  :: mLayerVolFracWat_summa_m3
   real(rkind)  :: scalarSnowSublimation_summa_m3
   real(rkind)  :: scalarSfcMeltPond_summa_m3
   real(rkind)  :: scalarGroundEvaporation_summa_m3
   real(rkind)  :: scalarExfiltration_summa_m3
   real(rkind)  :: mLayerBaseflow_summa_m3
   real(rkind)  :: scalarSoilDrainage_summa_m3
   real(rkind)  :: mLayerTranspire_summa_m3
   real(rkind)  :: scalarAquiferBaseflow_summa_m3
   real(rkind)  :: scalarAquiferRecharge_summa_m3
   real(rkind)  :: scalarAquiferStorage_summa_m3
   real(rkind)  :: scalarAquiferTranspire_summa_m3

   summaVars: associate(&
      timeStruct => summa1_struc%timeStruct, &
      fluxStruct => summa1_struc%fluxStruct, &
      nGRU       => summa1_struc%nGRU)

   ! Build simulation time array
   simtime(1) = timeStruct%var(iLookTIME%iyyy)
   simtime(2) = timeStruct%var(iLookTIME%im)
   simtime(3) = timeStruct%var(iLookTIME%id)
   simtime(4) = timeStruct%var(iLookTIME%ih)
   simtime(5) = timeStruct%var(iLookTIME%imin)

   hru_index = 0
   iy_r = 1  ! SUMMA has no y-direction
   iy_s = 1

   ! Loop through all GRUs and HRUs
   do iGRU = 1, nGRU
      do iHRU = 1, gru_struc(iGRU)%hruCount
         hru_index = hru_index + 1

         ! Skip water bodies
         if (summa1_struc%typeStruct%gru(iGRU)%hru(iHRU)%var(iLookTYPE%vegTypeIndex) == isWater) cycle

         ! ====================================================================
         ! Associate SUMMA variables
         ! ====================================================================
         DomainVars: associate( &
            hru_area_m2 => summa1_struc%attrStruct%gru(iGRU)%hru(iHRU)%var(iLookATTR%HRUarea) &
         )

         PrecipVars: associate( &
            scalarRainfall_summa_kg_m2_s        => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarRainfall)%dat(1), &
            scalarSnowfall_summa_kg_m2_s        => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSnowfall)%dat(1), &
            scalarThroughfallRain_summa_kg_m2_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallRain)%dat(1), &
            scalarThroughfallSnow_summa_kg_m2_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarThroughfallSnow)%dat(1) &
         )

         CanopyVars: associate( &
            scalarCanopyWat_summa_kg_m2             => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarCanopyWat)%dat(1), &
            scalarCanopySnowUnloading_summa_kg_m2_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySnowUnloading)%dat(1), &
            scalarCanopyLiqDrainage_summa_kg_m2_s   => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyLiqDrainage)%dat(1), &
            scalarCanopyTranspiration_summa_kg_m2_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyTranspiration)%dat(1), &
            scalarCanopyEvaporation_summa_kg_m2_s   => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopyEvaporation)%dat(1), &
            scalarCanopySublimation_summa_kg_m2_s   => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarCanopySublimation)%dat(1) &
         )

         RunoffVars: associate(&
            scalarSurfaceRunoff_m_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSurfaceRunoff)%dat(1), &
            scalarInfiltration_m_s  => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarInfiltration)%dat(1) &
         )

         Snow_SoilVars: associate(&
            current_nSnow                       => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSnow)%dat(1), &
            current_nSoil                       => summa1_struc%indxStruct%gru(iGRU)%hru(iHRU)%var(iLookINDEX%nSoil)%dat(1), &
            nSnow                               => gru_struc(iGRU)%hruInfo(iHRU)%nSnow, &
            nSoil                               => gru_struc(iGRU)%hruInfo(iHRU)%nSoil, &
            mLayerDepth_summa_m                 => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerDepth)%dat(:), &
            mLayerVolFracWat_summa_frac         => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%mLayerVolFracWat)%dat(:), &
            scalarSnowSublimation_summa_kg_m2_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSnowSublimation)%dat(1), &
            scalarSfcMeltPond_kg_m2             => summa1_struc%progStruct%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarSfcMeltPond)%dat(1), &
            iLayerLiqFluxSnow_summa_m_s         => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%iLayerLiqFluxSnow)%dat(:), &
            scalarGroundEvaporation_summa_kg_m2_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarGroundEvaporation)%dat(1), &
            iLayerLiqFluxSoil_summa_m_s         => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%iLayerLiqFluxSoil)%dat(:), &
            scalarExfiltration_summa_m_s        => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarExfiltration)%dat(1), &
            mLayerBaseflow_summa_m_s            => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerBaseflow)%dat(:), &
            scalarSoilDrainage_summa_m_s        => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarSoilDrainage)%dat(1), &
            mLayerTranspire_summa_m_s           => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%mLayerTranspire)%dat(:) &
         )

         AquiferVars: associate(&
            scalarAquiferStorage_summa_m    => progStruct_timestep_start%gru(iGRU)%hru(iHRU)%var(iLookPROG%scalarAquiferStorage)%dat(1), &
            scalarAquiferRecharge_summa_m_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferRecharge)%dat(1), &
            scalarAquiferBaseflow_summa_m_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferBaseflow)%dat(1), &
            scalarAquiferTranspire_summa_m_s => fluxStruct%gru(iGRU)%hru(iHRU)%var(iLookFLUX%scalarAquiferTranspire)%dat(1) &
         )

         ! ====================================================================
         ! Convert SUMMA units to OpenWQ units
         ! Volume: m3, Time: seconds
         ! kg/m2 -> m3: multiply by (hru_area_m2 / iden_water)
         ! kg/m2/s -> m3/timestep: multiply by (hru_area_m2 * data_step / iden_water)
         ! m/s -> m3/timestep: multiply by (hru_area_m2 * data_step)
         ! ====================================================================

         ! Precipitation
         scalarRainfall_summa_m3        = scalarRainfall_summa_kg_m2_s        * hru_area_m2 * data_step / iden_water
         scalarSnowfall_summa_m3        = scalarSnowfall_summa_kg_m2_s        * hru_area_m2 * data_step / iden_water
         scalarThroughfallRain_summa_m3 = scalarThroughfallRain_summa_kg_m2_s * hru_area_m2 * data_step / iden_water
         scalarThroughfallSnow_summa_m3 = scalarThroughfallSnow_summa_kg_m2_s * hru_area_m2 * data_step / iden_water

         ! Canopy
         canopyStorWat_kg_m3                = scalarCanopyWat_summa_kg_m2             * hru_area_m2 / iden_water
         scalarCanopySnowUnloading_summa_m3 = scalarCanopySnowUnloading_summa_kg_m2_s * hru_area_m2 * data_step / iden_water
         scalarCanopyLiqDrainage_summa_m3   = scalarCanopyLiqDrainage_summa_kg_m2_s   * hru_area_m2 * data_step / iden_water
         scalarCanopyTranspiration_summa_m3 = scalarCanopyTranspiration_summa_kg_m2_s * hru_area_m2 * data_step / iden_water
         scalarCanopyEvaporation_summa_m3   = scalarCanopyEvaporation_summa_kg_m2_s   * hru_area_m2 * data_step / iden_water
         scalarCanopySublimation_summa_m3   = scalarCanopySublimation_summa_kg_m2_s   * hru_area_m2 * data_step / iden_water

         ! Runoff
         scalarSurfaceRunoff_summa_m3 = scalarSurfaceRunoff_m_s * hru_area_m2 * data_step
         scalarInfiltration_summa_m3  = scalarInfiltration_m_s  * hru_area_m2 * data_step

         ! Snow/Soil (unlayered)
         scalarSnowSublimation_summa_m3   = scalarSnowSublimation_summa_kg_m2_s   * hru_area_m2 * data_step / iden_water
         scalarGroundEvaporation_summa_m3 = scalarGroundEvaporation_summa_kg_m2_s * hru_area_m2 * data_step / iden_water
         scalarSfcMeltPond_summa_m3       = scalarSfcMeltPond_kg_m2               * hru_area_m2 / iden_water
         scalarExfiltration_summa_m3      = scalarExfiltration_summa_m_s          * hru_area_m2 * data_step
         scalarSoilDrainage_summa_m3      = scalarSoilDrainage_summa_m_s          * hru_area_m2 * data_step

         ! Aquifer
         scalarAquiferStorage_summa_m3   = scalarAquiferStorage_summa_m    * hru_area_m2
         scalarAquiferRecharge_summa_m3  = scalarAquiferRecharge_summa_m_s * hru_area_m2 * data_step
         scalarAquiferBaseflow_summa_m3  = scalarAquiferBaseflow_summa_m_s * hru_area_m2 * data_step
         scalarAquiferTranspire_summa_m3 = scalarAquiferTranspire_summa_m_s * hru_area_m2 * data_step

         ! Initialize runoff volume tracker
         scalarRunoffVol_m3 = 0._rkind

         ! ====================================================================
         ! 1. CANOPY FLUXES
         ! ====================================================================
         if (scalarCanopyWat_summa_kg_m2 /= realMissing) then

            ! 1.1 Precipitation -> Canopy (intercepted portion)
            OpenWQindex_r = canopy_index_openwq
            iz_r = 1
            wflux_s2r = (scalarRainfall_summa_m3 - scalarThroughfallRain_summa_m3) &
                      + (scalarSnowfall_summa_m3 - scalarThroughfallSnow_summa_m3)
            err = openwq_obj%openwq_run_space_in( &
               simtime, 'PRECIP', &
               OpenWQindex_r, hru_index, iy_r, iz_r, &
               wflux_s2r)

            ! 1.2 Canopy -> Snow/Runoff (drainage + unloading)
            wflux_s2r = scalarCanopySnowUnloading_summa_m3 + scalarCanopyLiqDrainage_summa_m3
            OpenWQindex_s = canopy_index_openwq
            iz_s = 1
            wmass_source = canopyStorWat_kg_m3
            if (current_nSnow > 0) then
               OpenWQindex_r = snow_index_openwq
               iz_r = 1
            else
               OpenWQindex_r = runoff_index_openwq
               iz_r = 1
               scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r
            end if
            err = openwq_obj%openwq_run_space( &
               simtime, &
               OpenWQindex_s, hru_index, iy_s, iz_s, &
               OpenWQindex_r, hru_index, iy_r, iz_r, &
               wflux_s2r, wmass_source)

            ! 1.3 Canopy evaporation/transpiration/sublimation: NO MASS TRANSPORT
            ! Water leaves but chemicals stay, causing concentration enrichment

         endif

         ! ====================================================================
         ! 2. SNOW / RUNOFF FLUXES
         ! ====================================================================

         ! 2.1 Throughfall -> Snow/Runoff
         wflux_s2r = scalarThroughfallRain_summa_m3 + scalarThroughfallSnow_summa_m3
         if (current_nSnow > 0) then
            OpenWQindex_r = snow_index_openwq
            iz_r = 1
         else
            OpenWQindex_r = runoff_index_openwq
            iz_r = 1
            scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r
         end if
         err = openwq_obj%openwq_run_space_in( &
            simtime, 'PRECIP', &
            OpenWQindex_r, hru_index, iy_r, iz_r, &
            wflux_s2r)

         ! Snow layer fluxes (only if snow present)
         if (current_nSnow > 0) then

            ! 2.2 Snow sublimation: NO MASS TRANSPORT

            ! 2.3 Snow internal fluxes (between layers)
            do iLayer = 1, nSnow - 1
               OpenWQindex_s = snow_index_openwq
               iz_s = iLayer
               mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer) * hru_area_m2 * mLayerDepth_summa_m(iLayer)
               wmass_source = mLayerVolFracWat_summa_m3
               OpenWQindex_r = snow_index_openwq
               iz_r = iLayer + 1
               mLayerLiqFluxSnow_summa_m3 = iLayerLiqFluxSnow_summa_m_s(iLayer) * hru_area_m2 * data_step
               wflux_s2r = mLayerLiqFluxSnow_summa_m3
               err = openwq_obj%openwq_run_space( &
                  simtime, &
                  OpenWQindex_s, hru_index, iy_s, iz_s, &
                  OpenWQindex_r, hru_index, iy_r, iz_r, &
                  wflux_s2r, wmass_source)
            end do

            ! 2.4 Snow drainage (bottom layer) -> Runoff
            mLayerLiqFluxSnow_summa_m3 = iLayerLiqFluxSnow_summa_m_s(nSnow) * hru_area_m2 * data_step
            wflux_s2r = mLayerLiqFluxSnow_summa_m3
            OpenWQindex_s = snow_index_openwq
            iz_s = iLayer
            mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow) * hru_area_m2 * mLayerDepth_summa_m(nSnow)
            wmass_source = mLayerVolFracWat_summa_m3
            OpenWQindex_r = runoff_index_openwq
            iz_r = 1
            scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r
            err = openwq_obj%openwq_run_space( &
               simtime, &
               OpenWQindex_s, hru_index, iy_s, iz_s, &
               OpenWQindex_r, hru_index, iy_r, iz_r, &
               wflux_s2r, wmass_source)
         end if

         ! 2.5 Surface melt pond -> Runoff (snow without layer)
         if (nSnow > 0) then
            wflux_s2r = scalarSfcMeltPond_summa_m3
            OpenWQindex_s = snow_index_openwq
            iz_s = 1
            mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow) * hru_area_m2 * mLayerDepth_summa_m(nSnow)
            wmass_source = mLayerVolFracWat_summa_m3
            OpenWQindex_r = runoff_index_openwq
            iz_r = 1
            scalarRunoffVol_m3 = scalarRunoffVol_m3 + wflux_s2r
            err = openwq_obj%openwq_run_space( &
               simtime, &
               OpenWQindex_s, hru_index, iy_s, iz_s, &
               OpenWQindex_r, hru_index, iy_r, iz_r, &
               wflux_s2r, wmass_source)
         endif

         ! ====================================================================
         ! 3. RUNOFF FLUXES
         ! ====================================================================

         ! 3.1 Runoff -> Soil (infiltration)
         wflux_s2r = scalarInfiltration_summa_m3
         OpenWQindex_s = runoff_index_openwq
         iz_s = 1
         wmass_source = scalarRunoffVol_m3
         OpenWQindex_r = soil_index_openwq
         iz_r = 1
         scalarRunoffVol_m3 = scalarRunoffVol_m3 - wflux_s2r
         err = openwq_obj%openwq_run_space( &
            simtime, &
            OpenWQindex_s, hru_index, iy_s, iz_s, &
            OpenWQindex_r, hru_index, iy_r, iz_r, &
            wflux_s2r, wmass_source)

         ! 3.2 Runoff -> OUT (surface runoff leaving system)
         OpenWQindex_s = runoff_index_openwq
         iz_s = 1
         wmass_source = scalarRunoffVol_m3
         OpenWQindex_r = -1
         iz_r = -1
         wflux_s2r = scalarRunoffVol_m3
         err = openwq_obj%openwq_run_space( &
            simtime, &
            OpenWQindex_s, hru_index, iy_s, iz_s, &
            OpenWQindex_r, hru_index, iy_r, iz_r, &
            wflux_s2r, wmass_source)

         ! ====================================================================
         ! 4. SOIL FLUXES
         ! ====================================================================

         ! 4.1 Ground evaporation: NO MASS TRANSPORT

         ! 4.2 Exfiltration (upper soil -> OUT)
         OpenWQindex_s = soil_index_openwq
         iz_s = 1
         mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSnow + 1) * hru_area_m2 * mLayerDepth_summa_m(nSnow + 1)
         wmass_source = mLayerVolFracWat_summa_m3
         OpenWQindex_r = -1
         iz_r = -1
         wflux_s2r = scalarExfiltration_summa_m3
         err = openwq_obj%openwq_run_space( &
            simtime, &
            OpenWQindex_s, hru_index, iy_s, iz_s, &
            OpenWQindex_r, hru_index, iy_r, iz_r, &
            wflux_s2r, wmass_source)

         ! 4.3 Layer baseflow (each soil layer -> OUT)
         do iLayer = 1, nSoil
            OpenWQindex_s = soil_index_openwq
            iz_s = iLayer
            mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer + nSnow) * hru_area_m2 * mLayerDepth_summa_m(iLayer + nSnow)
            wmass_source = mLayerVolFracWat_summa_m3
            OpenWQindex_r = -1
            iz_r = -1
            mLayerBaseflow_summa_m3 = mLayerBaseflow_summa_m_s(iLayer) * hru_area_m2 * data_step
            if (iLayer == 1) then
               mLayerBaseflow_summa_m3 = mLayerBaseflow_summa_m3 - scalarExfiltration_summa_m3
            endif
            wflux_s2r = mLayerBaseflow_summa_m3
            err = openwq_obj%openwq_run_space( &
               simtime, &
               OpenWQindex_s, hru_index, iy_s, iz_s, &
               OpenWQindex_r, hru_index, iy_r, iz_r, &
               wflux_s2r, wmass_source)
         end do

         ! 4.4 Soil transpiration: NO MASS TRANSPORT

         ! 4.5 Soil internal fluxes (between layers)
         do iLayer = 1, nSoil - 1
            OpenWQindex_s = soil_index_openwq
            iz_s = iLayer
            mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(iLayer + nSnow) * hru_area_m2 * mLayerDepth_summa_m(iLayer + nSnow)
            wmass_source = mLayerVolFracWat_summa_m3
            OpenWQindex_r = soil_index_openwq
            iz_r = iLayer + 1
            iLayerLiqFluxSoil_summa_m3 = iLayerLiqFluxSoil_summa_m_s(iLayer) * hru_area_m2 * data_step
            wflux_s2r = iLayerLiqFluxSoil_summa_m3
            err = openwq_obj%openwq_run_space( &
               simtime, &
               OpenWQindex_s, hru_index, iy_s, iz_s, &
               OpenWQindex_r, hru_index, iy_r, iz_r, &
               wflux_s2r, wmass_source)
         end do

         ! 4.6 Soil drainage -> Aquifer
         OpenWQindex_s = soil_index_openwq
         iz_s = nSoil
         mLayerVolFracWat_summa_m3 = mLayerVolFracWat_summa_frac(nSoil) * hru_area_m2 * mLayerDepth_summa_m(nSoil)
         wmass_source = mLayerVolFracWat_summa_m3
         OpenWQindex_r = aquifer_index_openwq
         iz_r = 1
         wflux_s2r = scalarSoilDrainage_summa_m3
         err = openwq_obj%openwq_run_space( &
            simtime, &
            OpenWQindex_s, hru_index, iy_s, iz_s, &
            OpenWQindex_r, hru_index, iy_r, iz_r, &
            wflux_s2r, wmass_source)

         ! ====================================================================
         ! 5. AQUIFER FLUXES
         ! ====================================================================

         ! 5.1 Aquifer -> OUT (baseflow)
         OpenWQindex_s = aquifer_index_openwq
         iz_s = 1
         wmass_source = scalarAquiferStorage_summa_m3
         OpenWQindex_r = -1
         iz_r = -1
         wflux_s2r = scalarAquiferBaseflow_summa_m3
         err = openwq_obj%openwq_run_space( &
            simtime, &
            OpenWQindex_s, hru_index, iy_s, iz_s, &
            OpenWQindex_r, hru_index, iy_r, iz_r, &
            wflux_s2r, wmass_source)

         ! 5.2 Aquifer transpiration: NO MASS TRANSPORT

         end associate AquiferVars
         end associate Snow_SoilVars
         end associate RunoffVars
         end associate CanopyVars
         end associate PrecipVars
         end associate DomainVars

      end do
   end do

   end associate summaVars
end subroutine openwq_run_space_step


! ==============================================================================
! openwq_run_time_end: Finalize timestep
! ==============================================================================
! Called at the end of each timestep to solve equations and write outputs.
!
! Parameters:
!   summa1_struc - SUMMA master data structure (input)
! ==============================================================================
subroutine openwq_run_time_end(summa1_struc)
   USE summa_type, only: summa1_type_dec
   USE var_lookup, only: iLookTIME
   implicit none

   ! Input
   type(summa1_type_dec), intent(in) :: summa1_struc

   ! Local variables
   integer(i4b) :: simtime(5)
   integer(i4b) :: err

   summaVars: associate(&
      timeStruct => summa1_struc%timeStruct &
   )

   simtime(1) = timeStruct%var(iLookTIME%iyyy)
   simtime(2) = timeStruct%var(iLookTIME%im)
   simtime(3) = timeStruct%var(iLookTIME%id)
   simtime(4) = timeStruct%var(iLookTIME%ih)
   simtime(5) = timeStruct%var(iLookTIME%imin)

   err = openwq_obj%openwq_run_time_end(simtime)

   end associate summaVars
end subroutine openwq_run_time_end

end module summa_openwq
