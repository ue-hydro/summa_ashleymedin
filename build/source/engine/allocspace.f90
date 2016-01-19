! SUMMA - Structure for Unifying Multiple Modeling Alternatives
! Copyright (C) 2014-2015 NCAR/RAL
!
! This file is part of SUMMA
!
! For more information see: http://www.ral.ucar.edu/projects/summa
!
! This program is free software: you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation, either version 3 of the License, or
! (at your option) any later version.
!
! This program is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.

module allocspace_module
USE nrtype
implicit none
private
public::allocate_gru_struc
public::init_metad
public::alloc_stim
public::alloc_time
public::alloc_forc
public::alloc_attr
public::alloc_type
public::alloc_mpar
public::alloc_mvar
public::alloc_indx
public::alloc_bpar
public::alloc_bvar
! define missing values
integer(i4b),parameter :: missingInteger=-9999
real(dp),parameter     :: missingDouble=-9999._dp
contains


! ************************************************************************************************
 ! public subroutine read_gru_map: read in gru_hru spatial domain mask
 ! ************************************************************************************************
 subroutine allocate_gru_struc(nGRU,nHRU,err,message)
 USE netcdf
 USE nr_utility_module,only:arth
 ! provide access to subroutines
 USE netcdf_util_module,only:nc_file_open             ! open netCDF file
 ! provide access to data
 USE summaFileManager,only:SETNGS_PATH             ! path for metadata files
 USE summaFileManager,only:LOCAL_ATTRIBUTES        ! file containing information on local attributes
 USE data_struc,only: index_map                    ! relating different indexing system
 USE data_struc,only: gru_struc                    ! gru-hru structures
 
 implicit none
 ! define output
 integer(i4b),intent(out)             :: nGRU               ! number of grouped response units
 integer(i4b),intent(out)             :: nHRU               ! number of hydrologic response units
 integer(i4b),intent(out)             :: err                ! error code
 character(*),intent(out)             :: message            ! error message
 ! define general variables
 character(len=256)                   :: cmessage           ! error message for downwind routine
 character(LEN=256)                   :: infile             ! input filename
 
 ! define variables for NetCDF file operation
 integer(i4b)                         :: mode               ! netCDF file open mode
 integer(i4b)                         :: ncid               ! integer variables for NetCDF IDs
 integer(i4b)                         :: varid              ! variable id from netcdf file
 integer(i4b)                         :: hruDimID           ! integer variables for NetCDF IDs
 integer(i4b)                         :: gruDimID           ! integer variables for NetCDF IDs
 integer(i4b),allocatable             :: hru_id(:)          ! unique id of hru over entire domain
 integer(i4b),allocatable             :: gru_id(:)          ! unique ids of GRUs
 integer(i4b),allocatable             :: hru2gru_id(:)      ! unique GRU ids at each HRU
      
 ! define local variables
 integer(i4b)                         :: hruCount           ! number of hrus in a gru
 integer(i4b)                         :: iGRU               ! index of a GRU and HRU

 ! Start procedure here
 err=0; message="allocate_gru_hru_map/"

 ! check that gru_struc structure is initialized
 if(associated(gru_struc))then; deallocate(gru_struc); endif

 ! **********************************************************************************************
 ! (1) open hru_attributes file and allocate memory for data structure according to data file
 ! **********************************************************************************************
 ! build filename
 infile = trim(SETNGS_PATH)//trim(LOCAL_ATTRIBUTES)
 ! open file
 mode=nf90_NoWrite
 call nc_file_open(trim(infile), mode, ncid, err, cmessage)
 if(err/=0)then; message=trim(message)//trim(cmessage); return; endif

 ! get gru_ix dimension length
 err = nf90_inq_dimid(ncid, "ngru", gruDimID)
 err = nf90_inquire_dimension(ncid, gruDimID, len = nGRU)
 if(err/=0)then; message=trim(message)//'problem reading ngru'; return; endif

 ! get hru_dim dimension length (ie., the number of hrus in entire domain)
 err = nf90_inq_dimid(ncid, "nhru", hruDimID)
 err = nf90_inquire_dimension(ncid, hruDimID, len = nHRU)
 if(err/=0)then; message=trim(message)//'problem reading nhru'; return; endif
 
 allocate(gru_struc(nGRU), index_map(nHRU), stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for mapping structures'; return; endif
 
 allocate(gru_id(nGRU),hru_id(nHRU),hru2gru_id(nHRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//'problem allocating space for zLocalAttributes gru-hru correspondence vectors'; return; endif

 ! read gru_id from netcdf file
 err = nf90_inq_varid(ncid, "gruId", varid)
 err = nf90_get_var(ncid,varid,gru_id)

 ! read the HRU to GRU mapping
 ! (get the variable id)
 err = nf90_inq_varid(ncid, "hru2gruId", varid)
 if(err/=0)then; message=trim(message)//'problem inquiring hru2gruId'; return; endif
 ! (read the data)
 err = nf90_get_var(ncid,varid,hru2gru_id)
 if(err/=0)then; message=trim(message)//'problem reading hru2gruId'; return; endif

 ! allocate mapping array
 do iGRU=1,nGRU
  hruCount = count(hru2gru_id==gru_id(iGRU))
  gru_struc(iGRU)%hruCount = hruCount
  allocate(gru_struc(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//'problem allocating space for hru in gru_struc'; return; endif
 enddo

 ! close the HRU_ATTRIBUTES netCDF file
 err = nf90_close(ncid)
 if(err/=0)then; err=20; message=trim(message)//'error closing zLocalAttributes file'; return; endif
 end subroutine allocate_gru_struc


 ! ************************************************************************************************
 ! public subroutine init_metad: initialize metadata structures
 ! ************************************************************************************************
 subroutine init_metad(err,message)
 ! used to initialize the metadata structures
 USE var_lookup,only:maxvarTime,maxvarForc,maxvarAttr,maxvarType    ! maximum number variables in each data structure
 USE var_lookup,only:maxvarMpar,maxvarMvar,maxvarIndx               ! maximum number variables in each data structure
 USE var_lookup,only:maxvarBpar,maxvarBvar                          ! maximum number variables in each data structure
 USE data_struc,only:time_meta,forc_meta,attr_meta,type_meta        ! metadata structures
 USE data_struc,only:mpar_meta,mvar_meta,indx_meta                  ! metadata structures
 USE data_struc,only:bpar_meta,bvar_meta                            ! metadata structures
 implicit none
 ! declare variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! initialize errors
 err=0; message="init_model/"
 ! ensure metadata structures are deallocated
 if (associated(time_meta)) deallocate(time_meta)
 if (associated(forc_meta)) deallocate(forc_meta)
 if (associated(attr_meta)) deallocate(attr_meta)
 if (associated(type_meta)) deallocate(type_meta)
 if (associated(mpar_meta)) deallocate(mpar_meta)
 if (associated(mvar_meta)) deallocate(mvar_meta)
 if (associated(indx_meta)) deallocate(indx_meta)
 if (associated(bpar_meta)) deallocate(bpar_meta)
 if (associated(bvar_meta)) deallocate(bvar_meta)
 ! allocate metadata structures
 allocate(time_meta(maxvarTime),forc_meta(maxvarForc),attr_meta(maxvarAttr),type_meta(maxvarType),&
          mpar_meta(maxvarMpar),mvar_meta(maxvarMvar),indx_meta(maxvarIndx),&
          bpar_meta(maxvarBpar),bvar_meta(maxvarBvar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateMetadata"; return; endif
 end subroutine init_metad


 ! ************************************************************************************************
 ! public subroutine alloc_stim: initialize data structures for scalar time structures
 ! ************************************************************************************************
 subroutine alloc_stim(datastr,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:var_i,time_meta                 ! data structures
 implicit none
 ! dummy variables
 type(var_i),intent(out),pointer      :: datastr     ! data structure to allocate
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! initialize errors
 err=0; message="alloc_stim/"
 ! check that the metadata structure is allocated
 if(.not.associated(time_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(datastr)) deallocate(datastr)
 allocate(datastr,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 allocate(datastr%var(size(time_meta)),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! set values to missing
 datastr%var(:) = missingInteger
 end subroutine alloc_stim

 ! ************************************************************************************************
 ! public subroutine alloc_time: initialize data structures for time structures
 ! ************************************************************************************************
 subroutine alloc_time(nGRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:time_gru,time_meta              ! data structures
 USE data_struc,only:gru_struc            ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nGRU        ! number of GRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iGRU        ! loop through GRUs
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b)                         :: hruCount    ! number of HRUs within a GRU
 ! initialize errors
 err=0; message="alloc_time/"
 ! check that the metadata structure is allocated
 if(.not.associated(time_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(time_gru)) deallocate(time_gru)

 allocate(time_gru(nGRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 
 nVar = size(time_meta)
 do iGRU=1,nGRU 
  hruCount = gru_struc(iGRU)%hruCount
  allocate(time_gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  
  do iHRU=1,hruCount
   allocate(time_gru(iGRU)%hru(iHRU)%var(nVar),stat=err)
   if(err/=0)then; err=20; message=trim(message)//"problemAllocateData3rdLevel"; return; endif
   ! fill data with missing values
   time_gru(iGRU)%hru(iHRU)%var(:) = missingInteger
  enddo ! end of iHRU loop
 enddo ! end of iGRU loop
 end subroutine alloc_time


 ! ************************************************************************************************
 ! public subroutine alloc_forc: initialize data structures for model forcing data
 ! ************************************************************************************************
 subroutine alloc_forc(nGRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:forc_gru,forc_meta             ! data structures
 USE data_struc,only:gru_struc            ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nGRU        ! number of GRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iGRU        ! loop through GRUs
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b)                         :: hruCount    ! number of HRUs within a GRU
 ! initialize errors
 err=0; message="alloc_forc/"
 ! check that the metadata structure is allocated
 if(.not.associated(forc_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(forc_gru)) deallocate(forc_gru)
 allocate(forc_gru(nGRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(forc_meta)
 do iGRU=1,nGRU
  hruCount = gru_struc(iGRU)%hruCount
  allocate(forc_gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  
  do iHRU=1,hruCount
   allocate(forc_gru(iGRU)%hru(iHRU)%var(nVar),stat=err)
   if(err/=0)then; err=20; message=trim(message)//"problemAllocateData3rdLevel"; return; endif
   ! fill data with missing values
   forc_gru(iGRU)%hru(iHRU)%var(:) = missingDouble
  enddo ! end of iHRU loop
 enddo ! end of iGRU loop
 
 end subroutine alloc_forc


 ! ************************************************************************************************
 ! public subroutine alloc_attr: initialize data structures for local attributes
 ! ************************************************************************************************
 subroutine alloc_attr(nGRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:attr_meta,attr_gru   ! data structures
 USE data_struc,only:gru_struc            ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nGRU        ! number of GRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iGRU        ! loop through GRUs
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b)                         :: hruCount    ! number of HRUs within a GRU
 ! initialize errors
 err=0; message="alloc_attr/"
 ! check that the metadata structure is allocated
 if(.not.associated(attr_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(attr_gru)) deallocate(attr_gru)
 allocate(attr_gru(nGRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! initialize second level data structure
 nVar = size(attr_meta)
 do iGRU=1,nGRU
  hruCount = gru_struc(iGRU)%hruCount
  allocate(attr_gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif

  do iHRU=1,hruCount
   allocate(attr_gru(iGRU)%hru(iHRU)%var(nVar),stat=err)
   if(err/=0)then; err=20; message=trim(message)//"problemAllocateData3ndLevel"; return; endif
   ! fill data with missing values
   attr_gru(iGRU)%hru(iHRU)%var(:) = missingDouble
  end do  ! end of iHRU loop
 end do   ! end of iGRU loop
 end subroutine alloc_attr


 ! *************************************************************************************************
 ! public subroutine alloc_type: initialize data structures for local classification of veg, soil, etc.
 ! *************************************************************************************************
 subroutine alloc_type(nGRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:type_gru,type_meta             ! data structures
 USE data_struc,only:gru_struc                      ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nGRU        ! number of GRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iGRU        ! loop through HRUs
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b)                         :: hruCount    ! number of HRUs within a GRU
 ! initialize errors
 err=0; message="f-alloc_type/"
 ! check that the metadata structure is allocated
 if(.not.associated(type_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
! initialize top-level data structure
 if(associated(type_gru)) deallocate(type_gru)

 allocate(type_gru(nGRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 
 nVar = size(type_meta)
 do iGRU=1,nGRU 
  hruCount = gru_struc(iGRU)%hruCount
  allocate(type_gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  
  do iHRU=1,hruCount
   allocate(type_gru(iGRU)%hru(iHRU)%var(nVar),stat=err)
   if(err/=0)then; err=20; message=trim(message)//"problemAllocateData3rdLevel"; return; endif
   ! fill data with missing values
   type_gru(iGRU)%hru(iHRU)%var(:) = missingInteger
  enddo ! end of iHRU loop
 enddo ! end of iGRU loop
 end subroutine alloc_type


 ! *************************************************************************************************
 ! public subroutine alloc_mpar: initialize data structures for model parameters
 ! *************************************************************************************************
 subroutine alloc_mpar(nGRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:mpar_gru,mpar_meta             ! data structures
 USE data_struc,only:gru_struc                      ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nGRU        ! number of GRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iGRU        ! loop through GRUs
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nPar        ! number of parameters
 integer(i4b)                         :: hruCount    ! number of HRUs within a GRU
 ! initialize errors
 err=0; message="f-alloc_mpar/"
 ! check that the metadata structure is allocated
 if(.not.associated(mpar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
! initialize top-level data structure
 if(associated(mpar_gru)) deallocate(mpar_gru)

 allocate(mpar_gru(nGRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 
 nPar = size(mpar_meta)
 do iGRU=1,nGRU 
  hruCount = gru_struc(iGRU)%hruCount
  allocate(mpar_gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  
  do iHRU=1,hruCount
   allocate(mpar_gru(iGRU)%hru(iHRU)%var(nPar),stat=err)
   if(err/=0)then; err=20; message=trim(message)//"problemAllocateData3rdLevel"; return; endif
   ! fill data with missing values
   mpar_gru(iGRU)%hru(iHRU)%var(:) = missingDouble
  enddo ! end of iHRU loop
 enddo ! end of iGRU loop
 end subroutine alloc_mpar


 ! *************************************************************************************************
 ! public subroutine alloc_mvar: initialize data structures for model variables
 ! *************************************************************************************************
 subroutine alloc_mvar(nGRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:mvar_gru,mvar_meta             ! data structures
 USE data_struc,only:gru_struc                      ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nGRU        ! number of GRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iGRU        ! loop through GRUs
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b)                         :: hruCount    ! number of HRUs within a GRU
 ! initialize errors
 err=0; message="alloc_mvar/"
 ! check that the metadata structure is allocated
 if(.not.associated(mvar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(mvar_gru)) deallocate(mvar_gru)

 allocate(mvar_gru(nGRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 
 nVar = size(mvar_meta)
 do iGRU=1,nGRU 
  hruCount = gru_struc(iGRU)%hruCount
  allocate(mvar_gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  
  do iHRU=1,hruCount
   allocate(mvar_gru(iGRU)%hru(iHRU)%var(nVar),stat=err)
   if(err/=0)then; err=20; message=trim(message)//"problemAllocateData3rdLevel"; return; endif
  enddo ! end of iHRU loop
 enddo ! end of iGRU loop
 end subroutine alloc_mvar


 ! *************************************************************************************************
 ! public subroutine alloc_indx: initialize structure components for model indices
 ! *************************************************************************************************
 subroutine alloc_indx(nGRU,err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:indx_gru,indx_meta             ! data structures
 USE data_struc,only:gru_struc                      ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(in)              :: nGRU        ! number of GRUs
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iGRU        ! loop through HRUs
 integer(i4b)                         :: iHRU        ! loop through HRUs
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b)                         :: hruCount    ! number of HRUs within a GRU
 ! initialize errors
 err=0; message="alloc_indx/"
 ! check that the metadata structure is allocated
 if(.not.associated(indx_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(indx_gru)) deallocate(indx_gru)

 allocate(indx_gru(nGRU),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 
 nVar = size(indx_meta)
 do iGRU=1,nGRU 
  hruCount = gru_struc(iGRU)%hruCount
  allocate(indx_gru(iGRU)%hru(hruCount),stat=err)
  if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
  
  do iHRU=1,hruCount
   allocate(indx_gru(iGRU)%hru(iHRU)%var(nVar),stat=err)
   if(err/=0)then; err=20; message=trim(message)//"problemAllocateData3rdLevel"; return; endif
  enddo ! end of iHRU loop
 enddo ! end of iGRU loop
 end subroutine alloc_indx


 ! *************************************************************************************************
 ! public subroutine alloc_bpar: initialize data structures for basin-average model parameters
 ! *************************************************************************************************
 subroutine alloc_bpar(err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:bpar_data,bpar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: nPar        ! number of parameters
 ! initialize errors
 err=0; message="alloc_bpar/"
 ! check that the metadata structure is allocated
 if(.not.associated(bpar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(bpar_data)) deallocate(bpar_data)
 allocate(bpar_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! get the number of parameters
 nPar = size(bpar_meta)
 ! initialize second level data structure
 allocate(bpar_data%var(nPar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! set values to missing
 bpar_data%var(:) = missingDouble
 end subroutine alloc_bpar


 ! *************************************************************************************************
 ! public subroutine alloc_bvar: initialize data structures for basin-average model variables
 ! *************************************************************************************************
 subroutine alloc_bvar(err,message)
 ! used to initialize structure components for model variables
 USE data_struc,only:bvar_data,bvar_meta             ! data structures
 implicit none
 ! dummy variables
 integer(i4b),intent(out)             :: err         ! error code
 character(*),intent(out)             :: message     ! error message
 ! local variables
 integer(i4b)                         :: iVar        ! index of variables
 integer(i4b)                         :: nVar        ! number of variables
 integer(i4b),parameter               :: nTimeDelay=2000 ! number of elements in the time delay histogram
 ! initialize errors
 err=0; message="alloc_bvar/"
 ! check that the metadata structure is allocated
 if(.not.associated(bvar_meta))then
  err=10; message=trim(message)//"metadataNotInitialized"; return
 endif
 ! initialize top-level data structure
 if(associated(bvar_data)) deallocate(bvar_data)
 allocate(bvar_data,stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateDataTopLevel"; return; endif
 ! get the number of parameters
 nVar = size(bvar_meta)
 ! initialize second level data structure
 allocate(bvar_data%var(nVar),stat=err)
 if(err/=0)then; err=20; message=trim(message)//"problemAllocateData2ndLevel"; return; endif
 ! initialize third-level data structures
 do iVar=1,nVar
  select case(bvar_meta(ivar)%vartype)
   case('scalarv'); allocate(bvar_data%var(ivar)%dat(1),stat=err)
   case('routing'); allocate(bvar_data%var(ivar)%dat(nTimeDelay),stat=err)
   case default
    err=40; message=trim(message)//"unknownVariableType[name='"//trim(bvar_meta(ivar)%varname)//"'; &
                                   &type='"//trim(bvar_meta(ivar)%vartype)//"']"; return
  endselect
  bvar_data%var(ivar)%dat(:) = missingDouble
 end do ! (looping through model variables)
 end subroutine alloc_bvar


end module allocspace_module
