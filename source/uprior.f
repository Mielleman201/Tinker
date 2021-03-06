c
c
c     ###################################################
c     ##  COPYRIGHT (C)  2011  by  Jay William Ponder  ##
c     ##              All Rights Reserved              ##
c     ###################################################
c
c     #############################################################
c     ##                                                         ##
c     ##  module uprior  --  previous values of induced dipoles  ##
c     ##                                                         ##
c     #############################################################
c
c
c     maxpred   maximum number of predictor induced dipoles to save
c
c     nualt     number of sets of prior induced dipoles in storage
c     maxualt   number of sets of induced dipoles needed for predictor
c     gear      coefficients for Gear predictor binomial method
c     aspc      coefficients for always stable predictor-corrector
c     bpred     coefficients for induced dipole predictor polynomial
c     bpredp    coefficients for predictor polynomial in energy field
c     bpreds    coefficients for predictor for PB/GK solvation
c     bpredps   coefficients for predictor in PB/GK energy field
c     udalt     prior values for induced dipoles at each site
c     upalt     prior values for induced dipoles in energy field
c     usalt     prior values for induced dipoles for PB/GK solvation
c     upsalt    prior values for induced dipoles in PB/GK energy field
c     use_pred  flag to control use of induced dipole prediction
c     polpred   type of predictor polynomial (GEAR, ASPC or LSQR)
c
c
      module uprior
      implicit none
      integer maxpred
      parameter (maxpred=17)
      integer nualt
      integer maxualt
      real*8 gear(maxpred)
      real*8 aspc(maxpred)
      real*8 bpred(maxpred)
      real*8 bpredp(maxpred)
      real*8 bpreds(maxpred)
      real*8 bpredps(maxpred)
      real*8, allocatable :: udalt(:,:,:)
      real*8, allocatable :: upalt(:,:,:)
      real*8, allocatable :: usalt(:,:,:)
      real*8, allocatable :: upsalt(:,:,:)
      logical use_pred
      character*4 polpred
      save
      end
