      module qmmm
      implicit none
      integer qmatoms
      integer, allocatable :: qmlist(:)
      real*8, allocatable :: qmforces(:,:)
      character*20 qmfif
      character*20 qmcof
      save
      end