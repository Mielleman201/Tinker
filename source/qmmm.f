      module qmmm
c     qmatoms = the total number of qm atoms
c     qmlist = list with atomids of the qm atoms
c     qmforces = 3xqmatoms array with the forces for the qm atoms
c     qmmc = 3 array with the centre of mass of the qm molecule (for the grid)
c     qmfif = qm force input file
c     qmcof = qm coords output file
      implicit none
      integer qmatoms
      integer, allocatable :: qmlist(:)
      real*8, allocatable :: qmforces(:,:)
      real*8, allocatable :: qmmc(:)
      character*20 qmfif
      character*20 qmcof
      save
      end