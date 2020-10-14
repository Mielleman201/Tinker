      subroutine qmmminit()
        use qmmm
        use mdstuf
        use usage
        use iounit
        use output
        use inform
        use potent
        use limits

        if (.not. allocated(qmlist))  allocate (qmlist(qmatoms))

c
c       Write relevent data tot the console
c      
        if (debug) then
            write (iout,10)  qmlist
   10       format (/,' QM Atom list: ',(*(3x, i5.1)))
        end if

        if (verbose) then
            write(iout, 20) qmfif
   20       format (/, 'Forces will be read from:', A)
            write(iout, 30) qmcof
   30       format (/, 'Coords will be written to:', A)
        end if

        write(iout, 70) use_charge,use_chgdpl,use_dipole,
     &             use_mpole,use_polar,use_chgtrn,use_rxnfld
   70   format (/, 'Flags', (*(3x, L1)))

           write(iout, 80) use_lights,use_clist
   80   format (/, 'Flags', (*(3x, L1)))
          

      end

      subroutine qmmmwritecoords()
        use atoms
        use qmmm
        use mdstuf
        use usage
        use units
        use iounit
        use output
        use inform
        integer i, j, k
        integer freeunit, icof
        logical exist
        character*240 fstr
        character*240 file

        icof = freeunit ()
        file = qmcof
        inquire (file=file,exist=exist)
        if (exist) then
           open (unit=icof,file=file,status='old')
           rewind (unit=icof)
        else
           open (unit=icof,file=file,status='new')
        end if

c        fstr = '('' Current Timestep :'')'
c        write (icof,fstr(1:23))
        fstr = '(3f26.16)'
        do i = 1, nuse 
            k = iuse(i)
            if ( ANY(qmlist==k) ) then
                write (icof,fstr(1:9))  x(k)/bohr,y(k)/bohr,z(k)/bohr
            end if
        end do

        close (unit=icof)
        return
      end

      subroutine qmmmreadforces()
        use atoms
        use qmmm
        use mdstuf
        use usage
        use units
        use iounit
        use output
        use inform
        integer i, j, k
        integer freeunit, ifif
        logical exist
        character*240 file
        character*240 record

        if (.not. allocated(qmforces))  allocate (qmforces(3, qmatoms))
        if (.not. allocated(qmerrors))  allocate (qmerrors(3, qmatoms))

        ifif = freeunit ()
        file = qmfif
        inquire (file=file,exist=exist)
        if (exist) then
           open (unit=ifif,file=file,status='old')
           rewind (unit=ifif)
        else
            write (iout,40)
   40       format (/,' readforces  --  Unable to Find the forces',
     &                 ' Restart File')
            call qmmmwritecoords()
            call fatal
        end if

        do i = 1, qmatoms
            qmerrors(1, i) = 0.0d0
            qmerrors(2, i) = 0.0d0
            qmerrors(3, i) = 0.0d0

            read (ifif,50)  record
   50       format (a240)
            read (record,*,err=60,end=60)  qmforces(1, i),
     &            qmforces(2, i),qmforces(3, i),
     &            qmerrors(1,i),qmerrors(2,i),qmerrors(3,i)
            qmforces(1, i) = qmforces(1, i)*hartree/bohr
            qmforces(2, i) = qmforces(2, i)*hartree/bohr
            qmforces(3, i) = qmforces(3, i)*hartree/bohr
            qmerrors(1, i) = qmerrors(1, i)*hartree*ekcal/bohr
            qmerrors(2, i) = qmerrors(2, i)*hartree*ekcal/bohr
            qmerrors(3, i) = qmerrors(3, i)*hartree*ekcal/bohr

        end do
   60   continue

        close (unit=ifif)
        return
      end


      subroutine qmcentreofmass()
        use atoms
        use atomid
        use qmmm

        integer i,j
        real*8 summass
        real*8, allocatable :: sumcoord(:)

        if (.not. allocated(qmmc)) allocate (qmmc(3))
        allocate(sumcoord(3))

        summass = 0d0
        sumcoord = 0d0

        do i = 1, size(qmlist)
            j = qmlist(i)
            summass = summass + mass(j)
            sumcoord(1) = x(j)*mass(j)
            sumcoord(2) = y(j)*mass(j)
            sumcoord(3) = z(j)*mass(j)
        end do

        qmmc(1) = sumcoord(1)/summass
        qmmc(2) = sumcoord(2)/summass
        qmmc(3) = sumcoord(3)/summass

        deallocate(sumcoord)
      end

      subroutine qmmmakegrid()
        use atoms
        use atomid
        use qmmm

        integer i,j,k
        real*8 maxx,maxy,maxz
        integer cubesx,cubesy,cubesz
        real*8 margin
        real*8 cubesize

        margin = 1.1d0
c       cubesize in Angstrom
c       later set this in keyfile
        cubesize = 0.001d0

        maxx = 0.0d0
        maxy = 0.0d0
        maxz = 0.0d0

        call qmcentreofmass()

        do i = 1, size(qmlist)
           j = qmlist(i)

           maxx = max(maxx, abs(x(j) - qmmc(1)))
           maxy = max(maxy, abs(y(j) - qmmc(2)))
           maxz = max(maxz, abs(z(j) - qmmc(3)))
        end do

        maxx = maxx * margin
        maxy = maxy * margin
        maxz = maxz * margin

        cubesx = maxx/cubesize
        cubesy = maxy/cubesize
        cubesz = maxz/cubesize

        if (.not. allocated(qmcoordgrid)) allocate (qmcoordgrid(3,
     &          2*cubesx+1,2*cubesy+1,2*cubesz+1))

        do i = 1, 2*cubesx+1 
          do j = 1, 2*cubesy+1
             do k = 1, 2*cubesz+1
                qmcoordgrid(1,i,j,k) = (i - cubesx - 1) * cubesize
                qmcoordgrid(2,i,j,k) = (j - cubesy - 1) * cubesize
                qmcoordgrid(3,i,j,k) = (k - cubesz - 1) * cubesize
             end do
          end do
        end do

      end