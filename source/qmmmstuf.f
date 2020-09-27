      subroutine qmmminit()
        use qmmm
        use mdstuf
        use usage
        use iounit
        use output
        use inform

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

      end

      subroutine qmmmwritecoords()
        use atoms
        use qmmm
        use mdstuf
        use usage
        use qmmm
        use units
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
        use qmmm
        use units
        integer i, j, k
        integer freeunit, ifif
        logical exist
        character*240 file
        character*240 record

        if (.not. allocated(qmforces))  allocate (qmforces(3, qmatoms))

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
            call fatal
        end if

        do i = 1, qmatoms
            read (ifif,50)  record
   50       format (a240)
            read (record,*,err=60,end=60)  qmforces(1, i),
     &            qmforces(2, i),qmforces(3, i)
        end do

   60   continue
        close (unit=ifif)
        return
      end