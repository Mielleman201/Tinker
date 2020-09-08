c
c
c     ##############################################################
c     ##  COPYRIGHT (C) 1998 by Rohit Pappu & Jay William Ponder  ##
c     ##                   All Rights Reserved                    ##
c     ##############################################################
c
c     ##############################################################
c     ##                                                          ##
c     ##  subroutine sdstep  --  Verlet stochastic dynamics step  ##
c     ##                                                          ##
c     ##############################################################
c
c
c     "sdstep" performs a single stochastic dynamics time step
c     via the velocity Verlet integration algorithm
c
c     literature references:
c
c     M. P. Allen, "Brownian Dynamics Simulation of a Chemical
c     Reaction in Solution", Molecular Physics, 40, 1073-1087 (1980)
c
c     F. Guarnieri and W. C. Still, "A Rapidly Convergent Simulation
c     Method: Mixed Monte Carlo/Stochastic Dynamics", Journal of
c     Computational Chemistry, 15, 1302-1310 (1994)
c
c
      subroutine sdstep (istep,dt)
      use atoms
      use atomid
      use freeze
      use mdstuf
      use moldyn
      use units
      use usage
      use virial
      implicit none
      integer i,j,k
      integer istep
      real*8 dt,term
      real*8 epot,etot
      real*8 eksum
      real*8 temp,pres
      real*8 vxx,vyy,vzz
      real*8 vyx,vzx,vzy
      real*8 ekin(3,3)
      real*8 stress(3,3)
      real*8, allocatable :: xold(:)
      real*8, allocatable :: yold(:)
      real*8, allocatable :: zold(:)
      real*8, allocatable :: pfric(:)
      real*8, allocatable :: vfric(:)
      real*8, allocatable :: afric(:)
      real*8, allocatable :: prand(:,:)
      real*8, allocatable :: vrand(:,:)
      real*8, allocatable :: derivs(:,:)
c
c
c     perform dynamic allocation of some local arrays
c
      allocate (xold(n))
      allocate (yold(n))
      allocate (zold(n))
      allocate (pfric(n))
      allocate (vfric(n))
      allocate (afric(n))
      allocate (prand(3,n))
      allocate (vrand(3,n))
      allocate (derivs(3,n))
c
c     get frictional and random terms for position and velocity
c
      call sdterm (istep,dt,pfric,vfric,afric,prand,vrand)
c
c     store the current atom positions, then find full-step
c     positions and half-step velocities via modified Verlet
c
      do i = 1, nuse
         k = iuse(i)
         xold(k) = x(k)
         yold(k) = y(k)
         zold(k) = z(k)
         x(k) = x(k) + v(1,k)*vfric(k) + a(1,k)*afric(k) + prand(1,k)
         y(k) = y(k) + v(2,k)*vfric(k) + a(2,k)*afric(k) + prand(2,k)
         z(k) = z(k) + v(3,k)*vfric(k) + a(3,k)*afric(k) + prand(3,k)
         do j = 1, 3
            v(j,k) = v(j,k)*pfric(k) + 0.5d0*a(j,k)*vfric(k)
         end do
      end do
c
c     get constraint-corrected positions and half-step velocities
c
      if (use_rattle)  call rattle (dt,xold,yold,zold)
c
c     get the potential energy and atomic forces
c
      call gradient (epot,derivs)
c
c     use Newton's second law to get the next accelerations;
c     find the full-step velocities using modified Verlet
c
      do i = 1, nuse
         k = iuse(i)
         do j = 1, 3
            a(j,k) = -ekcal * derivs(j,k) / mass(k)
            v(j,k) = v(j,k) + 0.5d0*a(j,k)*vfric(k) + vrand(j,k)
         end do
      end do
c
c     correct internal virial to account for frictional forces
c
      do i = 1, nuse
         k = iuse(i)
         term = vfric(k)/dt - 1.0d0
         vxx = term * x(k) * derivs(1,k)
         vyx = 0.5d0 * term * (y(k)*derivs(1,k)+x(k)*derivs(2,k))
         vzx = 0.5d0 * term * (z(k)*derivs(1,k)+x(k)*derivs(3,k))
         vyy = term * y(k) * derivs(2,k)
         vzy = 0.5d0 * term * (z(k)*derivs(2,k)+y(k)*derivs(3,k))
         vzz = term * z(k) * derivs(3,k)
         vir(1,1) = vir(1,1) + vxx
         vir(2,1) = vir(2,1) + vyx
         vir(3,1) = vir(3,1) + vzx
         vir(1,2) = vir(1,2) + vyx
         vir(2,2) = vir(2,2) + vyy
         vir(3,2) = vir(3,2) + vzy
         vir(1,3) = vir(1,3) + vzx
         vir(2,3) = vir(2,3) + vzy
         vir(3,3) = vir(3,3) + vzz
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (xold)
      deallocate (yold)
      deallocate (zold)
      deallocate (pfric)
      deallocate (vfric)
      deallocate (afric)
      deallocate (prand)
      deallocate (vrand)
      deallocate (derivs)
c
c     find the constraint-corrected full-step velocities
c
      if (use_rattle)  call rattle2 (dt)
c
c     compute and control the temperature and pressure
c
      call kinetic (eksum,ekin,temp)
      call pressure (dt,epot,ekin,temp,pres,stress)
c
c     total energy is sum of kinetic and potential energies
c
      etot = eksum + epot
c
c     compute statistics and save trajectory for this step
c
      call mdstat (istep,dt,etot,epot,eksum,temp,pres)
      call mdsave (istep,dt,epot,eksum)
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine sdterm  --  frictional and random SD terms  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sdterm" finds the frictional and random terms needed to
c     update positions and velocities during stochastic dynamics
c
c
      subroutine sdterm (istep,dt,pfric,vfric,afric,prand,vrand)
      use atoms
      use atomid
      use bath
      use stodyn
      use units
      use usage
      implicit none
      integer i,j,k
      integer istep
      real*8 dt,ktm
      real*8 gdt,egdt
      real*8 gdt2,gdt3
      real*8 gdt4,gdt5
      real*8 gdt6,gdt7
      real*8 gdt8,gdt9
      real*8 pterm,vterm
      real*8 pnorm,vnorm
      real*8 normal
      real*8 psig,vsig
      real*8 rho,rhoc
      real*8 pfric(*)
      real*8 vfric(*)
      real*8 afric(*)
      real*8 prand(3,*)
      real*8 vrand(3,*)
      logical first
      external normal
      save first
      data first  / .true. /
c
c
c     perform dynamic allocation of some global arrays
c
      if (first) then
         first = .false.
         if (.not. allocated(fgamma))  allocate (fgamma(n))
c
c     set the atomic friction coefficients to the global value
c
         do i = 1, n
            fgamma(i) = friction
         end do
      end if
c
c     set the value of the friction coefficient for each atom
c
      if (use_sdarea)  call sdarea (istep)
c
c     get the frictional and random terms for stochastic dynamics
c
      do i = 1, nuse
         k = iuse(i)
         gdt = fgamma(k) * dt / mass(k)
c
c     stochastic dynamics reduces to simple MD for zero friction
c
         if (gdt .le. 0.0d0) then
            pfric(k) = 1.0d0
            vfric(k) = dt
            afric(k) = 0.5d0 * dt * dt
            do j = 1, 3
               prand(j,k) = 0.0d0
               vrand(j,k) = 0.0d0
            end do
c
c     analytical expressions when friction coefficient is large
c
         else
            if (gdt .ge. 0.05d0) then
               egdt = exp(-gdt)
               pfric(k) = egdt
               vfric(k) = (1.0d0-egdt) * mass(k) / fgamma(k)
               afric(k) = (dt-vfric(k)) * mass(k) / fgamma(k)
               pterm = 2.0d0*gdt - 3.0d0 + (4.0d0-egdt)*egdt
               vterm = 1.0d0 - egdt**2
               rho = (1.0d0-egdt)**2 / sqrt(pterm*vterm)
c
c     use series expansions when friction coefficient is small
c
            else
               gdt2 = gdt * gdt
               gdt3 = gdt * gdt2
               gdt4 = gdt2 * gdt2
               gdt5 = gdt2 * gdt3
               gdt6 = gdt3 * gdt3
               gdt7 = gdt3 * gdt4
               gdt8 = gdt4 * gdt4
               gdt9 = gdt4 * gdt5
               afric(k) = (gdt2/2.0d0 - gdt3/6.0d0 + gdt4/24.0d0
     &                       - gdt5/120.0d0 + gdt6/720.0d0
     &                       - gdt7/5040.0d0 + gdt8/40320.0d0
     &                       - gdt9/362880.0d0) * mass(k)**2 / fgamma(k)**2
               vfric(k) = dt - fgamma(k)*afric(k) / mass(k)
               pfric(k) = 1.0d0 - fgamma(k)*vfric(k) / mass(k)
               pterm = 2.0d0*gdt3/3.0d0 - gdt4/2.0d0
     &                    + 7.0d0*gdt5/30.0d0 - gdt6/12.0d0
     &                    + 31.0d0*gdt7/1260.0d0 - gdt8/160.0d0
     &                    + 127.0d0*gdt9/90720.0d0
               vterm = 2.0d0*gdt - 2.0d0*gdt2 + 4.0d0*gdt3/3.0d0
     &                    - 2.0d0*gdt4/3.0d0 + 4.0d0*gdt5/15.0d0
     &                    - 4.0d0*gdt6/45.0d0 + 8.0d0*gdt7/315.0d0
     &                    - 2.0d0*gdt8/315.0d0 + 4.0d0*gdt9/2835.0d0
               rho = sqrt(3.0d0) * (0.5d0 - gdt/16.0d0
     &                    - 17.0d0*gdt2/1280.0d0
     &                    + 17.0d0*gdt3/6144.0d0
     &                    + 40967.0d0*gdt4/34406400.0d0
     &                    - 57203.0d0*gdt5/275251200.0d0
     &                    - 1429487.0d0*gdt6/13212057600.0d0
     &                    + 1877509.0d0*gdt7/105696460800.0d0)
            end if
c
c     compute random terms to thermostat the nonzero friction case
c
            ktm = boltzmann * kelvin / mass(k)
            psig = sqrt(ktm*pterm) * mass(k) / fgamma(k)
            vsig = sqrt(ktm*vterm)
            rhoc = sqrt(1.0d0 - rho*rho)
            do j = 1, 3
               pnorm = normal ()
               vnorm = normal ()
               prand(j,k) = psig * pnorm
               vrand(j,k) = vsig * (rho*pnorm+rhoc*vnorm)
            end do
         end if
      end do
      return
      end
c
c
c     #############################################################
c     ##                                                         ##
c     ##  subroutine sdarea  --  scale SD friction coefficients  ##
c     ##                                                         ##
c     #############################################################
c
c
c     "sdarea" optionally scales the atomic friction coefficient
c     of each atom based on its accessible surface area
c
c     literature reference:
c
c     S. Yun-Yi, W. Lu and W. F. van Gunsteren, "On the Approximation
c     of Solvent Effects on the Conformation and Dynamics of
c     Cyclosporin A by Stochastic Dynamics Simulation Techniques",
c     Molecular Simulation, 1, 369-383 (1988)
c
c
      subroutine sdarea (istep)
      use atoms
      use atomid
      use couple
      use kvdws
      use math
      use stodyn
      use usage
      implicit none
      integer i,k,istep
      integer resurf,modstep
      real*8 probe,ratio,area
      real*8, allocatable :: radius(:)
c
c
c     determine new friction coefficients every few SD steps
c
      resurf = 100
      modstep = mod(istep,resurf)
      if (modstep .ne. 1)  return
c
c     perform dynamic allocation of some local arrays
c
      allocate (radius(n))
c
c     set the atomic radii to estimates of sigma values
c
      probe = 0.0d0
      do i = 1, n
         radius(i) = rad(class(i)) / twosix
         if (radius(i) .ne. 0.0d0)  radius(i) = radius(i) + probe
      end do
c
c     scale atomic friction coefficients by accessible area
c
      do i = 1, nuse
         k = iuse(i)
         if (radius(k) .ne. 0.0d0) then
            call surfatom (k,area,radius)
            ratio = area / (4.0d0*pi*radius(k)**2)
            fgamma(k) = ratio * friction
         end if
      end do
c
c     monovalent atoms with zero radius get attached atom value
c
      do i = 1, nuse
         k = iuse(i)
         if (radius(k).eq.0.0d0 .and. n12(k).eq.1) then
            fgamma(k) = fgamma(i12(1,k))
         end if
      end do
c
c     perform deallocation of some local arrays
c
      deallocate (radius)
      return
      end
