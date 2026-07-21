!+++++++++++++++++++++++++++++++++++++++++++++++++++
! Exact FGR grid code
!
! Features:
! 1. Automatic numerical convergence
! 2. Exact, HT and short-time FGR rates
! 3. Marcus rate
! 4. OpenMP parallel implementation
! 5. Parameter-grid partitioning for HPC clusters
!++++++++++++++++++++++++++++++++++++++++++++++++++++

program exact_fgr_grid
use omp_lib
implicit none

! Physical constants (SI units)
real(8), parameter :: pi   = 3.141592653589793d0
real(8), parameter :: hbar = 1.054571817d-34      ! J s
real(8), parameter :: kb   = 1.380649d-23         ! J/K
real(8), parameter :: c    = 2.99792458d10        ! cm/s
real(8), parameter :: hc   = 2.d0*pi*c*hbar       ! J cm
complex(8), parameter :: ii = (0.d0,1.d0)

! Progress variables
integer :: total_points
integer :: icount
real(8) :: job_start, job_end
real(8) :: point_start, point_end
integer :: nthreads

! System parameters
real(8), parameter :: Temp  = 300.d0
real(8), parameter :: Vc_cm = 30.d0

! Parameter grid
integer, parameter :: Nlambda = 20
integer, parameter :: Neta    = 20
integer, parameter :: NDG     = 20
real(8), parameter :: lambda_min = 10.d0
real(8), parameter :: lambda_max = 2000.d0
real(8), parameter :: eta_min = 10.d0
real(8), parameter :: eta_max = 1000.d0
real(8), parameter :: DG_min = 0.d0
real(8), parameter :: DG_max = 2000.d0
! Grid variables
integer :: ilam
integer :: ieta
integer :: iDG
real(8) :: lambda_cm
real(8) :: eta_cm
real(8) :: DG_cm
! DG_start and DG_end can be modified to divide the parameter grid among multiple independent jobs.
integer, parameter :: DG_start = 1
integer, parameter :: DG_end   = NDG
! SI variables
real(8) :: lambda
real(8) :: eta
real(8) :: DG
real(8) :: Vc
real(8) :: beta

! Numerical variables
real(8) :: tmax
real(8) :: dt
real(8) :: wmax
real(8) :: dw
integer :: Nt
integer :: Nw

! Working variables
integer :: jw
real(8) :: w
real(8) :: Jbath
real(8) :: lambda_int
real(8) :: lam_err

! Rate constants
real(8) :: kExact
real(8) :: kHT
real(8) :: kST
real(8) :: kMarcus
real(8) :: ratioExact
real(8) :: ratioHT
real(8) :: ratioST
real(8) :: lnExact
real(8) :: lnHT
real(8) :: lnST
real(8) :: k_old
real(8) :: err

!=============
! Output
!=============
open(unit=10,file='validity_map.dat',status='replace')
open(unit=20,file='convergence.log',status='replace')

write(10,'(A)') &
'# lambda(cm-1) eta(cm-1) DG(cm-1) '// &
'Exact(ps-1) HT(ps-1) ST(ps-1) Marcus(ps-1) '// &
'Exact/Marcus HT/Marcus ST/Marcus '// &
'ln(Exact/Marcus) ln(HT/Marcus) ln(ST/Marcus)'

write(20,'(A)') &
'# lambda eta DG tmax(ps) dt(fs) wmax(cm-1) dw(cm-1) Nt Nw lambda_error point_time(s)'

! Constants independent of the grid
Vc   = Vc_cm*hc
beta = 1.d0/(kb*Temp)

nthreads = omp_get_max_threads()

print *
print *, 'OpenMP threads = ', nthreads
print *

job_start = omp_get_wtime()

write(20,'(A,I4)') 'OpenMP threads = ', nthreads
write(20,*)

total_points = Nlambda*Neta*(DG_end-DG_start+1)
icount = 0

! Begin parameter grid
do ilam=1,Nlambda

   if (Nlambda == 1) then

      lambda_cm = lambda_min

   else

      lambda_cm = lambda_min * &
                  (lambda_max/lambda_min)** &
                  (dble(ilam-1)/dble(Nlambda-1))

   end if

   lambda = lambda_cm*hc
   
   do ieta = 1, Neta

      if (Neta == 1) then

         eta_cm = eta_max

      else

         eta_cm = eta_min * &
                  (eta_max/eta_min)** &
                  (dble(ieta-1)/dble(Neta-1))

      end if

      eta = eta_cm*2.d0*pi*c

      do iDG = DG_start, DG_end

         if (NDG == 1) then

            DG_cm = DG_min

         else

            DG_cm = DG_min + &
                    (DG_max-DG_min) * &
                    dble(iDG-1)/dble(NDG-1)

         end if

         DG = DG_cm*hc
         
         point_start = omp_get_wtime()

         write(*,'(/A,I6,A,I6)') &
         'POINT ', icount+1,' / ',total_points

         write(*,'(A,F8.2,A,F8.2,A,F8.2)') &
         'lambda=',lambda_cm, &
         '   eta=',eta_cm, &
         '   DG=',DG_cm
         
         icount = icount + 1

         ! Initialize numerical parameters
         dt   = 1.d-15                      ! 1 fs
         tmax = 10.d-12                     ! 10 ps
         dw   = 1.d0*2.d0*pi*c              ! 1 cm-1
         wmax = 8000.d0*2.d0*pi*c           ! 8000 cm-1
         Nt = int(tmax/dt) + 1
         Nw = int(wmax/dw) + 1
         
         !==================================================
         ! Recover reorganization energy
         !==================================================

         do

            lambda_int = 0.d0

            do jw = 1, Nw

               w = (jw-1)*dw

               if (jw == 1) w = 1.d-30

               Jbath = 2.d0*lambda*eta*w/(w*w + eta*eta)

               if (jw == 1 .or. jw == Nw) then

                  lambda_int = lambda_int + &
                               0.5d0*Jbath*dw/(pi*w)

               else

                  lambda_int = lambda_int + &
                               Jbath*dw/(pi*w)

               end if

            end do

            lam_err = abs(lambda_int-lambda)/lambda
            
            if (lam_err < 5.d-2) exit

            wmax = 2.d0*wmax
            dw   = dw/2.d0

            Nw = int(wmax/dw) + 1

         end do
         
         !==================================================
         ! Convergence with respect to tmax
         !==================================================
         
         Nt = int(tmax/dt) + 1

         call compute_fgr(kExact,kHT,kST)

         k_old = kExact

         tmax = 2.d0*tmax

         do

            Nt = int(tmax/dt) + 1

            call compute_fgr(kExact,kHT,kST)
    
            err = abs(kExact-k_old)/k_old
            
            if (err < 1.d-1) exit

            k_old = kExact

            tmax = 2.d0*tmax

         end do
         
         !==================================================
         ! Convergence with respect to dt
         !==================================================
         
         k_old = kExact

         dt = dt/2.d0

         do

            Nt = int(tmax/dt) + 1

            call compute_fgr(kExact,kHT,kST)

            err = abs(kExact-k_old)/k_old
            
            if (err < 1.d-1) exit

            k_old = kExact

            dt = dt/2.d0

         end do
         
         !==================================================
         ! Convergence with respect to wmax
         !==================================================
         
         k_old = kExact

         wmax = 2.d0*wmax

         do

            Nw = int(wmax/dw) + 1

            call compute_fgr(kExact,kHT,kST)

            err = abs(kExact-k_old)/k_old
            
            if (err < 1.d-1) exit

            k_old = kExact

            wmax = 2.d0*wmax

         end do

         !==================================================
         ! Convergence with respect to dw
         !==================================================
         
         k_old = kExact

         dw = dw/2.d0

         do

            Nw = int(wmax/dw) + 1

            call compute_fgr(kExact,kHT,kST)

            err = abs(kExact-k_old)/k_old
            
            if (err < 1.d-1) exit

            k_old = kExact

            dw = dw/2.d0

         end do    
              
         !==================================================
         ! Marcus rate
         !==================================================

         kMarcus = (Vc*Vc/hbar) * &
                   sqrt(pi/(lambda*kb*Temp)) * &
                   ( exp(-(DG+lambda)**2/(4.d0*lambda*kb*Temp)) + &
                     exp(-(lambda-DG)**2/(4.d0*lambda*kb*Temp)) )

         !==================================================
         ! Ratios
         !==================================================

         ratioExact = kExact/kMarcus
         ratioHT    = kHT/kMarcus
         ratioST    = kST/kMarcus

         lnExact = log(ratioExact)
         lnHT    = log(ratioHT)
         lnST    = log(ratioST)

         !==================================================
         ! Write one row
         !==================================================
         
         point_end = omp_get_wtime()
         
         write(10,'(13ES18.8)') &
         lambda_cm, &
         eta_cm, &
         DG_cm, &
         kExact*1.d-12, &
         kHT*1.d-12, &
         kST*1.d-12, &
         kMarcus*1.d-12, &
         ratioExact, &
         ratioHT, &
         ratioST, &
         lnExact, &
         lnHT, &
         lnST     
         flush(10)
         
         write(20,'(3F12.4,2F12.5,2F12.1,2I10,2ES14.5)') &
         lambda_cm, &
         eta_cm, &
         DG_cm, &
         tmax*1.d12, &
         dt*1.d15, &
         wmax/(2.d0*pi*c), &
         dw/(2.d0*pi*c), &
         Nt, &
         Nw, &
         lam_err, &
         point_end-point_start
         flush(20)
         
         point_end = omp_get_wtime()

         write(*,'(A)') &
         '------------------------------------------------------------'

         write(*,'(A,I6,A,I6)') &
         'DONE  ',icount,' / ',total_points

         write(*,'(A,F8.2,A,F8.2,A,F8.2)') &
         'lambda=',lambda_cm, &
         '   eta=',eta_cm, &
         '   DG=',DG_cm

         write(*,'(A,F6.1,A,A,F6.2,A,A,F8.0,A,A,F6.2,A)') &
         'tmax=',tmax*1.d12,' ps   ', &
         'dt=',dt*1.d15,' fs   ', &
         'wmax=',wmax/(2.d0*pi*c),' cm-1   ', &
         'dw=',dw/(2.d0*pi*c),' cm-1'
         
         write(*,'(A,I8,A,I8,A,ES10.2)') &
         'Nt=',Nt,'   Nw=',Nw,'   lambda err=',lam_err

         write(*,'(A,F8.2,A)') &
         'Point time = ',point_end-point_start,' s'

         print *      
              
      end do

   end do

end do

job_end = omp_get_wtime()

print *
print *, '======================================'
print '(A,F12.3)', 'TOTAL JOB TIME (s) = ', job_end-job_start
print '(A,F12.3)', 'TOTAL JOB TIME (h) = ', &
                  (job_end-job_start)/3600.d0
close(10)

write(20,*)
write(20,'(A)') &
'=============================================='

write(20,'(A,F12.3)') &
'TOTAL JOB TIME (s) = ',job_end-job_start

write(20,'(A,F12.3)') &
'TOTAL JOB TIME (h) = ',(job_end-job_start)/3600.d0

close(20)

contains

subroutine compute_fgr(kExact,kHT,kST)

   implicit none

   real(8), intent(out) :: kExact
   real(8), intent(out) :: kHT
   real(8), intent(out) :: kST
   integer :: i
   integer :: jw
   real(8) :: t
   real(8) :: w
   real(8) :: Jbath
   real(8) :: coth_exact
   real(8) :: coth_HT
   real(8) :: gR_exact
   real(8) :: gI_exact
   real(8) :: gR_HT
   real(8) :: gI_HT
   real(8) :: gR_ST
   real(8) :: gI_ST
   complex(8) :: Ff_exact
   complex(8) :: Fb_exact
   complex(8) :: Ff_HT
   complex(8) :: Fb_HT
   complex(8) :: Ff_ST
   complex(8) :: Fb_ST
   complex(8) :: sumExactF
   complex(8) :: sumExactB
   complex(8) :: sumHTF
   complex(8) :: sumHTB
   complex(8) :: sumSTF
   complex(8) :: sumSTB
   
   ! Initialise
   sumExactF = (0.d0,0.d0)
   sumExactB = (0.d0,0.d0)
   sumHTF = (0.d0,0.d0)
   sumHTB = (0.d0,0.d0)
   sumSTF = (0.d0,0.d0)
   sumSTB = (0.d0,0.d0)

   !===========
   ! Time loop
   !===========

   !$OMP PARALLEL DO DEFAULT(shared) &
   !$OMP PRIVATE(i,jw,t,w,Jbath,coth_exact,coth_HT, &
   !$OMP gR_exact,gI_exact,gR_HT,gI_HT,gR_ST,gI_ST, &
   !$OMP Ff_exact,Fb_exact,Ff_HT,Fb_HT,Ff_ST,Fb_ST) &
   !$OMP REDUCTION(+:sumExactF,sumExactB,sumHTF,sumHTB,sumSTF,sumSTB)

   do i = 1, Nt
 
      t = (i-1)*dt

      gR_exact = 0.d0
      gI_exact = 0.d0

      gR_HT = 0.d0
      gI_HT = 0.d0

      gR_ST = 0.d0
      gI_ST = 0.d0

      !------------------
      ! Frequency loop
      !------------------
      do jw = 1, Nw

         w = (jw-1)*dw

         if (jw == 1) w = 1.d-30

         Jbath = 2.d0*lambda*eta*w/(w*w + eta*eta)

         coth_exact = 1.d0/tanh(beta*hbar*w/2.d0)

         coth_HT = 2.d0/(beta*hbar*w)

         if (jw == 1 .or. jw == Nw) then

            gR_exact = gR_exact + 0.5d0*Jbath*(1.d0-cos(w*t))* &
                       coth_exact*dw/(pi*hbar*w*w)

            gI_exact = gI_exact - 0.5d0*Jbath*sin(w*t)* &
                       dw/(pi*hbar*w*w)

            gR_HT = gR_HT + 0.5d0*Jbath*(1.d0-cos(w*t))* &
                    coth_HT*dw/(pi*hbar*w*w)

            gI_HT = gI_HT - 0.5d0*Jbath*sin(w*t)* &
                    dw/(pi*hbar*w*w)

            !-----------------------------------
            ! Short-time expansion
            !-----------------------------------
            gR_ST = gR_ST + 0.5d0*Jbath*(0.5d0*w*w*t*t) * &
                    coth_exact * dw/(pi*hbar*w*w)            
                    
            gI_ST = gI_ST - 0.5d0*Jbath*(w*t)* &
                    dw/(pi*hbar*w*w)

         else

            gR_exact = gR_exact + Jbath*(1.d0-cos(w*t))* &
                       coth_exact*dw/(pi*hbar*w*w)

            gI_exact = gI_exact - Jbath*sin(w*t)* &
                       dw/(pi*hbar*w*w)

            gR_HT = gR_HT + Jbath*(1.d0-cos(w*t))* &
                    coth_HT*dw/(pi*hbar*w*w)

            gI_HT = gI_HT - Jbath*sin(w*t)* &
                    dw/(pi*hbar*w*w)

            gR_ST = gR_ST + Jbath*(0.5d0*w*w*t*t) * &
                    coth_exact * dw/(pi*hbar*w*w)

            gI_ST = gI_ST - Jbath*(w*t)* &
                    dw/(pi*hbar*w*w)

         end if

      end do
      
      !------------------------------------------
      ! Forward and backward propagators
      !------------------------------------------

      Ff_exact = exp(-gR_exact - ii*(gI_exact + DG*t/hbar))
      Fb_exact = exp(-gR_exact - ii*(gI_exact - DG*t/hbar))
      Ff_HT = exp(-gR_HT - ii*(gI_HT + DG*t/hbar))
      Fb_HT = exp(-gR_HT - ii*(gI_HT - DG*t/hbar))
      Ff_ST = exp(-gR_ST - ii*(gI_ST + DG*t/hbar))
      Fb_ST = exp(-gR_ST - ii*(gI_ST - DG*t/hbar))

      !------------------------------------------
      ! Trapezoidal integration over time
      !------------------------------------------

      if (i == 1 .or. i == Nt) then

         sumExactF = sumExactF + 0.5d0*Ff_exact*dt
         sumExactB = sumExactB + 0.5d0*Fb_exact*dt

         sumHTF = sumHTF + 0.5d0*Ff_HT*dt
         sumHTB = sumHTB + 0.5d0*Fb_HT*dt

         sumSTF = sumSTF + 0.5d0*Ff_ST*dt
         sumSTB = sumSTB + 0.5d0*Fb_ST*dt

      else

         sumExactF = sumExactF + Ff_exact*dt
         sumExactB = sumExactB + Fb_exact*dt

         sumHTF = sumHTF + Ff_HT*dt
         sumHTB = sumHTB + Fb_HT*dt

         sumSTF = sumSTF + Ff_ST*dt
         sumSTB = sumSTB + Fb_ST*dt

      end if

   end do
   !$OMP END PARALLEL DO
   !==================================================
   ! Final rate constants
   !==================================================

   kExact = 2*(Vc*Vc/(hbar*hbar)) * &
            (real(sumExactF) + real(sumExactB))

   kHT = 2*(Vc*Vc/(hbar*hbar)) * &
         (real(sumHTF) + real(sumHTB))

   kST = 2*(Vc*Vc/(hbar*hbar)) * &
         (real(sumSTF) + real(sumSTB))

end subroutine compute_fgr     

end program exact_fgr_grid 
