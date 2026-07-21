program exact_fgr
use omp_lib
implicit none

real(8) :: job_start, job_end
real(8) :: t1, t2
integer :: nthreads

! constants 
real(8), parameter :: pi   = 3.141592653589793d0
real(8), parameter :: hbar = 1.054571817d-34      ! J.s
real(8), parameter :: kb   = 1.380649d-23         ! J/K
real(8), parameter :: c    = 2.99792458d10        ! cm/s
real(8), parameter :: hc   = 2.d0*pi*c*hbar       ! J.cm

! Input parameters
real(8), parameter :: lambda_cm = 10.d0      ! cm-1
real(8), parameter :: DG_cm     = 316.d0       ! cm-1
real(8), parameter :: Vc_cm     = 30.d0      ! cm-1
real(8), parameter :: eta_cm    = 113.d0    ! cm-1
real(8), parameter :: Temp      = 300.d0     ! K

! Convergence variables
real(8) :: tmax                  
real(8) :: dt                  
real(8) :: wmax
real(8) :: dw
complex(8), parameter :: ii = (0.d0,1.d0)

! Variables
real(8) :: lambda, DG, Vc, eta, beta
real(8) :: w, Jbath
real(8) :: kFGR, kMarcus, ratio
integer :: jw, Nt, Nw
real(8) :: lambda_int
real(8) :: k_old, err, lam_err

! Convert to SI units
lambda = lambda_cm*hc
DG     = DG_cm*hc
Vc     = Vc_cm*hc
eta    = eta_cm*2.d0*pi*c
beta = 1.d0/(kb*Temp)

dt   = 1.d-15
wmax = 8000.d0*2.d0*pi*c
dw   = 1.d0*2.d0*pi*c

Nw = int(wmax/dw) + 1

nthreads = omp_get_max_threads()

print *
print *, 'OpenMP threads = ', nthreads
print *

job_start = omp_get_wtime()

! print input values
print *
print *, 'lambda (J) = ', lambda
print *, 'DG (J)     = ', DG
print *, 'Vc (J)     = ', Vc
print *, 'eta (s-1)  = ', eta
print *, 'beta (J-1) = ', beta
print *, 'Nw         = ', Nw
print *
!==================================
! Initial lambda recovery check
!==================================
do

   lambda_int = 0.d0

   do jw = 1, Nw

      w = (jw-1)*dw
      if (jw == 1) w = 1.d-30

      Jbath = 2.d0*lambda*eta*w/(w*w+eta*eta)

      if (jw == 1 .or. jw == Nw) then
         lambda_int = lambda_int + 0.5d0*Jbath*dw/(pi*w)
      else
         lambda_int = lambda_int + Jbath*dw/(pi*w)
      end if

   end do

   lam_err = abs(lambda_int-lambda)/lambda

   print *
   print *, 'Input lambda      = ', lambda
   print *, 'Integrated lambda = ', lambda_int
   print *, 'Relative error (%) = ',100.d0*lam_err

   if (lam_err < 0.01d0) exit

   print *, 'Improving initial frequency grid...'

   wmax = 2.d0*wmax
   dw   = dw/2.d0
   Nw   = int(wmax/dw)+1

end do

print * 

!========================
! Convergence in tmax
!========================
tmax = 10.d-12
Nt   = int(tmax/dt) + 1

t1 = omp_get_wtime()
call compute_fgr(kFGR)
t2 = omp_get_wtime()

print *, 'compute_fgr time (s) = ', t2-t1

print '(A,F6.1,A,ES12.4)', &
      'tmax = ', tmax*1.d12, ' ps   k = ', kFGR*1.d-12

k_old = kFGR
tmax  = 2.d0*tmax

do

   Nt = int(tmax/dt) + 1

   call compute_fgr(kFGR)

   err = abs(kFGR-k_old)/k_old

   print '(A,F6.1,A,ES12.4,A,F7.2,A)', &
         'tmax = ', tmax*1.d12, ' ps   k = ', &
         kFGR*1.d-12, ' ps-1   err = ',100.d0*err,' %'

   if (err < 5.d-2) exit

   k_old = kFGR
   tmax  = 2.d0*tmax

end do

!====================
! Convergence in dt
!====================
k_old = kFGR
dt = dt/2.d0

do

   Nt = int(tmax/dt) + 1
   
   call compute_fgr(kFGR)
   
   err = abs(kFGR-k_old)/k_old

   print '(A,F8.4,A,ES12.4,A,F7.2,A)', &
         'dt = ', dt*1.d15, ' fs   k = ', &
         kFGR*1.d-12, ' ps-1   err = ',100.d0*err,' %'

   if (err < 5.d-2) exit

   k_old = kFGR
   dt    = dt/2.d0

end do
  
!=====================
! Convergence in wmax
!=====================
k_old = kFGR
wmax  = 2.d0*wmax
!wmax  = 8000.d0*2.d0*pi*c

do

   Nw = int(wmax/dw) + 1

   call compute_fgr(kFGR)
   
   err = abs(kFGR-k_old)/k_old

   print '(A,F8.0,A,ES12.4,A,F7.2,A)', &
         'wmax = ', wmax/(2.d0*pi*c), ' cm-1   k = ', &
         kFGR*1.d-12, ' ps-1   err = ',100.d0*err,' %'

   if (err < 5.d-2) exit

   k_old = kFGR
   wmax  = 2.d0*wmax

end do
   
!====================
! Convergence in dw
!====================
k_old = kFGR
dw    = dw/2.d0
!dw    = 1.d0*2.d0*pi*c

do

   Nw = int(wmax/dw) + 1

   call compute_fgr(kFGR)
   
   err = abs(kFGR-k_old)/k_old

   print '(A,F8.3,A,ES12.4,A,F7.2,A)', &
         'dw = ', dw/(2.d0*pi*c), ' cm-1   k = ', &
         kFGR*1.d-12, ' ps-1   err = ',100.d0*err,' %'

   if (err < 5.d-2) exit

   k_old = kFGR
   dw    = dw/2.d0

end do      

! Converged numerical parameters
print *
print *, '======================================'
print *, 'Converged Numerical Parameters'
print *, '======================================'
print *, 'tmax (ps)    = ', tmax*1.d12
print *, 'dt (fs)      = ', dt*1.d15
print *, 'wmax (cm-1)  = ', wmax/(2.d0*pi*c)
print *, 'dw (cm-1)    = ', dw/(2.d0*pi*c)
print *, 'Nt           = ', Nt
print *, 'Nw           = ', Nw
print *, '======================================'
print *

!=================
! Marcus rate
!=================
kMarcus = (Vc*Vc/hbar)*sqrt(pi/(lambda*kb*Temp))* &
          ( exp(-(DG+lambda)**2/(4.d0*lambda*kb*Temp)) + &
            exp(-(lambda-DG)**2/(4.d0*lambda*kb*Temp)) )
ratio = kFGR/kMarcus

! Result
print *
print *, 'Exact FGR (s-1)    = ', kFGR
print *, 'Marcus    (s-1)    = ', kMarcus
print *, 'Ratio              = ', ratio
print *
print *, 'Exact FGR (ps-1)   = ', kFGR*1.d-12
print *, 'Marcus    (ps-1)   = ', kMarcus*1.d-12
print *

job_end = omp_get_wtime()

print *
print *, '======================================'
print '(A,F12.3)', 'TOTAL JOB TIME (s) = ', job_end-job_start
print '(A,F12.3)', 'TOTAL JOB TIME (h) = ', (job_end-job_start)/3600.d0

contains

subroutine compute_fgr(kFGR)

   real(8), intent(out) :: kFGR

   integer :: i, jw
   real(8) :: t, w
   real(8) :: gR, gI, Jbath, coth
   complex(8) :: Ff, Fb, sumF, sumB

   sumF = (0.d0,0.d0)
   sumB = (0.d0,0.d0)

   !$OMP PARALLEL DO DEFAULT(shared) &
   !$OMP PRIVATE(i,jw,t,w,gR,gI,Jbath,coth,Ff,Fb) &
   !$OMP REDUCTION(+:sumF,sumB)

   do i = 1, Nt

      t = (i-1)*dt

      gR = 0.d0
      gI = 0.d0

      do jw = 1, Nw

         w = (jw-1)*dw
         if (jw == 1) w = 1.d-30

         Jbath = 2.d0*lambda*eta*w/(w*w + eta*eta)
         coth  = 1.d0/tanh(beta*hbar*w/2.d0)

         if (jw == 1 .or. jw == Nw) then
            gR = gR + 0.5d0*Jbath*(1.d0-cos(w*t))*coth*dw/(pi*hbar*w*w)
            gI = gI - 0.5d0*Jbath*sin(w*t)*dw/(pi*hbar*w*w)
         else
            gR = gR + Jbath*(1.d0-cos(w*t))*coth*dw/(pi*hbar*w*w)
            gI = gI - Jbath*sin(w*t)*dw/(pi*hbar*w*w)
         end if

      end do
      
!      if (mod(i,10000)==0) then
!         print *, t*1.d12, gR
!      endif
      
      Ff = exp(-gR - ii*(gI + DG*t/hbar))
      Fb = exp(-gR - ii*(gI - DG*t/hbar))

      if (i == 1 .or. i == Nt) then
         sumF = sumF + 0.5d0*Ff*dt
         sumB = sumB + 0.5d0*Fb*dt
      else
         sumF = sumF + Ff*dt
         sumB = sumB + Fb*dt
      end if

   end do
   !$OMP END PARALLEL DO
   kFGR = 2*(Vc*Vc/(hbar*hbar))*(real(sumF)+real(sumB))

end subroutine compute_fgr

end program exact_fgr

