program main_serial
implicit none
! 1D Transient head conduction equation. Mid plane symmetry, and convection BC on left.
! Explicit Scheme
  integer, parameter :: PREC = kind( 1.0d0 ) ! double precision
  integer, parameter :: xnodes=100       ! x grid points
  integer, parameter :: ynodes=xnodes    ! y grid points

  real(PREC) :: pi = 4 * atan( 1._PREC ) ! pi
  real(PREC) :: a=1.0                    ! x domain length
  real(PREC) :: b=1.0                    ! y domain length
  real(PREC) :: x                        ! x-location
  real(PREC) :: y                        ! y-location
  real(PREC) :: dx                       ! x grid spacing
  real(PREC) :: dt                       ! t grid spacing
  real(PREC) :: R                        ! Cell Reynolds Number
  real(PREC) :: CFL=0.3_PREC	                 ! CLF Number
  real       :: start, finish
  real(PREC) :: Tn(xnodes) = 1._PREC
  real(PREC) :: Tnp1(xnodes) = 1._PREC
  real(PREC) :: eps = epsilon(1._PREC)
  real(PREC) :: norm_np1=0.0
  real(PREC) :: norm_n=1.0
  real(PREC) :: resd
  real(PREC) :: h=1.0
  real(PREC) :: k=1.0
  real(PREC) :: T_inf=0.0
  integer :: un=11
  integer :: ierror
  integer :: i                           ! x spatial counter
  integer :: j                           ! y spatial counter
  integer :: max_iter=1001               ! maximum iteration count
  integer :: count=0
  open(unit=un,file="Tn_output.txt",status="replace",iostat=ierror)
!  if (ierror==1) then
!    write(*,*) "Error opening Tn output file."
!  else       
!    write(un,100) xnodes,',',ynodes
!  end if
100 format(I10,1A,I10)
  dx = a / (xnodes-1)
  dt =  dx*dx*CFL
  write(*,*) dt
  call cpu_time(start)
    
! Main loop
resd=1
! Initialization
do i=1,xnodes
  Tn(i)=1.0
  Tnp1(i)=1.0
end do
!do while (resd>1e-3.and.count<max_iter)
do while (count<max_iter)
  resd=0
  Tnp1(1)=Tn(1)+dt/dx/dx*(2*Tn(2)-2*Tn(1))
  Tnp1(xnodes)=(Tnp1(xnodes)+h*dx/k*T_inf)/(1+h*dx/k)

  count = count+1
  do i=2,xnodes-1
    Tnp1(i)=Tn(i)+dt/dx/dx*(Tn(i+1)-2*Tn(i)+Tn(i-1)) 
  end do  
  Tnp1(1)=Tn(1)+dt/dx/dx*(2*Tn(2)-2*Tn(1))
  Tnp1(xnodes)=(Tnp1(xnodes)+h*dx/k*T_inf)/(1+h*dx/k)
  do i=2,xnodes-1
    resd=sqrt(abs(Tn(i)**2-Tnp1(i)**2))+resd
  end do
  do i=1,xnodes
    Tn(i)=Tnp1(i)
  end do
  write(*,*) resd
end do

! write(*,*) resd
  call cpu_time(finish)
  print '("Time = ",f9.3," seconds.")',finish-start
 do i=1,xnodes
!   do j=1,ynodes
      write(un,101) Tnp1(i)
!   end do
 end do
101 format(f10.5,',')

close(un)
end program main_serial
