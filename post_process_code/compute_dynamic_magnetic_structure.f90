! code written by Shaozhi Li, 09/17/2022
! This code is used to compute dynamical magnetic structure
! Use makefile to compile this code. Please change the path of fftw library in makefile
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Please change the read-in path on line 140 and 216
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


Module parameters
Implicit none
integer, parameter:: Nx=16   ! Lattice size
integer, parameter:: Ny=16   ! Lattice size
integer, parameter:: Nz=16   ! Lattice size
integer, parameter:: Nkpoint=(Nx+Ny+2*Nz)/2    ! Momentum number to be computed
integer, parameter:: Ntotal=40000              ! Vampire total time step
integer, parameter:: Nsub=1000                 ! Total time step in each measurement block
integer, parameter:: Nst=10000                 ! Thermalized time
integer, parameter:: Nremain=Ntotal-Nst        ! The remaining time step except for the thermalized time
integer, parameter:: Nset=Nremain/Nsub         ! Number of measurement blocks, here, I set Nset= 30
double precision:: time_unit=1.0d-15           ! Time unit between two time step
double precision:: frequency_unit

double precision:: ax=2.78d0
double precision:: ay=2.78d0
double precision:: az=2.78d0

double precision, dimension(0:Nkpoint,1:3):: qlist    ! Store all momentum along the path I defined
integer, dimension(1:2*Nx*Ny*Nz,1:4):: atoms          ! Store atoms' sublattice index and coordinate
End module

program main
use parameters
implicit none
double complex, dimension(1:Nremain,0:NKpoint,1:3,1:2):: spinq
double complex, dimension(1:Nset,1:Nsub,0:NKpoint,1:3,1:2):: spinqw
double precision, dimension(1:Nset,1:Nsub,0:NKpoint):: Sqw
double precision, dimension(1:Nsub,0:NKpoint):: aveSqw, smoothSqw
double precision, dimension(0:NKpoint):: avemagq
integer ti, i, j, k, l
double precision pi, arg
double precision eye
character file1*100, file2*100, file3*100

write(*,*) "Please change the read-in path on line 140 and 216!"

! Compute the frequency unit based on the time unit
frequency_unit=1/time_unit*4.1357d-15*1000.0d0


pi=2.0d0*asin(1.0d0)
eye=dcmplx(0.0d0,1.0d0)

! Read in atoms' coordinates
call build_atoms()


! Read in spin configuration m(r,t) after thermalized time
! Transfer m(r,t) to m(q,t)
do ti=1, Nremain
   i=ti+Nst-1
   call build_spin(i,spinq(ti,:,:,:))
enddo



! Fourier transfer m(q,t) to m(q,w)
k=0
do i=1,Nset
   k=(i-1)*Nsub
   do j=0, NKpoint 
      do l=1,3
         call fft_t2w(spinq(k+1:k+Nsub,j,l,1),spinqw(i,:,j,l,1),Nsub)
         call fft_t2w(spinq(k+1:k+nsub,j,l,2),spinqw(i,:,j,l,2),Nsub)
      enddo
   enddo
enddo




! Compute spin-spin correlation
! S_{a,b} (q,w) = m_{a}(q,w) * m_{b}(-q, -w) * e^{-iq*(r_a- r_b)}
! It is necessary to add the phase factor e^{-iq*(r_a- r_b)} due to the two sublattice in one unit cell
! m(-q,-w)=conjugate(m(q,w))
Sqw=0.0d0
do i=1, Nset
 do ti=1,Nsub
   do j=0, NKpoint
     arg=pi*(qlist(j,1)+qlist(j,2)+qlist(j,2))
     do l=1,3  
       Sqw(i,ti,j)=Sqw(i,ti,j)+real(exp(-eye*arg)*spinqw(i,ti,j,l,1)*conjg(spinqw(i,ti,j,l,2))) + &
                   real(exp(eye*arg)*spinqw(i,ti,j,l,2)*conjg(spinqw(i,ti,j,l,1)))+&
                   real(spinqw(i,ti,j,l,1)*conjg(spinqw(i,ti,j,l,1)))+&
                   real(spinqw(i,ti,j,l,2)*conjg(spinqw(i,ti,j,l,2)))
     enddo
   enddo
 enddo
enddo




! calculate the average value
aveSqw=0.0d0
do i=1,Nset
   aveSqw(:,:)=aveSqw(:,:)+Sqw(i,:,:)/dble(Nset)
enddo



! Output data
file2="spectral.dat"
open(file=file2,unit=58,action='write')
do ti=1,Nsub/2
   write(58,fmt='(f18.10,33f20.7)') (ti-1)/dble(Nsub)*frequency_unit, (aveSqw(ti,i),i=0,NKpoint)
enddo
close(58)


end program main

!------------------------------------------------
! Read in spin configuration and Fourier transfer to the momentum space
! Only specific momentum path is stored to reduce the storage
! The path of the read-in file is given in this subroutine.

subroutine build_spin(time,spinq)
use parameters
implicit none
integer time
double complex, dimension(0:NKpoint,1:3,1:2):: spinq
double precision, dimension(1:3):: spins
double precision, dimension(0:Nx-1,0:Ny-1,0:Nz-1,1:3):: spinA, spinB
double complex, dimension(0:Nx-1,0:Ny-1,0:Nz-1,1:3):: spinAq, spinBq
character file1*100, line1*100, ctmp*4
integer i, j, k, num, j1, j2

! Read in spin
write(unit=file1,fmt='("../T=0.01/spins-",i8.8,".data")') time
open(file=file1,unit=55,action='read')
read(55,*) num
spinA=0.0d0
do i=1,num
   read(55,*) spins(1), spins(2), spins(3)
   if(atoms(i,4).eq.1) then
       spinA(atoms(i,1),atoms(i,2),atoms(i,3),:)=spins
   else
      spinB(atoms(i,1),atoms(i,2),atoms(i,3),:)=spins
   endif
enddo

! Fourier transfer to the momentum space
do i=1,3
   call fft_r2q(spinA(:,:,:,i),spinAq(:,:,:,i),Nx,Ny,Nz)
   call fft_r2q(spinB(:,:,:,i),spinBq(:,:,:,i),Nx,Ny,Nz)
enddo


! Store data along the specific momentum path
! Here, I use the path along Gamma -- X --- M --- Z --- Gamma

! \Gamma -- X
spinq(0:Nx/2,:,1)=spinAq(0:Nx/2,0,0,:)
spinq(0:Nx/2,:,2)=spinBq(0:Nx/2,0,0,:)
do i=0, Nx/2
   qlist(i,1)=i/dble(Nx)
   qlist(i,2)=0.0d0
   qlist(i,3)=0.0d0
enddo

! X--M
spinq(Nx/2+1:Nx/2+Ny/2,:,1)=spinAq(Nx/2,1:Ny/2,0,:)
spinq(Nx/2+1:Nx/2+Ny/2,:,2)=spinBq(Nx/2,1:Ny/2,0,:)

do i=1, Ny/2
   qlist(Nx/2+i,1)=0.5d0
   qlist(Nx/2+i,2)=i/dble(Ny)
   qlist(Nx/2+i,3)=0.0d0
enddo

! M--Z
spinq(Nx/2+Ny/2+1:Nx/2+Ny/2+Nz/2,:,1)=spinAq(Nx/2,Ny/2,1:Nz/2,:)
spinq(Nx/2+Ny/2+1:Nx/2+Ny/2+Nz/2,:,2)=spinBq(Nx/2,Ny/2,1:Nz/2,:)

do i=1, Nz/2
   qlist(Nx/2+ny/2+i,1)=0.5d0
   qlist(Nx/2+ny/2+i,2)=0.5d0
   qlist(Nx/2+Ny/2+i,3)=i/dble(Nz)
enddo

! Z -- Gamma
k=1
do i=Nz/2-1,0,-1
   spinq(Nx/2+Ny/2+Nz/2+k,:,1)=spinAq(i,i,i,:)
   spinq(Nx/2+Ny/2+Nz/2+k,:,2)=spinBq(i,i,i,:)

   qlist(Nx/2+ny/2+Nz/2+k,1)=i/dble(Nx)
   qlist(Nx/2+ny/2+Nz/2+k,2)=i/dble(Ny)
   qlist(Nx/2+Ny/2+nz/2+k,3)=i/dble(Nz)
   k=k+1
enddo

close(55)
end subroutine

!-------------------------------------------
! Read in atoms coordinate
subroutine build_atoms()
use parameters
implicit none
integer i, j, k, num, j1, j2
character file1*100, line1*100, ctmp*4
double precision r1, r2, r3

file1="../T=0.01/atoms-coords.data"
open(file=file1,unit=55,action='read')
read(55,fmt=*) num
print*, "atom number =", num

atoms=0
do i=1, num
   read(55,fmt=*) j1,j2, r1, r2, r3
   
   atoms(i,1)=int(r1/ax)
   atoms(i,2)=int(r2/ay)
   atoms(i,3)=int(r3/az)

   if(abs(r1-atoms(i,1)*ax)>0.1.or.abs(r2-atoms(i,2))>0.1.or.abs(r3-atoms(i,3))>0.1) then
       atoms(i,4)=2
   else
       atoms(i,4)=1
   endif

enddo
close(55)
end subroutine

!---------------------------------------
! Use fifth to transfer real to momentum space
subroutine fft_r2q(f1,f2,Nx,Ny,Nz)
use, intrinsic:: iso_c_binding
implicit none
include 'fftw3.f03'
integer Nx, Ny, Nz
double precision, dimension(0:Nx-1,0:Ny-1,0:Nz-1):: f1
double complex, dimension(0:Nx-1,0:Ny-1,0:Nz-1):: f2
complex(C_double_complex), dimension(0:Nx-1,0:Ny-1,0:Nz-1):: in, out
type(C_PTR):: plan
integer iorb1, iorb2, kx, ky
double precision arg, pi
double complex eye

plan=fftw_plan_dft_3d(Nz,Ny,Nx,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
in=f1
call FFTW_EXECUTE_DFT(plan,in,out)
f2=out/dsqrt(1.0d0*Nx*Ny*Nz)

call fftw_destroy_plan(plan)
end subroutine

!-----------------------------------------------------
! Use fftw to transfer time to frequency space
subroutine fft_t2w(f1,f2,Nt)
use, intrinsic:: iso_c_binding
implicit none
include 'fftw3.f03'
integer Nt
double precision ave
double complex, dimension(0:Nt-1):: f1,f2
complex(C_double_complex), dimension(0:Nt-1):: in, out
type(C_PTR):: plan
integer i, iorb1, iorb2, kx, ky
double precision arg, pi
double complex eye
double precision delta

in=f1
plan=fftw_plan_dft_1d(Nt,in,out,FFTW_FORWARD,FFTW_ESTIMATE)
call FFTW_EXECUTE_DFT(plan,in,out)
f2=out/dsqrt(1.0d0*Nt)

call fftw_destroy_plan(plan)
end subroutine
