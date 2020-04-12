
!  DGV_kMaxwl_tools.f90
!
! Alex 09/12/16
!
! This modle contains useful procedures for working with K-Maxwell formulation 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!1

module DGV_kMaxwl_tools
use nrtype ! contains kind parameters (DP), (DP), (I4B) etc. 

implicit none

real (DP), parameter, private :: pi25DT = 3.141592653589793238462643d0
!
integer, parameter :: max_num_mxwls = 4 ! the maximum number of maxwellians that can be used in the model
! 
contains



!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMaxwl_SpatialOper
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMaxwl_RightSide0D(fmxwls,right_side)



use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R,mol_diam,L_inf,N_inf

use DGV_distributions_mod 
use DGV_dgvtools_mod

real (DP), dimension (:), intent(in) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), intent(out) :: right_side ! the time step for                       

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP):: grad_u, grad_v, grad_w, au,av,aw, g,h,ldist,diam,maxf 
                      ! scrap variables grad_u = to keep value of the gradient of the solution
                      ! au- to keep the approximation.
real (DP), dimension (:), allocatable :: max_approx ! array to store local extremuma (three 
                      ! DP to keep the velocity point and one for the local maximum
                      ! the length is 4K where K is the number of maxwellians that make up the 
                      ! solution. Format: velocity point, then max value5sz
real (DP), dimension (:,:), allocatable :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
real (DP), dimension (:), allocatable :: sample_u,sample_v,sample_w ! array to store local extremuma (three 
real (DP), dimension (:), allocatable :: coll_int ! array to keep values of the collision integral at sample velocity points                       
                      
integer :: loc_alloc_stat, nn, i,j, iter
real (DP) :: theta=5.0D-3 ! parameter determining the size of the gradient step. 
integer (I4B) :: max_iter = 1000 ! the maximum allowed number of step in gradient accent.
integer (I4B) :: sample_size_factor = 10! this coefficient determines how many samples is 
                        ! generated in the randomly generate velocity
                        ! the size of the sample is sample_size_factor*                                           
!!!!!!!!!!!!!! Debug
!!!
real (DP), dimension (:), allocatable :: coll_int1, right_side1 ! array to keep values of the collision integral at sample velocity points                       
!!!            END debug
!!!!!!!!!!!!!!!!!!!                                           
                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create velocity samples that are distributed according to the density given by the 
! solution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Step zero -- find maximum of the solution. We will solve this problem by running the algorithms of 
! gradient assent on the solution. We will make K runs of the algorithms, where K is the number of maxwellians
! that make up the solution. Then we select the best run and this will be approximation to max of the solution
!!!!!!!!!!!!!!!!!!!!!!
diam = abs(max(w_R-w_L,v_R-v_L,u_R-u_L)) ! rough measure of the size of the domain. 
nn = size(fmxwls,1)/5 ! this is the number of Maxwellians that make up the solution. 
allocate (max_approx(nn*4),sample_u(nn*5*sample_size_factor), &
 sample_v(nn*5*sample_size_factor),sample_w(nn*5*sample_size_factor),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D: Allocation error for (max_approx)"
 stop
endif 
! we repeat the gradient assent nn times where nn is the number of maxwellians
! that make up the solution.   
do i=1,nn
   ! set up the initial approximation to be equal to the bulk velocity of one of the 
   ! maxwellians from those that make up the solution
   au=fmxwls((i-1)*5+3)
   av=fmxwls((i-1)*5+4)
   aw=fmxwls((i-1)*5+5)
   ! we will now find the next local maximum  
   ldist = 1.0_DP ! this constant keeps track of the size of the last step. 
   iter = 0 ! reset iteration count 
   ! next we perform the gradient accent: evaluate the gradient and move in the direction of the 
   ! gradient. This is not 
   do while ((iter <= max_iter) .and. (ldist > diam*1.0D-5))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! now we will calculate the gradient of the solution
    grad_u=0
    grad_v=0
    grad_w=0 ! nullify the gradient before computing
    ! compute the gradient by taking derivative of each Maxwellian and adding it up
    do j=1,nn
     h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
     h = -h/fmxwls((j-1)*5+2)
     g = 2.0_DP*fmxwls((j-1)*5+1)/fmxwls((j-1)*5+2)*exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
     !!!!
     grad_u = grad_u - g*(au-fmxwls((j-1)*5+3))
     grad_v = grad_v - g*(av-fmxwls((j-1)*5+4))
     grad_w = grad_w - g*(aw-fmxwls((j-1)*5+5))
    end do
    ! The gradient is computed. Next we 
    ! caclulate the next approximation
    !
    h=sqrt(grad_u**2+grad_v**2+grad_w**2)
    if (h > 1.0_DP) then 
    !
     au=au+theta*grad_u/h
     av=av+theta*grad_v/h
     aw=aw+theta*grad_w/h
     ldist = theta ! keep track of the magnitude of the gradient step
    else 
    !
     au=au+theta*grad_u
     av=av+theta*grad_v
     aw=aw+theta*grad_w
     ldist = theta*sqrt(grad_u**2 + grad_v**2 + grad_w**2) ! keep track of the magnitude of the gradient step
    end if
    iter=iter+1 ! count the iterations
    !
   end do  
   ! we have computed the local minimum for maxwellian #i. We will now check if the min converged, or we exeeded the number of iterations: 
   if (iter > max_iter) then 
    print *, "kMaxwl_SpatialOper0D: max iteration exceeded in gradient assent. Maxwellian #", i
   end if
   ! record the last found local maximum
   max_approx((i-1)*4+1) = au
   max_approx((i-1)*4+2) = av
   max_approx((i-1)*4+3) = aw
   !!!! Evaluate the solution at the point of maximum !!!! 
   g = 0
   do j=1,nn
     g = g + maxwelveldist(fmxwls((j-1)*5+2),fmxwls((j-1)*5+3),fmxwls((j-1)*5+4),&
                   fmxwls((j-1)*5+5), fmxwls((j-1)*5+1),au,av,aw)
   end do 
   max_approx((i-1)*4+4) = g
   ! Recorded the solution at the point of local maximum.
end do 
! next we will find absolute maximum: 
maxf = max_approx(4)
do i=2,nn
 if (max_approx((i-1)*4+4)> maxf) then 
  maxf = max_approx((i-1)*4+4)
 end if  
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will generate a random sample of velocities that are disctributed according to the 
! velocity distribution dencity function.
! The approach to generate the density is acceptance-rejection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call random_seed
i=0
sample_u=0
sample_v=0
sample_w=0
do while (i < (nn*5*sample_size_factor))
 ! create a random uniformly distributed velocity point
 call random_number(au)
 au=(au-0.5)*(u_R-u_L)+(u_R+u_L)/2.0_DP
 call random_number(av)
 av=(av-0.5)*(v_R-v_L)+(v_R+v_L)/2.0_DP
 call random_number(aw)
 aw=(aw-0.5)*(w_R-w_L)+(w_R+w_L)/2.0_DP
 ! accept or reject the point using the distribution density given by the solution
 ! evaluate the solution at the new point,
 !!!! Evaluate the solution at the point of maximum !!!! 
 g = 0
 do j=1,nn
   g = g + maxwelveldist(fmxwls((j-1)*5+2),fmxwls((j-1)*5+3),fmxwls((j-1)*5+4),& 
              fmxwls((j-1)*5+5),fmxwls((j-1)*5 + 1),au,av,aw)
 end do 
 ! draw another random number
 call random_number(h)
 if (h < g/maxf) then ! acceptance-rejection 
  i=i+1
  sample_u(i)=au
  sample_v(i)=av
  sample_w(i)=aw
 end if 
 ! Recorded the solution at the point of local maximum.
 !!!!!!!!!!!
end do 
! sample velocity points are selected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the matrix that connects derivative of the parameters to the derivative 
! of the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (deriv_matrix(nn*5*sample_size_factor,nn*5), coll_int(nn*5*sample_size_factor), &
                stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D: Allocation error for (deriv_matrix)"
 stop
endif 
deriv_matrix = 0
do i=1,nn*5*sample_size_factor
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 do j=1,nn
   h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
   h = (-1)*h/fmxwls((j-1)*5+2)
   g = exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
   deriv_matrix(i,(j-1)*5+1)=g ! partial derivative of density
   deriv_matrix(i,(j-1)*5+2)=g*( (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2 + & 
       (aw-fmxwls((j-1)*5+5))**2/fmxwls((j-1)*5+2) - 1.5_DP )/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of temperature 
   deriv_matrix(i,(j-1)*5+3)=g*(au-fmxwls((j-1)*5+3))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{u} 
   deriv_matrix(i,(j-1)*5+4)=g*(av-fmxwls((j-1)*5+4))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{v} 
   deriv_matrix(i,(j-1)*5+5)=g*(aw-fmxwls((j-1)*5+5))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{w} 
 end do     
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The matrix hhas been evaluated. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we need to compute the vector of the solution derivative which is the value of the right hand side 
! at the sample points.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!! DEBUG
!right_side =4
!right_side(1)=3
!right_side(9)=2
!allocate (coll_int1(1:nn*5*sample_size_factor),right_side1(1:nn*5), stat=loc_alloc_stat)


!!! end Debug

!coll_int=matmul(deriv_matrix,right_side)

!call EvalCollArryS1Pts_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)
!coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

call EvalCollArryDecompS1PtsOMP_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)
coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the least squares problem to determine the derivaitve of the parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! solve linear least squares to determine values of the right side of the ODE 
 right_side = linleastsqrs(deriv_matrix,coll_int) 
 !right_side1 = linleastsqrs(deriv_matrix,coll_int1) 

 ! now right_side contain the values of the derivatives of the macroparameters in the 
 ! solution
  deallocate (deriv_matrix,coll_int,max_approx,sample_u,sample_v,sample_w)
 ! Debug
 !deallocate (coll_int1) 
 !end debug
 end subroutine kMaxwl_RightSide0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMaxwl_SpatialOper_MPI
!
!
! This is a modification of the above subroutine for MPI parallelization
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMaxwl_RightSide0D_MPI(fmxwls,right_side)



use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R,mol_diam,L_inf,N_inf

use DGV_distributions_mod 
use DGV_dgvtools_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

real (DP), dimension (:), pointer, intent(out) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), intent(out) :: right_side ! the time step for                       

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP):: grad_u, grad_v, grad_w, au,av,aw, g,h,ldist,diam,maxf 
                      ! scrap variables grad_u = to keep value of the gradient of the solution
                      ! au- to keep the approximation.
real (DP), dimension (:), allocatable :: max_approx ! array to store local extremuma (three 
                      ! DP to keep the velocity point and one for the local maximum
                      ! the length is 4K where K is the number of maxwellians that make up the 
                      ! solution. Format: velocity point, then max value5sz
real (DP), dimension (:,:), allocatable :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
real (DP), dimension (:), allocatable :: sample_u,sample_v,sample_w ! array to store local extremuma (three 
real (DP), dimension (:), allocatable :: coll_int ! array to keep values of the collision integral at sample velocity points                       
                      
integer :: loc_alloc_stat, nn, i,j, iter
real (DP):: theta=5.0D-3 ! parameter determining the size of the gradient step. 
integer (I4B) :: max_iter = 1000 ! the maximum allowed number of step in gradient accent.
integer (I4B) :: sample_size_factor = 10! this coefficient determines how many samples is 
                        ! generated in the randomly generate velocity
                        ! the size of the sample is sample_size_factor*                                           
!!!!!!!!!!!!!! Debug
!!!
real (DP), dimension (:), allocatable :: coll_int1, right_side1 ! array to keep values of the collision integral at sample velocity points                       
!!!            END debug
!!!!!!!!!!!!!!!!!!!                                           
                                
!!!!!!! MPI Variables
!
integer :: irank,ierr         ! scrap variable to keep the rank of the processor
integer :: mpi_dgv_action     ! scrap variable to keep the action code
integer, dimension(1:1) :: ibuff
real (DP), dimension (:), allocatable :: drbuffer  ! this is a scrap array to send double reals to other processors
!
!!!!!!!!!!!!!!!!!!!                                

!!!!!! Possible check !!!!!!!!!!!!!!!!!!!!!!!!! NOT critical to have 
!! call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
!!if (irank>0) then 
!! print *, "kMaxwl_RightSide0D_MPI: Error. Subroutine is started on processor with rank>0. Stop"
!! stop
!!end if
!!!!!!! End Possible check !!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The evaluation of the collision operator is distributed between the MPI processors.
! This subroutine should not be envoked on processor with Irank>0. Only on the processor with 
! irank ==0. To perform one step in time, the subroutine will prepare a sample, and then pass the sample 
! to all processors so that the processors could evaluate the collision opperator on the sample. 
! After that the results are combined on the master processor and the right side is computed.
                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create velocity samples that are distributed according to the density given by the 
! solution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Step zero -- find maximum of the solution. We will solve this problem by running the algorithms of 
! gradient assent on the solution. We will make K runs of the algorithms, where K is the number of maxwellians
! that make up the solution. Then we select the best run and this will be approximation to max of the solution
!!!!!!!!!!!!!!!!!!!!!!
diam = abs(max(w_R-w_L,v_R-v_L,u_R-u_L)) ! rough measure of the size of the domain. 
nn = size(fmxwls,1)/5 ! this is the number of Maxwellians that make up the solution. 
allocate (max_approx(nn*4),sample_u(nn*5*sample_size_factor), &
 sample_v(nn*5*sample_size_factor),sample_w(nn*5*sample_size_factor),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (max_approx)"
 stop
endif 
! we repeat the gradient assent nn times where nn is the number of maxwellians
! that make up the solution.   
do i=1,nn
   ! set up the initial approximation to be equal to the bulk velocity of one of the 
   ! maxwellians from those that make up the solution
   au=fmxwls((i-1)*5+3)
   av=fmxwls((i-1)*5+4)
   aw=fmxwls((i-1)*5+5)
   ! we will now find the next local maximum  
   ldist = 1.0_DP ! this constant keeps track of the size of the last step. 
   iter = 0 ! reset iteration count 
   ! next we perform the gradient accent: evaluate the gradient and move in the direction of the 
   ! gradient. This is not 
   do while ((iter <= max_iter) .and. (ldist > diam*1.0D-5))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! now we will calculate the gradient of the solution
    grad_u=0
    grad_v=0
    grad_w=0 ! nullify the gradient before computing
    ! compute the gradient by taking derivative of each Maxwellian and adding it up
    do j=1,nn
     h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
     h = -h/fmxwls((j-1)*5+2)
     g = 2.0_DP*fmxwls((j-1)*5+1)/fmxwls((j-1)*5+2)*exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
     !!!!
     grad_u = grad_u - g*(au-fmxwls((j-1)*5+3))
     grad_v = grad_v - g*(av-fmxwls((j-1)*5+4))
     grad_w = grad_w - g*(aw-fmxwls((j-1)*5+5))
    end do
    ! The gradient is computed. Next we 
    ! caclulate the next approximation
    !
    h=sqrt(grad_u**2+grad_v**2+grad_w**2)
    if (h > 1.0_DP) then 
    !
     au=au+theta*grad_u/h
     av=av+theta*grad_v/h
     aw=aw+theta*grad_w/h
     ldist = theta ! keep track of the magnitude of the gradient step
    else 
    !
     au=au+theta*grad_u
     av=av+theta*grad_v
     aw=aw+theta*grad_w
     ldist = theta*sqrt(grad_u**2 + grad_v**2 + grad_w**2) ! keep track of the magnitude of the gradient step
    end if
    iter=iter+1 ! count the iterations
    !
   end do  
   ! we have computed the local minimum for maxwellian #i. We will now check if the min converged, or we exeeded the number of iterations: 
   if (iter > max_iter) then 
    print *, "kMaxwl_SpatialOper0D: max iteration exceeded in gradient assent. Maxwellian #", i
   end if
   ! record the last found local maximum
   max_approx((i-1)*4+1) = au
   max_approx((i-1)*4+2) = av
   max_approx((i-1)*4+3) = aw
   !!!! Evaluate the solution at the point of maximum !!!! 
   g = 0
   do j=1,nn
     g = g + maxwelveldist(fmxwls((j-1)*5+2),fmxwls((j-1)*5+3),fmxwls((j-1)*5+4),&
                   fmxwls((j-1)*5+5), fmxwls((j-1)*5+1),au,av,aw)
   end do 
   max_approx((i-1)*4+4) = g
   ! Recorded the solution at the point of local maximum.
end do 
! next we will find absolute maximum: 
maxf = max_approx(4)
do i=2,nn
 if (max_approx((i-1)*4+4)> maxf) then 
  maxf = max_approx((i-1)*4+4)
 end if  
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will generate a random sample of velocities that are disctributed according to the 
! velocity distribution dencity function.
! The approach to generate the density is acceptance-rejection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call random_seed
i=0
sample_u=0
sample_v=0
sample_w=0
do while (i < (nn*5*sample_size_factor))
 ! create a random uniformly distributed velocity point
 call random_number(au)
 au=(au-0.5)*(u_R-u_L)+(u_R+u_L)/2.0_DP
 call random_number(av)
 av=(av-0.5)*(v_R-v_L)+(v_R+v_L)/2.0_DP
 call random_number(aw)
 aw=(aw-0.5)*(w_R-w_L)+(w_R+w_L)/2.0_DP
 ! accept or reject the point using the distribution density given by the solution
 ! evaluate the solution at the new point,
 !!!! Evaluate the solution at the point of maximum !!!! 
 g = 0
 do j=1,nn
   g = g + maxwelveldist(fmxwls((j-1)*5+2),fmxwls((j-1)*5+3),fmxwls((j-1)*5+4),& 
              fmxwls((j-1)*5+5),fmxwls((j-1)*5 + 1),au,av,aw)
 end do 
 ! draw another random number
 call random_number(h)
 if (h < g/maxf) then ! acceptance-rejection 
  i=i+1
  sample_u(i)=au
  sample_v(i)=av
  sample_w(i)=aw
 end if 
 ! Recorded the solution at the point of local maximum.
 !!!!!!!!!!!
end do 
! sample velocity points are selected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the matrix that connects derivative of the parameters to the derivative 
! of the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (deriv_matrix(nn*5*sample_size_factor,nn*5), coll_int(nn*5*sample_size_factor), &
                stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (deriv_matrix)"
 stop
endif 
deriv_matrix = 0
do i=1,nn*5*sample_size_factor
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 do j=1,nn
   h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
   h = (-1)*h/fmxwls((j-1)*5+2)
   g = exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
   deriv_matrix(i,(j-1)*5+1)=g ! partial derivative of density
   deriv_matrix(i,(j-1)*5+2)=g*( (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2 + & 
       (aw-fmxwls((j-1)*5+5))**2/fmxwls((j-1)*5+2) - 1.5_DP )/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of temperature 
   deriv_matrix(i,(j-1)*5+3)=g*(au-fmxwls((j-1)*5+3))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{u} 
   deriv_matrix(i,(j-1)*5+4)=g*(av-fmxwls((j-1)*5+4))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{v} 
   deriv_matrix(i,(j-1)*5+5)=g*(aw-fmxwls((j-1)*5+5))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{w} 
 end do     
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The matrix hhas been evaluated. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we need to compute the vector of the solution derivative which is the value of the right hand side 
! at the sample points.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!! DEBUG
!right_side =4
!right_side(1)=3
!right_side(9)=2
!allocate (coll_int1(1:nn*5*sample_size_factor),right_side1(1:nn*5), stat=loc_alloc_stat)

!!! end Debug

!coll_int=matmul(deriv_matrix,right_side)
!call EvalCollArryS1Pts_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)
!coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI EVALUATION OF THE COLLISION OPERATOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! send out the action code to slave processors
ibuff(1)=402 ! integer2, request evaluation of the collision operator using full decomposition model
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
 stop
end if    
! now that the code is broadcasted, the slave processors will need to have a copy of the 
! solution and the points where to evaluate the collision operator. The solution is contained 
! in the array maxwells(:) where macroparameters of all Maxwellians are given. 
! the points where the solution needs to be evaluated is given in the arrays sample_u, sample_v, sample_w
! all this information needs to be passed to the slave processors. 
!
! we will pass it in four steps: (1)  - send nuymber of entries of fmxwls, then (2) send maxwells
! similarly, (3) send number of entries in sample_u,sample_v,sample_w, then send (4) send the arrays 
nn = size(fmxwls,1) ! this is the length of the array containing macropatrameters of Maxwellians that make up the solution.
ibuff(1)=nn
! (1) send the size of the solution array
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(fmxwls,1) from proc 0 returned error", ierr
 stop
end if    
! (2) send the solution array
call mpi_bcast (fmxwls,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of fmxwls from proc 0 returned error", ierr
 stop
end if    
! (3) send the number of samples
nn = size(sample_u,1)*3 ! this should be equal to nn*sample_size_factor*3
ibuff(1) = nn ! let us hope the number is < 32,767
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(sample_u,1)*3 from proc 0 returned error", ierr
 stop
end if    
! (4) send the sample points, 
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer)"
 stop
endif
nn = size(sample_u,1)
drbuffer(1:nn) = sample_u(1:nn)
drbuffer(nn+1:2*nn) = sample_v(1:nn)
drbuffer(2*nn+1:3*nn) = sample_w(1:nn)
nn=nn*3 ! total number of records to be transmitted
call mpi_bcast (drbuffer,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of sample_u/_v/_w from proc 0 returned error", ierr
 stop
end if    
deallocate(drbuffer)
!!!!! The solution and the points to evaluate the collision operator has been transmitted. 
!!!!! Collision operator will be computed on othe transmitted points and send back to the 
!!!!! main processor./ 
nn = size(sample_u,1)
coll_int = 0
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer) 1"
 stop
endif
drbuffer =0 ! Make a empty array 
call mpi_allreduce(drbuffer,coll_int,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "fly trap" subroutine
if (ierr /= 0 ) then 
 print *,"GetRelaxRates1Donecell_DGV_MPI: (0) slave processor", irank,&
          "MPI allreduce of coll_int returned error", ierr
 stop
end if
deallocate(drbuffer)
! the values of the collision operator in the velocity sample points has been obtained. 
coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the least squares problem to determine the derivaitve of the parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! solve linear least squares to determine values of the right side of the ODE 
 right_side = linleastsqrs(deriv_matrix,coll_int) 
 !right_side1 = linleastsqrs(deriv_matrix,coll_int1) 

 ! now right_side contain the values of the derivatives of the macroparameters in the 
 ! solution
  deallocate (deriv_matrix,coll_int,max_approx,sample_u,sample_v,sample_w)
 ! Debug
 !deallocate (coll_int1) 
 !end debug
 end subroutine kMaxwl_RightSide0D_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMaxwl_RightSideSplMxwl0D_MPI
!
!
! This is a modification of the above MPI subroutine 
! In this subroutine, a different mechanism is used to create the velocity sample:
! the sample is generated distributed as the local maxwellian. 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMaxwl_RightSideSplMxwl0D_MPI(fmxwls,right_side)



use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R,mol_diam,L_inf,N_inf

use DGV_distributions_mod 
use DGV_dgvtools_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

real (DP), dimension (:), pointer, intent(out) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), intent(out) :: right_side ! the time step for                       

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP):: grad_u, grad_v, grad_w, au,av,aw, g,h,ldist,diam,maxf 
                      ! scrap variables grad_u = to keep value of the gradient of the solution
                      ! au- to keep the approximation.
real (DP), dimension (:,:), allocatable :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
real (DP), dimension (:), allocatable :: sample_u,sample_v,sample_w ! array to store local extremuma (three 
real (DP), dimension (:), allocatable :: coll_int ! array to keep values of the collision integral at sample velocity points   
                      
integer :: loc_alloc_stat, nn, i,j, iter,zz
real (DP):: theta=5.0D-3 ! parameter determining the size of the gradient step. 
integer (I4B) :: max_iter = 1000 ! the maximum allowed number of step in gradient accent.
integer (I4B) :: sample_size_factor = 10! this coefficient determines how many samples is 
                        ! generated in the randomly generate velocity
                        ! the size of the sample is sample_size_factor*                                           
real (DP) :: ndens,ubar,vbar,wbar,tempr ! scrap variables to keep macroparameters of the solution
                        
!!!!!!!!!!!!!! Debug
!!!
real (DP), dimension (:), allocatable :: coll_int1, right_side1 ! array to keep values of the collision integral at sample velocity points                       
!!!            END debug
!!!!!!!!!!!!!!!!!!!                                           
                                
!!!!!!! MPI Variables
!
integer :: irank,ierr         ! scrap variable to keep the rank of the processor
integer :: mpi_dgv_action     ! scrap variable to keep the action code
integer, dimension(1:1) :: ibuff
real (DP), dimension (:), allocatable :: drbuffer  ! this is a scrap array to send double reals to other processors
!
!!!!!!!!!!!!!!!!!!!                                

!!!!!! Possible check !!!!!!!!!!!!!!!!!!!!!!!!! NOT critical to have 
!! call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
!!if (irank>0) then 
!! print *, "kMaxwl_RightSide0D_MPI: Error. Subroutine is started on processor with rank>0. Stop"
!! stop
!!end if
!!!!!!! End Possible check !!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The evaluation of the collision operator is distributed between the MPI processors.
! This subroutine should not be envoked on processor with Irank>0. Only on the processor with 
! irank ==0. To perform one step in time, the subroutine will prepare a sample, and then pass the sample 
! to all processors so that the processors could evaluate the collision opperator on the sample. 
! After that the results are combined on the master processor and the right side is computed.
                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create velocity samples that are distributed according to the density given by the 
! solution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Step zero -- find maximum of the solution. We will solve this problem by running the algorithms of 
! gradient assent on the solution. We will make K runs of the algorithms, where K is the number of maxwellians
! that make up the solution. Then we select the best run and this will be approximation to max of the solution
!!!!!!!!!!!!!!!!!!!!!!
diam = abs(max(w_R-w_L,v_R-v_L,u_R-u_L)) ! rough measure of the size of the domain. 
nn = size(fmxwls,1)/5 ! this is the number of Maxwellians that make up the solution. 
allocate (sample_u(nn*5*sample_size_factor), &
 sample_v(nn*5*sample_size_factor),sample_w(nn*5*sample_size_factor),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (sample_u)"
 stop
endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we need to compute the macroparameters of the solution:
!!!!!!!!!!!
! Compute the macroparamters before the update 
ndens = 0 
ubar = 0
vbar = 0
wbar = 0
tempr = 0
do zz=1,nn
 ndens = ndens + fmxwls((zz-1)*5+1)
 ubar = ubar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+3)
 vbar = vbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+4)
 wbar = wbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+5)
 tempr = tempr + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+2)+ &
           2.0_DP/3.0_DP*fmxwls((zz-1)*5+1)*( fmxwls((zz-1)*5+3)**2 + &
          fmxwls((zz-1)*5+4)**2 + fmxwls((zz-1)*5+5)**2 )
end do
ubar=ubar/ndens
vbar=vbar/ndens
wbar=wbar/ndens
! caclulated total density and bulk velocity. Now we have an additiona term for the temperature:
tempr = tempr/ndens - 2.0_DP/3.0_DP*(ubar**2 + vbar**2 + wbar**2) 
maxf = ndens/sqrt(pi25DT*tempr)/(pi25DT*tempr) ! this is the largest value of the Maxwellian density
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will generate a random sample of velocities that are disctributed according to the 
! local Mawellian distribution density function.
! The approach to generate the density is acceptance-rejection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call random_seed
i=0
sample_u=0
sample_v=0
sample_w=0
do while (i < (nn*5*sample_size_factor))
 ! create a random uniformly distributed velocity point
 call random_number(au)
 au=(au-0.5)*(u_R-u_L)+(u_R+u_L)/2.0_DP
 call random_number(av)
 av=(av-0.5)*(v_R-v_L)+(v_R+v_L)/2.0_DP
 call random_number(aw)
 aw=(aw-0.5)*(w_R-w_L)+(w_R+w_L)/2.0_DP
 ! accept or reject the point using the distribution density given by the solution
 ! evaluate the solution at the new point,
 !!!! Evaluate the local Maxwellian at the point (au,av,aw)!!!! 
 g = maxwelveldist(tempr,ubar,vbar,wbar,ndens,au,aw,aw) 
 ! draw another random number between 0 and 1
 call random_number(h)
 if (h < g/maxf) then ! acceptance-rejection 
  i=i+1
  sample_u(i)=au
  sample_v(i)=av
  sample_w(i)=aw
 end if 
 ! Recorded the solution at the point of local maximum.
 !!!!!!!!!!!
end do 
! sample velocity points are selected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the matrix that connects derivative of the parameters to the derivative 
! of the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (deriv_matrix(nn*5*sample_size_factor,nn*5), coll_int(nn*5*sample_size_factor), &
                stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (deriv_matrix)"
 stop
endif 
deriv_matrix = 0
do i=1,nn*5*sample_size_factor
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 do j=1,nn
   h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
   h = (-1)*h/fmxwls((j-1)*5+2)
   g = exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
   deriv_matrix(i,(j-1)*5+1)=g ! partial derivative of density
   deriv_matrix(i,(j-1)*5+2)=g*( (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2 + & 
       (aw-fmxwls((j-1)*5+5))**2/fmxwls((j-1)*5+2) - 1.5_DP )/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of temperature 
   deriv_matrix(i,(j-1)*5+3)=g*(au-fmxwls((j-1)*5+3))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{u} 
   deriv_matrix(i,(j-1)*5+4)=g*(av-fmxwls((j-1)*5+4))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{v} 
   deriv_matrix(i,(j-1)*5+5)=g*(aw-fmxwls((j-1)*5+5))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{w} 
 end do     
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The matrix hhas been evaluated. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we need to compute the vector of the solution derivative which is the value of the right hand side 
! at the sample points.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!! DEBUG
!right_side =4
!right_side(1)=3
!right_side(9)=2
!allocate (coll_int1(1:nn*5*sample_size_factor),right_side1(1:nn*5), stat=loc_alloc_stat)

!!! end Debug

!coll_int=matmul(deriv_matrix,right_side)
!call EvalCollArryS1Pts_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)
!coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI EVALUATION OF THE COLLISION OPERATOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! send out the action code to slave processors
ibuff(1)=402 ! integer2, request evaluation of the collision operator using full decomposition model
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
 stop
end if    
! now that the code is broadcasted, the slave processors will need to have a copy of the 
! solution and the points where to evaluate the collision operator. The solution is contained 
! in the array maxwells(:) where macroparameters of all Maxwellians are given. 
! the points where the solution needs to be evaluated is given in the arrays sample_u, sample_v, sample_w
! all this information needs to be passed to the slave processors. 
!
! we will pass it in four steps: (1)  - send nuymber of entries of fmxwls, then (2) send maxwells
! similarly, (3) send number of entries in sample_u,sample_v,sample_w, then send (4) send the arrays 
nn = size(fmxwls,1) ! this is the length of the array containing macropatrameters of Maxwellians that make up the solution.
ibuff(1)=nn
! (1) send the size of the solution array
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(fmxwls,1) from proc 0 returned error", ierr
 stop
end if    
! (2) send the solution array
call mpi_bcast (fmxwls,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of fmxwls from proc 0 returned error", ierr
 stop
end if    
! (3) send the number of samples
nn = size(sample_u,1)*3 ! this should be equal to nn*sample_size_factor*3
ibuff(1) = nn ! let us hope the number is < 32,767
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(sample_u,1)*3 from proc 0 returned error", ierr
 stop
end if    
! (4) send the sample points, 
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer)"
 stop
endif
nn = size(sample_u,1)
drbuffer(1:nn) = sample_u(1:nn)
drbuffer(nn+1:2*nn) = sample_v(1:nn)
drbuffer(2*nn+1:3*nn) = sample_w(1:nn)
nn=nn*3 ! total number of records to be transmitted
call mpi_bcast (drbuffer,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of sample_u/_v/_w from proc 0 returned error", ierr
 stop
end if    
deallocate(drbuffer)
!!!!! The solution and the points to evaluate the collision operator has been transmitted. 
!!!!! Collision operator will be computed on othe transmitted points and send back to the 
!!!!! main processor./ 
nn = size(sample_u,1)
coll_int = 0
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer) 1"
 stop
endif
drbuffer =0 ! Make a empty array 
call mpi_allreduce(drbuffer,coll_int,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "fly trap" subroutine
if (ierr /= 0 ) then 
 print *,"GetRelaxRates1Donecell_DGV_MPI: (0) slave processor", irank,&
          "MPI allreduce of coll_int returned error", ierr
 stop
end if
deallocate(drbuffer)
! the values of the collision operator in the velocity sample points has been obtained. 
coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the least squares problem to determine the derivaitve of the parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! solve linear least squares to determine values of the right side of the ODE 
 right_side = linleastsqrs(deriv_matrix,coll_int) 
 !right_side1 = linleastsqrs(deriv_matrix,coll_int1) 

 ! now right_side contain the values of the derivatives of the macroparameters in the 
 ! solution
  deallocate (deriv_matrix,coll_int,sample_u,sample_v,sample_w)
 ! Debug
 !deallocate (coll_int1) 
 !end debug
 end subroutine kMaxwl_RightSideSplMxwl0D_MPI


!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMaxwl_RightSideSplCrclT0D_MPI
!
!
! This is a modification of the above MPI subroutine 
! In this subroutine, a different mechanism is used to create the velocity sample:
! the sample is generated uniformly distributed in a cricle 
! with the center of bulk velocity and the radius proportional to Temperature 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMaxwl_RightSideSplCrclT0D_MPI(fmxwls,right_side)



use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R,mol_diam,L_inf,N_inf,&
                  kMaxwell_sample_size_factor, kMaxwell_radius_temp_factor

use DGV_distributions_mod 
use DGV_dgvtools_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

real (DP), dimension (:), pointer, intent(out) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), intent(out) :: right_side ! the time step for                       

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP):: grad_u, grad_v, grad_w, au,av,aw, g,h,ldist,diam,maxf 
                      ! scrap variables grad_u = to keep value of the gradient of the solution
                      ! au- to keep the approximation.
real (DP), dimension (:,:), allocatable :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
real (DP), dimension (:), allocatable :: sample_u,sample_v,sample_w ! array to store local extremuma (three 
real (DP), dimension (:), allocatable :: coll_int ! array to keep values of the collision integral at sample velocity points   
                      
integer :: loc_alloc_stat, nn, i,j, iter,zz
real (DP):: theta=5.0D-3 ! parameter determining the size of the gradient step. 
integer (I4B) :: max_iter = 1000 ! the maximum allowed number of step in gradient accent.
real (DP) :: ndens,ubar,vbar,wbar,tempr,dist2,circ_rad2 ! scrap variables to keep macroparameters of the solution
                        
!!!!!!!!!!!!!! Debug
!!!
real (DP), dimension (:), allocatable :: coll_int1, right_side1 ! array to keep values of the collision integral at sample velocity points                       
!!!            END debug
!!!!!!!!!!!!!!!!!!!                                           
                                
!!!!!!! MPI Variables
!
integer :: irank,ierr         ! scrap variable to keep the rank of the processor
integer :: mpi_dgv_action     ! scrap variable to keep the action code
integer, dimension(1:1) :: ibuff
real (DP), dimension (:), allocatable :: drbuffer  ! this is a scrap array to send double reals to other processors
!
!!!!!!!!!!!!!!!!!!!                                

!!!!!! Possible check !!!!!!!!!!!!!!!!!!!!!!!!! NOT critical to have 
!! call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
!!if (irank>0) then 
!! print *, "kMaxwl_RightSide0D_MPI: Error. Subroutine is started on processor with rank>0. Stop"
!! stop
!!end if
!!!!!!! End Possible check !!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The evaluation of the collision operator is distributed between the MPI processors.
! This subroutine should not be envoked on processor with Irank>0. Only on the processor with 
! irank ==0. To perform one step in time, the subroutine will prepare a sample, and then pass the sample 
! to all processors so that the processors could evaluate the collision opperator on the sample. 
! After that the results are combined on the master processor and the right side is computed.
                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create velocity samples that are distributed according to the density given by the 
! solution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Step zero -- find maximum of the solution. We will solve this problem by running the algorithms of 
! gradient assent on the solution. We will make K runs of the algorithms, where K is the number of maxwellians
! that make up the solution. Then we select the best run and this will be approximation to max of the solution
!!!!!!!!!!!!!!!!!!!!!!
diam = abs(max(w_R-w_L,v_R-v_L,u_R-u_L)) ! rough measure of the size of the domain. 
nn = size(fmxwls,1)/5 ! this is the number of Maxwellians that make up the solution. 
allocate (sample_u(nn*5*kMaxwell_sample_size_factor), &
 sample_v(nn*5*kMaxwell_sample_size_factor),sample_w(nn*5*kMaxwell_sample_size_factor),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (sample_u)"
 stop
endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we need to compute the macroparameters of the solution:
!!!!!!!!!!!
! Compute the macroparamters before the update 
ndens = 0 
ubar = 0
vbar = 0
wbar = 0
tempr = 0
do zz=1,nn
 ndens = ndens + fmxwls((zz-1)*5+1)
 ubar = ubar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+3)
 vbar = vbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+4)
 wbar = wbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+5)
 tempr = tempr + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+2)+ &
           2.0_DP/3.0_DP*fmxwls((zz-1)*5+1)*( fmxwls((zz-1)*5+3)**2 + &
          fmxwls((zz-1)*5+4)**2 + fmxwls((zz-1)*5+5)**2 )
end do
ubar=ubar/ndens
vbar=vbar/ndens
wbar=wbar/ndens
! caclulated total density and bulk velocity. Now we have an additiona term for the temperature:
tempr = tempr/ndens - 2.0_DP/3.0_DP*(ubar**2 + vbar**2 + wbar**2) 
circ_rad2 = tempr*kMaxwell_radius_temp_factor*kMaxwell_radius_temp_factor ! radius of the acceptance circle squared 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will generate a random sample of velocities that are uniformly samples in a cirle with 
! center at the bulk velocity of the solution and radius sqrt(tempr)*radius_temp_factor 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call random_seed
i=0
sample_u=0
sample_v=0
sample_w=0
do while (i < (nn*5*kMaxwell_sample_size_factor))
 ! create a random uniformly distributed velocity point
 call random_number(au)
 au=(au-0.5)*(u_R-u_L)+(u_R+u_L)/2.0_DP
 call random_number(av)
 av=(av-0.5)*(v_R-v_L)+(v_R+v_L)/2.0_DP
 call random_number(aw)
 aw=(aw-0.5)*(w_R-w_L)+(w_R+w_L)/2.0_DP
 ! accept or reject the point by checking if the point is in the circle 
 ! with radius sqrt(tempr)*radius_temp_factor and center at (ubar,vbar,wbar) 
 dist2=(au-ubar)**2+(av-vbar)**2+(aw-wbar)**2 
 if (dist2 < circ_rad2) then ! acceptance-rejection 
  i=i+1
  sample_u(i)=au
  sample_v(i)=av
  sample_w(i)=aw
 end if 
 ! Recorded the solution at the point of local maximum.
 !!!!!!!!!!!
end do 
! sample velocity points are selected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the matrix that connects derivative of the parameters to the derivative 
! of the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (deriv_matrix(nn*5*kMaxwell_sample_size_factor,nn*5), coll_int(nn*5*kMaxwell_sample_size_factor), &
                stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (deriv_matrix)"
 stop
endif 
deriv_matrix = 0
do i=1,nn*5*kMaxwell_sample_size_factor
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 do j=1,nn
   h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
   h = (-1)*h/fmxwls((j-1)*5+2)
   g = exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
   deriv_matrix(i,(j-1)*5+1)=g ! partial derivative of density
   deriv_matrix(i,(j-1)*5+2)=g*( (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2 + & 
       (aw-fmxwls((j-1)*5+5))**2/fmxwls((j-1)*5+2) - 1.5_DP )/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of temperature 
   deriv_matrix(i,(j-1)*5+3)=g*(au-fmxwls((j-1)*5+3))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{u} 
   deriv_matrix(i,(j-1)*5+4)=g*(av-fmxwls((j-1)*5+4))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{v} 
   deriv_matrix(i,(j-1)*5+5)=g*(aw-fmxwls((j-1)*5+5))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative of \bar{w} 
 end do     
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The matrix hhas been evaluated. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we need to compute the vector of the solution derivative which is the value of the right hand side 
! at the sample points.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!! DEBUG
!right_side =4
!right_side(1)=3
!right_side(9)=2
!allocate (coll_int1(1:nn*5*kMaxwell_sample_size_factor),right_side1(1:nn*5), stat=loc_alloc_stat)

!!! end Debug

!coll_int=matmul(deriv_matrix,right_side)
!call EvalCollArryS1Pts_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)
!coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI EVALUATION OF THE COLLISION OPERATOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! send out the action code to slave processors
ibuff(1)=402 ! integer2, request evaluation of the collision operator using full decomposition model
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
 stop
end if    
! now that the code is broadcasted, the slave processors will need to have a copy of the 
! solution and the points where to evaluate the collision operator. The solution is contained 
! in the array maxwells(:) where macroparameters of all Maxwellians are given. 
! the points where the solution needs to be evaluated is given in the arrays sample_u, sample_v, sample_w
! all this information needs to be passed to the slave processors. 
!
! we will pass it in four steps: (1)  - send nuymber of entries of fmxwls, then (2) send maxwells
! similarly, (3) send number of entries in sample_u,sample_v,sample_w, then send (4) send the arrays 
nn = size(fmxwls,1) ! this is the length of the array containing macropatrameters of Maxwellians that make up the solution.
ibuff(1)=nn
! (1) send the size of the solution array
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(fmxwls,1) from proc 0 returned error", ierr
 stop
end if    
! (2) send the solution array
call mpi_bcast (fmxwls,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of fmxwls from proc 0 returned error", ierr
 stop
end if    
! (3) send the number of samples
nn = size(sample_u,1)*3 ! this should be equal to nn*kMaxwell_sample_size_factor*3
ibuff(1) = nn ! let us hope the number is < 32,767
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(sample_u,1)*3 from proc 0 returned error", ierr
 stop
end if    
! (4) send the sample points, 
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer)"
 stop
endif
nn = size(sample_u,1)
drbuffer(1:nn) = sample_u(1:nn)
drbuffer(nn+1:2*nn) = sample_v(1:nn)
drbuffer(2*nn+1:3*nn) = sample_w(1:nn)
nn=nn*3 ! total number of records to be transmitted
call mpi_bcast (drbuffer,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of sample_u/_v/_w from proc 0 returned error", ierr
 stop
end if    
deallocate(drbuffer)
!!!!! The solution and the points to evaluate the collision operator has been transmitted. 
!!!!! Collision operator will be computed on othe transmitted points and send back to the 
!!!!! main processor./ 
nn = size(sample_u,1)
coll_int = 0
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer) 1"
 stop
endif
drbuffer =0 ! Make a empty array 
call mpi_allreduce(drbuffer,coll_int,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "fly trap" subroutine
if (ierr /= 0 ) then 
 print *,"GetRelaxRates1Donecell_DGV_MPI: (0) slave processor", irank,&
          "MPI allreduce of coll_int returned error", ierr
 stop
end if
deallocate(drbuffer)
! the values of the collision operator in the velocity sample points has been obtained. 
coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the least squares problem to determine the derivaitve of the parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! solve linear least squares to determine values of the right side of the ODE 
 right_side = linleastsqrs(deriv_matrix,coll_int) 
 !right_side1 = linleastsqrs(deriv_matrix,coll_int1) 

 ! now right_side contain the values of the derivatives of the macroparameters in the 
 ! solution
  deallocate (deriv_matrix,coll_int,sample_u,sample_v,sample_w)
 ! Debug
 !deallocate (coll_int1) 
 !end debug
 end subroutine kMaxwl_RightSideSplCrclT0D_MPI

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMaxwl_RightSideCrcAddMxl0D_MPI
!
! THis is a modification of the above subroutine 
! a capability of adding an additional stream is added.
!
! This is a modification of the above MPI subroutine 
! In this subroutine, a different mechanism is used to create the velocity sample:
! the sample is generated uniformly distributed in a cricle 
! with the center of bulk velocity and the radius proportional to Temperature 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMaxwl_RightSideCrcAddMxl0D_MPI(fmxwls,right_side,num_mxwls)


use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R,mol_diam,L_inf,N_inf,&
                  kMaxwell_sample_size_factor, kMaxwell_radius_temp_factor

use DGV_distributions_mod 
use DGV_dgvtools_mod

!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     

real (DP), dimension (:), pointer, intent(out) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), pointer, intent(out) :: right_side ! the time step for     
integer (I4B), intent (out) :: num_mxwls  ! new number of maxwellians                 

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP):: grad_u, grad_v, grad_w, au,av,aw, g,h,ldist,diam,maxf 
                      ! scrap variables grad_u = to keep value of the gradient of the solution
                      ! au- to keep the approximation.
real (DP), dimension (:,:), allocatable :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
real (DP), dimension (:), allocatable :: sample_u,sample_v,sample_w, sample_p ! array to store sample velocities and their acceptance rejection rate 
real (DP), dimension (:), allocatable :: coll_int ! array to keep values of the collision integral at sample velocity points   
                      
integer :: loc_alloc_stat, nn, i,j, iter,zz
real (DP) :: ndens,ubar,vbar,wbar,tempr,dist2,circ_rad2 ! scrap variables to keep macroparameters of the solution
                        
!!!!!!!!!!!!!! Debug
!!!
real (DP), dimension (:), allocatable :: coll_int1, right_side1 ! array to keep values of the collision integral at sample velocity points                       
!!!            END debug
!!!!!!!!!!!!!!!!!!!                                           
                                
!!!!!!! MPI Variables
!
integer :: irank,ierr         ! scrap variable to keep the rank of the processor
integer :: mpi_dgv_action     ! scrap variable to keep the action code
integer, dimension(1:1) :: ibuff
real (DP), dimension (:), allocatable :: drbuffer  ! this is a scrap array to send double reals to other processors
!
!!!!!!!!!!!!!!!!!!!                                

!!!!!! Possible check !!!!!!!!!!!!!!!!!!!!!!!!! NOT critical to have 
!! call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...
!!if (irank>0) then 
!! print *, "kMaxwl_RightSide0D_MPI: Error. Subroutine is started on processor with rank>0. Stop"
!! stop
!!end if
!!!!!!! End Possible check !!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The evaluation of the collision operator is distributed between the MPI processors.
! This subroutine should not be envoked on processor with Irank>0. Only on the processor with 
! irank ==0. To perform one step in time, the subroutine will prepare a sample, and then pass the sample 
! to all processors so that the processors could evaluate the collision opperator on the sample. 
! After that the results are combined on the master processor and the right side is computed.
                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create velocity samples that are distributed according to the density given by the 
! solution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Step zero -- find maximum of the solution. We will solve this problem by running the algorithms of 
! gradient assent on the solution. We will make K runs of the algorithms, where K is the number of maxwellians
! that make up the solution. Then we select the best run and this will be approximation to max of the solution
!!!!!!!!!!!!!!!!!!!!!!
diam = abs(max(w_R-w_L,v_R-v_L,u_R-u_L)) ! rough measure of the size of the domain. 
nn = size(fmxwls,1)/5 ! this is the number of Maxwellians that make up the solution. should be the same as num_mxwls
allocate (sample_u(nn*5*kMaxwell_sample_size_factor), sample_p(nn*5*kMaxwell_sample_size_factor),&
 sample_v(nn*5*kMaxwell_sample_size_factor),sample_w(nn*5*kMaxwell_sample_size_factor),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_RightSideCrcAddMxl0D_MPI: Allocation error for (sample_u)"
 stop
endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we need to compute the macroparameters of the solution:
!!!!!!!!!!!
! Compute the macroparamters before the update 
ndens = 0 
ubar = 0
vbar = 0
wbar = 0
tempr = 0
do zz=1,nn
 ndens = ndens + fmxwls((zz-1)*5+1)
 ubar = ubar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+3)
 vbar = vbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+4)
 wbar = wbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+5)
 tempr = tempr + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+2)+ &
           2.0_DP/3.0_DP*fmxwls((zz-1)*5+1)*( fmxwls((zz-1)*5+3)**2 + &
          fmxwls((zz-1)*5+4)**2 + fmxwls((zz-1)*5+5)**2 )
end do
ubar=ubar/ndens
vbar=vbar/ndens
wbar=wbar/ndens
! caclulated total density and bulk velocity. Now we have an additiona term for the temperature:
tempr = tempr/ndens - 2.0_DP/3.0_DP*(ubar**2 + vbar**2 + wbar**2) 
circ_rad2 = tempr*kMaxwell_radius_temp_factor*kMaxwell_radius_temp_factor ! radius of the acceptance circle squared 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will generate a random sample of velocities that are uniformly samples in a cirle with 
! center at the bulk velocity of the solution and radius sqrt(tempr)*radius_temp_factor 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call random_seed
i=0
sample_u=0
sample_v=0
sample_w=0
do while (i < (nn*5*kMaxwell_sample_size_factor))
 ! create a random uniformly distributed velocity point
 call random_number(au)
 au=(au-0.5)*(u_R-u_L)+(u_R+u_L)/2.0_DP
 call random_number(av)
 av=(av-0.5)*(v_R-v_L)+(v_R+v_L)/2.0_DP
 call random_number(aw)
 aw=(aw-0.5)*(w_R-w_L)+(w_R+w_L)/2.0_DP
 ! accept or reject the point by checking if the point is in the circle 
 ! with radius sqrt(tempr)*radius_temp_factor and center at (ubar,vbar,wbar) 
 dist2=(au-ubar)**2+(av-vbar)**2+(aw-wbar)**2 
 if (dist2 < circ_rad2) then ! acceptance-rejection 
  i=i+1
  sample_p(i)=1.0_DP ! all points are accepted at rate 1
  sample_u(i)=au
  sample_v(i)=av
  sample_w(i)=aw
 end if 
 ! Recorded the solution at the point of local maximum.
 !!!!!!!!!!!
end do 
! sample velocity points are selected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the matrix that connects derivative of the parameters to the derivative 
! of the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (deriv_matrix(nn*5*kMaxwell_sample_size_factor,nn*5), coll_int(nn*5*kMaxwell_sample_size_factor), &
                stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (deriv_matrix)"
 stop
endif 
deriv_matrix = 0
do i=1,nn*5*kMaxwell_sample_size_factor
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 do j=1,nn
   h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
   h = (-1)*h/fmxwls((j-1)*5+2)
   g = exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
   deriv_matrix(i,(j-1)*5+1)=g ! partial derivative w.r. to density
   deriv_matrix(i,(j-1)*5+2)=g*( (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2 + & 
       (aw-fmxwls((j-1)*5+5))**2/fmxwls((j-1)*5+2) - 1.5_DP )/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to temperature 
   deriv_matrix(i,(j-1)*5+3)=g*(au-fmxwls((j-1)*5+3))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to \bar{u} 
   deriv_matrix(i,(j-1)*5+4)=g*(av-fmxwls((j-1)*5+4))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to \bar{v} 
   deriv_matrix(i,(j-1)*5+5)=g*(aw-fmxwls((j-1)*5+5))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to \bar{w} 
 end do     
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The derivative matrix hhas been evaluated. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! MPI EVALUATION OF THE COLLISION OPERATOR
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! send out the action code to slave processors
ibuff(1)=402 ! integer2, request evaluation of the collision operator using full decomposition model
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
 stop
end if    
! now that the code is broadcasted, the slave processors will need to have a copy of the 
! solution and the points where to evaluate the collision operator. The solution is contained 
! in the array maxwells(:) where macroparameters of all Maxwellians are given. 
! the points where the solution needs to be evaluated is given in the arrays sample_u, sample_v, sample_w
! all this information needs to be passed to the slave processors. 
!
! we will pass it in four steps: (1)  - send nuymber of entries of fmxwls, then (2) send maxwells
! similarly, (3) send number of entries in sample_u,sample_v,sample_w, then send (4) send the arrays 
nn = size(fmxwls,1) ! this is the length of the array containing macropatrameters of Maxwellians that make up the solution.
ibuff(1)=nn
! (1) send the size of the solution array
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(fmxwls,1) from proc 0 returned error", ierr
 stop
end if    
! (2) send the solution array
call mpi_bcast (fmxwls,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of fmxwls from proc 0 returned error", ierr
 stop
end if    
! (3) send the number of samples
nn = size(sample_u,1)*3 ! this should be equal to nn*sample_size_factor*3
ibuff(1) = nn ! let us hope the number is < 32,767
call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of size(sample_u,1)*3 from proc 0 returned error", ierr
 stop
end if    
! (4) send the sample points, 
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer)"
 stop
endif
nn = size(sample_u,1)
drbuffer(1:nn) = sample_u(1:nn)
drbuffer(nn+1:2*nn) = sample_v(1:nn)
drbuffer(2*nn+1:3*nn) = sample_w(1:nn)
nn=nn*3 ! total number of records to be transmitted
call mpi_bcast (drbuffer,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
if (ierr /= 0 ) then 
 print *,"kMaxwl_SpatialOper0D_MPI: master processor", &
          "MPI boradcast of sample_u/_v/_w from proc 0 returned error", ierr
 stop
end if    
deallocate(drbuffer)
!!!!! The solution and the points to evaluate the collision operator has been transmitted. 
!!!!! Collision operator will be computed on othe transmitted points and send back to the 
!!!!! main processor./ 
nn = size(sample_u,1)
coll_int = 0
allocate (drbuffer(nn), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer) 1"
 stop
endif
drbuffer =0 ! Make a empty array 
call mpi_allreduce(drbuffer,coll_int,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "fly trap" subroutine
if (ierr /= 0 ) then 
 print *,"GetRelaxRates1Donecell_DGV_MPI: (0) slave processor", irank,&
          "MPI allreduce of coll_int returned error", ierr
 stop
end if
deallocate(drbuffer)
! the values of the collision operator in the velocity sample points has been obtained. 
coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the least squares problem to determine the derivaitve of the parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! solve linear least squares to determine values of the right side of the ODE 
 right_side = linleastsqrs(deriv_matrix,coll_int) 
 ! now right_side contain the values of the derivatives of the macroparameters in the 
 !! Next we try to add an additional Maxwellian: This subroutine will attempt to add a new stream with zero density and 
 !! compute the updated (longer) right side. 
 call kMaxwl_AddMxwl0D_DGV(fmxwls,right_side,coll_int,deriv_matrix,num_mxwls,sample_u,sample_v,sample_w,sample_p)
 !! 
 deallocate (deriv_matrix,coll_int,sample_u,sample_v,sample_w,sample_p)
end subroutine kMaxwl_RightSideCrcAddMxl0D_MPI



!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMaxwl_SpatialOper
!
! This subroutine implements the right side in the kMaxwell model. This is an extension of the 
! above model. The new element in the extended model is a Maxwellian stream can be added to the existing 
! streams if the residual from the least squares problem "has too much mass". Macroparameters of the new  
!
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMaxwl_RightSideAddMxwls0D(fmxwls,right_side,num_mxwls,dt)

use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R,mol_diam,L_inf,N_inf

use DGV_distributions_mod 
use DGV_dgvtools_mod


real (DP), dimension (:), pointer, intent(out) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), pointer, intent(out) :: right_side ! the time  step for
integer (I4B), intent (out) :: num_mxwls    ! current number of the streams.                   
real (DP), intent (in) :: dt ! time step

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP):: grad_u, grad_v, grad_w,au,av,aw,g,h,ldist,diam,maxf 
                      ! scrap variables grad_u = to keep value of the gradient of the solution
                      ! au- to keep the approximation.
real (DP), dimension (:), allocatable :: max_approx ! array to store local extremuma (three 
                      ! DP to keep the velocity point and one for the local maximum
                      ! the length is 4K where K is the number of maxwellians that make up the 
                      ! solution. Format: velocity point, then max value5sz
real (DP), dimension (:,:), allocatable :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
real (DP), dimension (:), allocatable :: sample_u,sample_v,sample_w,f_sample ! arrays to sample points and avlues of the distribution 
real (DP), dimension (:), allocatable :: coll_int ! array to keep values of the collision integral at sample velocity points                       
                      
integer (I4B) :: loc_alloc_stat, nn, i,j, iter, ns
real (DP):: theta=5.0D-3 ! parameter determining the size of the gradient step. 
real (DP):: der_treshold = 1.0D+1 ! threshhold below which the new strem is not accepted. 
integer (I4B) :: max_iter = 1000 ! the maximum allowed number of step in gradient accent.
integer (I4B) :: sample_size_factor = 20! this coefficient determines how many samples is 
                        ! generated in the randomly generate velocity
                        ! the size of the sample is sample_size_factor*         
real (DP) :: dv, res_l1 ! scrap variable for volume element for stochastic integration
real (DP) :: dnpp1,Tpp1, upp1, vpp1, wpp1 ! scrap variables to hold macroparatmers of t5he new stream                                                           
logical :: accept_stream ! variable to keep decision on adding a new stream
real (DP), dimension (:), pointer :: scrap_arry_ptr ! scrap pointer to an allocatable array
!!!!!!!!!!!!!! Debug
!!!
real (DP), dimension (:), allocatable :: coll_int1, right_side1 ! array to keep values of the collision integral at sample velocity points                       
!!!            END debug
!!!!!!!!!!!!!!!!!!!                                           
                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create velocity samples that are distributed according to the density given by the 
! solution.
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Step zero -- find maximum of the solution. We will solve this problem by running the algorithms of 
! gradient assent on the solution. We will make K runs of the algorithms, where K is the number of maxwellians
! that make up the solution. Then we select the best run and this will be approximation to max of the solution
!!!!!!!!!!!!!!!!!!!!!!
diam = abs(max(w_R-w_L,v_R-v_L,u_R-u_L)) ! rough measure of the size of the domain. 
nn = size(fmxwls,1)/5 ! this is the number of Maxwellians that make up the solution. 
ns = nn*5*sample_size_factor ! total number of samples.
allocate (max_approx(nn*4),sample_u(ns),f_sample(ns),&
 sample_v(ns),sample_w(ns),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_RightSideAddMxwls0D: Allocation error for (max_approx, f_sample,sample_u/_v/_w)"
 stop
endif 
! we repeat the gradient assent nn times where nn is the number of maxwellians
! that make up the solution.   
do i=1,nn
   ! set up the initial approximation to be equal to the bulk velocity of one of the 
   ! maxwellians from those that make up the solution
   au=fmxwls((i-1)*5+3)
   av=fmxwls((i-1)*5+4)
   aw=fmxwls((i-1)*5+5)
   ! we will now find the next local maximum  
   ldist = 1.0_DP ! this constant keeps track of the size of the last step. 
   iter = 0 ! reset iteration count 
   ! next we perform the gradient accent: evaluate the gradient and move in the direction of the 
   ! gradient. This is not 
   do while ((iter <= max_iter) .and. (ldist > diam*1.0D-5))
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    ! now we will calculate the gradient of the solution
    grad_u=0
    grad_v=0
    grad_w=0 ! nullify the gradient before computing
    ! compute the gradient by taking derivative of each Maxwellian and adding it up
    do j=1,nn
     h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
     h = -h/fmxwls((j-1)*5+2)
     g = 2.0_DP*fmxwls((j-1)*5+1)/fmxwls((j-1)*5+2)*exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
     !!!!
     grad_u = grad_u - g*(au-fmxwls((j-1)*5+3))
     grad_v = grad_v - g*(av-fmxwls((j-1)*5+4))
     grad_w = grad_w - g*(aw-fmxwls((j-1)*5+5))
    end do
    ! The gradient is computed. Next we 
    ! caclulate the next approximation
    !
    h=sqrt(grad_u**2+grad_v**2+grad_w**2)
    if (h > 1.0_DP) then 
    !
     au=au+theta*grad_u/h
     av=av+theta*grad_v/h
     aw=aw+theta*grad_w/h
     ldist = theta ! keep track of the magnitude of the gradient step
    else 
    !
     au=au+theta*grad_u
     av=av+theta*grad_v
     aw=aw+theta*grad_w
     ldist = theta*sqrt(grad_u**2 + grad_v**2 + grad_w**2) ! keep track of the magnitude of the gradient step
    end if
    iter=iter+1 ! count the iterations
    !
   end do  
   ! we have computed the local minimum for maxwellian #i. We will now check if the min converged, or we exeeded the number of iterations: 
   if (iter > max_iter) then 
    print *, "kMaxwl_SpatialOper0D: max iteration exceeded in gradient assent. Maxwellian #", i
   end if
   ! record the last found local maximum
   max_approx((i-1)*4+1) = au
   max_approx((i-1)*4+2) = av
   max_approx((i-1)*4+3) = aw
   !!!! Evaluate the solution at the point of maximum !!!! 
   g = 0
   do j=1,nn
     g = g + maxwelveldist(fmxwls((j-1)*5+2),fmxwls((j-1)*5+3),fmxwls((j-1)*5+4),&
                   fmxwls((j-1)*5+5), fmxwls((j-1)*5+1),au,av,aw)
   end do 
   max_approx((i-1)*4+4) = g
   ! Recorded the solution at the point of local maximum.
end do 
! next we will find absolute maximum: 
maxf = max_approx(4)
do i=2,nn
 if (max_approx((i-1)*4+4)> maxf) then 
  maxf = max_approx((i-1)*4+4)
 end if  
end do
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will generate a random sample of velocities that are disctributed according to the 
! velocity distribution dencity function.
! The approach to generate the density is acceptance-rejection
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call random_seed
i=0
f_sample=0
sample_u=0
sample_v=0
sample_w=0
do while (i < ns)
 ! create a random uniformly distributed velocity point
 call random_number(au)
 au=(au-0.5)*(u_R-u_L)+(u_R+u_L)/2.0_DP
 call random_number(av)
 av=(av-0.5)*(v_R-v_L)+(v_R+v_L)/2.0_DP
 call random_number(aw)
 aw=(aw-0.5)*(w_R-w_L)+(w_R+w_L)/2.0_DP
 ! accept or reject the point using the distribution density given by the solution
 ! evaluate the solution at the new point,
 !!!! Evaluate the solution at the point of maximum !!!! 
 g = 0
 do j=1,nn
   g = g + maxwelveldist(fmxwls((j-1)*5+2),fmxwls((j-1)*5+3),fmxwls((j-1)*5+4),& 
              fmxwls((j-1)*5+5),fmxwls((j-1)*5 + 1),au,av,aw)
 end do 
 ! draw another random number
 call random_number(h)
 if (h < g/maxf) then ! acceptance-rejection 
  i=i+1
  f_sample(i)=g
  sample_u(i)=au
  sample_v(i)=av
  sample_w(i)=aw
 end if 
 ! Recorded the solution at the point of local maximum.
 !!!!!!!!!!!
end do 
! sample velocity points are selected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the matrix that connects derivative of the parameters to the derivative 
! of the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (deriv_matrix(ns,nn*5), coll_int(ns), &
                stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_SpatialOper0D: Allocation error for (deriv_matrix)"
 stop
endif 
deriv_matrix = 0
do i=1,ns
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 do j=1,nn
   h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
   h = (-1)*h/fmxwls((j-1)*5+2)
   g = exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
   deriv_matrix(i,(j-1)*5+1)=g ! partial derivative w.r. to density
   deriv_matrix(i,(j-1)*5+2)=g*( (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2 + & 
       (aw-fmxwls((j-1)*5+5))**2/fmxwls((j-1)*5+2) - 1.5_DP )/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative  w.r. to  temperature 
   deriv_matrix(i,(j-1)*5+3)=g*(au-fmxwls((j-1)*5+3))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative  w.r. to \bar{u} 
   deriv_matrix(i,(j-1)*5+4)=g*(av-fmxwls((j-1)*5+4))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative  w.r. to  \bar{v} 
   deriv_matrix(i,(j-1)*5+5)=g*(aw-fmxwls((j-1)*5+5))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative  w.r. to  \bar{w} 
 end do     
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The matrix has been evaluated. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we need to compute the vector of the solution derivative which is the value of the right hand side 
! at the sample points.
! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!! DEBUG
!right_side =4
!right_side(1)=3
!right_side(9)=2
!allocate (coll_int1(1:nn*5*sample_size_factor),right_side1(1:nn*5), stat=loc_alloc_stat)
!!!!!!!!!!!!!!!! End Debug

!coll_int=matmul(deriv_matrix,right_side)

call EvalCollArryDecompS1PtsOMP_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)
coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the least squares problem to determine the derivaitve of the parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! solve linear least squares to determine values of the right side of the ODE 
right_side = linleastsqrs(deriv_matrix,coll_int) 
! now right_side contain the values of the derivatives of the macroparameters in the 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! Now we will look at the residual and see if we need / can introduce an addition stream. 
coll_int = matmul(deriv_matrix,right_side) - coll_int ! now coll_int contains the residual 
res_l1 = sum(abs(coll_int))
!!! Attempt to find macroparameters of the new steam by using importance sampling
upp1 = u_L - .1 
vpp1 = v_L - .1 
wpp1 = w_L - .1 
Tpp1 = 0 
!!!
dnpp1 = sum(coll_int/f_sample)/Real(ns,DP) ! estimate derivative of the density of the new stream (dens=0 on this step)
!!!!!!!!!!!!!! DEBIG
           print *, "kMaxwl_RightSideAddMxwls0D: trying a stream: dn", dnpp1
!!!!!!!!!!!!! END DEBUG 
if (dnpp1 > der_treshold) then
 upp1 = sum(sample_u*coll_int/f_sample)/Real(ns,DP)/dnpp1
 vpp1 = sum(sample_v*coll_int/f_sample)/Real(ns,DP)/dnpp1
 wpp1 = sum(sample_w*coll_int/f_sample)/Real(ns,DP)/dnpp1
 Tpp1 = sum(((sample_u-upp1)**2+(sample_v-vpp1)**2+(sample_w-wpp1)**2)*coll_int/f_sample)&
           *2.0_DP/3.0_DP/Real(ns,DP)/dnpp1
           !!!!!!!!!!!!!! DEBIG
           print *, "kMaxwl_RightSideAddMxwls0D: trying a stream: T, u, v, w", Tpp1, upp1, vpp1, wpp1
           !!!!!!!!!!!!! END DEBUG 
end if
!!! A lot of sanity and other checks before accepting the new stream
accept_stream=.TRUE.
if ((dnpp1 < der_treshold) .or. (Tpp1 < 1.0D-5)) then !! Attention: hard limit coded! 
   accept_stream=.FALSE.
   print *, "kMaxwl_RightSideAddMxwls0D: dnpp1 or Tpp1 is small or negative - no new stream"
end if
if ((upp1 < u_L) .or. (upp1 > u_R)) then 
   accept_stream=.FALSE.
   print *, "kMaxwl_RightSideAddMxwls0D: upp1 is outside of the domain. no new stream"
end if
if ((vpp1 < v_L) .or. (vpp1 > v_R)) then 
   accept_stream=.FALSE.
   print *, "kMaxwl_RightSideAddMxwls0D: vpp1 is outside of the domain. no new stream"
end if
if ((wpp1 < w_L) .or. (wpp1 > w_R)) then 
   accept_stream=.FALSE.
   print *, "kMaxwl_RightSideAddMxwls0D: wpp1 is outside of the domain. no new stream"
end if
! check if maximum of the streams has been reached:
if (num_mxwls >= max_num_mxwls) then 
   accept_stream=.FALSE.
   print *, "kMaxwl_RightSideAddMxwls0D: num of streams at maximum. no new stream"
end if
!!! Last check --- did residual actualy became smaller
coll_int = coll_int - maxwelveldist(Tpp1,upp1,vpp1,wpp1,dnpp1,sample_u,sample_v,sample_w)
if (res_l1 > 8.0_DP*sum(abs(coll_int))) then 
   accept_stream=.FALSE.
   print *, "kMaxwl_RightSideAddMxwls0D: stream does not reduce residual. no new stream"
end if
!!! END OF checks to create new stream.
if (accept_stream) then  
  !!! Introduce an additional stream !!!!
  num_mxwls = num_mxwls+1
  ! create a memory block that contains p+1 streams
  allocate (scrap_arry_ptr((nn+1)*5), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then
   print *, "kMaxwl_RightSideAddMxwls0D: 1 Allocation error for (scrap_arry_ptr)"
   stop
  endif 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !Debug 
  ! print *, "nn", nn
  ! print *, "fmxwls", fmxwls
  ! end debug
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  ! copying macroparameters in the new/larger memory block 
  scrap_arry_ptr(1:nn*5) = fmxwls(:)
  scrap_arry_ptr(nn*5+1) = 0
  scrap_arry_ptr(nn*5+2) = Tpp1
  scrap_arry_ptr(nn*5+3) = upp1
  scrap_arry_ptr(nn*5+4) = vpp1
  scrap_arry_ptr(nn*5+5) = wpp1
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Debug 
  ! print *, "scrap", scrap_arry_ptr
  ! end debug
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  !!!!!!!!!!!!!!! 
  print *, "kMaxwl_RightSideAddMxwls0D: adding a stream: dn, T, u, v, w", dnpp1, Tpp1, upp1, vpp1, wpp1
  !!!!!!!!!!!!!!
  ! end writing to the larger memory block that contains p+1 streams
  ! next need to deallocate the old one
  deallocate (fmxwls)
  fmxwls => scrap_arry_ptr
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Debug 
  ! print *, "fmxwls", fmxwls
  ! end debug
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! make the scrap pointer point to nothing
  nullify (scrap_arry_ptr)
  !!!!!!!!!!!!!!! Updated fmxwls
  ! create a memory block that contains p+1 streams
  allocate (scrap_arry_ptr((nn+1)*5), stat=loc_alloc_stat)
  if (loc_alloc_stat >0) then
   print *, "kMaxwl_RightSideAddMxwls0D: 2 Allocation error for (scrap_arry_prt)"
   stop
  endif 
  ! copying macroparameters in the new/larger memory block 
  scrap_arry_ptr(1:nn*5) = right_side
  scrap_arry_ptr(nn*5+1) = dnpp1
  scrap_arry_ptr(nn*5+2) = 0
  scrap_arry_ptr(nn*5+3) = 0
  scrap_arry_ptr(nn*5+4) = 0
  scrap_arry_ptr(nn*5+5) = 0
  ! end writing to the larger memory block that contains p+1 streams
  ! next need to deallocate the old one
  deallocate (right_side)
  right_side => scrap_arry_ptr
  ! make the scrap pointer point to nothing
  nullify (scrap_arry_ptr)
  !!!!!!!!!!!!!!! Updated right_side
end if 
! now right_side contain the values of the derivatives of the macroparameters in the 
! solution
 deallocate (deriv_matrix,coll_int,max_approx,sample_u,sample_v,sample_w)
!!!!!!!!!!!!!!!!
!! Debug
! deallocate (coll_int1) 
!end debug
!!!!!!!!!!!!!!!!
end subroutine kMaxwl_RightSideAddMxwls0D

!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMaxwl_RightSideCrcAddMxl0D1
!
! THis is a modification of kMaxwl_RightSideCrcAddMxl0D_MPI to make it a serial version and aslo to use ES-BGK collision operator 
! Rather than full Boltzmann
!
! This is a modification of the above MPI subroutine 
! In this subroutine, a different mechanism is used to create the velocity sample:
! the sample is generated uniformly distributed in a cricle 
! with the center of bulk velocity and the radius proportional to Temperature 
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMaxwl_RightSideCrcAddMxl0Dv1(fmxwls,right_side,num_mxwls)


use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R,mol_diam,L_inf,N_inf,&
                  kBoltzmann,C_inf,gasR,gas_viscosity,&
                  kMaxwell_sample_size_factor, kMaxwell_radius_temp_factor

use DGV_distributions_mod 
use DGV_dgvtools_mod

real (DP), dimension (:), pointer, intent(out) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), pointer, intent(out) :: right_side ! the time step for     
integer (I4B), intent (out) :: num_mxwls  ! new number of maxwellians                 

!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP):: diam, au,av,aw,h,g ! scrap variables 
real (DP), dimension (:,:), allocatable :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
real (DP), dimension (:), allocatable :: sample_u,sample_v,sample_w,sample_p ! array to store sample velocities and their acceptance rejection rate 
real (DP), dimension (:), allocatable :: coll_int ! array to keep values of the collision integral at sample velocity points   
                      
integer :: loc_alloc_stat, nn, i,j, iter,zz
real (DP) :: ndens,ubar,vbar,wbar,tempr,circ_rad2,dist2 ! scrap variables to keep macroparameters of the solution
                        
!!!!!!!!!!!!!! Debug
!!!
real (DP), dimension (:), allocatable :: coll_int1, right_side1 ! array to keep values of the collision integral at sample velocity points                       
!!!            END debug
!!!!!!!!!!!!!!!!!!!                                           
                                
!!!!!!! MPI Variables
!
integer :: irank,ierr         ! scrap variable to keep the rank of the processor
integer :: mpi_dgv_action     ! scrap variable to keep the action code
integer, dimension(1:1) :: ibuff
real (DP), dimension (:), allocatable :: drbuffer  ! this is a scrap array to send double reals to other processors
!
!!!!!!!!!!!!!!!!!!!                                

                                           
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Create velocity samples that are distributed in a circle with the center at the bulk velocity 
! and radius proportional to square root of temperature
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! create arrays to store velocity samples
!!!!!!!!!!!!!!!!!!!!!!
diam = abs(max(w_R-w_L,v_R-v_L,u_R-u_L)) ! rough measure of the size of the domain. 
nn = size(fmxwls,1)/5 ! this is the number of Maxwellians that make up the solution. should be the same as num_mxwls
allocate (sample_u(nn*5*kMaxwell_sample_size_factor), sample_p(nn*5*kMaxwell_sample_size_factor), &
 sample_v(nn*5*kMaxwell_sample_size_factor),sample_w(nn*5*kMaxwell_sample_size_factor),stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_RightSideCrcAddMxl0Dv1: Allocation error for (sample_u)"
 stop
endif 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First, we need to compute the macroparameters of the solution:
!!!!!!!!!!!
! Compute the macroparamters before the update 
ndens = 0 
ubar = 0
vbar = 0
wbar = 0
tempr = 0
do zz=1,nn
 ndens = ndens + fmxwls((zz-1)*5+1)
 ubar = ubar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+3)
 vbar = vbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+4)
 wbar = wbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+5)
 tempr = tempr + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+2)+ &
           2.0_DP/3.0_DP*fmxwls((zz-1)*5+1)*( fmxwls((zz-1)*5+3)**2 + &
          fmxwls((zz-1)*5+4)**2 + fmxwls((zz-1)*5+5)**2 )
end do
ubar=ubar/ndens
vbar=vbar/ndens
wbar=wbar/ndens
! caclulated total density and bulk velocity. Now we have an additiona term for the temperature:
tempr = tempr/ndens - 2.0_DP/3.0_DP*(ubar**2 + vbar**2 + wbar**2) 
circ_rad2 = tempr*kMaxwell_radius_temp_factor*kMaxwell_radius_temp_factor ! radius of the acceptance circle squared 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next we will generate a random sample of velocities that are uniformly samples in a cirle with 
! center at the bulk velocity of the solution and radius sqrt(tempr)*radius_temp_factor 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
call random_seed
i=0
sample_u=0
sample_v=0
sample_w=0
sample_p=0
do while (i < (nn*5*kMaxwell_sample_size_factor))
 ! create a random uniformly distributed velocity point
 call random_number(au)
 au=(au-0.5)*(u_R-u_L)+(u_R+u_L)/2.0_DP
 call random_number(av)
 av=(av-0.5)*(v_R-v_L)+(v_R+v_L)/2.0_DP
 call random_number(aw)
 aw=(aw-0.5)*(w_R-w_L)+(w_R+w_L)/2.0_DP
 ! accept or reject the point by checking if the point is in the circle 
 ! with radius sqrt(tempr)*radius_temp_factor and center at (ubar,vbar,wbar) 
 dist2=(au-ubar)**2+(av-vbar)**2+(aw-wbar)**2 
 if (dist2 < circ_rad2) then ! acceptance-rejection 
  i=i+1
  sample_p(i)=1.0_DP
  sample_u(i)=au
  sample_v(i)=av
  sample_w(i)=aw
 end if 
 ! Recorded the solution at the point of local maximum.
 !!!!!!!!!!!
end do 
! sample velocity points are selected.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluate the matrix that connects derivative of the parameters to the derivative 
! of the solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
allocate (deriv_matrix(nn*5*kMaxwell_sample_size_factor,nn*5), coll_int(nn*5*kMaxwell_sample_size_factor), &
                stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_RightSideCrcAddMxl0Dv1: Allocation error for (deriv_matrix)"
 stop
endif 
deriv_matrix = 0
do i=1,nn*5*kMaxwell_sample_size_factor
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 do j=1,nn
   h = (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2+(aw-fmxwls((j-1)*5+5))**2  
   h = (-1)*h/fmxwls((j-1)*5+2)
   g = exp(h)/(( pi25DT*fmxwls((j-1)*5+2) )**1.5_DP)
   deriv_matrix(i,(j-1)*5+1)=g ! partial derivative w.r. to density
   deriv_matrix(i,(j-1)*5+2)=g*( (au-fmxwls((j-1)*5+3))**2+(av-fmxwls((j-1)*5+4))**2 + & 
       (aw-fmxwls((j-1)*5+5))**2/fmxwls((j-1)*5+2) - 1.5_DP )/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to temperature 
   deriv_matrix(i,(j-1)*5+3)=g*(au-fmxwls((j-1)*5+3))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to \bar{u} 
   deriv_matrix(i,(j-1)*5+4)=g*(av-fmxwls((j-1)*5+4))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to \bar{v} 
   deriv_matrix(i,(j-1)*5+5)=g*(aw-fmxwls((j-1)*5+5))*2.0_DP/fmxwls((j-1)*5+2)*fmxwls((j-1)*5+1) ! partial derivative w.r. to \bar{w} 
 end do     
end do 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! The derivative matrix hhas been evaluated. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Evaluation of the collision operator in the velocity sample points 

! uncomment next two lines if what to use Boltzmann collision operator
!! call EvalCollArryDecompS1PtsOMP_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)
!! coll_int = coll_int*(mol_diam/L_inf)**2*N_inf

! uncomment the next three lines is want to use ES operator in the kMaxwell model
call kMxl_EvalColESBGKPtsOMP_DGV(fmxwls,coll_int,sample_u,sample_v,sample_w) ! UNCOMMENT IF USING THE ES MODEL
coll_int = coll_int*(kBoltzmann*N_inf)*C_inf/(2*gasR)/L_inf**2/gas_viscosity

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Solve the least squares problem to determine the derivaitve of the parameters.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! solve linear least squares to determine values of the right side of the ODE 
 right_side = linleastsqrs(deriv_matrix,coll_int) 
 ! now right_side contain the values of the derivatives of the macroparameters in the 
 !! Next we try to add an additional Maxwellian: This subroutine will attempt to add a new stream with zero density and 
 !! compute the updated (longer) right side. 
 call kMaxwl_AddMxwl0D_DGV(fmxwls,right_side,coll_int,deriv_matrix,num_mxwls,sample_u,sample_v,sample_w,sample_p)
 !! 
 deallocate (deriv_matrix,coll_int,sample_u,sample_v,sample_w)
end subroutine kMaxwl_RightSideCrcAddMxl0Dv1


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!! kMaxwl_AddMxwl0D_DGV 
!!
!!
!! This subroutine adds an additional stream to the list of approximating maxwellians 
!! Information about the residual of the least squeares problem is used to add a new stream
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine kMaxwl_AddMxwl0D_DGV(fmxwls,right_side,coll_int,deriv_matrix,num_mxwls,sample_u,sample_v,sample_w,sample_p)

use DGV_commvar, only: u_L,u_R,v_L,v_R,w_L,w_R

use DGV_distributions_mod 
use DGV_dgvtools_mod

real (DP), dimension (:), pointer, intent(out) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocity
real (DP), dimension (:), pointer, intent(out) :: right_side ! the time step for
real (DP), dimension (:), intent (in) :: coll_int ! values of the collision integral in the sample points
real (DP), dimension (:), intent (in) :: sample_u,sample_v,sample_w,sample_p ! arrays that contain sample points,sample_p(i) is the acceptance rejection rate  
                      ! of the point sample_u(i),sample_v(i),sample_w(i) in the sample.
real (DP), dimension (:,:), intent (in) :: deriv_matrix ! matrix of the derivative of the folution with respect to the macroscopic parameters                        
					! the columns correspond to the derivatives of the macroscopic values in the following sequence:
					! density, temperature then three components of the bulk velocity for each stream
integer (I4B), intent (out) :: num_mxwls    ! current number of the streams.    
real (DP), dimension (:), pointer :: scrap_arry_ptr ! scrap pointer to an allocatable array
integer :: loc_alloc_stat
!!!!!!!!!!!!!!!!!!!!!!!
! Bunch of constants that determine the behaviour of the subroutine:
!
!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP) :: res_l2_tres = 1.0E-2   ! threshold for attempting reducing the residual				
real (DP) :: der_treshold = 1.0D+1  ! sensitivity threshold for the derivative  of the density
real (DP) :: allowed_temp_max = 0.8_DP ! maximum allowed temperature for the new stream
real (DP) :: allowed_temp_min = 0.01_DP ! maximum allowed temperature for the new stream
real (DP) :: res_improv_factor = 2.0_DP ! resudual improvement threshold. 

!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), allocatable :: res,res_temp,phihat ! array to sample values of the additional derivative column
real (DP), dimension (:,:), allocatable :: der_hat ! temp array to hold enlarged derivative
real (DP) :: res_l2,res_temp_l2 ! scrap variables to keep l-2 norm of the residual
real (DP) :: dnpp1,Tpp1,upp1,vpp1,wpp1 ! 
integer (I4B) :: i,j,ns ! scrap counters 
logical :: accept_stream !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
ns=size(sample_u,1)
allocate (res(1:ns), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then
 print *, "kMaxwl_AddMxwl0D_DGV: Allocation error for (res)"   
 stop
endif 
!! Now we will look at the residual and see if we need / can introduce an addition stream. 
res = matmul(deriv_matrix,right_side) - coll_int ! now res contains the residual 
res_l2 = sqrt(sum(res**2)/sum(coll_int**2)) ! relative l-2 norm of the residual
!!! Only try to add a new stream if residual is larger than the treshhold...
if (res_l2 > res_l2_tres) then 
  !!! We will attempt to add a new stream 
  !!! First step is to attempt to find macroparameters of the new steam by using importance sampling
  upp1 = u_L - .1 
  vpp1 = v_L - .1 
  wpp1 = w_L - .1 
  Tpp1 = 0 
  !!!
  dnpp1 = sum(res/sample_p)/Real(ns,DP) ! estimate derivative of the density of the new stream (dens=0 on this step)
  if (dnpp1 > der_treshold) then
   upp1 = sum(sample_u*coll_int/sample_p)/Real(ns,DP)/dnpp1
   vpp1 = sum(sample_v*coll_int/sample_p)/Real(ns,DP)/dnpp1
   wpp1 = sum(sample_w*coll_int/sample_p)/Real(ns,DP)/dnpp1
   Tpp1 = sum(((sample_u-upp1)**2+(sample_v-vpp1)**2+(sample_w-wpp1)**2)*coll_int/sample_p)&
           *2.0_DP/3.0_DP/Real(ns,DP)/dnpp1
           !!!!!!!!!!!!!! DEBUG
           print *, "kMaxwl_AddMxwl0D_DGV: trying a stream: n_delta'=", dnpp1
           print *, "kMaxwl_AddMxwl0D_DGV: trying a stream: T, u, v, w", Tpp1, upp1, vpp1, wpp1
           !!!!!!!!!!!!! END DEBUG 
   !!! A lot of sanity and other checks before accepting the new stream
   accept_stream=.TRUE.
   if ((Tpp1 > allowed_temp_max) .or. (Tpp1 < allowed_temp_min)) then 
    accept_stream=.FALSE.
    print *, "kMaxwl_AddMxwl0D_DGV: Tpp1 is outside bounds - no new stream. Tpp1=", Tpp1 
   end if
   if ((upp1 < u_L) .or. (upp1 > u_R)) then 
    accept_stream=.FALSE.
    print *, "kMaxwl_AddMxwl0D_DGV: upp1 is outside of the domain. no new stream", upp1
   end if
   if ((vpp1 < v_L) .or. (vpp1 > v_R)) then 
    accept_stream=.FALSE.
    print *, "kMaxwl_AddMxwl0D_DGV: vpp1 is outside of the domain. no new stream", vpp1
   end if
   if ((wpp1 < w_L) .or. (wpp1 > w_R)) then 
    accept_stream=.FALSE.
    print *, "kMaxwl_AddMxwl0D_DGV: wpp1 is outside of the domain. no new stream", wpp1
   end if
   ! check if maximum of the streams has been reached:
   if (num_mxwls >= max_num_mxwls) then 
    accept_stream=.FALSE.
    print *, "kMaxwl_AddMxwl0D_DGV: num of streams at maximum. no new stream", num_mxwls
   end if
   !!!!!!!!!!!!!!!!!!!!!!!! 
   !!! Last check --- did the residual actualy become smaller?
   !!! we need to re-solve the least squares... 
   !!! We only check this if previous cheks did not fail
   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
   if (accept_stream) then 
    allocate (der_hat(1:ns,1:size(deriv_matrix,2)+1), phihat(1:size(deriv_matrix,2)+1),&
      res_temp(1:ns), stat=loc_alloc_stat)
    if (loc_alloc_stat >0) then
     print *, "kMaxwl_AddMxwl0D_DGV: Allocation error for (der_hat,phihat,res_temp)"   
     stop
    endif 
    do i=1,size(deriv_matrix,2)
     do j=1,ns
      der_hat(j,i)=deriv_matrix(j,i)
     end do 
    end do   
    do j=1,ns
     der_hat(j,size(deriv_matrix,2)+1) = & 
        maxwelveldist(Tpp1,upp1,vpp1,wpp1,1.0_DP,sample_u(j),sample_v(j),sample_w(j)) 
    end do 
    ! derivative is ready
    ! reset the vector of derivatives:
    phihat = 0
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! Solve the least squares problem to determine the derivaitve of the parameters.
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    ! solve linear least squares to determine values of the right side of the ODE 
    phihat = linleastsqrs(der_hat,coll_int) 
    ! now phihat contains the values of the derivatives of the macroparameters in the 
    ! evaluate the new residual
    res_temp = matmul(der_hat,phihat) - coll_int ! now res contains the residual 
    res_temp_l2 = sqrt(sum(res_temp**2)/sum(coll_int**2)) ! relative l-2 norm of the residual
    if (res_improv_factor*res_temp_l2 > res_l2) then 
     accept_stream=.FALSE.
     print *, "kMaxwl_AddMxwl0D_DGV:  stream does not reduce residual. no new stream"
    end if
    !!! END OF checks to create new stream.
    if (accept_stream) then  
     !!! Introduce an additional stream !!!!
     num_mxwls = num_mxwls+1
     ! create a memory block that contains p+1 streams
     allocate (scrap_arry_ptr(num_mxwls*5), stat=loc_alloc_stat)
     if (loc_alloc_stat >0) then
      print *, "kMaxwl_AddMxwl0D_DGV: Allocation error for (scrap_arry_ptr)"
      stop
     endif 
     ! copying macroparameters in the new/larger memory block 
     ns=num_mxwls-1
     scrap_arry_ptr(1:ns*5) = fmxwls(:)
     scrap_arry_ptr(ns*5+1) = 0
     scrap_arry_ptr(ns*5+2) = Tpp1
     scrap_arry_ptr(ns*5+3) = upp1
     scrap_arry_ptr(ns*5+4) = vpp1
     scrap_arry_ptr(ns*5+5) = wpp1
     !!!!!!!!!!!!!!! 
     print *, "kMaxwl_AddMxwl0D_DGV: adding a stream: 0, T, u, v, w", Tpp1, upp1, vpp1, wpp1
     !!!!!!!!!!!!!!
     ! end writing to the larger memory block that contains p+1 streams
     ! next need to deallocate the old one
     deallocate (fmxwls)
     fmxwls => scrap_arry_ptr
     ! make the scrap pointer point to nothing
     nullify (scrap_arry_ptr)
     !!!!!!!!!!!!!!! Updated fmxwls
     ! create a memory block that contains p+1 streams
     allocate (scrap_arry_ptr((ns+1)*5), stat=loc_alloc_stat)
     if (loc_alloc_stat >0) then
      print *, "kMaxwl_AddMxwl0D_DGV:  2 Allocation error for (scrap_arry_prt)"
      stop
     endif 
     ! copying macroparameters in the new/larger memory block 
     scrap_arry_ptr(1:ns*5+1) = phihat
     scrap_arry_ptr(ns*5+2) = 0
     scrap_arry_ptr(ns*5+3) = 0
     scrap_arry_ptr(ns*5+4) = 0
     scrap_arry_ptr(ns*5+5) = 0
     ! end writing to the larger memory block that contains p+1 streams
     ! next need to deallocate the old one
     deallocate (right_side)
     right_side => scrap_arry_ptr
     ! make the scrap pointer point to nothing
     nullify (scrap_arry_ptr)
     !!!!!!!!!!!!!!! END Update right_side
    end if  
    deallocate (der_hat,phihat,res_temp)
   end if  
  end if  
end if
deallocate (res)
! all done
end subroutine kMaxwl_AddMxwl0D_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollArryS1Pts_DGV(sample_u,sample_v,sample_w,coll_int)
!
! This subroutine evaluates collision integral at a number of velocity points 
! using the Alekseenko-Eswar nodal-DG velocity formulation. 
! 
! The evaluation is done using A operator on a uniform grid in the velocity space. 
! 
! Each velocity point is placed on the on a uniform velocity grid associated with a pre-computed Kernel of 
! the nodal-DG formulation. 
!
! For each velocity point a grid shift vector is computed that corresponds to a shift from the celnter of the velocity 
! cell where the velocity point is located to the center of the velocity cells for which the collision operator is computed. 
! 
! This shit vector is used to evaluate the collision operator in that particular velocity:
!
! in the current setting, only s=1 piece-wise constant vectors will work properly. However, the method can be generalized to 
! an arbitrary order of nodal-DG discretization
!  
! Recall some properties of the collision kernel A: 
! For each velocity A assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!
!!!!!!!!!!!!!

subroutine EvalCollArryS1Pts_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,cells_gou,cells_gov,cells_gow,&
                   cells_lu,cells_lv,cells_lw,&  
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_u,nodes_v,nodes_w,nodes_gwts, &
                   u_R,u_L,v_R,v_L,w_R,w_L

use DGV_distributions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: sample_u,sample_v,sample_w ! the values of the velocity where the collision integral needs to be evaluated. 
real (DP), dimension (:), intent (out) :: coll_int ! the value of the collision operator for each velocity point.
real (DP), dimension (:), intent (in) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: i,j,jj ! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: capphi ! scrap variables. ! capphi keeps the number of records in A for the basis function phi
integer :: nn ! scrap variable to keep the number of Maxwellians 
real (DP) :: au, av, aw, du, dv, dw ! scrap variables to keep the velocity cell. 
real (DP) :: opA_lu,opA_lv,opA_lw
real (DP) :: shift_u,shift_v,shift_w ! scrap varibles to keep shift in the velocity variable 
integer (I4B) :: cell_opA! cell_bf scrap variable to keep the number of the cell where operator A is computed 
real (DP) :: f_xi,f_xi1 ! scrap variable to keep value of the solution at a velocity point
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! a quick check: dofc >1 the subroutine will return incorrect value
if (dofc>1) then 
 print *, "EvalCollArryS1Pts_DGV: Subroutine can only be applied to piece-wse constant DG formulation. Exit"
 stop
end if
! determine the number of the records in A (for s=1 there is only one basis function)
capphi = A_capphi(nodes_phican(1))
! IMPORTANT: We also assume that there is only one grid ! should use the parameters compatible with operator A 
pgcu=grids_cap_u(1)-1 ! number of cells in u 
pgcv=grids_cap_v(1)-1 ! number of cells in v
pgcw=grids_cap_w(1)-1 ! number of cells in w
! 
du=(u_R-u_L)/Real(pgcu,DP) ! cell size 
dv=(v_R-v_L)/Real(pgcv,DP) ! 
dw=(w_R-w_L)/Real(pgcw,DP) !
! determine the number of the cell where function to which operator A corresponds is located. for s=1 there will be 
! only one basis function. For s>1 there will be several, but thwy still will be on the same cell
cell_opA = nodes_pcell(A_phi(1)) ! This is a memory operation, perhaps can evaluate directly faster 
opA_lu = cells_lu(cell_opA) ! left boundary point of the canonical cell in u variable 
opA_lv = cells_lv(cell_opA) ! left boundary point of the canonical cell in v variable 
opA_lw = cells_lw(cell_opA) ! left boundary point of the canonical cell in w variable 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main loop over the arrays of sample velocity cells 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nn=size(fmxwls,1)/5
do i=1,size(sample_u,1)
 coll_int(i) = 0 ! nullify the result
 au = sample_u(i)
 av = sample_v(i)
 aw = sample_w(i)
 ! Now we need to figure shift from the cell_opA for each component of the velocity. 
 ! figuring shift in component u
 shift_u = u_L-opA_lu !
 j = 0
 if (au >= u_L) then  ! need to check where the point is located. 
  do while (au >= u_L+(j+1)*du)
   j = j+1
   shift_u=shift_u+du
   !!!!!!!!!!!!!!!!!!!!!!
   ! sanity check 
   if (j>1000) then
    print *, "EvalCollArryS1Pts_DGV: Stuck in endless loop 1. Exit"
    stop 
   end if
   ! end sanity check  
   !!!!!!!!!!!!!!!!!!!!!!
  end do
 else   
  do while (au > u_L-j*du)
   j = j+1
   shift_u=shift_u-du
   !!!!!!!!!!!!!!!!!!!!!!
   ! sanity check 
   if (j>1000) then
    print *, "EvalCollArryS1Pts_DGV: Stuck in endless loop 1. Exit"
    stop 
   end if
   ! end sanity check  
   !!!!!!!!!!!!!!!!!!!!!!
  end do
 end if
 ! shift in u is computed
 ! figuring shift in component v
 shift_v = v_L-opA_lv !
 j = 0
 if (av >= v_L) then  ! need to check where the point is located. 
  do while (av >= v_L+(j+1)*dv)
   j = j+1
   shift_v=shift_v+dv
   !!!!!!!!!!!!!!!!!!!!!!
   ! sanity check 
   if (j>1000) then
    print *, "EvalCollArryS1Pts_DGV: Stuck in endless loop 3. Exit"
    stop 
   end if
   ! end sanity check  
   !!!!!!!!!!!!!!!!!!!!!!
  end do
 else   
  do while (av > v_L-j*dv)
   j = j+1
   shift_v=shift_v-dv
   !!!!!!!!!!!!!!!!!!!!!!
   ! sanity check 
   if (j>1000) then
    print *, "EvalCollArryS1Pts_DGV: Stuck in endless loop 4. Exit"
    stop 
   end if
   ! end sanity check  
   !!!!!!!!!!!!!!!!!!!!!!
  end do
 end if
 ! shift in v is computed   
 ! figuring shift in component w
 shift_w = w_L-opA_lw !
 j = 0
 if (aw >= w_L) then  ! need to check where the point is located. 
  do while (aw >= w_L+(j+1)*dw)
   j = j+1
   shift_w=shift_w+dw
   !!!!!!!!!!!!!!!!!!!!!!
   ! sanity check 
   if (j>1000) then
    print *, "EvalCollArryS1Pts_DGV: Stuck in endless loop 5. Exit"
    stop 
   end if
   ! end sanity check  
   !!!!!!!!!!!!!!!!!!!!!!
  end do
 else   
  do while (aw > w_L-j*dw)
   j = j+1
   shift_w=shift_w-dw
   !!!!!!!!!!!!!!!!!!!!!!
   ! sanity check 
   if (j>1000) then
    print *, "EvalCollArryS1Pts_DGV: Stuck in endless loop 6. Exit"
    stop 
   end if
   ! end sanity check  
   !!!!!!!!!!!!!!!!!!!!!!
  end do
 end if
 ! shift in w is computed   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now that the shifts in each of the three variables are computed, we can evaluate the collision operator at the point
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  do j=1,capphi
    ! compute the values of solution at first velocity point
    au=nodes_u(A_xi(j))+shift_u
    av=nodes_v(A_xi(j))+shift_v
    aw=nodes_w(A_xi(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi = 0
    do jj=1,nn
     f_xi = f_xi + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
         fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    ! compute the value of the solution at second velocity point
    au=nodes_u(A_xi1(j))+shift_u
    av=nodes_v(A_xi1(j))+shift_v
    aw=nodes_w(A_xi1(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi1 = 0
    do jj=1,nn
     f_xi1 = f_xi1 + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
       fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    ! add a term to the collision integral 
    coll_int(i)=coll_int(i) + 2*f_xi*f_xi1*A(j)   
  end do 
  coll_int(i)=coll_int(i)/nodes_gwts(1) ! all weights are the same for s=1 otherwise use coorepsonding node on the canonical cell
 ! the value of the collision integral is computed for this velocity 
end do  ! End of the main loop in velocity points
!
end subroutine EvalCollArryS1Pts_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollArryDecompS1Pts_DGV(sample_u,sample_v,sample_w,coll_int)
!
! This subroutine is an analog of the above subroutine, except it is evaluated using decomposition f=f_M+Df algorithm
! 
!
! This subroutine evaluates collision integral at a number of velocity points 
! using the Alekseenko-Eswar nodal-DG velocity formulation. 
! 
! The evaluation is done using A operator on a uniform grid in the velocity space. 
! 
! Each velocity point is placed on the on a uniform velocity grid associated with a pre-computed Kernel of 
! the nodal-DG formulation. 
!
! For each velocity point a grid shift vector is computed that corresponds to a shift from the celnter of the velocity 
! cell where the velocity point is located to the center of the velocity cells for which the collision operator is computed. 
! 
! This shit vector is used to evaluate the collision operator in that particular velocity:
!
! in the current setting, only s=1 piece-wise constant vectors will work properly. However, the method can be generalized to 
! an arbitrary order of nodal-DG discretization
!  
! Recall some properties of the collision kernel A: 
! For each velocity A assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!
!!!!!!!!!!!!!

subroutine EvalCollArryDecompS1Pts_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,cells_gou,cells_gov,cells_gow,&
                   cells_lu,cells_lv,cells_lw,&  
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_u,nodes_v,nodes_w,nodes_gwts, &
                   u_R,u_L,v_R,v_L,w_R,w_L

use DGV_distributions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: sample_u,sample_v,sample_w ! the values of the velocity where the collision integral needs to be evaluated. 
real (DP), dimension (:), intent (out) :: coll_int ! the value of the collision operator for each velocity point.
real (DP), dimension (:), intent (in) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: i,j,jj ! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: capphi ! scrap variables. ! capphi keeps the number of records in A for the basis function phi
integer :: nn ! scrap variable to keep the number of Maxwellians 
real (DP) :: au, av, aw, du, dv, dw ! scrap variables to keep the velocity cell. 
real (DP) :: opA_cu,opA_cv,opA_cw
real (DP) :: shift_u,shift_v,shift_w ! scrap varibles to keep shift in the velocity variable 
integer (I4B) :: cell_opA! cell_bf scrap variable to keep the number of the cell where operator A is computed 
real (DP) :: temp,ubar,vbar,wbar,ndens !
real (DP) :: f_xi,f_xi1,fm_xi,fm_xi1 ! scrap variables to keep value of the solution and the local maxwellian at a velocity point
integer :: zz ! scrap integer
!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! a quick check: dofc >1 the subroutine will return incorrect value
if (dofc>1) then 
 print *, "EvalCollArryS1Pts_DGV: Subroutine can only be applied to piece-wse constant DG formulation. Exit"
 stop
end if
! determine the number of the records in A (for s=1 there is only one basis function)
capphi = A_capphi(nodes_phican(1))
! IMPORTANT: We also assume that there is only one grid ! should use the parameters compatible with operator A 
pgcu=grids_cap_u(1)-1 ! number of cells in u 
pgcv=grids_cap_v(1)-1 ! number of cells in v
pgcw=grids_cap_w(1)-1 ! number of cells in w
! 
du=(u_R-u_L)/Real(pgcu,DP) ! cell size 
dv=(v_R-v_L)/Real(pgcv,DP) ! 
dw=(w_R-w_L)/Real(pgcw,DP) !
! determine the number of the cell where function to which operator A corresponds is located. for s=1 there will be 
! only one basis function. For s>1 there will be several, but thwy still will be on the same cell
cell_opA = nodes_pcell(A_phi(1)) ! This is a memory operation, perhaps can evaluate directly faster 
opA_cu = cells_lu(cell_opA)+du/2.0_DP ! center point of the canonical cell in u variable 
opA_cv = cells_lv(cell_opA)+dv/2.0_DP ! center point of the canonical cell in v variable 
opA_cw = cells_lw(cell_opA)+dw/2.0_DP ! center point of the canonical cell in w variable 
! First we will evaluate paramters of the Maxwellian that is the sum of all the streams:
ndens = 0 
ubar = 0
vbar = 0
wbar = 0
temp = 0
nn=size(fmxwls,1)/5
do zz=1,nn
    ndens=ndens+fmxwls((zz-1)*5+1)
    ubar = ubar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+3)
    vbar = vbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+4)
    wbar = wbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+5)
    temp = temp + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+2)+ &
       2.0_DP/3.0_DP*fmxwls((zz-1)*5+1)*( fmxwls((zz-1)*5+3)**2 + &
       fmxwls((zz-1)*5+4)**2 + fmxwls((zz-1)*5+5)**2 )
end do
ubar=ubar/ndens
vbar=vbar/ndens
wbar=wbar/ndens
! caclulated total density and bulk velocity. Now we have an additiona term for the temperature:
temp = temp/ndens - 2.0_DP/3.0_DP*(ubar**2 + vbar**2 + wbar**2) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main loop over the arrays of sample velocity cells 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,size(sample_u,1)
 coll_int(i) = 0 ! nullify the result
 ! Now we need to figure shift from the cell_opA for each component of the velocity. 
 shift_u = sample_u(i) - opA_cu !
 shift_v = sample_v(i) - opA_cv !
 shift_w = sample_w(i) - opA_cw !
 ! shift is computed   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now that the shifts in each of the three variables are computed, we can evaluate the collision operator at the point
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !!! !!!!
 do j=1,capphi
    ! compute the values of solution at first velocity point
    au=nodes_u(A_xi(j))+shift_u
    av=nodes_v(A_xi(j))+shift_v
    aw=nodes_w(A_xi(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi = 0
    do jj=1,nn
     f_xi = f_xi + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
         fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    fm_xi = maxwelveldist(temp,ubar,vbar,wbar,ndens,au,av,aw)
    ! compute the value of the solution at second velocity point
    au=nodes_u(A_xi1(j))+shift_u
    av=nodes_v(A_xi1(j))+shift_v
    aw=nodes_w(A_xi1(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi1 = 0
    do jj=1,nn
     f_xi1 = f_xi1 + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
       fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    fm_xi1=maxwelveldist(temp,ubar,vbar,wbar,ndens,au,av,aw)
    ! add a term to the collision integral 
    coll_int(i) = coll_int(i) + (f_xi+fm_xi)*(f_xi1-fm_xi1)*A(j)+(f_xi-fm_xi)*(f_xi1+fm_xi1)*A(j)   
  end do 
  coll_int(i) = coll_int(i)/nodes_gwts(1) ! currently only works for s=1 models -- need to change for s>1./
 ! the value of the collision integral is computed for this velocity 
end do  ! End of the main loop in velocity points
!
end subroutine EvalCollArryDecompS1Pts_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollArryDecompS1PtsOMP_DGV(sample_u,sample_v,sample_w,coll_int)
!
! This subrtoutine is a copy of the above subroutine except it has openMP parallelization of the main loop.
! 
!
! This subroutine evaluates collision integral at a number of velocity points 
! using the Alekseenko-Eswar nodal-DG velocity formulation. 
! 
! The evaluation is done using A operator on a uniform grid in the velocity space. 
! 
! Each velocity point is placed on the on a uniform velocity grid associated with a pre-computed Kernel of 
! the nodal-DG formulation. 
!
! For each velocity point a grid shift vector is computed that corresponds to a shift from the celnter of the velocity 
! cell where the velocity point is located to the center of the velocity cells for which the collision operator is computed. 
! 
! This shit vector is used to evaluate the collision operator in that particular velocity:
!
! in the current setting, only s=1 piece-wise constant vectors will work properly. However, the method can be generalized to 
! an arbitrary order of nodal-DG discretization
!  
! Recall some properties of the collision kernel A: 
! For each velocity A assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!
!!!!!!!!!!!!!

subroutine EvalCollArryDecompS1PtsOMP_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,cells_gou,cells_gov,cells_gow,&
                   cells_lu,cells_lv,cells_lw,&  
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_u,nodes_v,nodes_w,nodes_gwts, &
                   u_R,u_L,v_R,v_L,w_R,w_L,Num_OMP_threads

use DGV_distributions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: sample_u,sample_v,sample_w ! the values of the velocity where the collision integral needs to be evaluated. 
real (DP), dimension (:), intent (out) :: coll_int ! the value of the collision operator for each velocity point.
real (DP), dimension (:), intent (in) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: i,j,jj ! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: capphi ! scrap variables. ! capphi keeps the number of records in A for the basis function phi
integer :: nn ! scrap variable to keep the number of Maxwellians 
real (DP) :: au, av, aw, du, dv, dw ! scrap variables to keep the velocity cell. 
real (DP) :: opA_cu,opA_cv,opA_cw
real (DP) :: shift_u,shift_v,shift_w ! scrap varibles to keep shift in the velocity variable 
integer (I4B) :: cell_opA! cell_bf scrap variable to keep the number of the cell where operator A is computed 
real (DP) :: temp,ubar,vbar,wbar,ndens !
real (DP) :: f_xi,f_xi1,fm_xi,fm_xi1 ! scrap variables to keep value of the solution and the local maxwellian at a velocity point
integer :: zz ! scrap integer
integer :: iiii ! a test variable to play with OpenMP runtime functions calls
!!!!!!!!!!!!!!!!!!!!!!!
interface 
 function omp_get_thread_num() result (y)
   integer :: y 
 end function omp_get_thread_num
 function omp_get_num_threads() result (y)
  integer :: y 
 end function omp_get_num_threads 
 function omp_get_num_procs() result (y)
  integer :: y 
 end function omp_get_num_procs
 function omp_get_stack_size () result (y)
  use nrtype
  integer (I2B) :: y
 end function omp_get_stack_size
end interface  
!!!!!!!!!!!!!!!!!
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! a quick check: dofc >1 the subroutine will return incorrect value
if (dofc>1) then 
 print *, "EvalCollArryS1Pts_DGV: Subroutine can only be applied to piece-wse constant DG formulation. Exit"
 stop
end if
! determine the number of the records in A (for s=1 there is only one basis function)
capphi = A_capphi(nodes_phican(1))
! IMPORTANT: We also assume that there is only one grid ! should use the parameters compatible with operator A 
pgcu=grids_cap_u(1)-1 ! number of cells in u 
pgcv=grids_cap_v(1)-1 ! number of cells in v
pgcw=grids_cap_w(1)-1 ! number of cells in w
! 
du=(u_R-u_L)/Real(pgcu,DP) ! cell size 
dv=(v_R-v_L)/Real(pgcv,DP) ! 
dw=(w_R-w_L)/Real(pgcw,DP) !
! determine the number of the cell where function to which operator A corresponds is located. for s=1 there will be 
! only one basis function. For s>1 there will be several, but thwy still will be on the same cell
cell_opA = nodes_pcell(A_phi(1)) ! This is a memory operation, perhaps can evaluate directly faster 
opA_cu = cells_lu(cell_opA)+du/2.0_DP ! center point of the canonical cell in u variable 
opA_cv = cells_lv(cell_opA)+dv/2.0_DP ! center point of the canonical cell in v variable 
opA_cw = cells_lw(cell_opA)+dw/2.0_DP ! center point of the canonical cell in w variable 
! First we will evaluate paramters of the Maxwellian that is the sum of all the streams:
ndens = 0 
ubar = 0
vbar = 0
wbar = 0
temp = 0
nn=size(fmxwls,1)/5
do zz=1,nn
    ndens=ndens+fmxwls((zz-1)*5+1)
    ubar = ubar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+3)
    vbar = vbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+4)
    wbar = wbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+5)
    temp = temp + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+2)+ &
       2.0_DP/3.0_DP*fmxwls((zz-1)*5+1)*( fmxwls((zz-1)*5+3)**2 + &
       fmxwls((zz-1)*5+4)**2 + fmxwls((zz-1)*5+5)**2 )
end do
ubar=ubar/ndens
vbar=vbar/ndens
wbar=wbar/ndens
! caclulated total density and bulk velocity. Now we have an additiona term for the temperature:
temp = temp/ndens - 2.0_DP/3.0_DP*(ubar**2 + vbar**2 + wbar**2) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main loop over the arrays of sample velocity cells 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,size(sample_u,1)
 coll_int(i) = 0 ! nullify the result
 ! Now we need to figure shift from the cell_opA for each component of the velocity. 
 shift_u = sample_u(i) - opA_cu !
 shift_v = sample_v(i) - opA_cv !
 shift_w = sample_w(i) - opA_cw !
 ! shift is computed   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now that the shifts in each of the three variables are computed, we can evaluate the collision operator at the point
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 !!! !!!!
! iiii=omp_get_num_procs
! OpenMP set the number of threads: 

call omp_set_num_threads(Num_OMP_threads)
!$OMP PARALLEL DO PRIVATE(j,au,av,aw,f_xi,jj,fm_xi,f_xi1,fm_xi1,iiii)  &
!$OMP  NUM_THREADS(Num_OMP_threads) SCHEDULE(DYNAMIC, 50)  
 !!! !!!!
 do j=1,capphi
    ! compute the values of solution at first velocity point
    au=nodes_u(A_xi(j))+shift_u
    av=nodes_v(A_xi(j))+shift_v
    aw=nodes_w(A_xi(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi = 0
    do jj=1,nn
     f_xi = f_xi + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
         fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    fm_xi = maxwelveldist(temp,ubar,vbar,wbar,ndens,au,av,aw)
    ! compute the value of the solution at second velocity point
    au=nodes_u(A_xi1(j))+shift_u
    av=nodes_v(A_xi1(j))+shift_v
    aw=nodes_w(A_xi1(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi1 = 0
    do jj=1,nn
     f_xi1 = f_xi1 + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
       fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    fm_xi1=maxwelveldist(temp,ubar,vbar,wbar,ndens,au,av,aw)
    ! add a term to the collision integral 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This assignment only done one processor at a time !!!
    !$omp critical   
    coll_int(i) = coll_int(i) + (f_xi + fm_xi)*(f_xi1-fm_xi1)*A(j)+(f_xi-fm_xi)*(f_xi1+fm_xi1)*A(j)   
    !$omp end critical
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   
    ! add this two lines to track progress and Print a hello message to Check if parallel stuff worked... !!!!
    ! iiii = omp_get_thread_num()
    ! print *, "j=", j, "Thread", iiii, "i=%i8", i    
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
  end do  ! End of parallel loop in j
  !!!!!!!!!!
  coll_int(i) = coll_int(i)/nodes_gwts(1) ! currently only works for s=1 models -- need to change for s>1./
  ! the value of the collision integral is computed for this velocity 
end do  ! End of the main loop in velocity points
!
end subroutine EvalCollArryDecompS1PtsOMP_DGV

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! EvalCollArryDecompS1Pts_MPI_DGV(sample_u,sample_v,sample_w,coll_int)
!
! This subrtoutine is a copy of the above subroutine except it was adopted to work with MPI 
! parallelization. In the supported implementation, the operator A is distributed between 
! the processors and the collision operator is evaluated at each velocity point in the sample_u,_v,_w 
! list. Evaluation of the collision oprator is distributed among MPI processors. The idea is that 
! each of the processors will compute just a portion of the collision operator for each velocity point
! then an MPI gather operation is performedto combine partial results on the master processor with Irank =0.
! 
!
! This subroutine evaluates collision integral at a number of velocity points 
! using the Alekseenko-Eswar nodal-DG velocity formulation. 
! 
! The evaluation is done using A operator on a uniform grid in the velocity space. 
! 
! Each velocity point is placed on the on a uniform velocity grid associated with a pre-computed Kernel of 
! the nodal-DG formulation. 
!
! For each velocity point a grid shift vector is computed that corresponds to a shift from the celnter of the velocity 
! cell where the velocity point is located to the center of the velocity cells for which the collision operator is computed. 
! 
! This shit vector is used to evaluate the collision operator in that particular velocity:
!
! in the current setting, only s=1 piece-wise constant vectors will work properly. However, the method can be generalized to 
! an arbitrary order of nodal-DG discretization
!  
! Recall some properties of the collision kernel A: 
! For each velocity A assumes a uniform velocity grid and a periodic structure of the 
! basis functions. Operator A is only evaluated for one basis function. The values 
! of A for other basis functions are obtained by a transposition from operator A on 
! the canonical cell. (The canonical cell should be selected at the center of the mesh 
! or close to the center. The value of A_phi can be used to determine which basis function 
! can be used a s the canonical. However, it is best to maintain proper record on how Akor was computed.)
! 
! The subroutine expects that the following specialized variables are prepared: 
!  nodes_phican
!
!!!!!!!!!!!!!

subroutine EvalCollArryDecompS1Pts_MPI_DGV(sample_u,sample_v,sample_w,coll_int,fmxwls)

use DGV_commvar, only: A,A_capphi,A_xi,A_xi1,A_phi,cells_gou,cells_gov,cells_gow,&
                   cells_lu,cells_lv,cells_lw,&  
                   grids_cap_u,grids_cap_v,grids_cap_w,nodes_Ashift,nodes_phican,&
                   nodes_dui,nodes_dvi,nodes_dwi,nodes_pcell,cells_ugi,cells_vgi,&
                   cells_wgi,nodes_u,nodes_v,nodes_w,nodes_gwts, &
                   u_R,u_L,v_R,v_L,w_R,w_L

use DGV_distributions_mod

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: sample_u,sample_v,sample_w ! the values of the velocity where the collision integral needs to be evaluated. 
real (DP), dimension (:), intent (out) :: coll_int ! the value of the collision operator for each velocity point.
real (DP), dimension (:), intent (in) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer (I4B) :: gou,gov,gow,dofc ! are the numbers of nodal points in each cell in directions u,v and w. Dofc=all degrees o f freedom on one cell
integer (I4B) :: i,j,jj ! scrap counters
integer (I4B) :: pgcu,pgcv,pgcw ! number of cells on the grid (we assume that there is only one uniform grid
integer (I4B) :: capphi ! scrap variables. ! capphi keeps the number of records in A for the basis function phi
integer :: nn ! scrap variable to keep the number of Maxwellians 
real (DP) :: au, av, aw, du, dv, dw ! scrap variables to keep the velocity cell. 
real (DP) :: opA_cu,opA_cv,opA_cw
real (DP) :: shift_u,shift_v,shift_w ! scrap varibles to keep shift in the velocity variable 
integer (I4B) :: cell_opA! cell_bf scrap variable to keep the number of the cell where operator A is computed 
real (DP) :: temp,ubar,vbar,wbar,ndens !
real (DP) :: f_xi,f_xi1,fm_xi,fm_xi1 ! scrap variables to keep value of the solution and the local maxwellian at a velocity point
integer :: zz ! scrap integer
integer :: iiii ! a test variable to play with OpenMP runtime functions calls
!
! IMPORTANT: We assume that all cells are identical, in particular that they have the same number of nodes in each dimension
!
gou=cells_gou(1)
gov=cells_gov(1)
gow=cells_gow(1)
dofc=gou*gov*gow
! a quick check: dofc >1 the subroutine will return incorrect value
if (dofc>1) then 
 print *, "EvalCollArryDecompS1Pts_MPI_DGV: Subroutine can only be applied to piece-wse constant DG formulation. Exit"
 stop
end if
! determine the number of the records in A (for s=1 there is only one basis function)
capphi = A_capphi(nodes_phican(1))
! IMPORTANT: We also assume that there is only one grid ! should use the parameters compatible with operator A 
pgcu=grids_cap_u(1)-1 ! number of cells in u 
pgcv=grids_cap_v(1)-1 ! number of cells in v
pgcw=grids_cap_w(1)-1 ! number of cells in w
! 
du=(u_R-u_L)/Real(pgcu,DP) ! cell size 
dv=(v_R-v_L)/Real(pgcv,DP) ! 
dw=(w_R-w_L)/Real(pgcw,DP) !
! determine the number of the cell where function to which operator A corresponds is located. for s=1 there will be 
! only one basis function. For s>1 there will be several, but thwy still will be on the same cell
cell_opA = nodes_pcell(A_phi(1)) ! This is a memory operation, perhaps can evaluate directly faster 
opA_cu = cells_lu(cell_opA)+du/2.0_DP ! center point of the canonical cell in u variable 
opA_cv = cells_lv(cell_opA)+dv/2.0_DP ! center point of the canonical cell in v variable 
opA_cw = cells_lw(cell_opA)+dw/2.0_DP ! center point of the canonical cell in w variable 
! First we will evaluate paramters of the Maxwellian that is the sum of all the streams:
ndens = 0 
ubar = 0
vbar = 0
wbar = 0
temp = 0
nn=size(fmxwls,1)/5
do zz=1,nn
    ndens=ndens+fmxwls((zz-1)*5+1)
    ubar = ubar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+3)
    vbar = vbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+4)
    wbar = wbar + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+5)
    temp = temp + fmxwls((zz-1)*5+1)*fmxwls((zz-1)*5+2)+ &
       2.0_DP/3.0_DP*fmxwls((zz-1)*5+1)*( fmxwls((zz-1)*5+3)**2 + &
       fmxwls((zz-1)*5+4)**2 + fmxwls((zz-1)*5+5)**2 )
end do
ubar=ubar/ndens
vbar=vbar/ndens
wbar=wbar/ndens
! caclulated total density and bulk velocity. Now we have an additiona term for the temperature:
temp = temp/ndens - 2.0_DP/3.0_DP*(ubar**2 + vbar**2 + wbar**2) 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! main loop over the arrays of sample velocity cells 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
do i=1,size(sample_u,1)
 coll_int(i) = 0 ! nullify the result
 ! Now we need to figure shift from the cell_opA for each component of the velocity. 
 shift_u = sample_u(i) - opA_cu !
 shift_v = sample_v(i) - opA_cv !
 shift_w = sample_w(i) - opA_cw !
 ! shift is computed   
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 ! Now that the shifts in each of the three variables are computed, we can evaluate the collision operator at the point
 !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 do j=1,capphi
    ! compute the values of solution at first velocity point
    au=nodes_u(A_xi(j))+shift_u
    av=nodes_v(A_xi(j))+shift_v
    aw=nodes_w(A_xi(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi = 0
    do jj=1,nn
     f_xi = f_xi + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
         fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    fm_xi = maxwelveldist(temp,ubar,vbar,wbar,ndens,au,av,aw)
    ! compute the value of the solution at second velocity point
    au=nodes_u(A_xi1(j))+shift_u
    av=nodes_v(A_xi1(j))+shift_v
    aw=nodes_w(A_xi1(j))+shift_w
    !!!!! evaluate solution at the point au, av, aw !!!!
    f_xi1 = 0
    do jj=1,nn
     f_xi1 = f_xi1 + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
       fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),au,av,aw)
    end do 
    fm_xi1=maxwelveldist(temp,ubar,vbar,wbar,ndens,au,av,aw)
    ! add a term to the collision integral 
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! This assignment only done one processor at a time !!!
    coll_int(i) = coll_int(i) + (f_xi + fm_xi)*(f_xi1-fm_xi1)*A(j)+(f_xi-fm_xi)*(f_xi1+fm_xi1)*A(j)   
 end do  ! End of loop in j
 !!!!!!!!!!
 coll_int(i) = coll_int(i)/nodes_gwts(1) ! currently only works for s=1 models -- need to change for s>1./
 ! the value of the collision integral is computed for this velocity 
end do  ! End of the main loop in velocity points
!
end subroutine EvalCollArryDecompS1Pts_MPI_DGV


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! kMxl_EvalColESBGKPtsOMP_DGV(fmxwls,coll_int,sample_u,sample_v,sample_w) 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Right hand side calculation for the ES-BGK distribution and collision operator
! Here, the RHS = nu * (f0 - f)
! where
! nu = collision frequency
! f0 = ESBGK distribution
! f = solution
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine kMxl_EvalColESBGKPtsOMP_DGV(fmxwls,coll_int,sample_u,sample_v,sample_w) 

use DGV_distributions_mod
use DGV_dgvtools_mod
use DGV_commvar, only: nodes_gwts, nodes_u, nodes_v, nodes_w, &
				   alpha, gas_viscosity, gas_T_reference, gas_alpha, C_inf, gasR

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), intent (in) :: sample_u,sample_v,sample_w ! the values of the velocity where the collision integral needs to be evaluated. 
real (DP), dimension (:), intent (out) :: coll_int ! the value of the collision operator for each velocity point.
real (DP), dimension (:), intent (in) :: fmxwls ! arrays that contains macroparameters of the 
                      ! maxwellians that make up the solution
                      ! Format: density, temperature, then three components of bulk velocit
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real (DP), dimension (1:size(nodes_u,1)) :: f ! scrap variable, to keep values of the assembles solution on the integration grid

real (DP) :: u_0, v_0, w_0 ! bulk velocities
real (DP) :: n ! density
real (DP) :: T ! temperature
real (DP) :: Pressure
real (DP) :: nu ! this is the collision frequency term
real (DP), dimension (:), allocatable :: f0 ! Distribution function

real (DP), parameter :: kBoltzmann = 1.3806503D-23
real (DP), dimension (3,3) :: Tensor, TensorInv ! both the tensor and the inverse tensor for the ES-BGK
real (DP) :: Determinant ! the determinant of the tensor
integer :: loc_alloc_stat

integer (I4B) :: ni,jj ! scrap integer
!!!!!!!!!!!!!!!!!!
!! We need to compute a few moments of the solution to define the model. We realize that all integrals can be computed analytically and 
!! in the future we will replace the integrals with analytic formulas. However, at the time being, we use a nodal-DG grid to compute the moments numerically
!! 
!!!!!!!!!!!!!!!!!!
ni=size(fmxwls,1)/5
do jj=1,ni
 f = f + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
         fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),nodes_u,nodes_v,nodes_w)
end do 
!!!!!!!!
! compute the macroparameters for use in the tensor computation
!!!!!!!!
! density (number density)
n=sum(f*nodes_gwts)
!!!!!!!!
! momentum 
u_0=sum(f*nodes_gwts*nodes_u)/n
v_0=sum(f*nodes_gwts*nodes_v)/n
w_0=sum(f*nodes_gwts*nodes_w)/n
!!!!!!!!
! temperature
T = sum(f*nodes_gwts*((nodes_u-u_0)**2+(nodes_v-v_0)**2+(nodes_w-w_0)**2))/n/3.0_DP*2.0_DP ! dimensionless temperature
!!

! Here we evaluate the stress tensor to use in the ES model
Tensor(1,1) = sum(f*nodes_gwts*(nodes_u-u_0)**2)/n*alpha*2.0_DP
Tensor(2,2) = sum(f*nodes_gwts*(nodes_v-v_0)**2)/n*alpha*2.0_DP
Tensor(3,3) = sum(f*nodes_gwts*(nodes_w-w_0)**2)/n*alpha*2.0_DP
Tensor(1,2) = sum(f*nodes_gwts*(nodes_u-u_0)*(nodes_v-v_0))/n*alpha*2.0_DP
Tensor(1,3) = sum(f*nodes_gwts*(nodes_u-u_0)*(nodes_w-w_0))/n*alpha*2.0_DP
Tensor(2,3) = sum(f*nodes_gwts*(nodes_v-v_0)*(nodes_w-w_0))/n*alpha*2.0_DP
Tensor(2,1) = Tensor(1,2)
Tensor(3,1) = Tensor(1,3)
Tensor(3,2) = Tensor(2,3)

Tensor(1,1) = Tensor(1,1) + (1-alpha)*T
Tensor(2,2) = Tensor(2,2) + (1-alpha)*T
Tensor(3,3) = Tensor(3,3) + (1-alpha)*T

! next we compute the inverse of the ES stress tensor

TensorInv = inv(Tensor)
Determinant = DetTensor(Tensor)

! TensorInv can now be used to compute the values of the ES operator

allocate (f0(1:size(sample_u,1)), stat=loc_alloc_stat)
if (loc_alloc_stat >0) then 
 print *, "kMxl_EvalColESBGKPtsOMP_DGV: Allocation error for variables (f0)"
 stop
end if
!!! Evaluate the solution on the sample points
f0=0 ! clear the array
do jj=1,ni
 f0 = f0 + maxwelveldist(fmxwls((jj-1)*5+2),fmxwls((jj-1)*5+3),fmxwls((jj-1)*5+4),&
         fmxwls((jj-1)*5+5),fmxwls((jj-1)*5+1),sample_u,sample_v,sample_w)
end do 
!!!!!!!
f0 = ESBGK_f0(TensorInv,Determinant,n,u_0,v_0,w_0,sample_u,sample_v,sample_w)-f0 ! f0 now contains the difference between ES and the solution
! now to evaluate the collision requency term
Pressure = n*T ! dimensionless Pressure is computed here
nu = Pressure/((1-alpha))*((gas_T_reference/C_inf**2*2.0d0*gasR)/T)**gas_alpha ! final dimensionless nu ! have gas_T_reference be dimensionless?
coll_int = nu*f0

deallocate (f0)
!
end subroutine kMxl_EvalColESBGKPtsOMP_DGV



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Reduce_Maxwellians
!
!
! When solution is represented as sum of maxwellian streams, this 
! subroutine will check if streams has different bulk velocities and temperature
! If two streams have bulk velocity and temperature that are the same or close, then
! these two streams are combined into one.  
!
! NOPT DONE YET: Also, Maxwellians with low density will be absorbed in maxwellians with high density
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

subroutine Reduce_Maxwellians(maxwells_out, n_out, right_side)

use DGV_distributions_mod

real (DP), dimension (:), pointer, intent(out) :: maxwells_out ! The maxwellians that will we reviewed
integer (I4B), intent (out) :: n_out ! number of maxwellians to send out
real (DP), dimension (:), pointer :: right_side ! The updates array of maxwellians after review and possible reduction
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
real (DP), dimension (:), pointer :: maxwells_tmp ! The updates array of maxwellians after review and possible reduction
real (DP) :: temper_thld=1.0D-4 ! relative error in temperature low which the temeperature is considered the same
real (DP) :: veloci_thld=1.0D-4 ! relative error in bulk velocity below which the bulk velocity is considered the same
real (DP) :: density=1.0D-8 ! relative error in density, after which the stream is absorbed 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: i,j,k, nn ! scrap counters and variables 
real (DP) :: n_a,T_a,u_a,v_a,w_a, n_b,T_b,u_b,v_b,w_b ! scrap variables to keep macroparameters of two streams.
real (DP) :: n,T,u,v,w ! scrap variables to keep macroparameters of the combined stream.
integer :: loc_alloc_stat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
nn = size(maxwells_out,1)/5 ! this is the number of streams that came in

! review bulk velocity and temperature:
i=1
do while (i<nn) ! it could be that the stream was eliminated by now.... nn keep the number of remaining streams
  n_a=maxwells_out((i-1)*5+1)
  T_a=maxwells_out((i-1)*5+2)
  u_a=maxwells_out((i-1)*5+3)
  v_a=maxwells_out((i-1)*5+4)
  w_a=maxwells_out((i-1)*5+5)
  ! got the parameters of stream A
  ! now need to compare it to remaining streams:
  j=i+1
  do while (j<= nn) ! it could be that the streams are eliminated and the actual number of streams is less than size of the array/5
     n_b=maxwells_out((j-1)*5+1)
     T_b=maxwells_out((j-1)*5+2)
     u_b=maxwells_out((j-1)*5+3)
     v_b=maxwells_out((j-1)*5+4)
     w_b=maxwells_out((j-1)*5+5)
     ! got the parameters of stream B
     ! now compare it to those of stream a:
     if ((abs(T_b-T_a)<temper_thld) .and. &
        ((u_a-u_b)**2+(v_a-v_b)**2+(w_a-w_b)**2 < veloci_thld**2)) then 
      ! we have our first coinsiding streams! Combining the streams 
      n=n_a+n_b
      u=(n_a*u_a+n_b*u_b)/n
      v=(n_a*v_a+n_b*v_b)/n
      w=(n_a*w_a+n_b*w_b)/n
      T=(T_a*n_a+T_b*n_b)/n+((u_a**2+v_a**2+w_a**2)*n_a+(u_b**2+v_b**2+w_b**2)*n_b)/n/3.0_DP*2.0_DP - &
        (u**2+v**2+w**2)/3.0_DP*2.0_DP
      ! there are the parameters of the combined stream.
      n_a=n
      maxwells_out((i-1)*5+1)=n
      T_a=T
      maxwells_out((i-1)*5+2)=T
      u_a=u
      maxwells_out((i-1)*5+3)=u
      v_a=v
      maxwells_out((i-1)*5+4)=v
      w_a=w
      maxwells_out((i-1)*5+5)=w
      ! the stream is updated. 
      ! now we need to remove the stream that was combined with the first stream
      do k=j,nn-1 ! push down the steams to eliminate the one that was combined with the first stream
       maxwells_out((k-1)*5+1)=maxwells_out(k*5+1)
       maxwells_out((k-1)*5+2)=maxwells_out(k*5+2)
       maxwells_out((k-1)*5+3)=maxwells_out(k*5+3)
       maxwells_out((k-1)*5+4)=maxwells_out(k*5+4)
       maxwells_out((k-1)*5+5)=maxwells_out(k*5+5)
      end do
      maxwells_out((nn-1)*5+1:(nn-1)*5+5)=0
      !
      nn = nn-1  ! reduce the total number of streams in nn
      ! do not advance j --- the streams were pushed down  
     else 
      j=j+1 ! advance j -- go to the next stream B  
     end if  
  end do 
  i=i+1  ! advance i -- go to the next stream A
end do 
!!!!!!
! nn has the new number of streams
! streams were reviewed
!!!!!!
if (n_out > nn) then
 print *,"Reduce_Maxwellians: Reducing the number of Maxwellians. old number", & 
                              n_out, "new number", nn   
 !        
 allocate (maxwells_tmp(1:nn*5), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "Reduce_Maxwellians: Allocation error for variables (maxwells_tmp)"
 end if   
 maxwells_tmp(:)=maxwells_out(1:nn*5) ! moved the new streams to the new memory location
 deallocate(maxwells_out) ! free the memory
 maxwells_out => maxwells_tmp ! make maxwells_out to point to the memory location that contains the reduced streams
 nullify (maxwells_tmp)
 ! also need to change the size of the right_side array: 
 deallocate(right_side) ! release the memory for the two arrays and th
 allocate (right_side(1:nn*5), stat=loc_alloc_stat)
 if (loc_alloc_stat >0) then 
  print *, "0DKmxwl: Allocation error for variables (num_streams,right_side)"
 end if 
 ! update the number of streams:
 n_out=nn
end if !
end subroutine Reduce_Maxwellians

subroutine EvalMaxwlsOnGrid(nodes_u,nodes_v,nodes_w,fmxwls,fval)

use DGV_distributions_mod

real (DP), dimension (:), intent (in) :: nodes_u, nodes_v, nodes_w ! arrays with points in veolcity space
real (DP), dimension (:), intent (in) :: fmxwls ! the arrays with macroparamters of maxwellians. (n,T,u,v,w) 
real (DP), dimension (:), intent (out) :: fval  ! the value of the solution computed on the grid. 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
integer :: j
 !!!!! evaluate solution at the point au, av, aw !!!!
    fval = 0
    do j=1,size(fmxwls,1)/5
     fval = fval + maxwelveldist(fmxwls((j-1)*5+2),fmxwls((j-1)*5+3),fmxwls((j-1)*5+4),&
         fmxwls((j-1)*5+5),fmxwls((j-1)*5+1),nodes_u,nodes_v,nodes_w)
    end do 

end subroutine EvalMaxwlsOnGrid

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WritekMxwlSolEU
! 
! This subroutine writes the solution of k-Maxwell model obtained using Euler time step,
! the subroutine will save the solution, the current time and the time step 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine WritekMxwlSolEU(fmxwls, t, dt)
!                   
use DGV_readwrite
!
intrinsic Trim
!
real (DP), dimension (:), intent(in) :: fmxwls ! the array of solution 
real (DP), intent (in) :: t, dt ! solution time and the time step
!
character (len=15) :: suff ! the string to hold the suffix
character (len=132) :: file_name ! the variable to store the file name 

! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
write (suff, "(F14.10)") t 
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_kMxwlsol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="UNFORMATTED", access="SEQUENTIAL")
write (15) t,dt ! record the current time and the used dt
write (15) size(fmxwls,1) ! record the length of the array 

! Next we save the solution .... 
write (15) fmxwls
! end save 
close (15)
!
end subroutine WritekMxwlSolEU 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! ReadkMxwlSolEU
! 
! This subroutine reads the solution of k-Maxwell model obtained using Euler time step,
! the subroutine will read the solution, the current time and the time step from the hard drive
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!1

subroutine ReadkMxwlSolEU(maxwells,suff, t, dt)
!
use DGV_readwrite
!                   
intrinsic Trim
!
real (DP), dimension (:), pointer :: maxwells ! pointer to the array of macroparamters.
character (len=15), intent (in) :: suff ! the string to hold the suffix
real (DP), intent (out) :: t, dt ! solution time and the time step
!
character (len=132) :: file_name ! the variable to store the file name 
integer (I4B) :: m ! the length of the arrays
integer (I4B) :: loc_alloc_stat ! some dump varaible


! first, we prepare the file name to store the solution
call MakeBaseNameDGV(file_name)
file_name = trim(Adjustl(file_name))//"_time"//trim(Adjustl(suff))//"_kMxwlsol.dat"  ! this file will keep the array 
!
! now we open the file for record and save some stuff in it: 
open (15, file=file_name, position ="REWIND", action="READ", &
                   form="UNFORMATTED", access="SEQUENTIAL")
read (15) t,dt ! read the solution time and the used dt
read (15) m !  read the length of the array 
! We now need to prepare the storage for the solution.
allocate (maxwells(1:m), stat=loc_alloc_stat)
    !
  if (loc_alloc_stat >0) then 
  print *, "ReadkMxwlSolEU: Allocation error for variable maxwells"
  close(15)
  stop
  end if
!
! Next we read  the solution .... 
read (15) maxwells
! end save 
close (15)
end subroutine ReadkMxwlSolEU 
 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! WriteRecMxwls_kMxsl0D (rec_mxwls)
!
! This subroutine will dump on the disk records of the kMxwl solution in time
! ATTN: Only work for 0D problems
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine WriteRecMxwls_kMxsl0D (rec_mxwls)
!
use DGV_readwrite
!
intrinsic Trim
!!!!!
real (DP), dimension (:,:), intent (in) :: rec_mxwls ! array that keeps recods of kMwls macroparameters with time steps.
!
character (len=132) :: file_name ! the variable to store the file name
character (len=750) :: text_line1  
integer (I4B) :: n,i,j ! some counters
character (len=20) :: suff  ! the string to hold the suffix
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! first, we prepare the file name to store the solution
call MakeBaseNameDGV2(file_name, "results/")
file_name = trim(Adjustl(file_name))//"_rec_kMxwl.txt"  ! this file will keep the array
!
! now we open the file in formatted write regime for record and save some stuff in it:
open (15, file=file_name, position ="REWIND", action="WRITE", &
                   form="FORMATTED")

write(15, "(A)") "This file contains the values of macroparameters of Maxwellians that approximate velocity", &
             " distribution function in kMxwls model "
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Prepare the line to write for te 

n = (size(rec_mxwls,2)-1)/5 ! max number of maxwellians 
text_line1 = " Time"
do i = 1,n
 write (suff, "(I2)") i
 text_line1 = trim(text_line1)// ", n" // trim(Adjustl(suff)) // ", " // "T" // trim(Adjustl(suff)) &
                 // ", " // "u" // trim(Adjustl(suff)) // ", " // "v" // trim(Adjustl(suff)) &
                 // ", " // "w" // trim(Adjustl(suff))  
end do 
write(15, "(A)") text_line1
! next we prepare a line to write the macroparameters

do i=1,size(rec_mxwls,1) 
 write (text_line1,"(F18.16)") rec_mxwls(i,1) ! time of the record
 text_line1 = trim(text_line1) 
 do j=1,n
  write (suff,"(F18.16)") rec_mxwls(i,(j-1)*5+2)
  text_line1 = trim(text_line1) // ", " // trim(Adjustl(suff)) 
  write (suff,"(F18.16)") rec_mxwls(i,(j-1)*5+3)
  text_line1 = trim(text_line1) // ", " // trim(Adjustl(suff))
  write (suff,"(F18.16)") rec_mxwls(i,(j-1)*5+4)
  text_line1 = trim(text_line1) // ", " // trim(Adjustl(suff))
  write (suff,"(F18.16)") rec_mxwls(i,(j-1)*5+5)
  text_line1 = trim(text_line1) // ", " // trim(Adjustl(suff)) 
  write (suff,"(F18.16)") rec_mxwls(i,(j-1)*5+6)
  text_line1 = trim(text_line1) // ", " // trim(Adjustl(suff)) 
 end do 
  !!! text_line1 = trim(text_line1) // "/"
 ! next we print array using nnn entries in a row
 write(15,"(A)") trim(text_line1)
end do
close (15)

end subroutine WriteRecMxwls_kMxsl0D




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FlyTrap0DKmxwl_DGV_MPI 
!
! This subroutine will be executed on idling processors. It contains pieces of MPI parallel algorithms that need to be 
! executed on processors with irank>0
!
! The work of this subroutine is determined by action codes. The subroutine is an infinite loop, 
! The first operation in the loop is a broadcast from irank=0. This broadcast will send an action code that 
! determines what the subroutine does next. The subroutine will perform action that are specified by the give action code and will return to the 
! expectation of broadcast.
!
! action codes that are currently implemented:
! 
! -777 exit code
! 401 -- evaluation of collision operator on the secondary mesh using full Boltzmann
! 402 -- evaluation of the collision operator on the seocondary mesh using decomposition f=fm+df
!
!!!!!!!!!!!!!
subroutine FlyTrap0DKmxwl_DGV_MPI

use DGV_commvar, only: nodes_uII,nodes_vII,nodes_wII
use DGV_distributions_mod
use DGV_collision_mod


!!!
include 'mpif.h' ! THIS IS TO ENABLE THE MPI ROUTINES
!!!     


!!!!!!!!!!!!!!!!!!!!!!!
integer (I4B) :: mm,nn ! scrap variables
real (DP), dimension(:), allocatable :: fII,fcolII,dpbuff,fcol1II,fmII
integer :: loc_alloc_stat
real (DP), dimension(:), allocatable :: fmxwls,sample_u,sample_v,sample_w,coll_int ! scrap variable to keep a copy fo the solution

!!!!!!!!!!!!!!!!!!!!!!!
integer :: irank,ierr          ! scrap variable to keep the rank of the processor
integer :: mpi_dgv_action ! scrap variable to keep the action code
integer, dimension(1:1) :: ibuff
!!!!!!!!!!!!!!!!!!!!!!
call mpi_comm_rank(MPI_COMM_WORLD,irank,ierr) ! check what processor the program is on...

mpi_dgv_action=0
do while (mpi_dgv_action /= -777)
 ibuff=0  
 call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
 if (ierr /= 0 ) then 
  print *,"FlyTrap0DKmxwl_DGV_MPI:  slave processor", irank,&
         "MPI boradcast of mpi_dgv_action from proc 0 returned error", ierr
  stop
 end if
 mpi_dgv_action=ibuff(1) 
 !!!!!!!!!!!!!!!!! MAIN FORK
 select case (mpi_dgv_action)
  case (402) ! this action code corresponds to evaluation of the collision operator using decomposition formulation
   ! now that the code was broadcasted, the slave processors will need to recieve a copy of the 
   ! solution and the points where to evaluate the collision operator. The solution is contained 
   ! in the array fmxwls(:) where macroparameters of all Maxwellians are given. 
   ! the points where the solution needs to be evaluated is given in the arrays sample_u, sample_v, sample_w 
   ! all this information needs to be passed to the slave processors. 
   !
   ! we will receive it in four steps: (1)  - send nuymber of entries of fmxwls, then (2) send maxwells
   ! similarly, (3) send number of entries in sample_u,sample_v,sample_w, then send (4) send the arrays 
   ! (1) receive the size of the solution array
   call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   if (ierr /= 0 ) then 
    print *,"FlyTrap0DKmxwl_DGV_MPI: (0) slave processor", irank, &
           "MPI broadcast of size(fmxwls,1) from proc 0 returned error", ierr
    stop
   end if
   nn = ibuff(1)
   allocate (fmxwls(nn), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then
    print *, "FlyTrap0DKmxwl_DGV_MPI: Allocation error for (drbuffer) (0) slave processor", irank
    stop
   endif
   ! (2) receive the solution array
   call mpi_bcast (fmxwls,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   if (ierr /= 0 ) then 
    print *,"FlyTrap0DKmxwl_DGV_MPI: slave processor", irank, &
           "MPI boradcast of fmxwls from proc 0 returned error", ierr
    stop
   end if    
   ! (3) receive the number of samples
   call mpi_bcast (ibuff,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
   if (ierr /= 0 ) then 
    print *,"FlyTrap0DKmxwl_DGV_MPI: slave processor", irank, &
           "MPI boradcast of size(sample_u,1)*3 from proc 0 returned error", ierr
    stop
   end if    
   nn = ibuff(1) ! this is the number of samples*3
   ! (4) receive  the sample points, 
   allocate (dpbuff(nn),sample_u(nn/3),sample_v(nn/3),sample_w(nn/3),stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then
    print *, "FlyTrap0DKmxwl_DGV_MPI: Allocation error for (drbuffer) (0) slave processor", irank
    stop
   endif
   call mpi_bcast (dpbuff,nn,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
   if (ierr /= 0 ) then 
    print *,"kMaxwl_SpatialOper0D_MPI: slave processor", irank, & 
            "MPI boradcast of sample_u/_v/_w from proc 0 returned error", ierr
    stop
   end if    
   nn = nn/3 
   sample_u(1:nn) = dpbuff(1:nn)
   sample_v(1:nn) = dpbuff(nn+1:2*nn)
   sample_w(1:nn) = dpbuff(2*nn+1:3*nn)
   deallocate(dpbuff)
   !!!!! The solution and the points to evaluate the collision operator has been transmitted. 
   !!!!! Collision operator will be now be computed on othe transmitted points and send back to the 
   !!!!! main processor./ 
   nn = size(sample_u,1)
   allocate (dpbuff(nn),coll_int(nn), stat=loc_alloc_stat)
   if (loc_alloc_stat >0) then
    print *, "kMaxwl_SpatialOper0D_MPI: Allocation error for (drbuffer) 1"
    stop
   endif
   coll_int = 0 
   !!!!!!! Evaluation of the collision operator. 
   call EvalCollArryDecompS1Pts_MPI_DGV (sample_u,sample_v,sample_w,coll_int,fmxwls)  
   !!!!!!!  
   deallocate(sample_u,sample_v,sample_w,fmxwls)
   dpbuff = 0 ! Make a empty array 
   call mpi_allreduce(coll_int,dpbuff,nn,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierr)  ! this will be matched in the "fly trap" subroutine
   if (ierr /= 0 ) then 
    print *,"GetRelaxRates1Donecell_DGV_MPI: (0) slave processor", irank,&
           "MPI allreduce of coll_int returned error", ierr
    stop
   end if
   deallocate(dpbuff,coll_int)
   ! done here
  case (-777)
   print *,"FlyTrap0DKmxwl_DGV_MPI: slave processor", irank,&
              "end of work signal recieved. Exit now" 
  case default
   print *,"FlyTrap0DKmxwl_DGV_MPI: slave processor", irank,&
         "Exit Error: unknown action case, mpi_dgv_action=", mpi_dgv_action
           stop
 end select   
end do 
end subroutine FlyTrap0DKmxwl_DGV_MPI


end module DGV_kMaxwl_tools