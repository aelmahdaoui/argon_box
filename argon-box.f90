!argon gas in a box simulation, molecular dynamics.
!compile with: gfortran argon-box.f90 $(pkg-config --cflags --libs plplotd-f95)

!the cubic geometry sides of length = L
!initial positions initialized according to fcc lattice structure	
!number of fcc cells per cartesian dimension = Ncell, 
!number of particles is N, (4 particles per cube)		

program argon_box
	use plplot
    implicit none
	
	integer, parameter :: N_cell_dim = 2, N_cell = N_cell_dim**3, N_part = N_cell*4
	real(8), parameter :: L_side = 1d0, T = 90000000, m = 1d0
	real(8), parameter :: s = 1d0, e = 1d0, r_cut = L_side, dt = 1d-7 ! lennard jones potential
	real(8), parameter :: t_stop = 1d0
	real(8), parameter :: Kb = 1 	!constants	
	integer, parameter :: N_avSteps = 100 ! #steps used for ensemble average
	integer :: i,j,k,l,n, step!iteration variables
	real(8):: xs(2) !two random numbers
	real(8) :: time, r, r_vec(3), F(3), dF(3), U(N_part), v_2(N_part), H, Ekin, P, F_R(N_avSteps), Temp
	real(8), dimension(1:3, 1:N_part) :: pos, vel 	! particle system arrays	
	
	
	!initial positions of all particles according to an fcc lattice structure
	print *, "initializing initial particle positions"
	n = 0
	do i = 1,N_cell_dim
	do j = 1,N_cell_dim
	do k = 1,N_cell_dim
		do l = 1,4
			n = n + 1
			pos(:,n) = L_side/N_cell_dim*((/ i-1, j-1, k-1 /) + fcc_cell(l))
!			print *, "particle:", n, "/",  N_part, "position:", pos(:,n)
		end do
	end do
	end do
	end do
		
	!initial particles velocities according to the maxwell distribution	
	! in the maxwell boltzman distribution each velocity component is normally distributed:
	! f(v) = sqrt(m/(2*PI*Kb*T)) * exp(-(v**2)/2 *m/(*Kb*T))
	! sigma**2 = Kb*T/m, and zero mean	
	print *, "initializing initial particle velocities"
	call init_random_seed
	do n = 1,N_part
		do i = 1,3
			CALL RANDOM_NUMBER(xs(1))
			CALL RANDOM_NUMBER(xs(2))						
			vel(i,n) = sqrt(Kb*T/m) * box_muller(xs)
		end do
!		print *,"particle:", n, "/",  N_part, "velocity:", VEL(:,n)
	end do

! linear aproximation: x = x + v*dt, dv= F/m*dt, read verlets paper to find out whether higher order terms are needed.
! also the order in whidh dx and dv applied to the system might be significant.
!time evolution for particles in lennard jones potential: U = 4*e*(s/r)**12-(s/r)**6,
!Fij = -du/dx = -du/dr*dr/dx = e*(48*s**12/r**13 - 6*s**6/r**7) * x/r, 
!r = sqrt(x**2+y**2+z**2)	
print *, "calculating time evolution"
time = 0d0
step = 0
call plot_init(0d0, L_side,0d0, L_side,0d0, L_side) 

!Time integration
do while (time < t_stop)
	call plot_points(pos)
	time = time + dt	
	step = step + 1	
	
	F_R(1:(N_avSteps-1)) = F_R(2:N_avSteps) !calculation of ensamble average as time average in virial theorem, pressure calculation
	F_R(N_avSteps) = 0
	
		! Force calculation	
		do n = 1,N_part 
		U(n) = 0 !potential
		F = 0
		do i = 	1,N_part !integrate over all particles inside box except i = n					
			do j = -1, 1 
			do k = -1, 1 !periodic boundary condition for potential forces..
			do l = -1, 1				
				if (n/=i)then !(.not.(n=i .and. (/j,k,l/)=(/0,0,0/)))	
					r_vec = (/pos(1,n)-pos(1,i), pos(2,n)-pos(2,i), pos(3,n)-pos(3,i)/) + L_side*(/j,k,l/)
					r = sqrt(dot_product(r_vec, r_vec))  
					if (r<r_cut) then
						! uitwerken op papier virial theorem.. potentiaal kracht in termen van elkaar..
						dF = e*(48*s**12/r**14 - 6*s**6/r**8) * r_vec
						F = F + dF 		
						U(n) = U(n) + 5d-1*(4*e*(s/r)**12-(s/r)**6)
						F_R(N_avSteps) = F_R(N_avSteps) + dot_product(r_vec, dF)	 !ensemble average
					end if
				end if		
			end do 
			end do 
			end do
		end do			
		vel(:,n) = vel(:,n) + F/m*dt	
		v_2(n) = dot_product(vel(:,n),vel(:,n))
!		print *,"particle:", n, "/",  N_part, "velocity:", (VEL(i,n), i=1,3)
	end do
	
	! Postion calculation
	do n = 1,N_part 
		pos(:,n) = pos(:,n) + vel(:,n)*dt			
		do i = 1,3 !implements periodic boundary conditions
			if (pos(i,n) < 0d0) then 
				pos(i,n) = pos(i,n) + L_side
			else if (pos(i,n) > L_side) then
				pos(i,n) = pos(i,n) - L_side 
			end if
		end do		
!		print *, "particle:", n, "/",  N_part, "position:", pos(:,n)
	end do	
	
	Ekin = sum(m/2*V_2)
	H = sum(U) + Ekin
	Temp = 2*Kb/3*Ekin/N_part
	P = N_part/(L_side**3)*Kb*T*(1 + 1/(6*Kb*T*N_part)* sum(F_R)/N_avSteps) 	! + correction cuttoff,.

	
	print *, step, "t = ", time, "H =", H,  "T =", Temp, "P =", P,  "vel1 =", vel(1,1), "pos1 =", pos(1,1)	
		
end do	

call plot_end

contains
	
	function fcc_cell(i) result(output)
		implicit none
		integer, intent(in) :: i
		real(8) :: output(3)
		! face centered cubic unit cell with basis particle positions:
		real(8), dimension(1:3), parameter :: &
			fcc_part1 = (/0d0,  0d0,  0d0/),  &
			fcc_part2 = (/0d0,  5d-1, 5d-1/), &
			fcc_part3 = (/5d-1, 0d0,  5d-1/), &
			fcc_part4 = (/5d-1, 5d-1, 0d0/)		
		real(8), dimension(1:3,1:4), parameter :: &
		R_cell = reshape( (/fcc_part1, fcc_part2, fcc_part3, fcc_part4/), (/3,4/))
		output = R_cell(:,i)
		
	!	print *, "face centered cubic unit cell"
	!	do i = 1,3		!		print *, (R_cell(i,j), j=1,4)
	!	end do
		
	end function fcc_cell
	
	! adapted from ICCP coding-notes
	function box_muller(xs) result(zs)
	  implicit none
	  real(8) :: pi = 4*atan(1d0), zs !(2)
	  real(8), intent(in) :: xs(2)
	  
	  zs = sqrt(-2d0*log(xs(1)))*cos(2*pi*xs(2))
	  !zs(2) = sqrt(-2d0*log(xs(1)))*sin(2*pi*xs(2))
	end function box_muller
	
	! copied from ICCP coding-notes
	subroutine init_random_seed()
	  implicit none
	  integer, allocatable :: seed(:)
	  integer :: i, n, un, istat, dt(8), pid, t(2), s
	  integer(8) :: count, tms

	  call random_seed(size = n)
	  allocate(seed(n))
	  open(newunit=un, file="/dev/urandom", access="stream",&
	       form="unformatted", action="read", status="old", &
	       iostat=istat)
	  if (istat == 0) then
	     read(un) seed
	     close(un)
	  else
	     call system_clock(count)
	     if (count /= 0) then
	        t = transfer(count, t)
	     else
	        call date_and_time(values=dt)
	        tms = (dt(1) - 1970)*365_8 * 24 * 60 * 60 * 1000 &
	             + dt(2) * 31_8 * 24 * 60 * 60 * 1000 &
	             + dt(3) * 24 * 60 * 60 * 60 * 1000 &
	             + dt(5) * 60 * 60 * 1000 &
	             + dt(6) * 60 * 1000 + dt(7) * 1000 &
	             + dt(8)
	        t = transfer(tms, t)
	     end if
	     s = ieor(t(1), t(2))
	     pid = getpid() + 1099279 ! Add a prime
	     s = ieor(s, pid)
	     if (n >= 3) then
	        seed(1) = t(1) + 36269
	        seed(2) = t(2) + 72551
	        seed(3) = pid
	        if (n > 3) then
	           seed(4:) = s + 37 * (/ (i, i = 0, n - 4) /)
	        end if
	     else
	        seed = s + 37 * (/ (i, i = 0, n - 1 ) /)
	     end if
	  end if
	  call random_seed(put=seed)
	end subroutine init_random_seed
	
	!adapted from coding notes
	subroutine plot_init(xmin,xmax, ymin,ymax,zmin,zmax)
		implicit none
		real(8), intent(in) :: xmin,xmax,ymin,ymax,zmin,zmax
		
	    call plsdev("xcairo")
		call plinit()
	    !call plparseopts(PL_PARSE_FULL)
		
		call pladv(0)
		call plvpor(0d0, 1d0, 0d0, 1d0)
		call plwind(-1d0, 1d0, -2d0 / 3, 4d0 / 3)
		call plw3d(1d0, 1d0, 1d0, xmin, xmax, ymin, ymax, &
		             zmin, zmax, 45d0, -45d0)
			
	end subroutine plot_init
	
	subroutine plot_end
		call plspause(.false.)
		call plend()
	end subroutine 
	
	!copied from coding notes
    subroutine plot_points(xyz)
	implicit none
		
      real(8), intent(in) :: xyz(:, :)

      call plclear()
      call plcol0(1)
      call plbox3("bnstu", "x", 0d0, 0, "bnstu", "y", &
                    0d0, 0, "bcnmstuv", "z", 0d0, 0)
      call plcol0(2)
      call plpoin3(xyz(1, :), xyz(2, :), xyz(3, :), 4)
      call plflush()
    end subroutine 

end program

