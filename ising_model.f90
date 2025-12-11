!====================================================================================================================================================
! Name        : ising_model.f90
! Author      : Pedro de Dios
! Description : This code computes observables (<E>, C_ν, <|M|>, Xv) of the 2D ising model via Markov-Monte Carlo method using the Metropolis algorithm.
!====================================================================================================================================================







module ising_utils  !periodical boundary conditions on the 2D lattice 

    
    implicit none
    contains
    
    integer function shiftup(num, NS)
        implicit none
        integer, intent(in) :: num, NS
        shiftup = num + 1
        if (num == NS) shiftup = 1
    end function shiftup

    integer function shiftdown(num, NS)
        implicit none
        integer, intent(in) :: num, NS
        shiftdown = num - 1
        if (num == 1) shiftdown = NS
    end function shiftdown



end module ising_utils








subroutine metropolis_ising(temperature, XS, YS, sweeps, thermalization, measuring_gap)

    use ising_utils 

    implicit none
    integer, intent(in) :: XS, YS, sweeps, thermalization, measuring_gap
    real(8), intent(in) :: temperature
    real(8) :: energy_density, energy_error, specific_heat
    real(8) :: abs_mag, magnetization_error, mag_susceptibility

    integer :: i,j,k,m,up_x,down_x,up_y,down_y,count, VOL
    real(8) :: beta, delta,H,H_prima,upsidedown,p,t,r, near_neighbours
    real(8) :: energy, magnetization
    real(8) :: energy_measure, magnetization_measure
    real(8) :: energy_square, magnetization_square
    real(8) :: energy_avg, energy_sqr_avg
    real(8) :: magnetization_avg, magnetization_sqr_avg
    real(8) :: energy_std, magnetization_std
    real(8) :: energy_error_local, magnetization_error_local
    real(8), allocatable :: lattice(:,:), total_energy(:), total_magnetization(:)


    allocate(lattice(XS,YS))
    allocate(total_energy(sweeps))
    allocate(total_magnetization(sweeps))


    beta = 1.0d0 / temperature
    VOL = XS*YS

   

    ! Initialize lattice (hot start)
    do i=1,XS
        do j=1,YS
            call random_number(r)
            if(r >= 0.5d0) then
                lattice(i,j) = 1.0d0
            else
                lattice(i,j) = -1.0d0
            endif
        end do
    end do

    ! --- Metropolis sweeps ---
    do m=1,sweeps
        do i=1,XS
            do j=1,YS

                up_x = shiftup(i, XS)
                down_x = shiftdown(i, XS)
                up_y = shiftup(j, YS)
                down_y = shiftdown(j, YS)

                upsidedown = -lattice(i,j)

                near_neighbours = lattice(i,up_y)+lattice(i,down_y)+lattice(up_x,j)+lattice(down_x,j)

                H = -lattice(i,j)*near_neighbours
                H_prima = -upsidedown*near_neighbours

                delta = H_prima - H

                if(H_prima < H) then
                    lattice(i,j) = upsidedown
                else
                    p = exp(-beta*delta)
                    call random_number(t)
                    if(t < p) lattice(i,j) = upsidedown
                endif

            end do !end i
        end do !end j

        ! Compute energy and magnetization after sweep
        energy = 0.0d0
        magnetization = 0.0d0

        do i=1,XS
            do j=1,YS

                up_x = shiftup(i, XS)
                down_x = shiftdown(i, XS)
                up_y = shiftup(j, YS)
                down_y = shiftdown(j, YS)

                energy = energy - lattice(i,j)*(lattice(i,up_y)+lattice(i,down_y)+lattice(up_x,j)+lattice(down_x,j))
                magnetization = magnetization + lattice(i,j)

            end do
        end do

        total_energy(m) = energy
        total_magnetization(m) = magnetization


        write(10,*) m,	beta*energy  
        write(11,*) m,	magnetization  

    end do  !end sweeps


    ! --- Measurements ---
    count = 0
    energy_measure = 0.0d0
    magnetization_measure = 0.0d0
    energy_square = 0.0d0
    magnetization_square = 0.0d0

    do k = thermalization, sweeps

        if(mod(k, measuring_gap) == 0) then

            count = count + 1
            energy_measure = energy_measure + total_energy(k)
            magnetization_measure = magnetization_measure + abs(total_magnetization(k))
            energy_square = energy_square + total_energy(k)**2
            magnetization_square = magnetization_square + total_magnetization(k)**2


        end if
    end do


    energy_avg = energy_measure/count
    energy_sqr_avg = energy_square/count
    specific_heat = (beta**2/VOL)*(energy_sqr_avg - energy_avg**2)
    energy_density = energy_avg/VOL

    magnetization_avg = magnetization_measure/count
    magnetization_sqr_avg = magnetization_square/count
    mag_susceptibility = (magnetization_sqr_avg - magnetization_avg**2)/VOL
    abs_mag = magnetization_avg/VOL


    ! Standard deviation (error)
    energy_std = 0.0d0
    magnetization_std = 0.0d0
    do m=thermalization,sweeps
        if(mod(m, measuring_gap) == 0) then
            energy_std = energy_std + (total_energy(m)/VOL - energy_density)**2
            magnetization_std = magnetization_std + (total_magnetization(m)/VOL - abs_mag)**2
        end if
    end do

    energy_error = sqrt(energy_std/count)
    magnetization_error = sqrt(magnetization_std/count)



    write(12,*) temperature, energy_density, energy_error, specific_heat
    write(13,*) temperature, abs_mag, magnetization_error, mag_susceptibility

    deallocate(lattice)
    deallocate(total_energy)
    deallocate(total_magnetization)



end subroutine metropolis_ising


















program run_ising  ! Run the Ising model (metropolis_ising subrotuine) for different temperature values 

    use ising_utils        
    implicit none

    ! Parameters
    integer :: XS, YS, sweeps, thermalization, measuring_gap, T_steps
    real(8) :: T_i, T_f, delta_T, temperature
    integer :: n

    ! --- Set free parameters ---
    XS = 12  !X dimension
    YS = 12   !Y dimension

    sweeps = 40000          !# of sweeps
    thermalization = 10000   !# of sweeps before measuring observables (aprox 40k )
    measuring_gap = 15      ! space between consecutive measuraments after thermalization 



    T_i = 1.0d0            ! Initial temperature
    T_f = 5.0d0            ! Final temperature
    delta_T = 0.2d0        ! Increase in temperature

    T_steps = int ((T_f - T_i) / delta_T)



   open(unit=10,file='energy', status='replace')  ! divide data in this file in segments of size = #sweeps to observe thermalization for different temperatures
   open(unit=11,file='magnetization', status='replace') ! same for this file with magnetization data across sweeps and temperatures
  
   open(unit=12,file='energy_measurament', status='replace')  ! Measuring data of energy (critical point at T_c)!
   open(unit=13,file='magnetization_measurament', status='replace') ! Measuring data of magnetization (critical point at T_c)!
 

    ! --- Loop over temperatures ---
    do n = 0, T_steps

        temperature = T_i + n*delta_T
        print *, 'Running temperature =', temperature

        call metropolis_ising(temperature, XS, YS, sweeps, thermalization, measuring_gap)


    end do


end program run_ising










