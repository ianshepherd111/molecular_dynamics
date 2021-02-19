program lj

  implicit none

  integer, parameter :: dp = SELECTED_REAL_KIND(15,307)
  integer, parameter :: large = SELECTED_INT_KIND(18)
  real (dp), parameter :: boltz_const = 1
  real (dp), parameter  :: pi=3.1415926535897932384626433
  integer, parameter :: dimen = 3

  !in/out file parameter no.s
  integer, parameter :: screen_unit = 6
  integer, parameter :: unit = 20, unit2 = 21, unit3 = 22, unit4 = 23 , unit5 = 24, unit_final_coords_txt = 35

  integer :: m, i, j, k, N, moves, out_no, io

  real (dp), allocatable, dimension (:,:) :: coord_positions, force_i, v, old_force_i

  !program running time variables
  real (dp) :: start_time, stop_time

  real (dp) :: sigma
  real (dp) :: epsilon
  real (dp) :: mass
  real (dp) :: time_step
  real (dp) :: potential_energy, kinetic_energy, target_temperature, real_temperature
  real (dp) :: box_l, r_cut, r_cut_sq
  real (dp) :: anderson_freq
  real (dp) :: density
  real (dp) :: pressure, pressure_component, pressure_component_mean
  real (dp) :: P_corr, lj_pot_corr, lj_trunc, lj_trunc_div_4_epsilon
  real (dp) :: potential_energy_per_N



  !path to save files to
  character*(*), parameter :: directory = 'output'
  character*(*), parameter :: path_out = './' // directory // '/'

  !command to make the directory
  character*(*), parameter :: makedirectory = 'mkdir ' // directory

  ! filename/path of the input/control files
  character*(*), parameter :: directory_in = 'input'
  character*(*), parameter :: path_in = './' // directory_in // '/'

  character(len=100) :: param_file, input_coord_file
  input_coord_file = path_in//'input.txt'
  param_file = path_in//'params.txt'




  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !     Program
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!



  !compile with fortran2008 standard for support for the execute_command_line
  call execute_command_line (makedirectory)

  call cpu_time (start_time)

  call read_param_file ( param_file, screen_unit, moves, time_step, target_temperature, sigma, &
  epsilon, mass, out_no, density, r_cut, anderson_freq )



  allocate ( coord_positions (dimen,1) )

  !read-in coord file
  open (unit, file=input_coord_file, status='old')
      !n = no. of lines - start at -1 beacuse of while loop structure below
      N = -1
      io = 0
      do while ( io == 0 )
          N = N + 1
          read(unit,*,iostat=io) coord_positions
      end do

  deallocate ( coord_positions )
  allocate ( coord_positions (dimen,N), v (dimen,N), force_i (dimen,N), old_force_i (dimen,N) )

  rewind(unit)

      read (unit,*,iostat=io)   ( (   coord_positions(i,j),   i=1,dimen   ),   j=1,N  )

  close(unit)


  !output initial coords/params
  open (unit,file=path_out//'input_coords.xyz', status='unknown')

    write (unit, *) N
    write (unit,*) ''

    do j = 1, N
        write (unit,*)   'Ar',   coord_positions(:,j)
    end do

  close (unit)

  box_l = (  (real (N, dp))/ (density) ) **(1.0_dp/3.0_dp)
  r_cut_sq = r_cut**2



  call pressure_and_lj_pot_correction (pi, sigma, epsilon, r_cut, density, P_corr, lj_pot_corr, lj_trunc )

lj_trunc_div_4_epsilon = lj_trunc / (epsilon * 4.0_dp)

  call init_random_seed()


  open (unit,file=path_out//'input_params.txt', status='unknown')
      write (unit, *) 'N', N
      write (unit, *) 'dimen', dimen
      write (unit, *) 'moves', moves
      write (unit, *) 'time_step', time_step
      write (unit, *) 'out_no', out_no
      write (unit, *) 'sigma', sigma
      write (unit, *) 'epsilon', epsilon
      write (unit, *) 'mass', mass
      write (unit, *) 'temperature', target_temperature
      write (unit, *) 'density', density
      write (unit, *) 'box size', box_l
      write (unit, *) 'cutoff distance', r_cut
      write (unit, *) 'pressure correction', P_corr
      write (unit, *) 'truncation error', lj_trunc

  close (unit)

  call initiate_velocities ( v, boltz_const, mass, target_temperature)

  call pbc (coord_positions, box_l)

  ! call force_calc ( coord_positions, r_cut_sq, sigma, epsilon, N, dimen, force_i, potential_energy, box_l, pressure_component, &
  ! lj_trunc_div_4_epsilon )

  call force_calc ( coord_positions, sigma, epsilon, N, dimen, force_i, potential_energy, box_l, pressure_component, &
  lj_trunc_div_4_epsilon )

  ! pressure = (   (boltz_const * target_temperature * N) + (pressure_component / dimen)   )    /   (box_l**3)
  ! write (screen_unit,*) 'initial pressure', pressure

  potential_energy_per_N = 0_dp
  pressure_component_mean = 0_dp


  !write out initial values
  open (unit,file=path_out//'input_data.dat', status='unknown')
      write (unit,*) 'Force'

      do j = 1, N
          write (unit,*)   force_i(:,j)
      end do
      write (unit,*) 'Potential energy'
      write (unit,*) potential_energy
  close (unit)


  !open files to write out to
  open (4,file=path_out//'traj.xyz', position='append', status='replace')
  open (3,file=path_out//'force.dat', position='append', status='replace')
  open (11,file=path_out//'energy.dat', position='append', status='replace')
  open (12,file=path_out//'temperature.dat', position='append', status='replace')

  write (11, '(a,25x,a,25x,a)' ) 'energy', 'kinetic', 'potential'

  !start of MD loops
  do m = 1, moves

    ! Velocity Verlet Algorithm
    coord_positions = coord_positions + (v * time_step) + ( ( force_i * (time_step**2)) / (mass + mass) )
    old_force_i = force_i

    call pbc (coord_positions, box_l)

    call force_calc ( coord_positions, sigma, epsilon, N, dimen, force_i, potential_energy, box_l, pressure_component, &
    lj_trunc_div_4_epsilon )
    !
    ! call force_calc ( coord_positions, r_cut_sq, sigma, epsilon, N, dimen, force_i, potential_energy, box_l, pressure_component, &
    ! lj_trunc_div_4_epsilon )

    v = v + (     (  ( force_i+old_force_i ) * time_step ) / (mass + mass)    )

    call anderson (v, boltz_const, mass, target_temperature, anderson_freq, time_step, m)


    kinetic_energy = real (0, dp)

    do i = 1, N
        kinetic_energy = kinetic_energy + ( sum (  (v(:,i) )**2  )   )*mass*(real (0.5, dp))
    end do

    real_temperature = (kinetic_energy*(real (2, dp))) / ( (real (3*N, dp))*boltz_const)

    potential_energy_per_N = potential_energy_per_N + potential_energy

    pressure_component_mean = pressure_component_mean +   ( (pressure_component - pressure_component_mean) / (real(m, dp))  )

    ! pressure = (   (boltz_const * target_temperature * N) + (pressure_component_mean / dimen)   )    /   (box_l**3)
    ! write (screen_unit,*) pressure



    ! File write outs
    if (modulo(m,out_no) == 0 ) then

        write (4, *) N
        write (4,*) m

        do j = 1, N
           write (4,*) 'Ar', coord_positions(:,j)
        end do

        write (3, *) N
        write (3,*) m

        do j = 1, N
            write (3,*) force_i(:,j)
        end do

        write (11, *) (kinetic_energy + potential_energy), kinetic_energy, potential_energy
        write (12, *) real_temperature

    end if


  end do
  !!!! End of Move loop

  !file housekeeping
  close (3)
  close (4)
  close (11)
  close (12)


  !convert pressure_component to pressure
  pressure = (   (boltz_const * target_temperature * N) + (pressure_component_mean / dimen)   )    /   (box_l**3)

  potential_energy_per_N = potential_energy_per_N / (real ((moves*N),dp) )


  open (40,file=path_out//'final_coords.xyz', status='unknown')

      write (40,*) N
      write (40,*) ''

      do j = 1, N
          write (40,*) 'Ar',   coord_positions(:,j)
      end do

  close (40)


  !create an output coordinate file that has the same formatting as the input.txt file
  open (unit_final_coords_txt,file=path_out//'final_coords.txt', status='unknown')

      do j = 1, N
          write (unit_final_coords_txt,*)  coord_positions(:,j)
      end do

  close (unit_final_coords_txt)



  deallocate (  coord_positions, v, force_i, old_force_i )

  write (screen_unit,'(/,a)') 'Results:'

  write (screen_unit,*) 'pressure', pressure
  write (screen_unit,*) 'P_corr', P_corr
  write (screen_unit,*) 'corrected pressure', pressure+P_corr

  write (screen_unit,*) 'potential_energy_per_N', potential_energy_per_N
  write (screen_unit,*) 'lj_pot_corr', lj_pot_corr
  write (screen_unit,*) 'final potential energy', potential_energy


  call cpu_time(stop_time)
  write (screen_unit,*) 'Time taken:', stop_time - start_time, "seconds"


  contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !
  !
  !     Subroutines and functions
  !
  !
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    subroutine read_param_file( param_file, screen_unit, moves, time_step, target_temperature, sigma, &
      epsilon, mass, out_no, density, r_cut, anderson_freq )
      ! new variable in needs: add to declaration, change int/real param no., a string & assignment at end of subroutine

      implicit none

      character(len=100), intent(in) :: param_file
      integer, intent (in) :: screen_unit

      ! Control file variables to read-in with this subroutine
      real (dp), intent (out) :: time_step, target_temperature, sigma, epsilon, mass, density, r_cut, anderson_freq
      integer, intent (out) :: moves, out_no

      ! Input related variables
      character(len=100) :: buffer, label
      integer :: pos, i, param_no, no_int_params, no_real_params
      integer, parameter :: unit5 = 15
      integer :: ios = 0
      integer :: iosfile = 0
      integer :: line = 0

      ! Timer related variables
      real(dp) :: start_time, stop_time

      integer, allocatable, dimension (:) :: chk
      character(len=100), allocatable, dimension (:) :: real_strings
      character(len=100), allocatable, dimension (:) :: int_strings
      real (dp), allocatable, dimension (:) :: reals
      integer, allocatable, dimension (:) :: ints

      ! Parameter type/no. variables
      no_int_params = 2
      no_real_params = 8

      param_no = no_int_params + no_real_params

      allocate( chk(param_no) )
      allocate( real_strings(no_real_params), reals(no_real_params) )
      allocate( int_strings(no_int_params), ints(no_int_params) )

      ! set "checking" array to 0 initially
      chk (:) = 0

      ! strings containing the names of variables
      int_strings(1) = 'moves'
      int_strings(2) = 'out_no'

      real_strings(1) = 'time_step'
      real_strings(2) = 'target_temperature'
      real_strings(3) = 'sigma'
      real_strings(4) = 'epsilon'
      real_strings(5) = 'mass'
      real_strings(6) = 'density'
      real_strings(7) = 'cutoff_dist'
      real_strings(8) = 'anderson_freq'


      open(unit5, file=param_file, status='old', iostat=iosfile)

      ! stop program if param file not found
      if (iosfile /= 0) then
        write (screen_unit,'(a,/,a)') 'parameter file not found', 'program terminating...'
        STOP
      end if

      ! start timer in case an error occurs or parameter file reading takes a very long time
      call cpu_time (start_time)
      call cpu_time (stop_time)

      !read first line into "buffer"
      read(unit5, '(A)', iostat=ios) buffer


      !do until reach end of file or find all variables (or an error occurs)
      do while (  (ios == 0)   .and.  ( sum(chk) < param_no ) .and. ( (stop_time-start_time) < 5_dp )   )

            !count line no.
            line = line + 1

            !remove blank space at the start of a line
            buffer = adjustl (buffer)

            !pos = place where whitepace starts (1 for blank lines)
            pos = index(buffer, ' ')

            !label excludes blank space and includes data's 'label'
            label = buffer(1:pos)

            !pos+1 to end -  holds the variable string and blank space
            buffer = buffer(pos+1:)

            !take no action if a line doesn't have the variable name listed in front of everything and then a space
            !e.g.:
            !N    7000    - correct
            !  N   7000   - correct
            !a   7000     - incorrect
            !  N7000      - incorrect

            !check line for integers
            do i = 1, no_int_params

              !if "label" matches a variable name, read value in "buffer" in
              if ( label == int_strings(i) ) then

                !get a value for the integer; if already have a value for the integer, don't read in another
                if (chk(i) == 0) then
                  read(buffer, *, iostat=ios)  ints (i)
                  if (ios == 0) then
                    chk(i) = 1
                  end if
                !otherwise set pos = 1 , so that the param counter below treats the line as blank
                else
                  pos = 1
                end if

              !change pos value to 1 for any non-variable lines
              else
                pos = 1
              end if
            end do

            !repeat for real values
            do i = 1, no_real_params

              if ( label == real_strings(i) ) then

                if (chk(i+no_int_params) == 0) then
                  read(buffer, *, iostat=ios)  reals (i)
                  if (ios == 0) then
                    chk(i+no_int_params) = 1
                  end if
                else
                  pos = 1
                end if
              else
                pos = 1
              end if
            end do

            !read next line into buffer
            read(unit5, '(A)', iostat=ios) buffer

            call cpu_time (stop_time)

      end do

      close (unit5)

      !timeout stop
      if  ( (stop_time-start_time) >= 5_dp  ) then
        write (screen_unit,'(a,/,a)') 'parameter file has taken over 5 seconds to read in', 'program terminating...'
        STOP
      end if


      !warning message if can't find/read-in all parameters
      if ( param_no /= sum(chk) ) then
        write (screen_unit,'(a,/,a,/,a,I0)') 'WARNING, could not read-in all parameters', &
        'Running calculation with unknown number', &
        'read parameter file down to line number  ', line
      else
        write (screen_unit,'(a,/)') 'all parameters read-in  -  see output/input_params for their values'
      end if

      !specific warning messages
      do i = 1, no_int_params
        if ( chk(i) /= 1) then
          write (screen_unit,'(a,3x,a)') 'Failed to read-in:', int_strings(i)
        end if
      end do

      do i = 1, no_real_params
        if ( chk(i+no_int_params) /= 1) then
          write (screen_unit,'(a,3x,a)') 'Failed to read-in a value for:', real_strings(i)
        end if
      end do

      ! assign the values read into the arrays to their respective variables
      !integer variables
      moves = ints(1)
      out_no = ints(2)

      !real variables
      time_step = reals(1)
      target_temperature = reals(2)
      sigma = reals(3)
      epsilon = reals(4)
      mass = reals(5)
      density = reals(6)
      r_cut = reals(7)
      anderson_freq = reals(8)

    end subroutine read_param_file



    subroutine init_random_seed()

      integer :: i, n, clock
      integer, dimension(:), allocatable :: seed

      call random_seed(size = n)
      allocate(seed(n))

      call system_clock(count=clock)

      seed = clock + 37 * (/ ((i - 1)**2, i = 1, n) /)
      call random_seed(put = seed)

      deallocate(seed)

    end subroutine init_random_seed



    subroutine pressure_and_lj_pot_correction (pi, sigmalj, epsilonlj, r_cut, density, P_corr, lj_pot_corr, lj_trunc )

      implicit none

      real(dp), intent (in) :: pi, sigmalj, epsilonlj, r_cut, density
      real(dp), intent (out) ::P_corr, lj_pot_corr, lj_trunc
      real (dp) :: sigma_div_rcut_3

      sigma_div_rcut_3 = (sigmalj/r_cut)**3

      ! Impulsive pressure correction
      P_corr = ( (sigma_div_rcut_3**3) - sigma_div_rcut_3 ) * (pi*epsilonlj*(sigmalj**3))  * (density**2) * (8.0_dp/3.0_dp)



      lj_pot_corr = (sigmalj**3)*pi*epsilonlj*density*(8.0_dp/3.0_dp) * &
       ( ( (sigma_div_rcut_3**3) / (3.0_dp) ) - (sigma_div_rcut_3) )



      ! ! tail corrections for pressure and energy
      ! P_corr = ( ( (2.0_dp/3.0_dp)*(sigma_div_rcut_3**3) ) - sigma_div_rcut_3 ) * &
      ! pi * epsilonlj * (sigmalj**3) * (density**2) * (16.0_dp/3.0_dp)
      !
      ! lj_pot_corr = ( ( (1.0_dp/3.0_dp)*(sigma_div_rcut_3**3) ) - sigma_div_rcut_3 ) * &
      ! pi * epsilonlj * (sigmalj**3) * density * (8.0_dp/3.0_dp)

      lj_trunc = 4*epsilonlj * ( (sigma_div_rcut_3**4) - (sigma_div_rcut_3**2) )


    end subroutine pressure_and_lj_pot_correction



    subroutine initiate_velocities ( v, boltz_const, mass, target_temp)

      real (dp), parameter  :: pi=3.1415926535897932384626433
      real (dp), intent (inout) :: v(:,:)
      real (dp), intent (in) :: mass, boltz_const, target_temp
      real (dp) :: uniform(3)
      real (dp) :: kinetic_energy, factor, real_temp, av_v
      integer :: i,j, n,m


      kinetic_energy = real (0, dp)

      m = size(v, 1)      !m = dimension no
      n = size(v, 2)      !n = No particles


      !Box-Muller - Uniform distribution to Normal distribution
      do i = 1,n

      	do j = 1,m
      		call random_number (uniform(j))
      	end do

      	v(1,i) = (  ( (real (-2, dp)) * log (uniform(1)) )**0.5 )* cos ((real (2, dp))*pi*uniform(2))
      	v(2,i) = (  ( (real (-2, dp)) * log (uniform(2)) )**0.5 )* cos ((real (2, dp))*pi*uniform(3))
      	v(3,i) = (  ( (real (-2, dp)) * log (uniform(3)) )**0.5 )* cos ((real (2, dp))*pi*uniform(1))

      end do

      !& multiply by 'factor' to get Boltzmann distribution from normal distribution

      factor = sqrt ( (boltz_const*target_temp)  /  (mass) )
      v(:,:) = factor*v(:,:)

      !correct velocities for drift
      do j = 1, m
      		av_v =	sum( v(j,:) ) / real (N, dp)		!use momentum when masses are non-uniform
          	v(j,:) = v(j,:) - av_v
      end do

      do i = 1, n
        	kinetic_energy = kinetic_energy + (    (  sum (v(:,i)**2)   )*mass*(real (0.5, dp))    )
      end do

      real_temp = ((real(2,dp))*kinetic_energy) / ( (real (3*n, dp))*boltz_const)

      !rescale velocities because wont sample whole boltz. distrib.
      v(:,:) = (sqrt (target_temp/real_temp) ) * v(:,:)

    end subroutine initiate_velocities



    subroutine anderson ( v, boltz_const, mass, target_temp, anderson_freq, time_step, m)

      real (dp), parameter  :: pi=3.1415926535897932384626433
      real (dp), intent (inout) :: v(:,:)
      real (dp), intent (in) :: mass, boltz_const, target_temp, anderson_freq, time_step
      real (dp) :: anderson_param, factor
      real (dp) :: uniform(6)

      integer :: i,j, n,m, s


      anderson_param = anderson_freq * time_step
      factor = sqrt ( (boltz_const*target_temp)  /  (mass) )


      !d = size(v, 1)      !d = dimension no
      n = size(v, 2)      !n = No particles


      do i = 1, n
          call random_number ( uniform(1) )
          if (uniform(1) > anderson_param ) then
            	call random_number (uniform(:))

              !Box-Muller -- Uniform distribution to Normal distribution
              v(1,i) = (  ( (real (-2, dp)) * log (uniform(1)) )**0.5 )* cos ((real (2, dp))*pi*uniform(2))
              v(2,i) = (  ( (real (-2, dp)) * log (uniform(3)) )**0.5 )* cos ((real (2, dp))*pi*uniform(4))
              v(3,i) = (  ( (real (-2, dp)) * log (uniform(5)) )**0.5 )* cos ((real (2, dp))*pi*uniform(6))
              !Uniform distribution to Boltzmann distribution
              v(:,i) = factor*v(:,i)
          end if
      end do

    end subroutine anderson



    subroutine pbc (r, L)

        real (dp), intent (inout) :: r(:, :)
        real (dp), intent (in) :: L
        real (dp) :: L_half


        integer :: n , m

        m = size(r, 1)      !m = dimension no
        n = size(r, 2)      !n = No particles

      L_half = (L/2.0_dp)

      ! check all particles coords/dimensions - move back inside of box - 0,0,0 centered
      ! do i = 1, m
      !     do j = 1, n
      !         if (r(i,j) > L_half ) then
      !             r(i,j) = r(i,j) - L
      !         else if ( r(i,j) <  - L_half ) then
      !             r(i,j) = r(i,j) + L
      !         end if
      !     end do
      ! end do

      !check all particles coords/dimensions - move back inside of box - 0,0,0 centred
      do i = 1, m
        do j = 1, n
                r(i,j) = r(i,j) - (     L * (nint( r(i,j) / L ))     )
        end do
      end do


    end subroutine pbc



    subroutine force_calc (r, sigma, epsilon, n, d, f, pot_energy, L, pressure_component, lj_trunc_div_4_epsilon)
    ! subroutine force_calc (r, r_cut_sq, sigma, epsilon, n, d, f, pot_energy, L, pressure_component, lj_trunc_div_4_epsilon)


      integer, intent (in) :: n, d
      real (dp), intent (in) :: r(:, :)
      ! real (dp), intent (in) :: sigma, epsilon, L, r_cut_sq, lj_trunc_div_4_epsilon
      real (dp), intent (in) :: sigma, epsilon, L, lj_trunc_div_4_epsilon
      real (dp), intent (out) :: f(:,:)
      real (dp), intent (out) :: pot_energy
      real (dp), intent (inout) :: pressure_component

      integer :: i,j, bin_no

      real (dp), dimension (d) :: del_dim

      real (dp) :: r_mod_sq, sigma_6_div_r_mod_6, force_pre, L_half


      pot_energy = 0.0_dp
      f (:,:) = 0.0_dp
      L_half = L/2.0_dp
      pressure_component = 0.0_dp

      !calculate force between every possible pair of particles (and the potential energy)
      do i = 1 , n-1
        do j = i+1 , n

          del_dim(:) = r (:,i) - r (:,j)

              !pbc conditions
              do k = 1, d
                  if (del_dim(k) > L_half ) then
                    del_dim(k) = del_dim(k) - L
                  end if
                  if (del_dim(k) < - L_half ) then
                    del_dim(k) = del_dim(k) + L
                  end if
              end do

              r_mod_sq = sum (del_dim**2)


              !cutoff condition
              ! if  (r_mod_sq  < r_cut_sq)  then

                    sigma_6_div_r_mod_6 = (sigma**6)/(r_mod_sq**3)

                    force_pre =  (sigma_6_div_r_mod_6 * (sigma_6_div_r_mod_6 - (real (0.5, dp)))) / (r_mod_sq)

                    pot_energy = (sigma_6_div_r_mod_6 * (sigma_6_div_r_mod_6 - 1.0_dp) ) + pot_energy
                    ! - lj_trunc_div_4_epsilon


                    !sum forces on i, take negative sum on j because of N3L.
                    f(:,i) = f(:,i) +  (force_pre * del_dim(:) )
                    f(:,j) = f(:,j) -  (force_pre * del_dim(:) )

                    !dot product of fij,rij summed over j>i
                    pressure_component = pressure_component + ( force_pre * ( dot_product( del_dim(:), del_dim(:)) )  )

              ! end if

        end do
      end do



      !multiply everything by a factor which wasn't included in each loop
      f = f * epsilon * 48.0_dp
      pot_energy = pot_energy * epsilon * 4.0_dp

      pressure_component = pressure_component * epsilon * 48.0_dp

    end subroutine force_calc



end program lj
