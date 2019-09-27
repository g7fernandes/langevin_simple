! Simple langevin dynamics with brownian motion, no interaction between partciles and custom temperature distribution
! conditions on the cold side will be the reference

module mod1
    !     use linkedlist
        integer, parameter:: dp=kind(0.d0)                   ! double precision
        ! vetor de strings
        type string
            character(len=:), allocatable :: str
            integer :: quanty                   
        end type string
end module mod1

module fisica

    contains 

    function St_slip(Kn0,T_loc, N) result(Cc)
        ! allen and raabe 1982
        use mod1 

        integer :: N
        real(dp) :: Kn0
        real(dp), dimension(:) :: T_loc
        real(dp), dimension(N) :: Cc

        Cc = 1 + (Kn0 * T_loc)*(1.257 + 0.4 * exp(-1.1 * Kn0 * T_loc))

    end function St_slip

    function densitycr(rhof,rho,T_loc,N) result(drho)
        ! corrige a densidade
        ! rho é adimensional
        use mod1

        real(dp) :: T_loc(:), rhof, rho
        real(dp), dimension(N) :: dp 

        drho = (rhof*rho - rhof/T)/(rhof*rho - rhof)

    end function densitycr


    subroutine comp_F(F,x,v,T_loc,St0,Pe0,Kn0,N,dt)
        use mod1

        real(dp), dimension(:,:), intent(inout) :: F
        real(dp), dimension(:,:), intent(in) :: x,v
        real(dp), dimension(:), intent(in) :: T_loc
        real(dp), intent(in) :: St0, Pe0, Kn0, dt
        integer :: N
        real(dp), dimension(N,2) :: xi 


        ! Gera numero aleatório e normaliza
        call random_seed()
        call random_number(xi)
        xi = xi*2 -1
        xi = xi / sqrt(x1(:,1)**2 + x1(:,2)**2)

        if (Pe0 > 0) then
            F = (1/(St/St_slip(Kn0,T_loc,N))) * (- v + GField + sqrt(6/(Pe0*dt)) )
        else

        end if

    end subroutine comp_F

    subroutine comp_v(F,v,dt)
        use mod1

        real(dp), dimension(:,:), intent(inout) :: v
        real(dp), dimension(:,:), intent(in) :: F
        real(dp), intent(in) :: dt 

        v = v + F*dt 

    end subroutine comp_v


    subroutine  comp_x(x,v,dimx,dimy,N,dt,wall)
        use mod1 
        
        real(dp), dimension(:,:), intent(inout) :: v,x
        real(dp), intent(in) :: dt,dimx,dimy
        integer, intent(in) :: N 
        character :: wall(4)

        x = x + v*dt 

        ! condições de contorno

        ! elastico 
        if (wall(3:4) == 'ee') then 
            do i = 1,N 
                if (x(i,1) > dimx) then 
                    x(i,1) = 2*dimx - x(i,1)
                    v(i,1) = -v(i,1)
                else if (x(i,1) < 0) then 
                    x(i,1) = - x(i,1)
                    v(i,1) = - v(i,1)
                end if
            end do 
        end if

        if (wall(1:2) == 'ee') then 
            do i = 1,N 
                if (x(i,2) > dimy) then 
                    x(i,2) = 2*dimx - x(i,2)
                    v(i,2) = - v(i,2)
                else if (x(i,1) < 0) then 
                    x(i,2) = - x(i,2)
                    v(i,2) = - v(i,2)
                end if
            end do 
        end if   

    end subroutine comp_x

end module fisica

program main

    use mod1 
    use m_config

    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use, intrinsic :: iso_fortran_env, only: real32

    implicit none

    ! Variables
    real(dp) :: Pe0, St0, Th, Tc, Kn, dt, t_fim, t
    real(dp),allocatable, dimension(:,:) :: v, x
    character :: wall(4)
    ! Auxiliary variables

    ! Leitura do arquivo de configuração

    call CFG_add(my_cfg, "global%N",1,&
        "number of particles")
    call CFG_add(my_cfg, "global%dt",0.0_dp,&
        "time step")
    call CFG_add(my_cfg, "global%t_fim",0.0_dp,&
        "simulation time") 
    call CFG_add(my_cfg, "global%nimpre",1,&
        "number of result files")
    call CFG_add(my_cfg, "global%dimX",1.0_dp,&
        "x-dimension of the region of calculus")  
    call CFG_add(my_cfg, "global%dimY",1.0_dp,&
        "y-dimension of the region of calculus")    
    call CFG_add(my_cfg,"global%Th",-1.0_dp, &
        "Thermostat hot wall")
    call CFG_add(my_cfg,"global%Tc",-1.0_dp, &
        "Thermostat cold wall") 
    call CFG_add(my_cfg, "global%St0",0.0_dp,&
        "Stokes number of reference")
    call CFG_add(my_cfg, "global%Pe0",0.0_dp,&
        "Peclet number of reference")
    call CFG_add(my_cfg, "global%Kn0",0.0_dp,&
        "Knuden number of reference")
    call CFG_add(my_cfg,"global%GField",(/0.0_dp, 0.0_dp/), &
        "Uniform Gravitational Field")
    call CFG_add(my_cfg,"global%x",(/0.0_dp, 0.0_dp/), &
        "position of the particles")
    call CFG_add(my_cfg,"global%wall",'eeee', &
        "wall's periodic vs elastic") 

    do while (t < t_fim)

    end do 
    

end program main