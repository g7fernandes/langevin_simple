! Simple langevin dynamics with brownian motion, no interaction between partciles and custom temperature distribution
! conditions on the cold side will be the reference

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
        real(dp), dimension(N) :: drho

        drho = (rho - rhof/T)/(rho - rhof)

    end function densitycr


    subroutine comp_F(F,GField,x,v,T_loc,St0,Pe0,Kn0,N,dt)
        use mod1

        real(dp), dimension(:,:), intent(inout) :: F
        real(dp), dimension(:,:), intent(in) :: x,v
        real(dp), dimension(:), intent(in) :: T_loc, Gfield
        real(dp), intent(in) :: St0, Pe0, Kn0, dt
        integer :: N
        real(dp), dimension(N,2) :: xi 


        ! Gera numero aleatório e normaliza
        call random_seed()
        call random_number(xi)
        xi = xi*2 -1
        xi(:,1) = xi(:,1) / sqrt(xi(:,1)**2 + xi(:,2)**2)
        xi(:,2) = xi(:,2) / sqrt(xi(:,1)**2 + xi(:,2)**2)

        if (Pe0 > 0) then
            F(:,1) = (1/(St0/St_slip(Kn0,T_loc,N))) * (- v(:,1) + GField(1) + sqrt(6*T_loc/(Pe0*dt))*xi(:,1))
            F(:,2) = (1/(St0/St_slip(Kn0,T_loc,N))) * (- v(:,2) + GField(2) + sqrt(6*T_loc/(Pe0*dt))*xi(:,2))
        else
            F(:,1) = (GField(1) + sqrt(6*T_loc/(Pe0*dt))*xi(:,1))/dt
            F(:,2) = (GField(2) + sqrt(6*T_loc/(Pe0*dt))*xi(:,2))/dt
        end if

    end subroutine comp_F

    subroutine comp_v(F,v,dt,Pe0)
        use mod1

        real(dp), dimension(:,:), intent(inout) :: v
        real(dp), dimension(:,:), intent(in) :: F
        real(dp), intent(in) :: dt,Pe0

        if (Pe0 > 0) then
            v = v + F*dt 
        else 
            v = F*dt
        end if 
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
        if (wall(3) == 'e') then 
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

        if (wall(1) == 'e') then 
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

        ! Periodico 

        if (wall(3) == 'p') then 
            do i = 1,N 
                if (x(i,1) > dimx) then 
                    x(i,1) = -dimx + x(i,1)
                else if (x(i,1) < 0) then 
                    x(i,1) = dimx + x(i,1)
                end if
            end do 
        end if

        if (wall(1) == 'p') then 
            do i = 1,N 
                if (x(i,2) > dimy) then 
                    x(i,2) = -dimy + x(i,2)
                else if (x(i,2) < 0) then 
                    x(i,2) = dimy + x(i,2)
                end if
            end do 
        end if   

    end subroutine comp_x

end module fisica

program main

    use mod1 
    use m_config
    use saida
    use fisica

    use, intrinsic :: ieee_arithmetic, only: IEEE_Value, IEEE_QUIET_NAN
    use, intrinsic :: iso_fortran_env, only: real32

    implicit none

    ! Variables
    real(dp) :: Pe0, St0, Th, Tc, Kn0, dt, t_fim, t, rhof, rho, dimx, dimy, Gfield(2),x0(2),start,finish
    real(dp),allocatable, dimension(:,:) :: v, x, F
    real(dp), allocatable, dimension(:) :: T_loc
    integer :: N, nimpre, i, ic1, cpu_countrate
    integer, allocatable :: interv(:)
    character :: wall(4)
    character(10) :: time 
    character(8) :: date
    type(CFG_t) :: my_cfg
    ! Auxiliary variables
    integer :: step, jprint
    

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
    call CFG_add(my_cfg, "global%rhof",0.0_dp,&
    "Knuden number of reference")
    call CFG_add(my_cfg,"global%GField",(/0.0_dp, 0.0_dp/), &
    "Uniform Gravitational Field")
    call CFG_add(my_cfg,"global%x",(/0.0_dp, 0.0_dp/), &
    "position of the particles")
    call CFG_add(my_cfg,"global%wall",'eeee', &
    "wall's periodic vs elastic") 
    
    call CFG_read_file(my_cfg,'settings.ini')
    
    call CFG_get(my_cfg, "global%N",N)
    call CFG_get(my_cfg, "global%dt",dt)
    call CFG_get(my_cfg, "global%t_fim", t_fim) 
    call CFG_get(my_cfg, "global%nimpre",nimpre)
    call CFG_get(my_cfg, "global%dimX",dimx)  
    call CFG_get(my_cfg, "global%dimY",dimy)    
    call CFG_get(my_cfg,"global%Th",Th)
    call CFG_get(my_cfg,"global%Tc",Tc) 
    call CFG_get(my_cfg, "global%St0",St0)
    call CFG_get(my_cfg, "global%Pe0",Pe0)
    call CFG_get(my_cfg, "global%Kn0",Kn0)
    call CFG_get(my_cfg, "global%Kn0",rhof)
    call CFG_get(my_cfg,"global%GField",Gfield)
    call CFG_get(my_cfg,"global%x",x0)
    call CFG_get(my_cfg,"global%wall",wall(4)) 

    allocate(interv(nimpre),T_loc(N),x(N,2),v(N,2),F(N,2))
    if (x0(1) < 0) then
        call random_seed()
        call random_number(x)
        x(:,1) = x(:,1)*dimx
        x(:,2) = x(:,2)*dimy
    else
        x(:,1) = x0(1) ! ou x = x0 funciona? 
        x(:,2) = x0(2) ! ou x = x0 funciona? 
    end if 
    v = 0 

    call system_clock(ic1,cpu_countrate)
    start = real(ic1,kind(0.d0))/real(cpu_countrate,kind(0.d0))

    interv =   (/((nint(t_fim/dt)/nimpre)*i,i=0,nimpre) /)

    jprint = 1
    step = 1

    call date_and_time(date,time)
    write(*,'("|| Program started at: ", a,":",a,":",a,2x, a,"/",a,"/",a, "||")') & 
            time(1:2),time(3:4),time(5:8),date(5:6),date(7:8),date(1:4)

    call system('mkdir temp') !pasta temporária para armazenar os resultados
    ! call linked2vec(malha,domx,domy,nxv,aux1)

    call vec2csv(x,N,2,'position',0,t,nimpre,start, .true.)
    call vec2csv(v,N,2,'velocity',0,t,nimpre,start, .true.)
    
    t = 0
    do while (t < t_fim)
        call comp_F(F,Gfield,x,v,T_loc,St0,Pe0,Kn0,N,dt)
        call comp_v(F,v,dt,Pe0)
        call comp_x(x,v,dimx,dimy,N,dt,wall) 

        T_loc = Tc + x(:,1)*((Th-Tc)/dimx)

        if (step == interv(jprint)) then 
            call vec2csv(v,N,2,'velocity',step,t,nimpre,start)
            call vec2csv(v,N,2,'postition',step,t,nimpre,start)
            jprint = jprint + 1
        end if 

        step = step + 1
        t = t + dt
    end do 


    call system_clock(ic1,cpu_countrate)
    finish = real(ic1)/real(cpu_countrate,kind(0.d0))
    print '("Time = ",f10.3," seconds.")',(finish-start)
    open(unit=22,file='settings.txt',status="old", position="append", action="write")
    write(22,'("#:Time to complete: ",f10.3," secounds.")') (finish-start)
    write(22,'("#:Execution date: ",a,"/",a,"/",a,2x,a,":",a,":",a)') & 
    date(5:6),date(7:8),date(1:4),time(1:2),time(3:4),time(5:8)
    close(22)
    

end program main