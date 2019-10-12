! Simple langevin dynamics with brownian motion, no interaction between partciles and custom temperature distribution
! conditions on the cold side will be the reference

module fisica

    contains 

    function St_slip(Kn0,T_loc, Tc, N, scorr) result(Cc)
        ! allen and raabe 1982
        use mod1 

        integer :: N
        real(dp) :: Kn0, Tc
        real(dp), dimension(:) :: T_loc
        real(dp), dimension(N) :: Cc
        logical :: scorr 

        if (scorr) then
            Cc = (1 + (Kn0 * T_loc)*(1.155 + 0.471 * exp(-0.596 /  (Kn0 * T_loc) ))) 
        else 
            Cc = 1
        end if

    end function St_slip

    function densitycr(rhof,T_loc,N) result(drho)
        ! corrige a densidade
        ! rho é adimensional
        use mod1

        real(dp) :: T_loc(:), rhof, rho
        real(dp), dimension(N) :: drho

        drho = (2 - rhof/T_loc)/(2 - rhof) 

    end function densitycr


    subroutine comp_F(F,GField,x,v,T_loc,Tc,St0,Pe0,Kn0,rhof,N,dt,scorr)
        use mod1

        real(dp), dimension(:,:), intent(inout) :: F
        real(dp), dimension(:,:), intent(in) :: x,v
        real(dp), dimension(:), intent(in) :: T_loc, Gfield
        real(dp), intent(in) :: St0, Pe0, Kn0, dt, rhof,Tc 
        integer :: N
        logical, INTENT(IN) :: scorr
        real(dp) :: xi(N,2),drho(N), Cc(N)


        ! Gera numero aleatório e normaliza
        call random_seed()
        call random_number(xi)
        xi = xi*2 -1
        xi(:,1) = xi(:,1) / sqrt(xi(:,1)**2 + xi(:,2)**2)
        xi(:,2) = xi(:,2) / sqrt(xi(:,1)**2 + xi(:,2)**2)

        ! read(*,*)

        drho = densitycr(rhof,T_loc,N)

        if (Pe0 > 0) then
            Cc = St_slip(Kn0,T_loc,Tc,N, scorr)

            ! F(:,1) = (1/(St0*drho/St_slip(Kn0,T_loc,Tc,N,scorr))) * (- v(:,1) + GField(1) + sqrt(6*T_loc/(Pe0*dt))*xi(:,1))
            ! F(:,2) = (1/(St0*drho/St_slip(Kn0,T_loc,Tc,N,scorr))) * (- v(:,2) + GField(2) + sqrt(6*T_loc/(Pe0*dt))*xi(:,2))

            F(:,1) = (1/(St0*drho)) * (- v(:,1)/Cc + GField(1) + sqrt(6*T_loc/(Cc*Pe0*dt))*xi(:,1))
            F(:,2) = (1/(St0*drho)) * (- v(:,2)/Cc + GField(2) + sqrt(6*T_loc/(Cc*Pe0*dt))*xi(:,2))

            ! u1x = (Cc(T_ax)/(St*Ddensity(T_ax)))*(-u0x + (6/(Pe*tau))^.5*ksi(1)*temper/Cc(T_ax));

            ! F(:,1) = (1/(St0)) * (- v(:,1) + GField(1) + sqrt(6*T_loc/(Pe0*dt))*xi(:,1))
            ! F(:,2) = (1/(St0)) * (- v(:,2) + GField(2) + sqrt(6*T_loc/(Pe0*dt))*xi(:,2))

            ! print*, "F =", F
            ! print*, "drho =", drho
            ! print*, "v =", v
            ! print*, "x =", x
            ! read(*,*)

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
            ! print*, "F*dt =", F, "x", dt, "=", F*dt
        else 
            v = F*dt
        end if 
    end subroutine comp_v


    subroutine  comp_x(x,v,dimx,dimy,N,dt,wall)
        use mod1 
        
        real(dp), dimension(:,:), intent(inout) :: v,x
        real(dp), intent(in) :: dt,dimx,dimy
        integer, intent(in) :: N 
        character, intent(in) :: wall(4)

        x = x + v*dt 
        ! print*, "v*dt =", v*dt
        ! read(*,*)

        do i = 1,N
            if (v(i,1)*dt > dimx/2 .or. v(i,2)*dt > dimy/2) then
                print*, "Partículas rápidas demais"
                call system("killall langevin")
            end if 
        end do
        ! print*, "v*dt", v*dt
        ! read(*,*)

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
                    x(i,2) = 2*dimy - x(i,2)
                    v(i,2) = - v(i,2)
                else if (x(i,2) < 0) then 
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
    real(dp) :: Pe0, St0, Th, Tc, Kn0, dt, t_fim, t, rhof, dimx, dimy, Gfield(2),x0(2),v0(2),start,finish
    real(dp),allocatable, dimension(:,:) :: v, x, F
    real(dp), allocatable, dimension(:) :: T_loc
    integer :: N, nimpre, i, ic1, cpu_countrate
    integer, allocatable :: interv(:)
    character(4) :: wall
    character(10) :: time 
    character(8) :: date
    type(CFG_t) :: my_cfg
    logical :: scorr 
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
    call CFG_add(my_cfg,"global%slip_correction", .True., &
        "Stokes slip correction")
    call CFG_add(my_cfg, "global%rhof",0.0_dp,&
    "sensitivity of the density to temperature")
    call CFG_add(my_cfg,"global%GField",(/0.0_dp, 0.0_dp/), &
    "Uniform Gravitational Field")
    call CFG_add(my_cfg,"global%wall",'eeee', &
        "wall's periodic vs elastic") 
    call CFG_add(my_cfg,"global%x",(/0.0_dp, 0.0_dp/), &
    "position of the particles")
    call CFG_add(my_cfg,"global%v",(/0.0_dp, 0.0_dp/), &
    "initial velocity of the particles")

    
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
    call CFG_get(my_cfg, "global%slip_correction", scorr)
    call CFG_get(my_cfg, "global%Kn0",rhof)
    call CFG_get(my_cfg,"global%GField",Gfield)
    call CFG_get(my_cfg,"global%x",x0)
    call CFG_get(my_cfg,"global%v",v0)
    call CFG_get(my_cfg,"global%wall",wall)

    call CFG_write(my_cfg, "settings.txt") 



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
    
    v(:,1) = v0(1)
    v(:,2) = v0(2)
    print*, "vini = ", v0

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

    t = 0
    call vec2csv(x,N,2,'position',0,t,nimpre,start, .true.)
    call vec2csv(v,N,2,'velocity',0,t,nimpre,start, .true.)
    
    print*, "Begin"
    ! print*, interv
    T_loc = Tc + x(:,1)*((Th-Tc)/dimx)
    do while (t < t_fim)
        call comp_F(F,Gfield,x,v,T_loc,Tc,St0,Pe0,Kn0,rhof,N,dt,scorr)
        call comp_v(F,v,dt,Pe0)
        call comp_x(x,v,dimx,dimy,N,dt,wall) 

        T_loc = Tc + x(:,1)*((Th-Tc)/dimx)

        if (jprint == 1) then
            print*, "Avg velocity:", sum(sqrt(v(:,1)**2 + v(:,2)**2))/N 
            print*, "Avg displacement:", dt*sum(sqrt(v(:,1)**2 + v(:,2)**2))/N 
        end if
        if (jprint == interv(step+1)) then 
            call vec2csv(v,N,2,'velocity',step,t,nimpre,start)
            call vec2csv(x,N,2,'position',step,t,nimpre,start)
            step = step + 1
        end if 

        jprint = jprint + 1
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
    
    call system("python csv2vtk_particles.py")

end program main