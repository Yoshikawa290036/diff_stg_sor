program main
    implicit none

    integer :: i,j
    integer :: nr,ntheta
    integer :: max_step,nstep
    double precision :: PI,Pe,peinv,sh,alpha,alphainv
    double precision :: lambda,acerr,minac,maxac
    double precision :: drmin,drmininv
    double precision :: dtheta,dthetainv
    double precision :: time,dt

    double precision,dimension(:),allocatable :: rs,thetas,dr3inv,drinv,dsinv
    double precision,dimension(:,:),allocatable :: ur,uw,c,adv,xs,ys,coef
    double precision :: rtmp,thetatmp
    integer :: nmkxyc,nmkstd
    character(32) fname

    PI = atan(1.0d0)*4.0d0
    read(*,*) Pe

! ================ system set ================ !
    max_step = 160000
    nr = 512
    ntheta = 256
    nmkxyc = 10000
    nmkstd = 100

    drmin = 0.001d0
    alpha = 1.02d0
    lambda = 3.0d0/4.0d0
! ================ system set ================ !

    peinv = 1.0d0/Pe
    alphainv = 1.0d0/alpha
    dtheta = PI/dble(ntheta)
    dthetainv = dble(ntheta)/PI
    drmininv = 1.0d0/drmin
    nstep = 0
    time = 0.0d0
    dt = 0.0d0
    minac = 1.d-8
    maxac = 1.d-10
    acerr = minac

    include'allocate.h'

    call cal_rs(nr,drmin,alpha,rs,dr3inv,drinv)
    call cal_thetas(ntheta,dtheta,thetas,dsinv)
    call cal_xs_ys(nr,ntheta,rs,thetas,xs,ys)

    write(*,'("max step                 ",1i9)') max_step
    write(*,'("nmkxyc                   ",1i9)') nmkxyc
    write(*,'("nmkstd                   ",1i9)') nmkstd
    write(*,*)
    write(*,'("Nr                       ",1i9)') nr
    write(*,'("Ntheta                   ",1i9)') ntheta
    write(*,'("Pe                       ",20e20.10)') Pe
    write(*,'("Rmax                     ",20e20.10)') rs(nr)
    write(*,'("drmin                    ",20e20.10)') drmin
    write(*,'("alpha                    ",20e20.10)') alpha
    write(*,'("lambda                   ",20e20.10)') lambda
    write(*,'("acerr                    ",20e20.10)') acerr

    write(*,*)
    write(*,*)

    ! do i = 0,nr-1
    !     rtmp = (rs(i)+rs(i+1))*0.5d0
    !     do j = 0,ntheta-1
    !         thetatmp = (thetas(j)+thetas(j+1))*0.5
    !         ! c(i,j)=(1/rs(nr-1))*((rs(i)+rs(i+1))*0.5-1.0)
    !         ! c(i,j) = 1.0d0/(rs(nr-1)-1)*(rs(nr-1)/rtmp -1)
    !         c(i,j)=1.0d0/rtmp
    !     enddo
    ! enddo


    call bndset(nr,ntheta,alpha,alphainv,c)
    include'mk_xyc.h'

    call cal_vel(nr,ntheta,lambda,xs,ys,rs,thetas,ur,uw)
    call cal_sh(nr,ntheta,drmininv,dtheta,thetas,c,sh)
    call cal_dt(pe,drmin,dtheta,dt)
    call cal_coef(nr,ntheta,peinv,dthetainv,dt,drinv,dr3inv,dsinv,rs,thetas,coef)

    open(12,file='tdtsh')
    write(12,'(20e20.10)') time,dt,sh

    do nstep = 1,max_step
        call bndset(nr,ntheta,alpha,alphainv,c)
        call cal_adv(nr,ntheta,peinv,rs,thetas,drinv,dr3inv,dsinv,dthetainv,ur,uw,c,adv)
        ! if(nstep<=20000) then
        !     acerr = (minac-maxac)/dble(20000)*dble(nstep)+maxac
        ! endif
        call sor2(nr,ntheta,acerr,peinv,dt,dthetainv,rs,thetas,dr3inv,drinv,dsinv,c,adv,coef)
        ! stop
        if(mod(nstep,nmkstd)==0) then
            call cal_sh(nr,ntheta,drmininv,dtheta,thetas,c,sh)
            time = dt*dble(nstep)
            write(12,'(20e20.10)') time,dt,sh
            write(*,*) '---------------------------------------'
            write(*,'("nstep        ",1i9.9)') nstep
            write(*,'("time         ",20e20.10)') time
            write(*,'("dt           ",20e20.10)') dt
            write(*,'("Sh           ",20e20.10)') Sh
            write(*,'("acerr        ",20e20.10)') acerr

        endif

        if(mod(nstep,nmkxyc)==0) then
            include'mk_xyc.h'
        endif
        call flush (6)

    enddo
    close(12)

endprogram main
