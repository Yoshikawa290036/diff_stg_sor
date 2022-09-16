subroutine cal_coef(nr,ntheta,peinv,dthetainv,dt,drinv,dr3inv,dsinv,rs,thetas,coef)
    implicit none
    integer :: nr,ntheta
    double precision :: peinv,dt,dthetainv
    double precision :: drinv(0:nr)
    double precision :: dr3inv(0:nr-1)
    double precision :: dsinv(0:ntheta-1)
    double precision :: rs(0:nr)
    double precision :: thetas(0:ntheta)
    double precision :: coef(0:nr-1,0:ntheta-1)

    integer :: i,j
    double precision :: beta,dtinv
    ! dtinv=1.0d0/dt

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do j = 0,ntheta-1
        do i = 0,nr-1
            coef(i,j) = 1.0d0/(1.0d0-2.0d0*peinv*dt*dr3inv(i)* &
                      & (-(rs(i)**2)*drinv(i)-(rs(i+1)**2)*drinv(i+1)+ &
                      & dsinv(j)*dthetainv*(-sin(thetas(j))-sin(thetas(j+1)))))
            ! coef(i,j) = 1.0d0/(-dtinv+(2.0d0*peinv*dr3inv(i)* &
            !           & (-(rs(i)**2)*drinv(i)-(rs(i+1)**2)*drinv(i+1)+ &
            !           & dsinv(j)*dthetainv*(-sin(thetas(j))-sin(thetas(j+1))))))
        enddo
    enddo
!$OMP  END PARALLEL DO

    ! do j = 0,ntheta-1
    !     do i = 0,nr-1
    !         write(*,*) coef(i,j)
    !     enddo
    ! enddo
    ! stop
endsubroutine cal_coef
