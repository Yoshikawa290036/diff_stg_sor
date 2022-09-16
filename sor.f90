subroutine sor(nr,ntheta,acerr,peinv,dt,dthetainv,rs,thetas,dr3inv,drinv,dsinv,c,adv,coef)
    implicit none
    integer :: nr,ntheta
    double precision :: peinv,dt,dthetainv
    double precision :: rs(0:nr)
    double precision :: thetas(0:ntheta)
    double precision :: dr3inv(0:nr-1)
    double precision :: drinv(0:nr)
    double precision :: dsinv(0:ntheta-1)
    double precision :: c(-1:nr,0:ntheta-1)
    double precision :: pc(-1:nr,0:ntheta-1)

    double precision :: adv(-1:nr,0:ntheta-1)
    double precision :: coef(0:nr-1,0:ntheta-1)
    double precision :: omega

    double precision :: err,diff,l2err,dr
    double precision :: dc(0:nr-1,0:ntheta-1)
    integer :: i,j,k,MAX,mid
    double precision :: cp1,cm1,aw,ae,as,an

    double precision :: acerr
    MAX = nr*ntheta*4
    ! acerr = 1.d-9
    omega = 0.8d0

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do j = 0,ntheta-1
        do i = -1,nr
            pc(i,j) = c(i,j)
        enddo
    enddo
!$OMP  END PARALLEL DO

    do k = 1,MAX
        err = 0.0d0

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,diff,dr,cm1,cp1,aw,ae,as,an)
        do j = 0,ntheta-1
            do i = 0,nr-1

                if(j.ge.1) then
                    cm1 = c(i,j-1)
                else
                    cm1 = c(i,j)
                endif
                if(j.le.ntheta-2) then
                    cp1 = c(i,j+1)
                else
                    cp1 = c(i,j)
                endif
                dr = rs(i+1)-rs(i)
                ! adv = dr3inv(i)*dsinv(j)*0.5d0*( &
                !     ((-c(i-1,j)+c(i,j))*ur(i,j)+(-c(i,j)+c(i+1,j))*ur(i+1,j)) &
                !     +((-cm1+c(i,j))*uw(i,j)+(-c(i,j)+cp1)*uw(i,j+1)) &
                !     )

                aw = dr3inv(i)*rs(i)**2*drinv(i)
                ae = dr3inv(i)*rs(i+1)**2*drinv(i+1)
                as = dr3inv(i)*dr*dsinv(j)*sin(thetas(j))*dthetainv
                an = dr3inv(i)*dr*dsinv(j)*sin(thetas(j+1))*dthetainv

                diff = 2.0d0*peinv*( &
                       +aw*(c(i-1,j)-c(i,j)) &
                       +ae*(c(i+1,j)-c(i,j)) &
                       +as*(cm1-c(i,j)) &
                       +an*(cp1-c(i,j)) &
                       )
                dc(i,j) = coef(i,j)*(pc(i,j)-c(i,j)+dt*(adv(i,j)-diff))
                pc(i,j) = pc(i,j)-omega*dc(i,j)
                err = err+dc(i,j)**2
            enddo
        enddo
!$OMP  END PARALLEL DO

        l2err = sqrt(err/dble(ntheta*nr))
        ! if(mod(k,100)==0) then
        !     write(*,*) k,l2err
        ! endif

        if(l2err<acerr) then
            ! write(*,*) 'SOR roop',k
            exit
        endif

    enddo

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
    do j = 0,ntheta-1
        do i = -1,nr
            c(i,j) = pc(i,j)
        enddo
    enddo
!$OMP  END PARALLEL DO

    ! stop
endsubroutine sor
