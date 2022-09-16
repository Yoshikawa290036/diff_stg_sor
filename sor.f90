
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
    double precision :: omega,dtinv

    double precision :: err,diff,l2err
    double precision :: dc(0:nr-1,0:ntheta-1)
    integer :: i,j,k,MAX,mid

    double precision :: acerr
    MAX = nr*ntheta*4
    ! acerr = 1.d-9
    omega=1.0d0
dtinv=1.0d0/dt

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
do j = 0, ntheta-1
    do i = -1, nr
        pc(i,j)=c(i,j)
    end do
end do
!$OMP  END PARALLEL DO

    do k = 1,MAX
        err = 0.0d0

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,diff) &
!$OMP& REDUCTION(+:err)
        do j = 1,ntheta-2
            do i = 0,nr-1
                diff = 2.0d0*peinv*dr3inv(i)* &
                    & ((-(rs(i)**2)*(-pc(i-1,j)+pc(i,j))*drinv(i)+(rs(i+1)**2)*(-pc(i,j)+pc(i+1,j))*drinv(i+1))+ &
                    & dsinv(j)*dthetainv*(-sin(thetas(j))*(-pc(i,j-1)+pc(i,j))+sin(thetas(j+1))*(-pc(i,j)+pc(i,j+1))))
                dc(i,j) = coef(i,j)*((pc(i,j)-c(i,j))*dtinv +adv(i,j) - diff)
! pc(i,j) = pc(i,j)+omega*dc(i,j)
                err = err+dc(i,j)**2
            enddo
        enddo
!$OMP  END PARALLEL DO


! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i,j,diff) &
! !$OMP& REDUCTION(+:err)
!         do i = 0,nr-1
!             j = 0
!             diff = 2.0d0*peinv*dr3inv(i)* &
!                 & ((-(rs(i)**2)*(-c(i-1,j)+c(i,j))*drinv(i)+(rs(i+1)**2)*(-c(i,j)+c(i+1,j))*drinv(i+1))+ &
!                 & dsinv(j)*dthetainv*(sin(thetas(j+1))*(-c(i,j)+c(i,j+1))))
!             dc(i,j) = coef(i,j)*((pc(i,j)-c(i,j))*dtinv +adv(i,j) - diff)
! ! c(i,j) = c(i,j)+omega*dc(i,j)
!             err = err+dc(i,j)**2
!             ! err = max(err,abs(dc(i,j)))

!             j = ntheta-1
!             diff = 2.0d0*peinv*dr3inv(i)* &
!                 & ((-(rs(i)**2)*(-c(i-1,j)+c(i,j))*drinv(i)+(rs(i+1)**2)*(-c(i,j)+c(i+1,j))*drinv(i+1))+ &
!                 & dsinv(j)*dthetainv*(-sin(thetas(j))*(-c(i,j-1)+c(i,j))))
!             dc(i,j) = coef(i,j)*((pc(i,j)-c(i,j))*dtinv +adv(i,j) - diff)
! ! c(i,j) = c(i,j)+omega*dc(i,j)
!             err = err+dc(i,j)**2
!             ! err = max(err,abs(dc(i,j)))
!         enddo
! !$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
        do j = 0,ntheta-1
            do i = 0,nr-1
                pc(i,j) = pc(i,j)-omega*dc(i,j)
            enddo
        enddo

!$OMP  END PARALLEL DO

        l2err = sqrt(err/dble(ntheta*nr))
        ! if(mod(k,10)==0) then
        !     write(*,*) k,l2err
        ! endif
            write(*,*) k,l2err

        if(l2err<acerr) then
            write(*,*) 'SOR roop',k
            exit
        endif

    enddo
    ! write(*,*) 'SOR roop end --------------'

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j)
do j = 0, ntheta-1
    do i = -1, nr
        c(i,j)=pc(i,j)
    end do
end do
!$OMP  END PARALLEL DO

    stop
endsubroutine sor
