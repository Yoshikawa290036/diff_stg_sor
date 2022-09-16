subroutine cal_adv(nr,ntheta,peinv,rs,thetas,drinv,dr3inv,dsinv,dthetainv,ur,uw,c,adv)
    implicit none
    integer :: nr,ntheta
    double precision :: peinv,dthetainv
    double precision :: rs(0:nr)
    double precision :: thetas(0:ntheta)
    double precision :: dr3inv(0:nr-1)
    double precision :: drinv(0:nr)
    double precision :: dsinv(0:ntheta-1)
    double precision :: ur(0:nr,0:ntheta-1)
    double precision :: uw(0:nr-1,0:ntheta)
    double precision :: c(-1:nr,0:ntheta-1)
    double precision :: adv(-1:nr,0:ntheta-1)
    double precision :: cp1,cm1,aw,ae,as,an

    integer :: i,j
    ! double precision :: adv,vis

! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i,j)
!     do j = 1,ntheta-2
!         do i = 0,nr-1
!             adv(i,j) = dr3inv(i)*dsinv(j)*0.5d0* &
!                 & ((-((c(i-1,j)+c(i,j))*ur(i,j))+((c(i,j)+c(i+1,j))*ur(i+1,j)))+ &
!                 &  (-((c(i,j-1)+c(i,j))*uw(i,j))+((c(i,j)+c(i,j+1))*uw(i,j+1))))

!             ! adv = dr3inv(i)*dsinv(j)* &
!             !     & ((c(i,j)*(-ur(i,j)+ur(i+1,j)-uw(i,j)+uw(i,j+1)))+ &
!             !     & 0.5d0*((-c(i-1,j)+c(i,j))*ur(i,j)+(-c(i,j)+c(i+1,j))*ur(i+1,j)+ &
!             !     & (-c(i,j-1)+c(i,j))*uw(i,j)+(-c(i,j)+c(i,j+1))*uw(i,j+1)))

!             ! vis = 2.0d0*peinv*dr3inv(i)* &
!             !     & ((-(rs(i)**2)*(-c(i-1,j)+c(i,j))*drinv(i)+(rs(i+1)**2)*(-c(i,j)+c(i+1,j))*drinv(i+1))+ &
!             !     & dsinv(j)*dthetainv*(-sin(thetas(j))*(-c(i,j-1)+c(i,j))+sin(thetas(j+1))*(-c(i,j)+c(i,j+1))))

!             ! advis(i,j) = -adv+vis

!         enddo
!     enddo
! !$OMP  END PARALLEL DO

! !$OMP  PARALLEL DO &
! !$OMP& SCHEDULE(static,1) &
! !$OMP& DEFAULT(SHARED) &
! !$OMP& PRIVATE(i,j)
!     do i = 0,nr-1
!         j = 0
!         ! adv = dr3inv(i)*dsinv(j)* &
!         !     & (-(0.5d0*(c(i-1,j)+c(i,j))*ur(i,j))+(0.5d0*(c(i,j)+c(i+1,j))*ur(i+1,j)))
!         ! vis = 2.0d0*peinv*dr3inv(i)* &
!         !     & (-(rs(i)**2)*(-c(i-1,j)+c(i,j))*drinv(i)+rs(i+1)**2*(-c(i,j)+c(i+1,j))*drinv(i+1))

!         adv(i,j) = dr3inv(i)*dsinv(j)* &
!             & ((-(0.5d0*(c(i-1,j)+c(i,j))*ur(i,j))+(0.5d0*(c(i,j)+c(i+1,j))*ur(i+1,j)))+ &
!             &  (-(c(i,j)*uw(i,j))+(0.5d0*(c(i,j)+c(i,j+1))*uw(i,j+1))))

!         ! vis = 2.0d0*peinv*dr3inv(i)* &
!         !     & ((-(rs(i)**2)*(-c(i-1,j)+c(i,j))*drinv(i)+rs(i+1)**2*(-c(i,j)+c(i+1,j))*drinv(i+1))+ &
!         !     & dsinv(j)*dthetainv*(sin(thetas(j+1))*(-c(i,j)+c(i,j+1))))
!         ! advis(i,j) = -adv+vis

!         j = ntheta-1
!         adv(i,j) = dr3inv(i)*dsinv(j)* &
!             & ((-(0.5d0*(c(i-1,j)+c(i,j))*ur(i,j))+(0.5d0*(c(i,j)+c(i+1,j))*ur(i+1,j)))+ &
!             &  (-(0.5d0*(c(i,j-1)+c(i,j))*uw(i,j))+(c(i,j)*uw(i,j+1))))

!         ! vis = 2.0d0*peinv*dr3inv(i)* &
!         !     & ((-(rs(i)**2)*(-c(i-1,j)+c(i,j))*drinv(i)+rs(i+1)**2*(-c(i,j)+c(i+1,j))*drinv(i+1))+ &
!         !     & dsinv(j)*dthetainv*(-sin(thetas(j))*(-c(i,j-1)+c(i,j))))

!         ! advis(i,j) = -adv+vis
!     enddo
! !$OMP  END PARALLEL DO

!$OMP  PARALLEL DO &
!$OMP& SCHEDULE(static,1) &
!$OMP& DEFAULT(SHARED) &
!$OMP& PRIVATE(i,j,cm1,cp1)
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

            adv(i,j) = dr3inv(i)*dsinv(j)*0.5d0*( &
                       ((-c(i-1,j)+c(i,j))*ur(i,j)+(-c(i,j)+c(i+1,j))*ur(i+1,j)) &
                       +((-cm1+c(i,j))*uw(i,j)+(-c(i,j)+cp1)*uw(i,j+1)) &
                       )

            ! aw = dr3inv(i)*rs(i)**2*drinv(i)
            ! ae = dr3inv(i)*rs(i+1)**2*drinv(i+1)
            ! as = dr3inv(i)/drinv(i)*dsinv(j)*sin(thetas(j))*dthetainv
            ! an = dr3inv(i)/drinv(i)*dsinv(j+1)*sin(thetas(j+1))*dthetainv

            ! vis = 2.0d0*peinv*( &
            !       +aw*(c(i-1,j)-c(i,j)) &
            !       +ae*(c(i+1,j)-c(i,j)) &
            !       +as*(cm1-c(i,j)) &
            !       +an*(cp1-c(i,j)) &
            !       )

            ! advis(i,j) = -adv+vis

        enddo
    enddo
!$OMP  END PARALLEL DO

endsubroutine cal_adv
