subroutine cal_dt(pe,drmin,dt)
    implicit none
    double precision :: pe,drmin,dt

    dt =drmin*0.1d0
    dt = min(dt,(drmin**2)*pe/4.0d0)
endsubroutine cal_dt
