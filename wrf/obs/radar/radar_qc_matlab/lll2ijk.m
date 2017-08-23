


SUBROUTINE lll2ijk(rlon,rlat,rlev,ri,rj,rk)
  IMPLICIT NONE
  REAL(r_size),INTENT(IN) :: rlon
  REAL(r_size),INTENT(IN) :: rlat
  REAL(r_size),INTENT(IN) :: rlev
  REAL(r_size),INTENT(OUT) :: ri
  REAL(r_size),INTENT(OUT) :: rj
  REAL(r_size),INTENT(OUT) :: rk

    ri=(rlon-LON(1,1))/DLON + 1.0d0
    rj=(rlat-LAT(1,1))/DLAT + 1.0d0    
    rk=(rlev-Z(1,1,1))/DZ   + 1.0d0
  

  RETURN
END SUBROUTINE lll2ijk
