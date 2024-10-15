! Compute derivative and intra-sigma levels after sigma levels (Z) have
! been set
subroutine set_sigma
      use mod_prec, only: SP
      use mod_sizes, only: mgl
      use mod_lims, only: KB, KBM1, KBM2
      use mod_hydrovars, only: Z, DZ, DZZ, ZZ, Z2DG
      implicit none

      integer :: I, K

      Do K = 1, KBM1
        ZZ (K) = 0.5_SP * (Z(K)+Z(K+1))
        DZ (K) = Z (K) - Z (K+1)
      End Do
        ZZ (KB) = 2.0_SP * ZZ (KBM1) - ZZ (KBM2)
      Do K = 1, KBM2
        DZZ (K) = ZZ (K) - ZZ (K+1)
      End Do
      DZZ (KBM1) = 0.0
      DZ (KB) = 0.0

      ! These spatially varying values may get overwritten in ncdio
      ! if SETSIGMA is on
      DO I=1,MGL
    		DO K=1,KB
    			Z2DG (I,K) = Z(K)
    		ENDDO
  		ENDDO
      ! The other spatially varying variables get set in pdomdec
end subroutine set_sigma
