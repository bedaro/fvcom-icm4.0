program co2sysbench

      use mod_prec, only: SP
      use mod_CO2SYS, only: co2sysconst_alloc, co2sysconstants, &
        co2sysconst_dealloc, which_K1K2, whoseKSO4
      use mod_sizes, only: mgl, ngl
      use mod_lims, only: kb, kbm1, mtloc, mloc, nloc, ntloc
      use mod_control, only: msr
      use mod_hydrovars, only: zz, d, zz2d, hydro_geom_alloc, &
        hydro_alloc, hydro_geom_dealloc, hydro_dealloc
      use mod_wqm, only: wqm_alloc, t, salt, talk, po4, wqm_dealloc
      implicit none

      integer, parameter :: n = 1000, seed = 12345

      integer :: i, j
      real :: start, finish

      msr = .True.
      kb = 10
      kbm1 = 9
      mgl = 100000
      mloc = 100000
      mtloc = 100000
      ngl = 1
      nloc = 1
      ntloc = 1

      call hydro_geom_alloc()
      call hydro_alloc()
      call wqm_alloc()
      call co2sysconst_alloc()

      call srand(seed)
      ! Construct a test grid with some depths and sigma layers
      do i = 1, mgl
        d(i) = rand() * 150
      end do
      zz(1) = -0.0158
      zz(2) = -0.0605
      zz(3) = -0.1269
      zz(4) = -0.2086
      zz(5) = -0.3033
      zz(6) = -0.4092
      zz(7) = -0.5252
      zz(8) = -0.6506
      zz(9) = -0.7847
      zz(10) = -0.9269
      do i = 0, mgl
        do j = 1, kbm1
          t(i, j) = 18*rand()
          salt(i, j) = 25+8*rand()
          talk(i, j) = 7700*rand()
          po4(i, j) = 0.6*rand()
        end do
        zz2d(i, kb) = zz(kb)
      end do

      which_K1K2 = 10
      whoseKSO4 = 1

      call cpu_time(start)

      call co2sysconstants()

      call cpu_time(finish)

      print '("Time = ",f6.3," seconds.")', finish - start

      call co2sysconst_dealloc()
      call wqm_dealloc()
      call hydro_dealloc()
      call hydro_geom_dealloc()

end program
