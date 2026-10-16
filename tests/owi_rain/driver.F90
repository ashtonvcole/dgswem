program owi_rain_driver
    !! Driver to exercise the OWI_RAIN module outside of DG-SWEM.
    !!
    !! Usage: owi_rain_driver DIR T0 TIMINC MASK MMHR
    !! DIR must contain fort.27, fort.425 (and fort.427), points.txt, and
    !! times.txt. points.txt holds a count followed by one "lon lat" per
    !! line. times.txt holds a count followed by one query time per line, in
    !! increasing order. For each time, one line is printed to stdout:
    !! "RESULT time value_1 ... value_n".

    use owi_rain, only: owi_rain_init, get_rain, owi_rain_final
    use sizes, only: sz

    implicit none

    character(len=512) :: arg !! A command line argument
    character(len=512) :: dir !! The input directory
    real(sz) :: t0 !! Reference start time
    real(sz) :: timinc !! Snapshot interval
    integer :: mask !! Nonzero to mask negative values
    integer :: mmhr !! Nonzero to convert from mm/hr to m/s
    integer :: unit !! A file unit
    integer :: npt !! Number of interpolation points
    integer :: ntime !! Number of query times
    integer :: i !! Loop variable
    real(sz), allocatable :: lons(:) !! Point longitudes
    real(sz), allocatable :: lats(:) !! Point latitudes
    real(sz), allocatable :: times(:) !! Query times
    real(sz), allocatable :: rain(:) !! Rainfall at the points

    call get_command_argument(1, dir)
    call get_command_argument(2, arg); read(arg, *) t0
    call get_command_argument(3, arg); read(arg, *) timinc
    call get_command_argument(4, arg); read(arg, *) mask
    call get_command_argument(5, arg); read(arg, *) mmhr

    open(newunit=unit, file=trim(dir) // "/points.txt", status="old", &
        action="read")
    read(unit, *) npt
    allocate(lons(npt), lats(npt), rain(npt))
    do i = 1, npt
        read(unit, *) lons(i), lats(i)
    end do
    close(unit)

    open(newunit=unit, file=trim(dir) // "/times.txt", status="old", &
        action="read")
    read(unit, *) ntime
    allocate(times(ntime))
    do i = 1, ntime
        read(unit, *) times(i)
    end do
    close(unit)

    ! The logging module writes to fort.16
    open(unit=16, file=trim(dir) // "/fort.16", status="replace", &
        action="write")

    call owi_rain_init(trim(dir), lons, lats, t0, timinc, mask /= 0, &
        mmhr /= 0, 0)
    do i = 1, ntime
        call get_rain(times(i), rain, 0)
        write(*, "(a, *(1x, es24.16))") "RESULT", times(i), rain
    end do
    call owi_rain_final()
end program owi_rain_driver
