module owi_rain
    !! Module to handle time-series mass sources and sinks on rectilinear
    !! longitude-latitude grids. For the shallow water equations, this
    !! corresponds to rainfall and evaporation.
    !!
    !! Added by Ashton Cole in October 2026. Based on work by Kendra Dresback
    !! at the University of Oklahoma in 2013.
    !!
    !! ## Functions
    !! - owi_rain_init: Call once at the beginning of the simulation to open
    !!   precipitation files and read metadata
    !! - get_rain: Call at each time step (or SSP-RK Euler sub-step) to get
    !!   spatially and temporally interpolated rainfall rates at points of
    !!   interest.
    !! - cycle_snapshots: Called internally to update the two snapshots
    !!   used for time interpolation.

    use owi_util, only: &
        get_owi_header, &
        get_owi_record_scalar, &
        apply_owi_interpolator, &
        get_owi_precipitation_parameters, &
        owi_grid_params_t, &
        owi_grid_t, &
        owi_grid_interpolator_t
    use logging, only: log_info, log_error
    use sizes, only: sz

    implicit none

    private
    public :: owi_rain_init, get_rain, owi_rain_final

    ! Committing sin of global variables, for now

    integer, parameter :: parameter_file_unit = 27
    !! Unit number for fort.27 parameter file.
    character(len=7), parameter :: parameter_file_name = "fort.27"
    !! File name for fort.27 parameter file.
    integer :: num_sets
    !! Number of rainfall data grids used. Accepted values are 1 (use fort.425)
    !! and 2 (use fort.425 and optional refined fort.427).
    integer :: snapshot_offset
    !! Number of blank snapshots to apply before reading time series data. If
    !! the value is negative, then snapshots will be skipped.
    real(sz), allocatable :: lons_interp(:)
    !! Longitudes of points of interest to interpolate at, e.g. mesh nodes,
    !! element centers, or Gauss points.
    real(sz), allocatable :: lats_interp(:)
    !! Latitudes of points of interest to interpolate at, e.g. mesh nodes,
    !! element centers, or Gauss points.
    real(sz) :: ta
    !! First time for time-series interpolation.
    real(sz) :: tb
    !! Second time for time-series interpolation.
    real(sz) :: snapshot_increment
    !! Fixed time increment between snapshots.
    real(sz), allocatable :: rain_ta(:)
    !! First rainfall value for time-series interpolation.
    real(sz), allocatable :: rain_tb(:)
    !! Second rainfall value for time-series interpolation.
    logical :: exclude_negatives
    !! Whether to mask negative values from the OWI data.
    logical :: mmhr2ms
    !! Whether to convert input from mm/hr to m/s.

    ! Variables for fort.425
    ! Coarse "basin" file
    integer, parameter :: basin_file_unit = 425
    !! Unit number for fort.425 parameter file.
    character(len=8), parameter :: basin_file_name = "fort.425"
    !! File name for fort.425 parameter file.
    ! type(owi_grid_interpolator_t) :: basin_interpolator_ta
    ! Not needed, spatial interpolation always happens at tb
    logical :: basin_file_eof
    !! Whether the basin file has reached EOF. This prevents additional reads,
    !! which would trigger a format error.
    type(owi_grid_interpolator_t) :: basin_interpolator_tb
    !! Basin spatial interpolation object at time tb.

    ! Variables for fort.427
    ! Optional refined "region" file
    ! Basin values are overwritten
    integer, parameter :: region_file_unit = 427
    !! Unit number for fort.427 parameter file.
    character(len=8), parameter :: region_file_name = "fort.427"
    !! File name for fort.425 parameter file.
    ! type(owi_grid_interpolator_t) :: region_interpolator_ta
    ! Not needed, spatial interpolation always happens at tb
    logical :: region_file_eof
    !! Whether the region file has reached EOF. This prevents additional reads,
    !! which would trigger a format error.
    type(owi_grid_interpolator_t) :: region_interpolator_tb
    !! Region spatial interpolation object at time tb.

    contains

    subroutine owi_rain_init(gblinputdir, lons, lats, t, timinc, &
        mask_negatives, use_mmhr, myproc)
        !! Initialize the OWI_RAIN module and its global variables.
        !!
        !! This means reading metadata and inital time series data from OWI-
        !! style precipitation data files fort.27, fort.425, and fort.427. File
        !! 27 is read and closed. Files 425 and 427 are left open.

        implicit none

        character(len=*), intent(in) :: gblinputdir
        !! The directory in which to look for precipitation files.
        real(sz), intent(in) :: lons(:)
        !! Longitudes of points of interest to interpolate at, e.g. mesh nodes,
        !! element centers, or Gauss points.
        real(sz), intent(in) :: lats(:)
        !! Latitudes of points of interest to interpolate at, e.g. mesh nodes,
        !! element centers, or Gauss points.
        real(sz), intent(in) :: t
        !! Reference start time for precipitation. Since using dates from the
        !! OWI file is not yet supported, this number is technically arbitrary.
        real(sz), intent(in) :: timinc
        !! The time increment to the next snapshot. Currently, using dates in
        !! the OWI file is not supported.
        logical, intent(in) :: mask_negatives
        !! Whether to mask negative values from the OWI data.
        logical, intent(in) :: use_mmhr
        !! Whether the file uses mm/hr as its units.
        integer, intent(in) :: myproc
        !! The MPI process number.

        integer :: iostat !! Error code for OWI file read
        integer :: i !! Snapshot loop variable
        integer :: npt !! Number of interpolation points
        character(len=256) :: msg !! Buffer for log messages

        npt = min(size(lons), size(lats))

        ! Initialize some module globals
        lons_interp = lons
        lats_interp = lats
        ta = t
        tb = t + timinc
        snapshot_increment = timinc
        allocate(rain_ta(npt))
        rain_ta = 0.0
        allocate(rain_tb(npt))
        rain_tb = 0.0
        exclude_negatives = mask_negatives
        mmhr2ms = use_mmhr


        ! Read precipitation metadata fort.27
        open( &
            unit=parameter_file_unit, &
            file=trim(gblinputdir) // "/" // parameter_file_name, &
            status="old", &
            action="read" &
        )
        call get_owi_precipitation_parameters( &
            unit=parameter_file_unit, &
            iostat=iostat, &
            num_sets=num_sets, &
            snapshot_offset=snapshot_offset &
        )
        if (iostat /= 0) then
            write(msg, "(a, i4)") "Error in module owi_rain: " // &
                "precipitation parameter read failed with code ", iostat
            call log_error(msg, myproc)
        end if
        close(parameter_file_unit)

        ! Check number of data sets
        ! If 1, use fort.425 only
        ! If 2, use fort.425 and optional refined fort.427
        ! Data from the refined regional grid overwrites the basin data
        ! Otherwise, throw an error
        if (num_sets /=1 .and. num_sets /= 2) then
            call log_error("Error in module owi_rain: num_sets" // &
                " must be 1 or 2", myproc)
        end if

        ! Record the timing assumption, since OWI snapshot dates are ignored
        ! Screen output on rank 0 only, since every rank calls this routine
        write(msg, "(a, f0.1, a, f0.1, a)") &
            "OWI rainfall: dates in the OWI files are ignored. Snapshots " // &
            "are assumed to be spaced ", snapshot_increment, &
            " s apart, with the first at t = ", ta, " s."
        call log_info(trim(msg), myproc=myproc, to_screen=myproc == 0)

        ! Read coarse basin file fort.425
        ! Get metadata
        ! Skip the necessary number of snapshots (if negative offset)
        ! Read in first two snapshots
        ! Take into account the number of blanks (if positive offset)
        ! Leave open for reading records
        open( &
            unit=basin_file_unit, &
            file=trim(gblinputdir) // "/" // basin_file_name, &
            status="old", &
            action="read" &
        )
        basin_file_eof = .false.
        call get_owi_header( &
            unit=basin_file_unit, &
            iostat=iostat &
        ) ! OWI dates not yet supported
        if (iostat /= 0) then
            write(msg, "(a, i4)") "Error in module owi_rain: " // &
                "basin header read failed with code ", iostat
            call log_error(msg, myproc)
        end if
        skip_basin_snapshots: do i = 1, -1 * snapshot_offset
            call get_owi_record_scalar( &
                unit=basin_file_unit, &
                iostat=iostat &
            )
            if (iostat < 0) then
                ! EOF is fine
                call log_info("Module owi_rain: basin file EOF " // &
                    "during skipping phase.", myproc=myproc, to_screen=myproc == 0)
                basin_file_eof = .true.
                exit skip_basin_snapshots
            else if (iostat > 0) then
                ! Format error
                call log_error("Error in module owi_rain: basin file " // &
                    "format error during skipping phase.", myproc)
            end if
        end do skip_basin_snapshots

        ! Optionally, read refined region file fort.427
        ! Get metadata
        ! Skip the necessary number of snapshots (if negative offset)
        ! Read in first two snapshots
        ! Take into account the number of blanks (if positive offset)
        ! Leave open for reading records
        if (num_sets == 2) then
            open( &
                unit=region_file_unit, &
                file=trim(gblinputdir) // "/" // region_file_name, &
                status="old", &
                action="read" &
            )
            region_file_eof = .false.
            call get_owi_header( &
                unit=region_file_unit, &
                iostat=iostat &
            ) ! OWI dates not yet supported
            if (iostat /= 0) then
                write(msg, "(a, i4)") "Error in module owi_rain: " // &
                    "region header read failed with code ", iostat
                call log_error(msg, myproc)
            end if
            skip_region_snapshots: do i = 1, -1 * snapshot_offset
                call get_owi_record_scalar( &
                    unit=region_file_unit, &
                    iostat=iostat &
                )
                if (iostat < 0) then
                    ! EOF is fine
                    call log_info("Module owi_rain: region file EOF " // &
                        "during skipping phase.", myproc=myproc, to_screen=myproc == 0)
                    region_file_eof = .true.
                    exit skip_region_snapshots
                else if (iostat > 0) then
                    ! Format error
                    call log_error("Error in module owi_rain: region file " // &
                        "format error during skipping phase.", myproc)
                end if
            end do skip_region_snapshots
        end if

        ! Read in the first two snapshots
        ! This initializes rain_ta, rain_tb, and interpolators
        ! ta and tb are already set, so don't update them
        call cycle_snapshots(advance_time=.false., myproc=myproc)
        call cycle_snapshots(advance_time=.false., myproc=myproc)
    end subroutine owi_rain_init

    subroutine get_rain(t, rain, myproc)
        !! Get rainfall values at certain points at a certain time step.
        !!
        !! Crucially, this function is intended to be called for monotonically
        !! increasing t, since past snapshots are discarded after they are read
        !! and used. If this function is called for t < ta, an error will be
        !! thrown, and the program will be terminated.

        implicit none

        real(sz), intent(in) :: t
        !! The current simulation time.
        real(sz), intent(out) :: rain(:)
        !! The rainfall values at points of interest.
        integer, intent(in) :: myproc
        !! The MPI process number.

        ! Local variables
        real(sz) :: a !! Normalized distance from ta

        ! Check if t is in range
        ! Assumption is that time steps are computed in temporal order
        ! This is not necessarily guaranteed for RK substages
        ! If before the current snapshots, error
        ! If after the current snapshots, read in the next record
        if (t < ta - 1.0e-6_sz) then
            call log_error("Error in module owi_rain: t < ta", myproc)
        else if (ta <= t .and. t <= tb) then
            ! Do nothing
            continue
        else
            ! Shift to the appropriate snapshot interval
            do while (t > tb)
                call cycle_snapshots(myproc=myproc)
            end do
        end if

        ! Set the time-interpolated value
        a = (t - ta) / (tb - ta)
        rain = (1 - a) * rain_ta + a  * rain_tb
    end subroutine get_rain

    subroutine cycle_snapshots(advance_time, myproc)
        !! Read in data and update the rainfall snapshots for time
        !! interpolation.
        !!
        !! This means that the data at tb is pushed to ta, and new data is
        !! read into tb.

        implicit none

        logical, optional, intent(in) :: advance_time
        !! Whether to advance the reference times. Generally this is true,
        !! except for when skipping snapshots at initialization.
        integer, intent(in) :: myproc
        !! The MPI process number.

        ! For defaulting
        logical :: the_advance_time !! Whether to advance the reference times

        ! Local variables
        type(owi_grid_params_t) :: basin_grid_params_tb !! Local parameters for new basin grid
        type(owi_grid_t) :: basin_grid_tb !! Local new basin grid
        real(sz), allocatable :: basin_data_tb(:, :) !! Local new basin data
        integer :: basin_iostat !! Status of basin file read
        integer(8) :: basin_dt !! Date of basin snapshot as YYYYMMDDHHMM, currently unused
        type(owi_grid_params_t) :: region_grid_params_tb !! Local parameters for new region grid
        type(owi_grid_t) :: region_grid_tb !! Local new region grid
        real(sz), allocatable :: region_data_tb(:, :) !! Local new region data
        integer :: region_iostat !! Status of region file read
        integer(8) :: region_dt !! Date of region snapshot as YYYYMMDDHHMM, currently unused

        if (present(advance_time)) then
            the_advance_time = advance_time
        else
            the_advance_time = .true.
        end if

        ! Send tb variables to ta
        if (the_advance_time) ta = tb
        rain_ta = rain_tb
        ! basin_interpolator_ta = basin_interpolator_tb ! Not used

        ! Read new data to tb variables
        ! Leave file-based timing unimplemented
        ! Only update tb interpolators if grids are not equal
        ! Otherwise, use the existing weights to update rain_tb
        if (the_advance_time) tb = ta + snapshot_increment

        if (snapshot_offset > 0) then
            ! Apply a blank snapshot and decrement
            rain_tb = 0.0
            snapshot_offset = snapshot_offset - 1
        else if (num_sets == 1) then
            ! Only use basin
            ! Read basin record
            call get_owi_record_scalar( &
                unit=basin_file_unit, &
                iostat=basin_iostat, &
                grid_params=basin_grid_params_tb, &
                dt=basin_dt, &
                data=basin_data_tb &
            )
            if (basin_file_eof .or. basin_iostat < 0) then
                ! Handle EOF with zeros
                call log_info("Module owi_rain: basin file EOF. " // &
                    "Inserting blank snapshot.", myproc=myproc, to_screen=myproc == 0)
                rain_tb = 0.0
                basin_file_eof = .true.
                return
            else if (basin_iostat > 0) then
                ! Format error
                call log_error("Error in module owi_rain: basin file " // &
                    "format error", myproc)
            end if
            ! Unit converion
            if (mmhr2ms) basin_data_tb = basin_data_tb / 3600000.0
            ! Check if weights need to be recomputed
            if (basin_grid_params_tb /= &
                basin_interpolator_tb%grid%grid_params &
            ) then
                call basin_grid_tb%init(basin_grid_params_tb)
                call basin_interpolator_tb%init( &
                    grid=basin_grid_tb, &
                    lons_interp=lons_interp, &
                    lats_interp=lats_interp &
                )
            end if
            call apply_owi_interpolator( &
                data_grid=basin_data_tb, &
                interpolator_grid=basin_interpolator_tb, &
                data_interp=rain_tb, & ! Output
                exclude_negatives=exclude_negatives &
            )
        else if (num_sets == 2) then
            ! Use both basin and region
            ! Read basin record
            call get_owi_record_scalar( &
                unit=basin_file_unit, &
                iostat=basin_iostat, &
                grid_params=basin_grid_params_tb, &
                dt=basin_dt, &
                data=basin_data_tb &
            )
            if (basin_file_eof .or. basin_iostat < 0) then
                ! Handle EOF with zeros
                ! If one file is not present, neither is used
                call log_info("Module owi_rain: basin file EOF. " // &
                    "Inserting blank snapshot.", myproc=myproc, to_screen=myproc == 0)
                rain_tb = 0.0
                basin_file_eof = .true.
                return
            else if (basin_iostat > 0) then
                ! Format error
                call log_error("Error in module owi_rain: basin file " // &
                    "format error", myproc)
            end if
            ! Unit conversion
            if (mmhr2ms) basin_data_tb = basin_data_tb / 3600000.0
            ! Read region record
            call get_owi_record_scalar( &
                unit=region_file_unit, &
                iostat=region_iostat, &
                grid_params=region_grid_params_tb, &
                dt=region_dt, &
                data=region_data_tb &
            )
            if (region_file_eof .or. region_iostat < 0) then
                ! Handle EOF with zeros
                ! If one file is not present, neither is used
                call log_info("Module owi_rain: region file EOF. " // &
                    "Inserting blank snapshot.", myproc=myproc, to_screen=myproc == 0)
                rain_tb = 0.0
                region_file_eof = .true.
                return
            else if (region_iostat > 0) then
                ! Format error
                call log_error("Error in module owi_rain: region file " // &
                    "format error", myproc)
            end if
            ! Unit conversion
            if (mmhr2ms) region_data_tb = region_data_tb / 3600000.0
            ! Check if weights need to be recomputed
            if (basin_grid_params_tb /= &
                basin_interpolator_tb%grid%grid_params &
            ) then
                call basin_grid_tb%init(basin_grid_params_tb)
                call basin_interpolator_tb%init( &
                    grid=basin_grid_tb, &
                    lons_interp=lons_interp, &
                    lats_interp=lats_interp &
                )
            end if
            if (region_grid_params_tb /= &
                region_interpolator_tb%grid%grid_params &
            ) then
                call region_grid_tb%init(region_grid_params_tb)
                call region_interpolator_tb%init( &
                    grid=region_grid_tb, &
                    lons_interp=lons_interp, &
                    lats_interp=lats_interp &
                )
            end if
            call apply_owi_interpolator( &
                data_grid=basin_data_tb, &
                interpolator_grid=basin_interpolator_tb, &
                data_interp=rain_tb, & ! Output
                data_subgrid=region_data_tb, &
                interpolator_subgrid=region_interpolator_tb, &
                exclude_negatives=exclude_negatives &
            )
        end if
    end subroutine cycle_snapshots

    subroutine owi_rain_final()
        !! Finalize the OWI_RAIN module and its global variables.
        !!
        !! This means closing the precipitation data files fort.425 and
        !! fort.427 and releasing the allocatable arrays. Call once at the end
        !! of the simulation. Calling it twice, or without a prior call to
        !! owi_rain_init, is harmless.

        implicit none

        logical :: is_open !! Whether a file unit is currently open
        type(owi_grid_interpolator_t) :: empty_interpolator !! Unallocated, used to reset

        ! Close the files left open by owi_rain_init
        inquire(unit=basin_file_unit, opened=is_open)
        if (is_open) close(basin_file_unit)
        inquire(unit=region_file_unit, opened=is_open)
        if (is_open) close(region_file_unit)

        ! Release the module arrays
        if (allocated(lons_interp)) deallocate(lons_interp)
        if (allocated(lats_interp)) deallocate(lats_interp)
        if (allocated(rain_ta)) deallocate(rain_ta)
        if (allocated(rain_tb)) deallocate(rain_tb)

        ! Assigning an unallocated object deallocates the components
        basin_interpolator_tb = empty_interpolator
        region_interpolator_tb = empty_interpolator
    end subroutine owi_rain_final


end module owi_rain
