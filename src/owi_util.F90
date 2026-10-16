module owi_util
    !! Module with routines to read and use OWI-format data files.
    !!
    !! Created by Ashton Cole in October 2026.

    use sizes, only: sz

    implicit none

    type :: owi_grid_params_t
        !! Structure to hold parameters defining a regular lon-lat grid.

        integer :: ilat = 0
        !! The number of horizontal parallels.
        integer :: ilon = 0
        !! The number of vertical meridians.
        real(sz) :: dx = 0.0
        !! The horizontal/longitudnal spacing of the data grid.
        real(sz) :: dy = 0.0
        !! The vertical/latitudnal spacing of the data grid.
        real(sz) :: swlat = 0.0
        !! The latitude of the southwest corner of the data grid.
        real(sz) :: swlon = 0.0
        !! The longitude of the southwest corner of the data grid.

        contains

        procedure, public :: equal => owi_grid_params_equal
        procedure, public :: not_equal => owi_grid_params_not_equal

        generic, public :: operator(==) => equal
        ! generic, public :: operator(.eq.) => equal
        generic, public :: operator(/=) => not_equal
        ! generic, public :: operator(.ne.) => not_equal
    end type owi_grid_params_t

    type :: owi_grid_t
        !! Structure to hold a regular lon-lat grid.

        type(owi_grid_params_t) :: grid_params
        !! The specifications used to construct the grid.
        real(sz), allocatable :: lons(:)
        !! Longitudes for the grid, from west to east.
        real(sz), allocatable :: lats(:)
        !! Latitudes for the grid, from south to north.

        contains

        procedure, public :: init => owi_grid_init
        procedure, public :: equal => owi_grid_equal
        procedure, public :: not_equal => owi_grid_not_equal

        generic, public :: operator(==) => equal
        ! generic, public :: operator(.eq.) => equal
        generic, public :: operator(/=) => not_equal
        ! generic, public :: operator(.ne.) => not_equal
    end type owi_grid_t

    type :: owi_grid_interpolator_t
        !! Structure to hold data for interpolation.

        type(owi_grid_t) :: grid
        !! The grid interpolated upon.
        real(sz), allocatable :: interp_coeff(:, :)
        !! Resultant interpolation weights for the point series, representing
        !! the interpolation point's position in an OWI box. Dimension is
        !! (num_points, 4), with counter-clockwise ordering SW, SE, NE, NW. If
        !! the point is in the subgrid region, these weights will be zero.
        integer, allocatable :: sw_i(:)
        !! First, horizontal, longitudnal index for the SW corner of the OWI
        !! grid box containing the interpolated point. By convention, the east
        !! corners are at i + 1.
        integer, allocatable :: sw_j(:)
        !! Second, vertical, latitudnal index for the SW corner of the OWI
        !! grid box containing the interpolated point. By convention, the north
        !! corners are at j + 1.

        contains

        procedure, public :: init => owi_grid_interpolator_init
    end type owi_grid_interpolator_t

    contains

    subroutine get_owi_header(unit, iostat, title, start_date, end_date)
        !! Read the header of an OWI-format data file and return the metadata.
        !!
        !! The format for a header is (a55,i10,5x,i10)
        !! - 55 characters for a descriptive title
        !! - 10 characters for an integer start date YYYYMMDDHH
        !! - 5 ignored characters
        !! - 10 characters for an integer end date YYYYMMDDHH
        !!
        !! This assumes that the file has already been opened and that the
        !! buffer points to the beginning of the file. After execution, the
        !! buffer points to the first data record in the file.

        implicit none

        integer, intent(in) :: unit
        !! The Fortran file unit number used to open the OWI data file.
        integer, intent(out) :: iostat
        !! The status of the file read operation. This marks if the read was
        !! successful, reached the end of the file, or encountered a format
        !! error.
        character(len=*), optional, intent(out) :: title
        !! The title associated with the data file. Ususally, this describes the
        !! data contents, e.g., "Large Grid Pressure at 10m".
        integer, optional, intent(out) :: start_date
        !! The start date for the data set.
        integer, optional, intent(out) :: end_date
        !! The end date for the data set.

        ! Dummy parameters
        character(len=55) :: the_title
        integer :: the_start_date
        integer :: the_end_date

        read(unit, "(a55,i10,5x,i10)", iostat=iostat) the_title, &
            the_start_date, the_end_date

        ! Return requested information
        if (present(title)) title = trim(the_title)
        if (present(start_date)) start_date = the_start_date
        if (present(end_date)) end_date = the_end_date
    end subroutine get_owi_header

    subroutine get_owi_record_scalar(unit, iostat, grid_params, dt, data)
        !! Read a single time slice of an OWI-format data file for a scalar
        !! and return the metadata.
        !!
        !! The format for a record is (5x, i4, 6x, i4, 3x, f6.0, 3x, f6.0, 6x,
        !! 8.0f, 6x, 8.0f, 3x, i12) and (8f10.0)
        !! - 5 ignored characters for "iLat="
        !! - 4 characters for an integer number of latitude parallels
        !! - 6 ignored characters for "iLong="
        !! - 4 characters for an integer number of longitude meridians
        !! - 3 ignored characters for "DX="
        !! - 6 characters for a float longitudnal meridian spacing
        !! - 3 ignored characters for "DY="
        !! - 6 characters for a float latitudnal parallel spacing
        !! - 6 ignored characters for "SWLat="
        !! - 8 characters for a float southwest corner latitude for the grid
        !! - 6 ignored characters for "SWLon="
        !! - 8 characters for a float southwest corner longitude for the grid
        !! - 3 ignored characters for "DT="
        !! - 12 characters for an integer value of the time in YYYYMMDDHHMM
        !! - A sequence of data values for iLat * iLon grid points.
        !!
        !! This assumes that the file has already been opened and that the
        !! buffer points to the beginning of a data record. After execution, the
        !! buffer points to the next data record in the file. If the buffer is
        !! misaligned or reaches the end of file prematurely, this will be
        !! reflected in the iostat parameter.

        implicit none

        integer, intent(in) :: unit
        !! The Fortran file unit number used to open the OWI data file.
        integer, intent(out) :: iostat
        !! The status of the file read operation. This marks if the read was
        !! successful, reached the end of the file, or encountered a format
        !! error.
        type(owi_grid_params_t), optional, intent(out) :: grid_params
        !! The parameters for the grid on which the snapshot is defined. This
        !! usually is, but is not necessarily, constant for all time.
        integer(8), optional, intent(out) :: dt
        !! The time of the snapshot.
        real(sz), allocatable, optional, intent(out) :: data(:, :)
        !! The grid data. Of dimension (ilon, ilat).

        ! Dummy parameters
        integer :: the_ilat
        integer :: the_ilon
        real(sz) :: the_dx
        real(sz) :: the_dy
        real(sz) :: the_swlat
        real(sz) :: the_swlon
        integer(8) :: the_dt
        real(sz) :: dummy_real

        ! Loop parameters
        integer :: i
        integer :: j

        ! Read the time step header
        read(unit, "(5x,i4,6x,i4,3x,f6.0,3x,f6.0,6x,f8.0,6x,f8.0,3x,i12)", &
            iostat=iostat) &
            the_ilat, the_ilon, the_dx, the_dy, the_swlat, the_swlon, the_dt
        if (iostat /= 0) return

        ! Read the data array
        if (present(data)) then ! Read into data
            if (allocated(data)) deallocate(data)
            allocate(data(the_ilon, the_ilat))
            data = 0.0
            read(unit, "(8f10.0)", iostat=iostat) ((data(i, j), &
                i = 1, the_ilon), j = 1, the_ilat)
        else ! Discard by reading into a dummy variable
            read(unit, "(8f10.0)", iostat=iostat) ((dummy_real, &
                i = 1, the_ilon), j = 1, the_ilat)
        end if

        ! Return requested information
        if (present(grid_params)) then
            grid_params%ilat = the_ilat
            grid_params%ilon = the_ilon
            grid_params%dx = the_dx
            grid_params%dy = the_dy
            grid_params%swlat = the_swlat
            grid_params%swlon = the_swlon
        end if
        if (present(dt)) dt = the_dt
        ! data already handled
    end subroutine get_owi_record_scalar

    pure subroutine apply_owi_interpolator(data_grid, interpolator_grid, &
        data_interp, data_subgrid, interpolator_subgrid, exclude_negatives)
        !! Apply interpolation coefficients for OWI data at longitude-latitude
        !! cells onto a series of longitudes and latitudes.
        !!
        !! The interpolated value at a given longitude and latitude is a
        !! weighted sum of the values at the four corners around that point.
        !! This routine provides the weights for each point, applied to those
        !! values.
        !!
        !! The general formula, ignoring region nesting, is as follows. The
        !! ordering for weights is counter-clockwise: SW, SE, NE, NW.
        !!
        !! data(i) = interp_coeff(i, 1) * data_grid(sw_i(i), sw_j(i)) +
        !! interp_coeff(i, 2) * data_grid(sw_i(i) + 1, sw_j(i)) +
        !! interp_coeff(i, 3) * data_grid(sw_i(i) + 1, sw_j(i) + 1) +
        !! interp_coeff(i, 4) * data_grid(sw_i(i), sw_j(i) + 1)
        !!
        !! This optionally permits a nested subgrid for higher-resolution data.
        !! At points in the subgrid region, the subgrid values overwrite the
        !! coarser ones.

        implicit none

        real(sz), intent(in) :: data_grid(:, :)
        !! Data values at each (lon, lat) point on the lon-lat grid, from
        !! southwest to northeast. Dimension is (num_lons, num_lats).
        class(owi_grid_interpolator_t), intent(in) :: interpolator_grid
        !! Interpolator object for the main grid data.
        real(sz), intent(out) :: data_interp(:)
        !! The interpolated data at points of interest. The array is assumed to
        !! have already been allocated for efficiency.
        real(sz), optional, intent(in) :: data_subgrid(:, :)
        !! Data values at each (lon, lat) point on the subgrid, from
        !! southwest to northeast. Dimension is (num_sublons, num_sublats).
        class(owi_grid_interpolator_t), optional, intent(in) :: interpolator_subgrid
        !! Optional interpolator object for the subgrid data.
        logical, optional, intent(in) :: exclude_negatives
        !! Whether to mask negative grid values to zero. Default is false.

        ! For defaulting
        logical :: the_exclude_negatives !! Whether to mask negatives

        ! Local variables
        real(sz), dimension(4) :: vals !! Interpolation values for a point
        real(sz), dimension(4) :: weights !! Interpolation weights for a point
        integer :: i !! For grid indices
        integer :: j !! For grid indices
        integer :: npt !! Number of interpolation points

        ! Loop indices
        integer :: pt !! Index of a point to interpolate upon

        npt = size(interpolator_grid%interp_coeff, dim=1)

        ! Determine whether to exclude negatives
        if (present(exclude_negatives)) then
            the_exclude_negatives = exclude_negatives
        else
            the_exclude_negatives = .false.
        end if

        ! First, apply the interpolation for points in the grid
        do pt = 1, npt
            ! Save indices locally for checks
            i = interpolator_grid%sw_i(pt)
            j = interpolator_grid%sw_j(pt)

            ! Check that the indices are valid to set data_vals
            ! The assumption is that i and j, if non-zero, are set in bounds
            ! Otherwise, default is 0
            if (i > 0 .and. j > 0) then
                vals(1) = data_grid(i, j) ! SW
                vals(2) = data_grid(i + 1, j) ! SE
                vals(3) = data_grid(i + 1, j + 1) ! NE
                vals(4) = data_grid(i, j + 1) ! NW
            else
                vals = 0.0
            end if

            ! Set interpolation weights
            weights(1) = interpolator_grid%interp_coeff(pt, 1)
            weights(2) = interpolator_grid%interp_coeff(pt, 2)
            weights(3) = interpolator_grid%interp_coeff(pt, 3)
            weights(4) = interpolator_grid%interp_coeff(pt, 4)

            ! Check if negative values are excluded
            ! In such a case, choose an appropriate resolution
            if (the_exclude_negatives) then
                where (vals < 0.0)
                    ! Recognize that it is a flag, not a true negative
                    ! Mask to 0
                    vals = 0.0
                end where
            end if

            ! Apply the weighted sum
            data_interp(pt) = sum(vals * weights)
        end do

        ! Optionally, apply the interpolation for points in the subgrid
        if (present(data_subgrid) .and. present(interpolator_subgrid)) then
            do pt = 1, npt
                ! Save the indices locally for checks
                i = interpolator_subgrid%sw_i(pt)
                j = interpolator_subgrid%sw_j(pt)

                ! Check that the indices are valid to set data_vals
                ! The assumption is that i and j, if non-zero, are set in bounds
                ! Otherwise, do not set the point value
                ! Cycle to the next iteration
                ! If it is not in the subgrid
                ! Either it is set by the main grid
                ! Or it is set to 0
                if (i > 0 .and. j > 0) then
                    vals(1) = data_subgrid(i, j) ! SW
                    vals(2) = data_subgrid(i + 1, j) ! SE
                    vals(3) = data_subgrid(i + 1, j + 1) ! NE
                    vals(4) = data_subgrid(i, j + 1) ! NW
                else
                    cycle ! Do not change the value
                end if

                ! Set interpolation weights
                weights(1) = interpolator_subgrid%interp_coeff(pt, 1)
                weights(2) = interpolator_subgrid%interp_coeff(pt, 2)
                weights(3) = interpolator_subgrid%interp_coeff(pt, 3)
                weights(4) = interpolator_subgrid%interp_coeff(pt, 4)

                ! Check if negative values are excluded
                ! In such a case, choose an appropriate resolution
                if (the_exclude_negatives) then
                    where (vals < 0.0)
                        ! Recognize that it is a flag, not a true negative
                        ! Mask to 0
                        vals = 0.0
                    end where
                end if

                ! Apply the weighted sum
                ! Overwrite any existing value
                data_interp(pt) = sum(vals * weights)
            end do
        end if
    end subroutine apply_owi_interpolator

    subroutine get_owi_precipitation_parameters(unit, iostat, num_sets, snapshot_offset)
        !! Read a precipitation parameter file and return the parameters.
        !!
        !! This assumes that the file has already been opened and that the
        !! buffer points to the beginning of the file. After execution, the
        !! buffer points to the end of the file.

        implicit none

        integer, intent(in) :: unit
        !! The Fortran file unit number used to open the OWI data file.
        integer, intent(out) :: iostat
        !! The status of the file read operation. This marks if the read was
        !! successful, reached the end of the file, or encountered a format
        !! error.
        integer, intent(out) :: num_sets
        !! The number of overlapping nested data files used.
        integer, intent(out) :: snapshot_offset
        !! Adjustment of the timing of the precipitation data. A positive number
        !! means that blank frames are inserted at the beginning of the run. A
        !! negative number means that frames are skipped.

        read(unit, *, iostat=iostat) num_sets
        if (iostat /= 0) then
            num_sets = 0
            snapshot_offset = 0
            return
        end if
        read(unit, *, iostat=iostat) snapshot_offset
        if (iostat /= 0) then
            snapshot_offset = 0
            return
        end if
    end subroutine get_owi_precipitation_parameters

    pure subroutine owi_grid_init(self, grid_params)
        !! Constructor for owi_grid_t.

        implicit none

        class(owi_grid_t), intent(inout) :: self
        !! The derived type.
        class(owi_grid_params_t), intent(in) :: grid_params
        !! A specification for a rectangular lon-lat grid.

        ! Loop variables
        integer :: i

        self%grid_params = grid_params
        if (allocated(self%lons)) deallocate(self%lons)
        allocate(self%lons(grid_params%ilon))
        self%lons = [(grid_params%swlon + (i - 1) * grid_params%dx, &
            i = 1, grid_params%ilon)]
        if (allocated(self%lats)) deallocate(self%lats)
        allocate(self%lats(grid_params%ilat))
        self%lats = [(grid_params%swlat + (i - 1) * grid_params%dy, &
            i = 1, grid_params%ilat)]
    end subroutine owi_grid_init

    pure subroutine owi_grid_interpolator_init(self, grid, lons_interp, &
        lats_interp)
        !! Constructor for owi_grid_interpolator_t.
        !!
        !! This computes interpolation coefficients for OWI data at
        !! longitude-latitude cells onto a series of longitudes and latitudes.
        !!
        !! The interpolated value at a given longitude and latitude is a
        !! weighted sum of the values at the four corners around that point.
        !! This routine provides the weights for each point, applied to those
        !! values.
        !!
        !! The general formula, ignoring region nesting, is as follows. The
        !! ordering for weights is counter-clockwise: SW, SE, NE, NW.
        !!
        !! data(i) = interp_coeff(i, 1) * data_grid(sw_i(i), sw_j(i)) +
        !! interp_coeff(i, 2) * data_grid(sw_i(i) + 1, sw_j(i)) +
        !! interp_coeff(i, 3) * data_grid(sw_i(i) + 1, sw_j(i) + 1) +
        !! interp_coeff(i, 4) * data_grid(sw_i(i), sw_j(i) + 1)
        !!
        !! This optionally permits a nested subgrid for higher-resolution data.
        !! At points in the subgrid region, the subgrid values overwrite the
        !! coarser ones.

        implicit none

        class(owi_grid_interpolator_t) , intent(inout) :: self
        !! The derived type.
        class(owi_grid_t), intent(in) :: grid
        !! The grid being used for interpolation.
        real(sz), intent(in) :: lons_interp(:)
        !! The longitudes of points to be interpolated at.
        real(sz), intent(in) :: lats_interp(:)
        !! The latitudes of points to be interpolated at.

        ! Local variables
        integer :: npt !! Number of interpolation points
        real(sz) :: x !! Unit square coordinate
        real(sz) :: y !! Unit square coordinate

        ! Loop variables
        integer :: pt !! Index of a point to interpolate upon
        integer :: i !! Longitude meridian index
        integer :: j !! Latitude parallel index

        npt = min(size(lons_interp), size(lats_interp))

        self%grid = grid
        if (allocated(self%sw_i)) deallocate(self%sw_i)
        allocate(self%sw_i(npt))
        self%sw_i = 0
        if (allocated(self%sw_j)) deallocate(self%sw_j)
        allocate(self%sw_j(npt))
        self%sw_j = 0
        if (allocated(self%interp_coeff)) deallocate(self%interp_coeff)
        allocate(self%interp_coeff(npt, 4))
        self%interp_coeff = 0.0

        ! For each point
        ! - Find the cell its in
        ! - Note the index
        ! - Compute interpolation weights
        ! - If not inside any, indices and weights should all be 0
        do pt = 1, npt
            ! First search for longitude/horizontal
            longitude_search: do i = 1, grid%grid_params%ilon - 1
                if ( &
                    grid%lons(i) <= lons_interp(pt) &
                    .and. lons_interp(pt) < grid%lons(i + 1) &
                ) then
                    self%sw_i(pt) = i
                    exit longitude_search ! No need to keep looking
                end if
            end do longitude_search

            ! In a separate loop, search for latitude/vertical
            latitude_search: do j = 1, grid%grid_params%ilat - 1
                if ( &
                    grid%lats(j) <= lats_interp(pt) &
                    .and. lats_interp(pt) < grid%lats(j + 1) &
                ) then
                    self%sw_j(pt) = j
                    exit latitude_search ! No need to keep looking
                end if
            end do latitude_search

            ! If one of the indices is zero, make everything zero
            ! Otherwise, compute the weights
            if (self%sw_i(pt) == 0 .or. self%sw_j(pt) == 0) then
                self%sw_i(pt) = 0
                self%sw_j(pt) = 0
                self%interp_coeff(pt, :) = 0.0
            else
                ! Normalize to a unit square
                x = (lons_interp(pt) - grid%lons(i)) &
                    / (grid%lons(i + 1) - grid%lons(i))
                y = (lats_interp(pt) - grid%lats(j)) &
                    / (grid%lats(j + 1) - grid%lats(j))

                ! Standard weights for bilinear interpolation
                self%interp_coeff(pt, 1) = (1 - x) * (1 - y) ! SW
                self%interp_coeff(pt, 2) = x * (1 - y) ! SE
                self%interp_coeff(pt, 3) = x * y ! NE
                self%interp_coeff(pt, 4) = (1 - x) * y ! NW
            end if
        end do
    end subroutine owi_grid_interpolator_init

    pure logical function owi_grid_params_equal(self, other) &
        result(res)
        !! Whether two grid parameter specifications are equal. This is true
        !! if all parameters are the same, and false otherwise.

        implicit none

        class(owi_grid_params_t), intent(in) :: self
        !! The first grid parameter specification.
        class(owi_grid_params_t), intent(in) :: other
        !! The second grid parameter specification.

        res = self%ilat == other%ilat &
            .and. self%ilon == other%ilon &
            .and. self%dx == other%dx &
            .and. self%dy == other%dy &
            .and. self%swlon == other%swlon &
            .and. self%swlat == other%swlat
    end function owi_grid_params_equal

    pure logical function owi_grid_params_not_equal(self, &
        other) result(res)
        !! Whether two grid parameter specifications are not equal. This is true
        !! if all parameters are the same, and false otherwise.

        implicit none

        class(owi_grid_params_t), intent(in) :: self
        !! The first grid parameter specification.
        class(owi_grid_params_t), intent(in) :: other
        !! The second grid parameter specification.

        res = .not. owi_grid_params_equal(self, other)
    end function owi_grid_params_not_equal

    pure logical function owi_grid_equal(self, other) result(res)
        !! Whether two grids are equal. This is true if the parameters used to
        !! construct them are the same, and false otherwise.

        implicit none

        class(owi_grid_t), intent(in) :: self
        !! The first grid.
        class(owi_grid_t), intent(in) :: other
        !! The second grid.

        res = self%grid_params == other%grid_params
    end function owi_grid_equal

    pure logical function owi_grid_not_equal(self, other) result(res)
        !! Whether two grids are not equal. This is true if the parameters used
        !! to construct them are not the same, and false otherwise.

        implicit none

        class(owi_grid_t), intent(in) :: self
        !! The first grid.
        class(owi_grid_t), intent(in) :: other
        !! The second grid.

        res = .not. owi_grid_equal(self, other)
    end function owi_grid_not_equal

end module owi_util
