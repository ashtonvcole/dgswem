#!/usr/bin/env python3
"""
Generate inputs for the draining square case.

A closed, flat-bottomed square of water (1 m deep, initially at rest) drains
under a uniform evaporative sink of -1e-7 m/s, supplied through an OWI-style
2 x 2 precipitation grid (fort.27 / fort.425). There are no open boundaries,
no wind, no tides, and no flux through the walls, so the exact solution is
spatially uniform:

    zeta(t) = -E * t,    E = 1e-7 m/s,

until the water depth reaches the wetting/drying threshold H0 = 0.1 m.

Run from this directory:
    python3 make_inputs.py
"""

from pathlib import Path

# Domain: 0.01 degree square near 29 N, 95 W (about 1.1 km side)
LON0, LAT0 = -95.0, 29.0
SIDE = 0.01
NX = NY = 5  # nodes per side, giving 4 x 4 cells and 32 triangles

BED_DEPTH = 1.0  # m, positive below datum; initial water depth at rest
EVAP = -1.0e-7  # m/s, negative means a sink (evaporation)

# OWI grid: 2 x 2 points, extended beyond the mesh so every node lies
# inside the single cell [swlon, swlon + dx) x [swlat, swlat + dy)
OWI_DX = OWI_DY = 0.04
OWI_SWLON = LON0 - 0.02
OWI_SWLAT = LAT0 - 0.02
N_DAYS = 120
OWI_INTERVAL = 86400  # s, one snapshot per day
N_SNAPSHOTS = N_DAYS + 2  # one more than the run needs, so tb is always valid

DT = 20.0  # s
NSPOOL = int(86400 / DT)  # one global output per day


def write_fort14(path: Path) -> None:
    lons = [LON0 + SIDE * i / (NX - 1) for i in range(NX)]
    lats = [LAT0 + SIDE * j / (NY - 1) for j in range(NY)]

    def node(i, j):
        return i * NY + j  # 0-based, matching the sloping-beach generator

    lines = ["Draining square, flat bed, uniform evaporative sink"]
    n_elem = (NX - 1) * (NY - 1) * 2
    n_node = NX * NY
    lines.append(f"{n_elem} {n_node}")
    for i in range(NX):
        for j in range(NY):
            lines.append(
                f"{node(i, j) + 1:>12d} {lons[i]:12.8f} {lats[j]:12.8f} "
                f"{BED_DEPTH:12.8f}"
            )

    e = 1
    for i in range(NX - 1):
        for j in range(NY - 1):
            a, b = node(i, j) + 1, node(i + 1, j) + 1
            c, d = node(i, j + 1) + 1, node(i + 1, j + 1) + 1
            # Counterclockwise triangles
            lines.append(f"{e:>12d} 3 {a:d} {b:d} {c:d}")
            lines.append(f"{e + 1:>12d} 3 {c:d} {b:d} {d:d}")
            e += 2

    # Perimeter loop, counterclockwise from the southwest corner
    perimeter = []
    perimeter += [node(0, j) for j in range(NY)]  # west side, south to north
    perimeter += [node(i, NY - 1) for i in range(1, NX)]  # north side
    perimeter += [node(NX - 1, j) for j in range(NY - 2, -1, -1)]  # east side
    perimeter += [node(i, 0) for i in range(NX - 2, 0, -1)]  # south side

    lines.append("0")  # open boundaries (none)
    lines.append("0")
    lines.append("1")  # land boundaries: one closed loop
    lines.append(f"{len(perimeter):d}")
    lines.append(f"{len(perimeter):d} 0")  # normal-flow (no flux) boundary
    lines += [f"{p + 1:d}" for p in perimeter]

    path.write_text("\n".join(lines) + "\n")


FORT15 = f"""\
Draining square                          ! 32 CHARACTER ALPHANUMERIC RUN DESCRIPTION
draining_square                          ! 24 CHARACTER ALPHANUMERIC RUN IDENTIFICATION
1 20.0 1 50 1000.0                       ! NFOVER - NONFATAL ERROR OVERRIDE OPTION
0                                        ! NABOUT - ABREVIATED OUTPUT OPTION PARAMETER
1000                                     ! NSCREEN - OUTPUT TO UNIT 6 PARAMETER
0                                        ! IHOT - HOT START OPTION PARAMETER (0 = cold start)
2                                        ! ICS - COORDINATE SYSTEM OPTION PARAMETER (lon/lat)
0                                        ! IM - MODEL RUN TYPE: 0 = 2DDI
1                                        ! NOLIBF - NONLINEAR BOTTOM FRICTION OPTION
2                                        ! NOLIFA - FINITE AMPLITUDE TERMS, wetting and drying on
1                                        ! NOLICA - OPTION TO INCLUDE CONVECTIVE ACCELERATION TERMS
1                                        ! NOLICAT - OPTION TO CONSIDER TIME DERIVATIVE OF CONV ACC TERMS
0                                        ! NWP - Number of nodal attributes
0                                        ! NCOR - VARIABLE CORIOLIS IN SPACE OPTION PARAMETER
0                                        ! NTIP - TIDAL POTENTIAL OPTION PARAMETER
0                                        ! NWS - WIND STRESS AND BAROMETRIC PRESSURE OPTION PARAMETER
0                                        ! NRAMP - RAMP FUNCTION OPTION
9.81                                     ! G - ACCELERATION DUE TO GRAVITY
0.05000                                  ! TAU0 - WEIGHTING FACTOR IN GWCE
{DT:f}                                   ! DT - TIME STEP (IN SECONDS)
0.00000                                  ! STATIM - STARTING SIMULATION TIME IN DAYS
0.00000                                  ! REFTIME - REFERENCE TIME (IN DAYS)
{N_DAYS:.5f}                             ! RNDAY - TOTAL LENGTH OF SIMULATION (IN DAYS)
1.00000                                  ! DRAMP - DURATION OF RAMP FUNCTION (IN DAYS)
0.800000 0.200000 0.000000               ! TIME WEIGHTING FACTORS FOR THE GWCE EQUATION
0.100000 2 10 0.010000                   ! H0, NODEDRYMIN, NODEWETMIN, VELMIN
{LON0 + SIDE / 2:.6f} {LAT0 + SIDE / 2:.6f}   ! SLAM0, SFEA0 - projection center (lon, lat)
0.002500                                 ! FFACTOR - 2DDI BOTTOM FRICTION COEFFICIENT
0.0                                      ! ESLM - HORIZONTAL EDDY VISCOSITY
0.000000                                 ! CORI - CONSTANT CORIOLIS COEFFICIENT
0                                        ! NTIF - NUMBER OF TIDAL POTENTIAL CONSTITUENTS
0                                        ! NBFR - NUMBER OF PERIODIC FORCING FREQUENCIES (none)
110                                      ! ANGINN - MINIMUM ANGLE FOR TANGENTIAL FLOW
0 0.000000 0.000000 0                    ! NOUTE, TOUTSE, TOUTFE, NSPOOLE - FORT 61 OPTIONS
0                                        ! NSTAE - NUMBER OF ELEVATION RECORDING STATIONS
0 0.000000 0.000000 0                    ! NOUTV, TOUTSV, TOUTFV, NSPOOLV - FORT 62 OPTIONS
0                                        ! NSTAV - NUMBER OF VELOCITY RECORDING STATIONS
1 0.000000 {float(N_DAYS):.6f} {NSPOOL:d}     ! NOUTGE, TOUTSGE, TOUTFGE, NSPOOLGE - GLOBAL ELEVATION OUTPUT (fort.63)
0 0.000000 0.000000 0                    ! NOUTGV, TOUTSGV, TOUTFGV, NSPOOLGV - GLOBAL VELOCITY OUTPUT
0                                        ! NHARF - NUMBER OF FREQUENCIES IN HARMONIC ANALYSIS
0.000000 0.000000 0 0.000000             ! THAS, THAF, NHAINC, FMV - HARMONIC ANALYSIS PARAMETERS
0 0 0 0                                  ! NHASE, NHASV, NHAGE, NHAGV
0 0                                      ! NHSTAR, NHSINC - HOT START FILE GENERATION PARAMETERS
1 0 1.000000000E-007 25                  ! ITITER, ISLDIA, CONVCR, ITMAX - ALGEBRAIC SOLUTION PARAMETERS
"""

FORT_DG = """\
! dgswem fort.dg in keyword format, draining square case

dgswe = 1
dg_to_cg = 1
modal_ic = 0

! p enrichment options (fixed linear p)
padapt = 0
pflag = 2
gflag = 1
dis_tol = 8
pl = 1
ph = 1
px = 1
slimit = 0.00005
plimit = 10
k = 1
ks = 0.5
L = 2
nelem = 3564104

! flux option
fluxtype = 1          ! 1=Roe

! timestepping options
rk_stage = 2
rk_order = 2

! hot start options
dghot = 0
dghotspool = 43200

! slopelimiting options
slopeflag = 5
weight = 1

! sediment options
sedflag = 0
porosity = 0.0001
sevdm = 0.00001
layers = 1
sed_equationX = (ZE_ROE+bed_ROE)**-1 *QX_ROE
sed_equationY = (ZE_ROE+bed_ROE)**-1 *QY_ROE

! chemistry options
rxn_rate = 1

! artificial diffusion options
artdif = 0
kappa = -1.0
s0 = 0.0
uniform_dif = 2.5e-6
tune_by_hand = 0

! precipitation options
rainfall = 2         ! OWI fort.27 / fort.425
rtiminc = 86400.0    ! snapshot interval (s)
rexneg = 0           ! 0 = keep negative values (sink); see owi_rain docs
rusemmhr = 0         ! 0 = file values are already in m/s
"""


def write_fort27(path: Path) -> None:
    # One basin grid, no region grid, no blank or skipped snapshots
    path.write_text("1\n0\n")


def write_fort425(path: Path) -> None:
    # Dates are not used by the reader (snapshots are assumed regular), so
    # they are written as zeros rather than invented timestamps
    title = "Uniform evaporation over draining square".ljust(55)
    lines = [f"{title}{0:10d}     {0:10d}"]
    for _ in range(N_SNAPSHOTS):
        dt = 0
        lines.append(
            f"iLat={2:4d}"
            f"iLong={2:4d}"
            f"DX={OWI_DX:6.2f}"
            f"DY={OWI_DY:6.2f}"
            f"SWLat={OWI_SWLAT:8.3f}"
            f"SWLon={OWI_SWLON:8.3f}"
            f"DT={dt:12d}"
        )
        # 2 x 2 grid, uniform value, 8 fields per line per the OWI format
        lines.append("".join(f"{EVAP:10.3E}" for _ in range(4)))
    path.write_text("\n".join(lines) + "\n")


def main() -> None:
    here = Path(__file__).resolve().parent
    write_fort14(here / "fort.14")
    (here / "fort.15").write_text(FORT15)
    (here / "fort.dg").write_text(FORT_DG)
    write_fort27(here / "fort.27")
    write_fort425(here / "fort.425")


if __name__ == "__main__":
    main()
