#!/bin/bash

gfortran -ffixed-line-length-132 -o owi_rain_driver \
    ../../src/sizes.F \
    ../../src/logging.F90 \
    ../../src/owi_util.F90 \
    ../../src/owi_rain.F90 \
    driver.F90
