cmake_minimum_required(VERSION 3.5)
project (Coarrays)
enable_language (Fortran)

if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE "Release")
endif()


if(CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    set(dialect "-ffree-form -std=f2008 -fimplicit-none -cpp -Wall -pedantic")
    set(bounds "-fbounds-check")
    set(realeaseopts "-O2  -pthread -funroll-all-loops")
    set(debugopts "-O0 -g3 -fbacktrace -fcheck=all -fdump-fortran-optimized")
elseif(CMAKE_Fortran_COMPILER_ID MATCHES "Intel")
    set(dialect "-std08  -module . -implicitnone -fpp -Wall -pedantic")
    set(bounds " -check bounds")
    set(realeaseopts "-funroll-all-loops -O2")
    set(debugopts "-O0 -g -traceback -check all -debug extended -debug-parameters all")
endif()

set(CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} ${dialect}")
set(CMAKE_Fortran_FLAGS_DEBUG "${CMAKE_Fortran_FLAGS_DEBUG} ${bounds} ${debugopts}")
set(CMAKE_Fortran_RELEASE_FLAGS "${CMAKE_Fortran_RELEASE_FLAGS} ${realeaseopts}")

set(SRC "${PROJECT_SOURCE_DIR}/src")
set(TESTS "${PROJECT_SOURCE_DIR}/tests")
set(mm_SRC
  ${SRC}/matrix_lib.F90
  ${SRC}/run_lib.F90
  )
set(mm_MAIN main.F90)

# --------------------------------------------------------------------
# Executables
# --------------------------------------------------------------------

add_executable("${PROJECT_NAME}" ${mm_MAIN} ${mm_SRC})