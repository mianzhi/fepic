CMAKE_MINIMUM_REQUIRED(VERSION 3.11)
PROJECT(fepic)
ENABLE_LANGUAGE(Fortran C)

# build FOSolverS in-place
INCLUDE(FetchContent)
FETCHCONTENT_DECLARE(
  fosolvers
  GIT_REPOSITORY https://github.com/mianzhi/fosolvers.git
  GIT_TAG master
)
FETCHCONTENT_MAKEAVAILABLE(fosolvers)

# compiler & linker arguments
SET(CMAKE_Fortran_COMPILER mpif90)
SET(CMAKE_Fortran_FLAGS "-Wall -Ofast -march=native -fopenmp")
SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_C_FLAGS "-Wall -Ofast -march=native -fopenmp -std=c11")
SET(CMAKE_EXE_LINKER_FLAGS "-Wall -Ofast -march=native -fopenmp")
SET(CMAKE_SHARED_LINKER_FLAGS "-Wall -Ofast -march=native -fopenmp")

# module, executable & lib output directories
SET(CMAKE_Fortran_MODULE_DIRECTORY ${PROJECT_BINARY_DIR}/include)
SET(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)

# include/library directories
INCLUDE_DIRECTORIES(
  ${PROJECT_BINARY_DIR}/include
)

# external packages
SET(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_SOURCE_DIR}/CMakeModules)
# BLAS
FIND_PACKAGE(BLAS)
# LAPACK
FIND_PACKAGE(LAPACK)
# SUNDIALS
FIND_PACKAGE(SUNDIALS)
IF(SUNDIALS_FOUND)
  INCLUDE_DIRECTORIES(${SUNDIALS_INCLUDE_DIRS})
ENDIF(SUNDIALS_FOUND)
# umfpack
FIND_PACKAGE(UMFPACK)
# matheval
FIND_PACKAGE(MATHEVAL)
IF(MATHEVAL_FOUND)
  INCLUDE_DIRECTORIES(${MATHEVAL_INCLUDE_DIR})
ENDIF(MATHEVAL_FOUND)
# METIS
FIND_PACKAGE(METIS)

# subdirectories
ADD_SUBDIRECTORY(src)
