PROJECT(EXCITON CXX C)
#set(CMAKE_BUILD_TYPE "RELEASE")

cmake_host_system_information(RESULT _host_name QUERY HOSTNAME)
##=====================================================================
if(${_host_name} MATCHES "kelvin2") #QUB KELVIN2 AMD MACHINES
##=====================================================================
#set(CMAKE_CXX_FLAGS "-O0 -g -ansi")
#set(CMAKE_CXX_FLAGS "-O3 -march=znver2 -ansi")
set(CMAKE_CXX_FLAGS "-O3 -march=znver2 -ansi") # -Wall -Wno-error=misleading-indentation")

##VALGRIND FLAGS:
##set(CMAKE_CXX_FLAGS "-O0 -g -fsanitize=memory -march=znver2 -ansi")
##set(CMAKE_CXX_FLAGS "-O0 -g -fsanitize=memory -fsanitize=thread -march=znver2 -ansi")
##set(CMAKE_CXX_FLAGS "-O0 -g -fsanitize=address -fno-omit-frame-pointer -march=znver2 -ansi")
##set(CMAKE_CXX_FLAGS "-O3 -fsanitize=address -fno-omit-frame-pointer -fsanitize=undefined -march=znver2 -ansi") #-Wall ##
##-fsanitize=bounds -fsanitize=undefined -fsanitize=leak -fsanitize=address dataflow,cfi,safe-stack leak
## -stdlib=libstdc++

##FLAGS FOR AMD LIBRARIES:
##  set(BLA_VENDOR "FLAME")
##  set(MATHLFLAGS "-lscalapack -lflang -lpgmath")
##FLAGS FOR INTEL MKL LIBRARY
  set(MATHLFLAGS "-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl  -I${MKLROOT}/include")

##=====================================================================
else(TRUE) # INTEL MACHINES
##=====================================================================
  set(CMAKE_CXX_FLAGS "-DUSEMKL -mkl=cluster -ansi -O3 -g")#-check-pointers=rw")
  set(MATHLFLAGS "-L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl  -I${MKLROOT}/include")
endif(${_host_name} MATCHES "kelvin2")
##=====================================================================
##=====================================================================

MESSAGE("-- _host_name variable is set to: " ${_host_name})
MESSAGE("CMAKE_CXX_FLAGS = " ${CMAKE_CXX_FLAGS} \n)
MESSAGE("CMAKE_EXE_LINKER_FLAGS = " ${CMAKE_EXE_LINKER_FLAGS} \n)

set(VERSION_MAJOR "1")
set(VERSION_MINOR "0")
set(VERSION_PATCH "0")

set(AUTHOR "TCD-QUB")
set(AUTHOR_DETAILS "Charles Patterson for EXCITON; QUB Green group for antimatter extensions")
set(DESCRIPTION "EXCITON21")

cmake_minimum_required(VERSION 2.6)
if(COMMAND cmake_policy)
  cmake_policy(SET CMP0004 NEW)
endif(COMMAND cmake_policy)


##########################
###source files and target name
set(src_dir source)
set(target_name exciton21_mol)
##########################

########################################
## place to put some useful stuff
set(CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/Modules")


#set(BLA_VENDOR "FLAME")
#set(BLA_VENDOR "OpenBLAS")
find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
##########################################
### Make sure the compiler can find include files from our project

find_package(MPI REQUIRED)
set(LIBS ${LIBS} ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})

INCLUDE_DIRECTORIES(${PROJECT_SOURCE_DIR}/include ${MPI_INCLUDE_PATH} ${PROJECT_SOURCE_DIR}/libxc-5.1.5/include)

##
##########################
## external libs

SET(CMAKE_FIND_LIBRARY_SUFFIXES .a)

##find_package(GSL REQUIRED)
##message("Using this linkage for gsl: ${GSL_LIBRARIES} with flags  ${CMAKE_GSL_CXX_FLAGS}\n\n")
##INCLUDE_DIRECTORIES(${GSL_INCLUDE_DIR})
##set(LIBS ${LIBS} ${CMAKE_GSL_CXX_FLAGS} ${GSL_LIBRARIES}
##${PROJECT_SOURCE_DIR}/libxc/lib/libxc.a)
set(LIBS ${LIBS} ${PROJECT_SOURCE_DIR}/libxc-5.1.5/lib/libxc.a)

##set(MKL_LINKER_FLAGS  -L${MKLROOT}/lib/intel64 -lmkl_scalapack_lp64 -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -lmkl_blacs_intelmpi_lp64 -lpthread -lm -ldl  -I${MKLROOT}/include)

set(LIBS ${LIBS} ${BLAS_LIBRARIES} ${LAPACK_LIBRARIES} ${MATHLFLAGS}) ##

message("\nAll libraries: ${LIBS} \n")

##########################

message("MPI_CXX_FOUND ${MPI_CXX_FOUND}\n\nMPI_CXX_COMPILER ${MPI_CXX_COMPILER}\n\nMPI_CXX_COMPILE_FLAGS ${MPI_CXX_COMPILE_FLAGS}\n\nMPI_CXX_INCLUDE_PATH ${MPI_CXX_INCLUDE_PATH}\n\nMPI_CXX_LINK_FLAGS ${MPI_CXX_LINK_FLAGS}\n\nMPI_CXX_LIBRARIES ${MPI_CXX_LIBRARIES}\n\n")


##set the output folder for libs and bins
set(LIBRARY_OUTPUT_PATH ${PROJECT_BINARY_DIR}/lib)
set(EXECUTABLE_OUTPUT_PATH ${PROJECT_BINARY_DIR}/bin)
###################


####################
#~ set_source_files_properties(${src_dir}/SCF2.f90
#~ 	PROPERTIES COMPILE_FLAGS /Od)
####################
##### adds the main target
add_executable(${target_name}
#COMMON SOURCE FILES
${src_dir}/TOOLS.cpp                            #              ${src_dir}/TOOLS.cpp
${src_dir}/ERRORS.cpp                           #              ${src_dir}/ERRORS.cpp
${src_dir}/HEADER.cpp                           #              ${src_dir}/HEADER.cpp
${src_dir}/SETUP_ATOMS.cpp                      #              ${src_dir}/SETUP_ATOMS.cpp
${src_dir}/SETUP_RECIPROCAL_LATTICE.cpp         #              ${src_dir}/SETUP_RECIPROCAL_LATTICE.cpp
${src_dir}/FOURIER_TRANSFORM.cpp                #              ${src_dir}/FOURIER_TRANSFORM.cpp
${src_dir}/IVANIC_RUEDENBERG.cpp                #              ${src_dir}/IVANIC_RUEDENBERG.cpp
${src_dir}/SYMMETRY_ADAPTATION.cpp              #              ${src_dir}/SYMMETRY_ADAPTATION.cpp
${src_dir}/LINEAR_ALGEBRA_UTIL.cpp              #              ${src_dir}/LINEAR_ALGEBRA_UTIL.cpp
${src_dir}/MATRIX_UTIL.cpp                      #              ${src_dir}/MATRIX_UTIL.cpp
${src_dir}/PARALLEL.cpp                         #              ${src_dir}/PARALLEL.cpp
${src_dir}/SCALAPACK.cpp                        #              ${src_dir}/SCALAPACK.cpp
${src_dir}/ROTATIONS_MOLECULE.cpp               #              ${src_dir}/ROTATIONS_MOLECULE.cpp
${src_dir}/SCF_ATOM.cpp                         #              ${src_dir}/SCF_ATOM.cpp
${src_dir}/SCF_MOLECULE.cpp                     #              ${src_dir}/SCF_MOLECULE.cpp
${src_dir}/GW_BSE_MOLECULE.cpp                  #              ${src_dir}/GW_BSE_MOLECULE.cpp
${src_dir}/DENSITY_MATRIX_MOLECULE.cpp          #              ${src_dir}/DENSITY_MATRIX_MOLECULE.cpp
${src_dir}/DENSITY_FITTING_MOLECULE.cpp         #              ${src_dir}/DENSITY_FITTING_MOLECULE.cpp
#${src_dir}/INTEGRALS1.cpp                      #              ${src_dir}/INTEGRALS1.cpp 
${src_dir}/INTEGRALS_2C_MOLECULE.cpp            #              ${src_dir}/INTEGRALS_2C_MOLECULE.cpp
${src_dir}/INTEGRALS_3C_MOLECULE.cpp            #              ${src_dir}/INTEGRALS_3C_MOLECULE.cpp
${src_dir}/INTEGRALS_4C_MOLECULE.cpp            #              ${src_dir}/INTEGRALS_4C_MOLECULE.cpp
${src_dir}/E_COEFFICIENTS.cpp                   #              ${src_dir}/E_COEFFICIENTS.cpp
${src_dir}/MCMURCHIE_DAVIDSON.cpp               #              ${src_dir}/MCMURCHIE_DAVIDSON.cpp
${src_dir}/CARTESIAN_TO_SH.cpp                  #              ${src_dir}/CARTESIAN_TO_SH.cpp
${src_dir}/RECURSION.cpp                        #              ${src_dir}/RECURSION.cpp
${src_dir}/INCOMPLETE_GAMMA.cpp                 #              ${src_dir}/INCOMPLETE_GAMMA.cpp
${src_dir}/OPTICAL_SPECTRUM_MOLECULE.cpp        #              ${src_dir}/OPTICAL_SPECTRUM_MOLECULE.cpp
${src_dir}/PRINT_MOLECULE.cpp                   #              ${src_dir}/PRINT_MOLECULE.cpp
${src_dir}/PAIRS_QUADS.cpp                      #              ${src_dir}/PAIRS_QUADS.cpp
#DISTINCT SOURCE FILES                          #
${src_dir}/SETUP_RUNTIME_MOLECULE.cpp           #              ${src_dir}/SETUP_RUNTIME_CRYSTAL.cpp
${src_dir}/SETUP_SYMMETRY_MOLECULE.cpp          #              
${src_dir}/BUILD_FOCK_MATRIX_MOLECULE.cpp       #              ${src_dir}/BUILD_FOCK_MATRIX.cpp
${src_dir}/ALLOCATE_MEMORY_MOLECULE.cpp         #              ${src_dir}/ALLOCATE_MEMORY.cpp
${src_dir}/PLOTTING_MOLECULE.cpp                #              ${src_dir}/PAIRS_QUADS.cpp
                                                #              ${src_dir}/DENSITY_FITTING_CRYSTAL.cpp
                                                #              ${src_dir}/SCF_CRYSTAL.cpp
                                                #              ${src_dir}/INTEGRALS_2C_CRYSTAL.cpp
                                                #              ${src_dir}/INTEGRALS_3C_CRYSTAL.cpp
                                                #              ${src_dir}/INTEGRALS_4C_CRYSTAL.cpp
						#              ${src_dir}/KPOINTS.cpp
                                                #              ${src_dir}/DENSITY_MATRIX_CRYSTAL.cpp
                                                #              ${src_dir}/TDHF_CRYSTAL.cpp
                                                #              ${src_dir}/SETUP_REAL_LATTICE.cpp
						#              ${src_dir}/OPTICAL_SPECTRUM.cpp
${src_dir}/INPUT_MOLECULE.cpp                   #              ${src_dir}/INPUT_ALL.cpp
#                                               #
#                                               #
#                                               #
#${src_dir}/SETUP_SYMMETRY.cpp                  #              ${src_dir}/INPUT_ALL.cpp
#${src_dir}/SETUP_REAL_LATTICE.cpp              #              #
#                                               #
#${src_dir}/PLOTTING.cpp                        #              ${src_dir}/SETUP_REAL_LATTICE.cpp
#${src_dir}/DFT.cpp                             #              ${src_dir}/SETUP_SYMMETRY.cpp
#${src_dir}/LEBEDEV_LAIKOV.cpp                  #              ${src_dir}/SETUP_CRYSTAL.cpp
#${src_dir}/INTEGRALS_THREE_CENTRE.cpp          #              ${src_dir}/CRYSTAL09.cpp
#${src_dir}/DIELECTRIC_FUNCTION.cpp                           #
#${src_dir}/INTEGRALS_TEST.cpp                                #${src_dir}/PLOTTING.cpp
)                                                             #${src_dir}/DFT.cpp
####################                                          #${src_dir}/LEBEDEV_LAIKOV.cpp
                                                              #${src_dir}/INTEGRALS_THREE_CENTRE.cpp
#################
## sets the linking
link_directories (${LIBRARY_OUTPUT_PATH})
target_link_libraries(${target_name} ${LIBS})
#target_link_libraries(${target_name} algebra ${LIBS})
#target_link_libraries(${target_name} algebra ${LIBS} ${MPI_CXX_LIBRARIES})
#OMPItarget_link_libraries(${target_name} algebra ${LIBS} -static-libgfortran -static-libgcc -static-libstdc++)
#target_link_libraries(${target_name} algebra ${LIBS} -nostartfiles -nostdlib -nodefaultlibs -static-libgfortran -static-libgcc -static-libstdc++)
#################

if(DOXYGEN)
  configure_file(${CMAKE_SOURCE_DIR}/Doxyfile.cmake Doxyfile)
  add_custom_target(docs
    ${DOXYGEN_EXECUTABLE} ${PROJECT_BINARY_DIR}/Doxyfile)
endif(DOXYGEN)
