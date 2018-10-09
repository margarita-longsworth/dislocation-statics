cmake_minimum_required(VERSION 3.9)
project(montecarlo LANGUAGES CXX)

### General Build System Setup #################################################

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# set install prefix
get_filename_component(PROJECT_ROOT
                       ${CMAKE_SOURCE_DIR}
                       DIRECTORY ABSOLUTE PATH)
set(CMAKE_INSTALL_PREFIX ${PROJECT_ROOT} CACHE ABSOLUTE PATH FORCE)

# set default build type
if (CMAKE_GENERATOR STREQUAL "Unix Makefiles" AND NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE RelWithDebInfo)
endif()

# check build type and override RELWITHDEBINFO flags
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g -O3")
string(TOUPPER ${CMAKE_BUILD_TYPE} _cmake_build_type_uppercase)


### Monte Carlo Library Component ##############################################

set(MONTECARLO_LIB MontecarloAPI-${CMAKE_CXX_COMPILER_ID})

add_library(${MONTECARLO_LIB}
        acceptanceCheck.cpp
        acceptanceCheck.h
        Cell.cpp
        Cell.h
        Domain.cpp
        Domain.h
        localMD.cpp
        localMD.h
        Particle.cpp
        Particle.h
        TrialMove.cpp
        TrialMove.h
        UtilityMethods.cpp
        UtilityMethods.h
)
target_compile_features(${MONTECARLO_LIB} PUBLIC cxx_std_11 cxx_nullptr)


### Monte Carlo Executable #####################################################

get_filename_component(CXX_COMPILER_NAME ${CMAKE_CXX_COMPILER} NAME)
set(MONTECARLO_EXECUTABLE montecarlo_${CMAKE_SYSTEM_PROCESSOR}_${CXX_COMPILER_NAME})

add_executable(${MONTECARLO_EXECUTABLE}
        IntVector3d.h
        Job.cpp
        Job.h
        Master.cpp
        Master.h
        MasterSlavePostTest.cpp
        MonteCarloGlobals.h
        MonteCarloHeaders.h
        montecarlo_imd.h
        Slave.cpp
        Slave.h
        Vector3d.h
        VectorRange.h
)
target_compile_features(${MONTECARLO_EXECUTABLE} PUBLIC cxx_std_11 cxx_nullptr)
target_link_libraries(${MONTECARLO_EXECUTABLE} PRIVATE ${MONTECARLO_LIB})


### ITAP MD library ############################################################

message(STATUS "Setting up auto-build of ITAP MD library")

set(IMD_SOURCE_DIR ${PROJECT_ROOT}/0K_src_a2951c0)

if (CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
    set(IMDSYS native-gcc)
elseif (CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
    set(IMDSYS jureca-icc)
endif()

add_custom_command(
    TARGET ${MONTECARLO_LIB}
    PRE_BUILD
    COMMAND mkdir
    ARGS -p ${PROJECT_ROOT}/lib
    COMMENT "Ensuring existence of ${PROJECT_ROOT}/lib/"
    )
add_custom_command(
    TARGET ${MONTECARLO_LIB}
    PRE_BUILD
    COMMAND make
    ARGS clean IMDSYS=${IMDSYS} imd_eam_glok_virtual_atoms_nblfast_mclib BIN_DIR=${PROJECT_ROOT}/lib
    WORKING_DIRECTORY ${IMD_SOURCE_DIR}
    DEPENDS ${IMD_SOURCE_DIR}/*.c ${IMD_SOURCE_DIR}/*.h
    COMMENT "Building/Updating ITAP MD library"
    )
target_link_libraries(${MONTECARLO_LIB} PUBLIC -L${PROJECT_ROOT}/lib -limd)
target_link_libraries(${MONTECARLO_EXECUTABLE} PUBLIC -L${PROJECT_ROOT}/lib -limd)

message(STATUS "Setting up auto-build of ITAP MD library - done")


### MPI ########################################################################

find_package(MPI REQUIRED)

if (MPI_CXX_FOUND)
    # Sanitize MPI specific flags comming from FindMPI.cmake
    string(STRIP "${MPI_CXX_COMPILE_FLAGS}" MPI_CXX_COMPILE_FLAGS)
    string(STRIP "${MPI_CXX_INCLUDE_PATH}" MPI_CXX_INCLUDE_PATH)
    string(STRIP "${MPI_CXX_LINK_FLAGS}" MPI_CXX_LINK_FLAGS)
    string(STRIP "${MPI_CXX_LIBRARIES}" MPI_CXX_LIBRARIES)

    # Add MPI specific flags to properties
    target_compile_options(${MONTECARLO_LIB} PRIVATE ${MPI_CXX_COMPILE_FLAGS})
    target_include_directories(${MONTECARLO_LIB} PRIVATE ${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(${MONTECARLO_LIB} PUBLIC ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})

    target_compile_options(${MONTECARLO_EXECUTABLE} PRIVATE ${MPI_CXX_COMPILE_FLAGS})
    target_include_directories(${MONTECARLO_EXECUTABLE} PRIVATE ${MPI_CXX_INCLUDE_PATH})
    target_link_libraries(${MONTECARLO_EXECUTABLE} PUBLIC ${MPI_CXX_LINK_FLAGS} ${MPI_CXX_LIBRARIES})
endif(MPI_CXX_FOUND)


### Finalize setup #############################################################

# set installation directories
install(TARGETS ${MONTECARLO_EXECUTABLE} ${MONTECARLO_LIB}
        RUNTIME DESTINATION bin
        LIBRARY DESTINATION lib
        ARCHIVE DESTINATION lib)

