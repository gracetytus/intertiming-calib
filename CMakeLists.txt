#----------------------------------------------------------------------------
# Setup the project
cmake_minimum_required(VERSION 3.13 FATAL_ERROR)
project(intertiming_calib)

#----------------------------------------------------------------------------

find_package(ROOT REQUIRED)
find_package(GAPS REQUIRED)
find_package(Boost 1.74 REQUIRED program_options timer log log_setup)

add_definitions(-DUSE_BOOST_PROGRAM_OPTIONS)
add_definitions(-DUSE_BOOST_TIMER)
add_definitions(-DUSE_BOOST_LOG)
add_definitions(-DBOOST_LOG_DYN_LINK)

include(${ROOT_USE_FILE})
set(CMAKE_EXPORT_COMPILE_COMMANDS ON) # only needed for working with my ide

#----------------------------------------------------------------------------
# Add the executable, and link it to the ROOT/GAPS Libs
#

set(LIB_DIRS ${ROOT_LIBRARIES}
                            ${Geant4_LIBRARIES}
                            ${GAPS_LIBRARIES}/libGAPSCommon.so
                            ${GAPS_LIBRARIES}/libGAPSAnalysis.so
                            ${GAPS_LIBRARIES}/libGAPSReconstruction.so
                            Boost::program_options
                            Boost::timer
                            Boost::log
                            Boost::log_setup)

add_executable(IntertimingCalib IntertimingCalib.cc)
target_include_directories(IntertimingCalib PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(IntertimingCalib ${LIB_DIRS})

add_executable(IntertimingCalibCombo IntertimingCalibCombo.cc)
target_include_directories(IntertimingCalibCombo PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(IntertimingCalibCombo ${LIB_DIRS})

add_executable(LeastSquaresIntertiming lsq_intertiming.cc)
target_include_directories(LeastSquaresIntertiming PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(LeastSquaresIntertiming ${LIB_DIRS})

add_executable(RBIntertiming rb_intertiming.cc)
target_include_directories(RBIntertiming PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(RBIntertiming ${LIB_DIRS})

add_executable(IntertimingShifted IntertimingCalib_stable_shifted.cc)
target_include_directories(IntertimingShifted PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(IntertimingShifted ${LIB_DIRS})

#add_executable(MakeHist make_hist.cc)
#target_include_directories(MakeHist PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
#target_link_libraries(MakeHist ${LIB_DIRS})
