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

add_executable(BigBoy glob_plot.cc)
target_include_directories(BigBoy PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(BigBoy ${LIB_DIRS})

add_executable(NFGlob no_fit_glob_plot.cc)
target_include_directories(NFGlob PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(NFGlob ${LIB_DIRS})

add_executable(doubleGauss double_gauss_shifted.cc)
target_include_directories(doubleGauss PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(doubleGauss ${LIB_DIRS})

add_executable(HBigBoy horiz_glob_plot.cc)
target_include_directories(HBigBoy PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(HBigBoy ${LIB_DIRS})

add_executable(VBigBoy vert_glob_plot.cc)
target_include_directories(VBigBoy PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(VBigBoy ${LIB_DIRS})

add_executable(TailEdep edeps_tails.cc)
target_include_directories(TailEdep PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(TailEdep ${LIB_DIRS})

add_executable(UCTailEdep uncut_edeps_tails.cc)
target_include_directories(UCTailEdep PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(UCTailEdep ${LIB_DIRS})

#add_executable(TailPID paddle_dist_tails.cc)
#target_include_directories(TailPID PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
#target_link_libraries(TailPID ${LIB_DIRS})

add_executable(IntertimingCutShifted IntertimingCutShifted.cc)
target_include_directories(IntertimingCutShifted PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(IntertimingCutShifted ${LIB_DIRS})

add_executable(SingleLady single_gauss_glob.cc)
target_include_directories(SingleLady PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(SingleLady ${LIB_DIRS})

add_executable(AllPanelGlob panel_all_glob_plot.cc)
target_include_directories(AllPanelGlob PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(AllPanelGlob ${LIB_DIRS})

add_executable(HVStackGlob horiz_vert_stacked_glob_hist.cc)
target_include_directories(HVStackGlob PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(HVStackGlob ${LIB_DIRS})

add_executable(P2P panel_to_panel.cc)
target_include_directories(P2P PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(P2P ${LIB_DIRS})

add_executable(P2PB panel_to_panel_beta.cc)
target_include_directories(P2PB PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(P2PB ${LIB_DIRS})

add_executable(P2P_offsets get_panel_to_panel_offsets.cc)
target_include_directories(P2P_offsets PRIVATE ${GAPS_INCLUDE_DIRS}/gaps)
target_link_libraries(P2P_offsets ${LIB_DIRS})
