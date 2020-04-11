cmake_minimum_required(VERSION 3.8)
project(catalyst_adaptor C CXX)

# Default to optimised build
if (NOT CMAKE_BUILD_TYPE OR CMAKE_BUILD_TYPE STREQUAL "")
  set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
endif()

# Adds option -DBUILD_TESTING
include(CTest)

set(USE_CATALYST ON CACHE BOOL "Link with Catalyst")
option(BUILD_SHARED_LIBS "Build shared library to simplify linking" ON)

if(USE_CATALYST)

  set(Adaptor_SRCS adaptor.cxx coprocessor.cxx vtk_pipeline.cxx)
  add_library(catalyst_adaptor ${Adaptor_SRCS})
  set_property(TARGET catalyst_adaptor PROPERTY CXX_STANDARD 11)

  find_package(ParaView REQUIRED)

  if (ParaView_VERSION GREATER_EQUAL 5.7)

    target_link_libraries(catalyst_adaptor
                          INTERFACE
                            VTK::PythonUsed
                          PRIVATE
                            ParaView::Catalyst
                            ParaView::PythonCatalyst
                            VTK::CommonDataModel)
  else()

    include("${PARAVIEW_USE_FILE}")
    target_link_libraries(catalyst_adaptor vtkPVPythonCatalyst)

  endif()

  add_definitions("-DUSE_CATALYST")
  if(NOT PARAVIEW_USE_MPI)
    message(SEND_ERROR "ParaView must be built with MPI enabled")
  endif()

  if(BUILD_TESTING)
    add_subdirectory(test)
  endif()

else()
  find_package(MPI REQUIRED)
  include_directories(${MPI_CXX_INCLUDE_PATH})
endif()

install(TARGETS catalyst_adaptor LIBRARY DESTINATION lib)
