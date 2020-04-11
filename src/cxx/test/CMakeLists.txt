# Include Catch CMake integration module
include("${CMAKE_CURRENT_SOURCE_DIR}/Catch2/cmake/Catch.cmake")

set(UnitTest_SRCS catch_main.cxx adaptor_test.cxx coprocessor_test.cxx vtk_pipeline_test.cxx)

add_executable (catalyst_adaptor_unittests ${UnitTest_SRCS})
set_property(TARGET catalyst_adaptor_unittests PROPERTY CXX_STANDARD 11)
target_include_directories(catalyst_adaptor_unittests PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/..)

if (ParaView_VERSION GREATER_EQUAL 5.7)

  # New build system creates module libraries which cannot be linked,
  # need to build additional shared library for testing
  set(lib_SRCS ../adaptor.cxx ../coprocessor.cxx ../vtk_pipeline.cxx)
  add_library(catalystAdaptorShared ${lib_SRCS})
  target_link_libraries(catalystAdaptorShared PRIVATE
    ParaView::Catalyst
    VTK::CommonDataModel)

  target_link_libraries(catalyst_adaptor_unittests PRIVATE
    catalystAdaptorShared
    ParaView::Catalyst
    ParaView::PythonCatalyst
    VTK::CommonDataModel)

else()

  target_link_libraries(catalyst_adaptor_unittests LINK_PRIVATE catalyst_adaptor ${VTK_LIBRARIES})
  include(vtkModuleMacros)
  include(vtkMPI)
  vtk_mpi_link(catalyst_adaptor_unittests)

endif()

catch_discover_tests(catalyst_adaptor_unittests)
