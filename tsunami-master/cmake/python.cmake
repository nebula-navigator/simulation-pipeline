find_package(SWIG 4.0 REQUIRED COMPONENTS python)
include(UseSWIG)

list(APPEND CMAKE_SWIG_FLAGS "-doxygen")

if(UNIX AND NOT APPLE)
    list(APPEND CMAKE_SWIG_FLAGS "-DSWIGWORDSIZE64")
endif()

# Find Python 3
find_package(Python3 REQUIRED COMPONENTS Interpreter Development NumPy)
list(APPEND CMAKE_SWIG_FLAGS "-globals" "-fastproxy" "-builtin" "-c++")

set_property(SOURCE swig/tsunami4py.i PROPERTY CPLUSPLUS ON)
set_property(SOURCE swig/tsunami4py.i PROPERTY SWIG_MODULE_NAME tsunami)
swig_add_library(tsunami
        TYPE MODULE
        LANGUAGE python
        OUTPUT_DIR "${CMAKE_BINARY_DIR}/python"
        OUTFILE_DIR "${CMAKE_CURRENT_BINARY_DIR}${CMAKE_FILES_DIRECTORY}"
        SOURCES swig/tsunami4py.i)

add_library(tsunamilibshared SHARED src/IO.cpp src/classification.cpp src/keplerutils.cpp src/stopcond.cpp src/kozailidov.cpp)
set_target_properties(tsunamilibshared PROPERTIES OUTPUT_NAME tsunamilib)

if(stopcond)
    message("-- Compiling the stopping conditions, see stopcond_test.py for example")
    target_compile_options(tsunamilibshared PUBLIC -DSTOPCOND)
endif(stopcond)

target_include_directories(tsunami PRIVATE lib ${Python3_INCLUDE_DIRS})
set_property(TARGET tsunami PROPERTY SWIG_USE_TARGET_INCLUDE_DIRECTORIES ON)
target_link_libraries(tsunami ${PYTHON_LIBRARIES} Python3::NumPy)

if(APPLE)
    set_target_properties(tsunami PROPERTIES
            SUFFIX ".so"
            #INSTALL_RPATH "@loader_path;@loader_path/../../${PYTHON_PROJECT}/.libs"
            )
    set_property(TARGET tsunami APPEND PROPERTY
            LINK_FLAGS "-flat_namespace -undefined suppress"
            )
elseif(UNIX)
#    set_target_properties(tsunami PROPERTIES
#            INSTALL_RPATH "$ORIGIN:$ORIGIN/../../${PYTHON_PROJECT}/.libs"
#            )
endif()

set_target_properties(tsunami PROPERTIES LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/python")

target_link_libraries(tsunami tsunamilibshared)

add_custom_target(copy-python-dir ALL COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_SOURCE_DIR}/src/python/ ${CMAKE_BINARY_DIR}/python/)

