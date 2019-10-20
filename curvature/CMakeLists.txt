get_filename_component(EXERCISENAME ${CMAKE_CURRENT_LIST_DIR} NAME)

file(GLOB_RECURSE HEADERS "*.h")
file(GLOB_RECURSE SHADERS "*.glsl")

#================================================================
# Curvature
#================================================================

add_executable(${EXERCISENAME}
    curvature.cpp
    viewer.cpp
    ${HEADERS}
    ${SHADERS})
target_link_libraries(${EXERCISENAME} surface_mesh)
# Lastly, additional libraries may have been built for you.  In addition to linking
# against NanoGUI, we need to link against those as well.
target_link_libraries(${EXERCISENAME} nanogui ${NANOGUI_EXTRA_LIBS})

#================================================================
# Unit Test
#================================================================

add_executable(curvature_test
    curvature_test.cpp
    viewer.cpp
    ${HEADERS}
    ${SHADERS})
target_link_libraries(curvature_test surface_mesh)
# Lastly, additional libraries may have been built for you.  In addition to linking
# against NanoGUI, we need to link against those as well.
target_link_libraries(curvature_test nanogui ${NANOGUI_EXTRA_LIBS})