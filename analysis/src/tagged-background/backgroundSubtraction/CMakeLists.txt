
cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(basic_subtraction basic_subtraction.cpp )
target_link_libraries(basic_subtraction ${ROOT_LIBRARIES} )

add_executable(testSubtraction testSubtraction.cpp )
target_link_libraries(testSubtraction ${ROOT_LIBRARIES} )

add_executable(create_plots create_plots.cpp )
target_link_libraries(create_plots ${ROOT_LIBRARIES} )
