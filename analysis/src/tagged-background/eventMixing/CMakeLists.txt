
cmake_minimum_required (VERSION 2.6)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/include )
include_directories(${INCLUDES})

add_executable(inclusive_skim inclusive_skim.cpp )
target_link_libraries(inclusive_skim ${ROOT_LIBRARIES} )

add_executable(neutron_skim neutron_skim.cpp )
target_link_libraries(neutron_skim ${ROOT_LIBRARIES} )


add_executable(event_mixing event_mixing.cpp )
target_link_libraries(event_mixing ${ROOT_LIBRARIES} )
