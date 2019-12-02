set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/calibration/include )

include_directories(${INCLUDES})

# Spectral fitting for each individual file for all PMTs:
add_executable(timewalk timewalk.cpp )
target_link_libraries(timewalk ${ROOT_LIBRARIES} Hipo4 Banks4 )

add_executable(test_parameters test_parameters.cpp )
target_link_libraries(test_parameters ${ROOT_LIBRARIES} Hipo4 Banks4 )

install( TARGETS timewalk
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/timewalkfit )
install( TARGETS test_parameters
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/timewalkfit )
