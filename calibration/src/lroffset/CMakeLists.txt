set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/dependencies/hipolib ${CMAKE_SOURCE_DIR}/calibration/include )

include_directories(${INCLUDES})

# Spectral fitting for each individual file for all PMTs:
add_executable(offsets offsets.cpp ../calibclass.cpp ../helpers.cpp )
target_link_libraries(offsets ${ROOT_LIBRARIES} Hipo4 Banks4 )

add_executable(test_offsets test_offsets.cpp ../calibclass.cpp ../helpers.cpp )
target_link_libraries(test_offsets ${ROOT_LIBRARIES} Hipo4 Banks4 )

install( TARGETS offsets
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/lroffset )
install( TARGETS test_offsets
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/lroffset )
