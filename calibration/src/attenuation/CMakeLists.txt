set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/calibration/include )

include_directories(${INCLUDES})

# Spectral fitting for each individual file for all PMTs:
add_executable(attenuation attenuation.cpp )
target_link_libraries(attenuation ${ROOT_LIBRARIES} Hipo4 Banks4 )

install( TARGETS attenuation
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/attenuation )
