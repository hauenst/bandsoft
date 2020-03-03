set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/dependencies/hipolib ${CMAKE_SOURCE_DIR}/calibration/include )

include_directories(${INCLUDES})

# Spectral fitting for each individual file for all PMTs:
add_executable(barbybar barbybar.cpp ../calibclass.cpp ../helpers.cpp )
target_link_libraries(barbybar ${ROOT_LIBRARIES} Hipo4 Banks4 )

#add_executable(layerbylayer layerbylayer.cpp )
#target_link_libraries(layerbylayer ${ROOT_LIBRARIES} Hipo4 Banks4 calibclass )

install( TARGETS barbybar
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/baroffset )

#install( TARGETS layerbylayer
#		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/baroffset )
