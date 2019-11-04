set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_CURRENT_SOURCE_DIR}/../include )

include_directories(${INCLUDES})
add_executable(singlebar singlebar.cpp )
target_link_libraries(singlebar ${ROOT_LIBRARIES} Hipo4 Banks4 )
install( TARGETS singlebar
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/sourcefit )
