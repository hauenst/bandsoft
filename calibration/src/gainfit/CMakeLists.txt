set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/calibration/include )

include_directories(${INCLUDES})

# Spectral fitting for each individual file for all PMTs:
add_executable(specFit specFit.cpp )
target_link_libraries(specFit ${ROOT_LIBRARIES} Hipo4 Banks4 )
install( TARGETS specFit
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/gainfit )

# Combine spectral fitting results for the HV setting for that run:
add_executable(combFiles combFiles.cpp )
target_link_libraries(combFiles ${ROOT_LIBRARIES} Hipo4 Banks4 )
install( TARGETS combFiles
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/gainfit )

# Do the gain fit on 5 files 
add_executable(gainFit gainFit.cpp )
target_link_libraries(gainFit ${ROOT_LIBRARIES} Hipo4 Banks4 )
install( TARGETS gainFit
		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/gainfit )

# Convert the result to 
#add_executable(compareSLC compareSLC.cpp )
#target_link_libraries(compareSLC ${ROOT_LIBRARIES} Hipo4 Banks4 )
#install( TARGETS compareSLC
#		DESTINATION ${CMAKE_SOURCE_DIR}/calibration/bin/gainfit )
