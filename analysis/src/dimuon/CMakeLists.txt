set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/dependencies/hipolib ${CMAKE_SOURCE_DIR}/analysis/include )

include_directories(${INCLUDES})
add_executable(dimuon dimuon.cpp)
target_link_libraries(dimuon ${ROOT_LIBRARIES} Hipo4 Banks4  -lpthread -lz -lm -ldl  )
install( TARGETS dimuon
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )

