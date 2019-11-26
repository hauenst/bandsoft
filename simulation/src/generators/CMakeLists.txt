
enable_language(Fortran)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/simulation/include )

include_directories(${INCLUDES})
add_executable(dis dis.cpp )
target_link_libraries(dis ${ROOT_LIBRARIES} Foam dis_lib dis_tree )
install( TARGETS dis
		DESTINATION ${CMAKE_SOURCE_DIR}/simulation/bin )
