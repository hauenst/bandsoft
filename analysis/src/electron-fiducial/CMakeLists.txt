find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/dependencies/hipolib ${CMAKE_SOURCE_DIR}/analysis/include ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})



add_executable(fiducials fiducials.cpp )
target_link_libraries(fiducials ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS fiducials
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )


add_executable(draw_fiducials draw_fiducials.cpp )
target_link_libraries(draw_fiducials ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS draw_fiducials
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )
