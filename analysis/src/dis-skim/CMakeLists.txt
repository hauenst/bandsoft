find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/analysis/include ${MYSQL_INCLUDE_DIR} )

include_directories(${INCLUDES})

add_executable(hipo2root hipo2root.cpp )
target_link_libraries(hipo2root ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS hipo2root
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )

