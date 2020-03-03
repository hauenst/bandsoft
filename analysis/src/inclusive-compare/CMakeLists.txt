find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/dependencies/hipolib ${CMAKE_SOURCE_DIR}/analysis/include ${MYSQL_INCLUDE_DIR})

include_directories(${INCLUDES})


add_executable(inclusive-skim inclusive-skim.cpp )
target_link_libraries(inclusive-skim ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS inclusive-skim
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )

add_executable(compare compare.cpp )
target_link_libraries(compare ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS compare
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )
