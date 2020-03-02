find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/analysis/include ${MYSQL_INCLUDE_DIR} )

include_directories(${INCLUDES})


add_executable(tagged-skim tagged-skim.cpp )
target_link_libraries(tagged-skim ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS tagged-skim
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )

add_executable(compare-tagged compare-tagged.cpp )
target_link_libraries(compare-tagged ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
install( TARGETS compare-tagged
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )
