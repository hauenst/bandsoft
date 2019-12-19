find_package(MySQL REQUIRED)
add_definitions(-DRCDB_MYSQL)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/analysis/include ${MYSQL_INCLUDE_DIR} )

include_directories(${INCLUDES})
add_executable(eepler eepler.cpp )
add_executable(eepnler eepnler.cpp )
target_link_libraries(eepler ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )
target_link_libraries(eepnler ${ROOT_LIBRARIES} Hipo4 Banks4  -lmysqlclient -lpthread -lz -lm -ldl -lssl -lcrypto )

install( TARGETS eepler
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )
install( TARGETS eepnler
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )

add_executable(makesimhists makesimhists.cpp )
target_link_libraries(makesimhists ${ROOT_LIBRARIES} )
install( TARGETS makesimhists
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )

add_executable(eepcompare eepcompare.cpp )
target_link_libraries(eepcompare ${ROOT_LIBRARIES} )
install( TARGETS eepcompare
		DESTINATION ${CMAKE_SOURCE_DIR}/analysis/bin )
