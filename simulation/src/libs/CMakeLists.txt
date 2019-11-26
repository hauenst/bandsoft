enable_language(Fortran)
SET(LIBNAME "dis_lib")

FILE(GLOB SRCS *.cpp)
file( GLOB HEADERS ./*.h )

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/simulation/include )
include_directories(${INCLUDES})


add_library(dis_lib SHARED crossdis.cpp crossincl.cpp deuteronwf.cpp alekhin.f F1F209_test.f ${HEADERS} )
TARGET_LINK_LIBRARIES(${LIBNAME} ${ROOT_LIBRARIES} )
INSTALL ( TARGETS ${LIBNAME} DESTINATION "${CMAKE_SOURCE_DIR}/simulation/lib")

# install the grid files so that the cross section code can access them
INSTALL (FILES grids/a09.dsfs_lNNC grids/a09.sfs_lNNC
	DESTINATION "${CMAKE_SOURCE_DIR}/simulation/bin/grids/")
