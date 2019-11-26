# find source and header files, set paths, includes...
SET(LIBNAME "dis_tree")
FILE(GLOB SRCS *.cpp)


FILE(GLOB HEADERS ${CMAKE_SOURCE_DIR}/simulation/include/gen_tree.h ${CMAKE_SOURCE_DIR}/simulation/include/band_tree.h ${CMAKE_SOURCE_DIR}/simulation/include/recon_tree.h ${CMAKE_SOURCE_DIR}/simulation/include/digit_tree.h ${CMAKE_SOURCE_DIR}/simulation/include/prop_tree.h)
set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include ${CMAKE_SOURCE_DIR}/simulation/include )
INCLUDE_DIRECTORIES(${INCLUDES})

# create root dictionary
SET(LINKDEF ${CMAKE_CURRENT_SOURCE_DIR}/LinkDef.h)
SET(DICTIONARY ${LIBNAME}Dict.cpp)
ROOT_GENERATE_DICTIONARY("${HEADERS}" "${LINKDEF}" "${DICTIONARY}" "${INCLUDES}")
LIST(APPEND SRCS ${DICTIONARY})

# build the library
ADD_LIBRARY(${LIBNAME} SHARED ${SRCS} ${HEADERS})
TARGET_LINK_LIBRARIES(${LIBNAME} ${ROOT_LIBRARIES} )
INSTALL ( TARGETS ${LIBNAME} DESTINATION "${CMAKE_SOURCE_DIR}/simulation/lib")
