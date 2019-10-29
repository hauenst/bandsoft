project(hipolib)

set(INCLUDES ${LZ4_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include )
include_directories(${INCLUDES})

add_library(Hipo4 
		bank.cpp dictionary.cpp event.cpp reader.cpp record.cpp recordbuilder.cpp utils.cpp wrapper.cpp writer.cpp )

target_link_libraries(Hipo4 ${LZ4_LIBRARY} )
