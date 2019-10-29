project(banklib)

set(INCLUDES ${ROOT_INCLUDE_DIR} ${CMAKE_SOURCE_DIR}/dependencies/include )
include_directories(${INCLUDES})

add_library(Banks4 
		BBand.cpp BCalorimeter.cpp BConfig.cpp BEvent.cpp BParticle.cpp BScaler.cpp BScintillator.cpp )
target_link_libraries(Banks4 Hipo4 )
