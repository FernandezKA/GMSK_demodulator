cmake_minimum_required(VERSION 3.25)

project (GMSK_demodulator)

set(CMAKE_BUILD_TYPE Debug)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wextra")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Weffc++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Woverloaded-virtual")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wctor-dtor-privacy")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wnon-virtual-dtor")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wold-style-cast")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Werror")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wconversion -Wsign-conversion")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Winit-self")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wunreachable-code")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -pedantic")


include_directories(Includes)
include_directories(Sources)

set(HEADERS
    ./Includes/GMSK_demodulator.hpp
    ./Includes/tests.hpp
)

set(EXECUTABLE 
    ./Sources/main.cpp
    ./Sources/GMSK_demodulator.cpp
    ./Sources/tests.cpp    
)

add_executable(
    main ${EXECUTABLE} ${HEADERS}
)
