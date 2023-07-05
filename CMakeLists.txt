cmake_minimum_required(VERSION 3.25)

project (GMSK_demodulator)

set(CMAKE_BUILD_TYPE Debug)

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