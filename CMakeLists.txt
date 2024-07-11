cmake_minimum_required(VERSION 3.0 FATAL_ERROR)
project(TestDrive)


add_executable(TestDrive src/main.cpp)

#add_subdirectory(TestDrive src/main.cpp)

#target_link_libraries(TestDrive fftw)



set(FFTW_PATH "D:/pract/TestDrive/fftw-lib/fftw")
include_directories(${FFTW_PATH})
set(FFTW_LIB "${FFTW_PATH}/libfftw3-3.lib")

# добавить FFTW_LIB в проект
target_link_libraries(TestDrive ${FFTW_LIB})

