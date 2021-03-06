cmake_minimum_required(VERSION 3.0)

project(libtfhe-tests)

include_directories(${GTEST_INCLUDE_DIRS})

set(GOOGLETEST_SOURCES
    arithmetic_test.cpp
    lwe_test.cpp
    polynomial_test.cpp
    tlwe_test.cpp
    lwekeyswitch_test.cpp
    bootstrapping_test.cpp
    )

foreach (FFT_PROCESSOR IN LISTS FFT_PROCESSORS) 

if (FFT_PROCESSOR STREQUAL "fftw")
set(RUNTIME_LIBS
    tfhe-fftw
    ${FFTW_LIBRARIES}
    )

else ()
set(RUNTIME_LIBS
    tfhe-${FFT_PROCESSOR}
    )

endif (FFT_PROCESSOR STREQUAL "fftw") 


add_executable(test-bootstrapping-${FFT_PROCESSOR} test-bootstrapping.cpp ${TFHE_HEADERS})
add_executable(test-bootstrapping-fft2-${FFT_PROCESSOR} test-bootstrapping-fft2.cpp ${TFHE_HEADERS})

add_executable(test-lwe-${FFT_PROCESSOR} test-lwe.cpp ${TFHE_HEADERS})

add_executable(test-tlwe-${FFT_PROCESSOR} test-tlwe.cpp ${TFHE_HEADERS})


target_link_libraries(test-bootstrapping-${FFT_PROCESSOR} ${RUNTIME_LIBS})
target_link_libraries(test-bootstrapping-fft2-${FFT_PROCESSOR} ${RUNTIME_LIBS})

target_link_libraries(test-lwe-${FFT_PROCESSOR} ${RUNTIME_LIBS})

target_link_libraries(test-tlwe-${FFT_PROCESSOR} ${RUNTIME_LIBS})


add_test(unittests-${FFT_PROCESSOR} unittests-${FFT_PROCESSOR})

endforeach (FFT_PROCESSOR IN LISTS FFT_PROCESSORS) 

