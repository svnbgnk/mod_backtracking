cmake_minimum_required (VERSION 3.2) # require at least version 3.2 of cmake
project (seqan3_infrastructure_tutorial CXX) # this is a c++-only project
set(SeqAn3_DIR "${CMAKE_SOURCE_DIR}/../seqan3/build_system")
find_package (SeqAn3 REQUIRED)

#set (CMAKE_CXX_FLAGS "-g -fopenmp -pthread")
#set (CMAKE_CXX_FLAGS "-pthread")

add_executable(test test.cpp)
add_executable(benchmark benchmark.cpp)  # build an executable from the hello.cpp source
#add_executable(factorial_test factorial_test.cpp)
target_link_libraries (test seqan3::seqan3)
target_link_libraries (benchmark seqan3::seqan3)

#enable_testing()
#add_executable(factorial_test factorial_test.cpp)
#add_test(factorial factorial_test)

#add_subdirectory (googletest)
#enable_testing()
#add_executable(factorial_test factorial_test.cpp)
#target_link_libraries(factorial_test gtest_main)
#add_test(factorial factorial_test)
