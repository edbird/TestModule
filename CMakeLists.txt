cmake_minimum_required(VERSION 3.3)
project(TestModule)

# tempoarily enable debugging symbols
set(CMAKE_BUILD_TYPE Debug)

find_package(Falaise REQUIRED)
#find_package(Caffe REQUIRED)

add_library(TestModule SHARED TestModule.h TestModule.cpp)

#include_directories(${Caffe_INCLUDE_DIRS})
#target_link_libraries(TestModule ${Caffe_LIBRARIES})
#set(Caffe_INCLUDE_DIRS "/opt/caffe/src;/usr/include;/opt/caffe/bin/include;/usr/include/opencv;/usr/include/atlas;/opt/caffe/include")
#set(Caffe_LIBRARIES "/opt/caffe/bin/lib/libcaffe.so")
set(Caffe_INCLUDE_DIRS "/usr/include/caffe")
set(Caffe_LIBRARIES "/usr/lib/libcaffe.so")
target_link_libraries(TestModule PUBLIC Falaise::FalaiseModule)
#target_link_libraries(TestModule FalaiseModule)

#install(TARGETS TestModule DESTINATION ${CMAKE_INSTALL_PLUGINDIR})
