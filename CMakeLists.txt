cmake_minimum_required(VERSION 3.3)
project(TestModule)

find_package(Falaise REQUIRED)
find_package(caffe REQUIRED)

add_library(TestModule SHARED TestModule.h TestModule.cpp)

target_link_libraries(TestModule PUBLIC Falaise::FalaiseModule)
#target_link_libraries(TestModule FalaiseModule)

#install(TARGETS TestModule DESTINATION ${CMAKE_INSTALL_PLUGINDIR})
