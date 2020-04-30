cmake_minimum_required(VERSION 3.15)
project(pre_processor)

set(CMAKE_CXX_STANDARD 14)

#add_definitions(-DDISABLE_PYTHON_BINDINGS_BUILD)

include_directories(c:/users/sadeep/appdata/local/programs/python/python36/include)
link_directories(c:/users/sadeep/appdata/local/programs/python/python36/libs)

add_executable(pre_processor main.cpp vectorizer.cpp  vectorizer.h )