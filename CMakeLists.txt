cmake_minimum_required(VERSION 3.15)
project(pre_processor)

set(CMAKE_CXX_STANDARD 14)

add_executable(pre_processor main.cpp vectorizer.cpp vectorizer.h)