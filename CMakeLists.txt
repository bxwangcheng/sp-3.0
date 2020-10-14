project(libsp++3)
cmake_minimum_required(VERSION 3.17)
FILE(GLOB SRC *.cpp *.h)
set(SOURCE
        ${SRC}
        )
ADD_LIBRARY( sp++3 SHARED ${SOURCE} )