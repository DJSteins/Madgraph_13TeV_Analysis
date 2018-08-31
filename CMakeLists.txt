cmake_minimum_required(VERSION 3.9)
project(madgraph_gambit)

set(CMAKE_CXX_STANDARD 14)

include_directories(
        "/home/dsteiner/usra/gambit/contrib/RestFrames"
        "/home/dsteiner/usra/gambit/Utils/include"
        "/home/dsteiner/usra/gambit/Logs/include"
        "/home/dsteiner/usra/gambit/contrib/heputils/include"
        "/home/dsteiner/usra/gambit/contrib/mcutils/include"
        "/home/dsteiner/usra/gambit/ColliderBit/include"
        "/usr/lib/gcc/x86_64-linux-gnu/7.3.0/include/"
        "/home/dsteiner/usra/root/include/")


add_executable(
        madgraph_gambit
        /home/dsteiner/usra/gambit/ColliderBit/src/analyses/mt2_bisect.cpp
        /home/dsteiner/usra/gambit/ColliderBit/src/Utils.cpp
        tev13_analysis.C)

target_compile_definitions(madgraph_gambit PUBLIC DANIELS_COMP)
target_link_libraries(madgraph_gambit /usr/lib/gcc/x86_64-redhat-linux/4.4.7/libgomp.so)