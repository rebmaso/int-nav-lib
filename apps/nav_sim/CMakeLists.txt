cmake_minimum_required(VERSION 3.7)

project(nav_sim)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake")

find_package(Eigen3 REQUIRED)
find_package(OpenCV 4 REQUIRED)
find_package(glog REQUIRED)
find_package(OsgEarth REQUIRED)
set(OpenGL_GL_PREFERENCE LEGACY)
find_package(OpenGL REQUIRED)
find_package(EGL REQUIRED)
find_package(OSG REQUIRED COMPONENTS osg osgDB osgGA osgUtil osgViewer)
find_package(Fontconfig REQUIRED)

find_package(intnavlib REQUIRED)

add_executable(nav_sim 
                    ${CMAKE_CURRENT_SOURCE_DIR}/src/nav_sim.cpp
                    )
 
target_link_libraries(nav_sim  
                        PRIVATE
                            glog::glog
                            ${OPENTHREADS_LIBRARY}
                            Eigen3::Eigen
                            intnavlib::intnavlib
                            ${OpenCV_LIBS}
                            ${EGL_LIBRARY} 
                            ${OSG_LIBRARY} 
                            ${OSGVIEWER_LIBRARY}
                            ${OSGDB_LIBRARY} 
                            ${OSGGA_LIBRARY} 
                            ${OSGUTIL_LIBRARY}
                            ${OSGEARTH_LIBRARY}
                        )

# manually include and link JSBSim
target_include_directories(nav_sim PRIVATE /usr/include/JSBSim)
target_link_libraries(nav_sim PRIVATE /usr/lib/libJSBSim.a)