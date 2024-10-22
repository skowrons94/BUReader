cmake_minimum_required(VERSION 3.0)
project(BUReader VERSION 1.0 LANGUAGES CXX)

set(CMAKE_INSTALL_DIR /opt/BUReader)

set(SRC_DIR ${CMAKE_CURRENT_SOURCE_DIR}/src)
set(HEADER_DIR ${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package(ROOT REQUIRED COMPONENTS RIO Net)

include(${ROOT_USE_FILE})

include_directories(
  ${HEADER_DIR}
  ${ROOT_INCLUDE_DIR}
  ${CAENSLC_DIR}
  ${CMAKE_CURRENT_BINARY_DIR}
)

add_executable(BUReader
  ${SRC_DIR}/main.cc
  ${SRC_DIR}/BUReader.cc
)

target_compile_features(BUReader PUBLIC cxx_std_14)

target_link_libraries(BUReader PRIVATE ${ROOT_LIBRARIES})

install(TARGETS BUReader
  DESTINATION ${CMAKE_INSTALL_DIR})
