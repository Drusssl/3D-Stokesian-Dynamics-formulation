cmake_minimum_required(VERSION 3.16)
project(surge_sd3d LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Eigen3 REQUIRED NO_MODULE)

add_library(surge_sd3d
  src/coefficients.cpp
  src/ewald.cpp
  src/near_field.cpp
  src/projection.cpp
  src/solver.cpp
  src/time_integrator.cpp
)

target_include_directories(surge_sd3d
  PUBLIC
    ${CMAKE_CURRENT_SOURCE_DIR}/include
    ${CMAKE_CURRENT_SOURCE_DIR}
    ${CMAKE_CURRENT_SOURCE_DIR}/..
)

target_link_libraries(surge_sd3d PUBLIC Eigen3::Eigen)
target_compile_options(surge_sd3d PRIVATE
  $<$<CXX_COMPILER_ID:Clang,AppleClang,GNU>:-Wall -Wextra -pedantic>
)

add_executable(two_sphere_trajectory examples/two_sphere_trajectory.cpp)
target_link_libraries(two_sphere_trajectory PRIVATE surge_sd3d)

add_executable(smoke_checks examples/smoke_checks.cpp)
target_link_libraries(smoke_checks PRIVATE surge_sd3d)
