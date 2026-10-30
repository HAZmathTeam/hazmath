# Checks that an installed HAZmath can be consumed with find_package(hazmath).
# HAZmath must be configured with -DHAZMATH_RELOCATABLE_INSTALL=ON -DCMAKE_INSTALL_PREFIX=<prefix>.
# Usage: cmake -S tests/cmake_consumer -B build-consumer -DCMAKE_PREFIX_PATH=<hazmath prefix>
cmake_minimum_required(VERSION 3.16)
project(hazmath_consumer C)

find_package(hazmath 0.1 CONFIG REQUIRED)

# hazmath::hazmath must carry all its link dependencies (e.g. libm)
add_executable(amg_laplace main.c)
target_link_libraries(amg_laplace PRIVATE hazmath::hazmath)

enable_testing()
add_test(NAME amg_laplace COMMAND amg_laplace)
