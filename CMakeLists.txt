# heyoka requires at least CMake 3.20.
cmake_minimum_required(VERSION 3.20.0)


set(CMAKE_CXX_COMPILER /opt/gcc-15/bin/g++)

# The name of our project.
project(sample_project)

include_directories(include)

# Look for an installation of heyoka in the system.
find_package(heyoka REQUIRED)

# Create an executable, and link it to the heyoka::heyoka imported target.
# This ensures that, in the compilation of 'main', heyoka's include
# dirs are added to the include path of the compiler and that heyoka's
# dependencies are transitively linked to 'main'.

add_executable(CR3BP_heyoka src/CR3BP_heyoka.cpp)
target_link_libraries(CR3BP_heyoka 
    heyoka::heyoka
    /home/jojo/local/lib/libssl.so.3
    /home/jojo/local/lib/libcrypto.so.3
    /home/jojo/.local/lib/libboost_serialization.so.1.85.0
    /home/jojo/.local/lib/libboost_charconv.so.1.85.0
    /home/jojo/.local/lib/libboost_regex.so.1.85.0
    /home/jojo/local/lib/libfmt.so.12
    /home/jojo/local/lib/libspdlog.so.1.16
    /home/jojo/local/lib/libtbb.so.12
)

add_executable(test-heyoka-cr3bp-impl src/test-heyoka-cr3bp-impl.cpp)
target_link_libraries(test-heyoka-cr3bp-impl 
    heyoka::heyoka
    /home/jojo/local/lib/libssl.so.3
    /home/jojo/local/lib/libcrypto.so.3
    /home/jojo/.local/lib/libboost_serialization.so.1.85.0
    /home/jojo/.local/lib/libboost_charconv.so.1.85.0
    /home/jojo/.local/lib/libboost_regex.so.1.85.0
    /home/jojo/local/lib/libfmt.so.12
    /home/jojo/local/lib/libspdlog.so.1.16
    /home/jojo/local/lib/libtbb.so.12
)