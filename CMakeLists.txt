# heyoka requires at least CMake 3.20.
cmake_minimum_required(VERSION 3.20.0)

set(CMAKE_CXX_COMPILER /opt/gcc-15/bin/g++)

project(sample_project)

include_directories(include /usr/include/eigen3)

find_package(heyoka REQUIRED)
add_executable(CR3BP_heyoka src/CR3BP_heyoka.cpp)
#[[
add_executable(CR3BP_rk src/CR3BP_rk.cpp)
add_executable(CR3BP_vgl src/CR3BP_vgl.cpp)
add_executable(test-heyoka-cr3bp-impl src/test-heyoka-cr3bp-impl.cpp)
add_executable(fehlberg67-test src/boost_odeint.cpp)
]]
add_executable(benchmark-bulirsch_stoer src/benchmark/bulirsch-stoer-benchmark.cpp)
add_executable(benchmark-fehlberg78 src/benchmark/fehlberg78-benchmark.cpp)
add_executable(benchmark-DP5 src/benchmark/DP5-benchmark.cpp)
add_executable(benchmark-heyoka src/benchmark/heyoka-benchmark.cpp)

set(heyoka_link_TARGETS
CR3BP_heyoka
#[[
    CR3BP_rk
    CR3BP_vgl
    test-heyoka-cr3bp-impl
    fehlberg67-test
    ]]
    benchmark-bulirsch_stoer
    benchmark-fehlberg78
    benchmark-DP5
    benchmark-heyoka
)

foreach(tgt IN LISTS heyoka_link_TARGETS)
    target_link_libraries(${tgt} PRIVATE
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
endforeach()