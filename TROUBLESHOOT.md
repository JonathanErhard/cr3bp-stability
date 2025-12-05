
# heyoka

building heyoka can be a bit of a journey, since the library relies on modern features of the library. If you installed outdated libraries (eg boost < 1_79) the library or the project including the library will not build.
For me, I am running Ubuntu 20.04 to support the correct version of ros for a different project and both the libraries and compiler available by default were not sufficient.

The process of building and installing a library from source is very streamlined.

## TBB

old versions of TBB have a slightly different file structure, where header files are not in include/oneapi/tbb but rather in include/tbb. In case the compiler complains it cannot find tbb header files, check your tbb include directories structur. To install a new version, run:

```bash
git clone https://github.com/oneapi-src/oneTBB.git
cd oneTBB
git checkout v2021.9.0
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=/usr/local ..
cmake --build . -- -j$(nproc)
sudo cmake --install .
```

## boost

You need boost > 1_79. to compile from source, use:

```bash
wget https://boostorg.jfrog.io/artifactory/main/release/1.84.0/source/boost_1_84_0.tar.gz
tar xf boost_1_84_0.tar.gz
cd boost_1_84_0
./bootstrap.sh --prefix=/usr/local
sudo ./b2 install
```

## g++

heyoka requires modern c++23 features. g++-13 can compile c++23, but the std library is missing some modern features. To compile heyoka, g++-14 or newer is necessary.
If you don't have an ancient os, you can just install it with your favorite package manager. In my case my os was designed to be used on stone tablets and punsh cards, so I had to build it myself using

```bash
wget https://ftp.gnu.org/gnu/gcc/gcc-15.1.0/gcc-15.1.0.tar.xz
tar xf gcc-15.1.0.tar.xz
cd gcc-15.1.0
mkdir build
cd build
../configure --prefix=/usr/local/gcc-15 --enable-languages=c,c++ --disable-multilib --program-suffix=15
make -j$(nproc)
sudo make install
```

to set g++15 as the compiler, use

```bash
export CC=/usr/local/gcc-15/bin/gcc15
export CXX=/usr/local/gcc-15/bin/g++15
```

to set it for your terminal, or add

```cmake
set(CMAKE_CXX_COMPILER /opt/gcc-15/bin/g++)
```

## select correct verion

in case you installed multiple versions of some library and cmake doest chose the correct one, set include dirs in the cmake with

```cmake
set(<lib_name>_ROOT_DIRS "/path/to/lib")
```

## compiling your own project

in case there are any linker issues, you can see which .so libraries your executable depends on with

```bash
ldd ./<executable-name>
```

if there are and lines

```bash
<lib-name>.so.X.Y.Z => not found
```

locate them using

```bash
locate <lib-name>.so.X.Y.Z 
```

and link them in the cmake

```cmake
target_link_libraries(test_heyoka
    PRIVATE
        heyoka::heyoka
        path/to/shared-object-lib.so
)
```

notice errors with shared object libraries can also occur at runtime, not just as a linker error.
