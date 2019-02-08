mkdir -p build #make dir if first build
cd build

cmake ../ -DCMAKE_INSTALL_PREFIX="/usr/local/" -DCMAKE_TOOLCHAIN_FILE=./scripts/toolchains/HiDEM-ubuntu.cmake -DCMAKE_BUILD_TYPE=Release

#If you want a debug build instead:
#cmake ../ -DCMAKE_INSTALL_PREFIX="/usr/local/" -DCMAKE_TOOLCHAIN_FILE=./scripts/toolchains/HiDEM-ubuntu.cmake -DCMAKE_BUILD_TYPE=Debug
make
sudo make install
