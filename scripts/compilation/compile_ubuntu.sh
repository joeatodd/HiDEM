mkdir -p build #make dir if first build
cd build

cmake ../ -DCMAKE_INSTALL_PREFIX="/usr/local/" -DCMAKE_TOOLCHAIN_FILE=./scripts/toolchains/HiDEM-ubuntu.cmake
make
sudo make install
