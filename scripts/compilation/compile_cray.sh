#module swap PrgEnv-gnu PrgEnv-cray
mkdir -p build #make dir if first build
cd build

cmake ../ -DCMAKE_INSTALL_PREFIX="." -DCMAKE_TOOLCHAIN_FILE=./scripts/toolchains/HiDEM-cray.cmake
make
make install
