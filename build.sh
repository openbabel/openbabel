
set -x

mkdir build
cd build
cmake -DCMAKE_INSTALL_PREFIX=~/usr ../
make -j 8
make install
