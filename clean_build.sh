#!/bin/bash

# run uninstall target
if [ -d "build" ]; then
	cd build
	sudo make uninstall
	cd ..
	sudo rm -rf build
fi

# manually remove rogue shared objects (not
# doing this step will not affect the build,
# but there may be runtime Python import issues)
#sudo rm -f /usr/local/lib/libct.so*

# build and install
mkdir build
cd build
cmake -DBUILD_DRIVERS_LIMEPCIE=ON ..
make -j$(nproc)
sudo make install
sudo ldconfig
cd ..
