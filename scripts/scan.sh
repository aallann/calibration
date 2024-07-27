#!/bin/bash

dependencies=(
    "build-essential"
    "git" 
    "cmake" 
    "libeigen3-dev" 
    "libboost-all-dev" 
    "valgrind" 
    "liblapack-dev" 
    "libopenblas-dev" 
    "liblua5.3-dev" 
    "python3" 
    "python3-pip" 
)

scanPackage() {
    package=$1

    if dpkg -s "$package" &> /dev/null; then
        echo "$package installed"

    else
        echo "$package is not installed."
    fi
}

# loop 
for dep in "${dependencies[@]}"; do
    scanPackage "$dep"
    echo "---------------------------------"
done

echo "Dependency check completed."