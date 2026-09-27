## Usage of the codes
We give a brief introduction on how to run our code on a linux platform.

1. Install Gurobi (our version is 9.1.2) and configure the required environment variables such as "GUROBI_HOME" and "LD_LIBRARY_PATH".

2. Open the file "SuperpolyBGL.cpp" to set the *cube_index* and *rounds*, and set the maximum number of optimal conditions you would like to search for as the input parameter to the function `MITM_framework.start()`.


3. Compile the source files with multi-threading support. This should generate an executable program. Let us say its name is "mitm". 

4. Type `./mitm` in the console to start the search for optimal conditions. After the search completes, you should see the optimal conditions output to the standard output.

Note that the header file "dynamic_bitset.hpp" used in the codes is from the C++ Boost Library, which can be downloaded from (https://www.boost.org/).

This is an example of the compilation command `g++  SuperpolyBGL.cpp deg.cpp -o mitm -std=c++17 -O2 -lm -lpthread -I/$GUROBI_HOME/include/   -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm`
