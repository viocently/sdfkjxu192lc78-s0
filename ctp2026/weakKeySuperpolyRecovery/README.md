## Usage of the codes
We give a brief introduction on how to run our code on a linux platform.

1. Install Gurobi (our version is 9.1.2) and configure the required environment variables such as "GUROBI_HOME" and "LD_LIBRARY_PATH".

2. Open the file "SuperpolyBGL.cpp" to set the *cube_index*, *rounds* and *bit conditions*.
   For example, if you would like regenerate the weak-key superpoly of $I_1$ listed in the paper, you can set the *cube_index* to $I_1$, *rounds* to $852$ and set *bit conditions* as
   "bit_conditions[58] = BooleanPolynomial(288, "0");bit_conditions[59] = BooleanPolynomial(288, "s32");bit_conditions[42] = BooleanPolynomial(288, "s40s41+s15")", which correspond to the bit conditions
   $k[58] = 0, k[59]+k[57]k[58]+k[32] = 0$ and $k[42]+k[40]k[41]+k[15] = 0$. 

4. Create three folders named "STATE", "LOG" and "TERM" in the current directory and compile the source files with multi-threading support. This should generate an executable program, let us say it is "mitm".

5. Type `./mitm` in the console to start the superpoly recovery. After the program completes, you shall see a file "superpoly.txt" in the folder "TERM", which contains the weak-key superpoly.


## Dependencies
Note that the header file "dynamic_bitset.hpp" used in the codes is from the C++ Boost Library, which can be downloaded from (https://www.boost.org/).

This is an example of the compilation command `g++  SuperpolyBGL.cpp deg.cpp -o mitm -std=c++17 -O2 -lm -lpthread -I/$GUROBI_HOME/include/   -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm`
