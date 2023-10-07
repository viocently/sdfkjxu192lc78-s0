# Reference codes for the paper "Massive Superpoly Recovery with a Meet-in-the-middle Framework -- Improved Cube Attacks on Trivium and Kreyvium"

## 1. Structure of the codes

The header file "BooleanPolynomial.h" defines the class for the representation and operations of Boolean polynomials

The header file "flag.h" defines the class for the representation and operations of flags.

The header file "listsOfPolynomials" defines the class for computing the contributions of core monomial trails.

The header file "cipher.h" defines the update functions for different ciphers and gives the concrete implementation of our CMP-based approach.

The header file "framework.h" gives the concrete implementation of the MITM framework. 

The cpp file "SuperpolyBGL.cpp" configures parameters for different ciphers and then calls the MITM framework to recover the superpoly.

## 2. Usage of the codes
We give a brief introduction on how to run our code on a linux platform.

1. Install Gurobi (our version is 9.5.1) and configure the required environment variables such as "GUROBI_HOME" and "LD_LIBRARY_PATH".

2. Open the file "SuperpolyBGL.cpp" to set the *cipher_name*, *cube_index* and *rounds*. 
We provides two options for outputting the superpoly via the variable *solver_mode*. 

   Setting *solver_mode* to *mode::OUTPUT_EXP* means once the MITM framework finishes, the concrete expression of the superpoly will be output to the folder "TERM" as a file named "superpoly.txt".

   Setting *solver_mode* to *mode::OUTPUT_FILE* means the program will save the contribution of each core monomial trail in an unexpanded form to the folder "TERM". This mode can reduce the memory usage during the running of the program. To recover the exact superpoly under this mode, you can set the variable *isAccurate* to *true* in the main function, then after the program teminates there will be a file named "superpoly.txt" in the folder "TERM" that contains the final superpoly; if you set *isAccurate* to *false*, there will also be a file named "superpoly.txt" in the folder "TERM", but this file contains the unexpanded contributions that appears odd-number times.

   Regardless of what *solver_mode* is set to, some information about the superpoly (such as the algebraic degree, the number of monomials appearing in the superpoly, etc.) will end up being output in the standard output.

Tips: due to the special memory usage structure of the linux system, if you encounter an out of the memory (OOM) problem when running the program, you can adjust the environment variable *GLIBC_TUNABLES=glibc.malloc.arena_max* in the file "exec.sh". The smaller the value of this variable, the lower the memory consumption, but with some loss of speed.

3. Create three folders named "STATE", "LOG" and "TERM" and type `sh exec.sh` in the console. This should generate an executable program "mitm". "LOG" contains log files; "STATE" stores the hash table $P$ after each expansion; "TERM" stores the information of the superpoly that has been extracted by the CMP-based approach.

4. Type `./mitm` in the console to start the superpoly recovery. While the program is running, the status of the program will be recorded in the log files.

## 3. Dependencies
Note that the header file "dynamic_bitset.hpp" used in the codes is from the C++ Boost Library, which can be downloaded from (https://www.boost.org/).
