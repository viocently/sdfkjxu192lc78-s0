# Reference codes for the paper "Revisiting Core Monomial Prediction: Graph-Based MITM Framework for Superpoly Recovery"

## 1. Structure of the codes

The header file "BooleanPolynomial.h" defines the class for the representation and operations of Boolean polynomials

The header file "flag.h" defines the class for the representation and operations of flags

The header file "cipher.h" defines the update functions for different ciphers and gives the concrete implementation of our CMP-based approach

The header file "framework.h" gives the concrete implementation of the MITM framework. 

The cpp file "SuperpolyBGL.cpp" configures parameters for different ciphers and then calls the MITM framework to recover the superpoly.

## 2. Usage of the codes
We give a brief introduction on how to run our code on a linux platform.

1. Install Gurobi (our version is 9.5.1) and configure the required environment variables such as "GUROBI_HOME" and "LD_LIBRARY_PATH".

2. Open the file "SuperpolyBGL.cpp" to set the *cipher_name*, *cube_index* and *rounds*. 
We provides two options for outputting the superpoly via the variable *solver_mode*. 

   Setting *solver_mode* to *mode::OUTPUT_EXP* means once the MITM framework finishes, the concrete expression of the superpoly will be output to the folder "TERM" as a file named "superpoly.txt".
   
   Setting *solver_mode* to *mode::OUTPUT_FILE* means only the necessary information for superpoly recovery is recorded, so that we can calculate the concrete expression of the superpoly from the "TERM" and "STATE" folders after the program terminates. This mode can reduce the memory usage during the running of the program.
   
   If you want to recover the expression of the superpoly from the folders "TERM" and "STATE" under *mode::OUTPUT_FILE*, you can uncommenting the line of code *MITM_framework.read_sols_and_output()* in the file "SuperpolyBGL.cpp", and the superpoly will also be output to the folder "TERM" as a file "superpoly.txt". However, this may cause the system to kill the process due to high memory usage if the superpoly is very complex (e.g., the superpoly of trivium).

3. Type `sh exec.sh` in the console, this should create three folders named "STATE", "LOG" and "TERM", and generate an executable program "mitm". "LOG" contains log files; "STATE" stores the hash table $P$ after each expansion; "TERM" 
stores the information of the superpoly that has been extracted by the CMP-based approach.

4. Type `./mitm` in the console to start the superpoly recovery. While the program is running, the status of the program will be recorded in the log files.

## 3. Dependencies
Note that the header file "dynamic_bitset.hpp" used in the codes is from the C++ Boost Library, which can be downloaded from (https://www.boost.org/).
