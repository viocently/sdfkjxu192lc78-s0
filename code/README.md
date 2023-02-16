# Reference codes for the paper "Revisiting Core Monomial Prediction: Graph-Based MITM Framework for Superpoly Recovery"

## 1. Structure of the codes

The header file "BooleanPolynomial.h" defines the class for the representation and operations of Boolean polynomials

The header file "flag.h" defines the class for the representation and operations of flags

The header file "cipher.h" defines the update functions for different ciphers and gives the concrete implementation of our CMP-based approach

The header file "framework.h" gives the concrete implementation of the MITM framework. 

The cpp file "SuperpolyBGL.cpp" configures parameters for different ciphers and then calls the MITM framework to recover the superpoly.

## 2. Usage of the codes
We give a brief introduction on how to run our code on a linux platform.

1. Install Gurobi and configure the required environment variables such as "GUROBI_HOME" and "LD_LIBRARY_PATH".

2. Open the file "SuperpolyBGL.cpp" to set the *cipher_name*, *cube_index* and *rounds*. 
We provides two options for outputting the superpoly via the variable *solver_mode*. 

2.1 Setting *solver_mode* to "mode::OUTPUT_FILE" means only the necessary information for superpoly recovery is recorded, so that we can calculate the concrete expression of the superpoly from the "TERM" and "STATE" folders after the program terminates. This mode can reduce the memory usage during the running of the program.

2.2 Setting *solver_mode* to "mode::OUTPUT_EXP" means once the MITM framework finishes, the concrete expression of the superpoly will be output to the folder "TERM" as a file named "superpoly.txt". 

3. Type `sh exec.sh` in the console, this should create three folders named "STATE", "LOG" and "TERM", and generate an executable program "mitm". "LOG" contains log files; "STATE" stores the hash table $P$ after each expansion; "TERM" 
stores the information of the superpoly that has been extracted by the CMP-based approach.

4. Type `./mitm` in the console to start the superpoly recovery. While the program is running, the status of the program will be recorded in the log files.
