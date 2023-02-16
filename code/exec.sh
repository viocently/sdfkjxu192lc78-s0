mkdir LOG
mkdir TERM
mkdir STATE
g++ deg.cpp SuperpolyBGL.cpp  -o mitm  -std=c++17 -O3 -lm -lpthread -I/$GUROBI_HOME/include/     -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi95 -lm



