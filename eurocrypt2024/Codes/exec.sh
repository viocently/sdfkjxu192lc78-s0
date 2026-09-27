export GLIBC_TUNABLES=glibc.malloc.arena_max=8
rm nohup.out
rm STATE/*
rm TERM/*
rm LOG/*
g++  SuperpolyBGL.cpp deg.cpp -o mitm -std=c++17 -O2 -lm -lpthread -I/$GUROBI_HOME/include/   -L/$GUROBI_HOME/lib -lgurobi_c++ -lgurobi91 -lm



