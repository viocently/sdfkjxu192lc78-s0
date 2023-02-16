# Reference codes for the paper "Revisiting Core Monomial Prediction: Graph-Based MITM Framework for Superpoly Recovery"

## 1. Structure of the codes

The header file "BooleanPolynomial.h" defines the class for the representation and operations of Boolean polynomials

The header file "flag.h" defines the class for the representation and operations of flags

The header file "cipher.h" defines the update functions for different ciphers and gives the concrete implementation of our CMP-based approach

The header file "framework.h" gives the concrete implementation of the MITM framework. 

The cpp file "SuperpolyBGL.cpp" configures parameters for different ciphers and then calls the MITM framework to recover the superpoly.

## 2. Usage of the codes
The codes are 
