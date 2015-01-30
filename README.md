cart
====

High-Order Cartesian Flow Solver
----
The following options for ilhs (Implicit Solver) are implemented:
   
0) LU-SGS (Lower-Upper Symmetric Gauss-Seidel)

1) ADI (Alternating Direction Implicit)

2) DD-DDI (Diagonally-Dominant ADI)

3) Gauss-Seidel Line Relaxation

4) Same as 3 but a more streamlined implementation

Upcoming:

5) Point Gauss-Seidel with special boundary treatment for concurrent	data transfer

Testing and demos:
-----------------
INSTALL_DIR is the directory where you told CMAKE to install the
executable. See file "BUILDING" for  example on specifying the
INSTALL directory. 

1. Accuracy check
   verifies accuracy of the discretized Navier-Stokes 
   against exact divergence using complex arithemetic using 
   method of manufactured solutions

   cd $INSTALL_DIR/tests/accuracy;
   ../../bin/nscheck.exe < verify.input

2. Linearization test
   verifies consistency of the linearization of the 2nd order inviscid terms (for now
   cd $INSTALL_DIR/tests/linearization
   ../../bin/linearization_test 

3. F90 executable demo

   cd $INSTALL_DIR/demo/example;
   ../../bin/cart.exe

   Runs the Taylor-Green problem for 10 steps

4. python execution demo

   cd $INSTALL_DIR/demo/python_example
   python test.py
    
   Runs Taylor-Green through python

5. Time-Spectral execution demo (requires libTS)

   cd $INSTALL_DIR/demo/time_spectral
   mpirun.lsf (or equivalent) python cartTS.py

   * First requires appropriate path information to libTS (default assumes same installation depth as cart)