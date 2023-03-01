# Examples

This directory contains two examples of how to use SCL.

* `01_finite_fields.cc` shows how SCL can be used to define and work with finite
  fields.
  
* `02_matrix_and_vector.cc` shows usage of matrices and vectors over finite
  fields.
  
* `03_secret_sharing.cc` shows how to work with both Shamir and additive secret
  sharing.
  
* `04_networking.cc` shows a how two parties can communicate.

* `05_discovery.cc` shows how a network config can be discovered by using a
  central server that receives ports, IDs and hostnames from other parties and
  then sends complete network information back to everyone.

## Running the examples

Build `libscl.so.0.3` as instructed on the front page. Put `libscl.so.0.3` into
the examples directory, and then run:

```
cmake . -B build
cd build
make
```

which will build an executable for each example.
