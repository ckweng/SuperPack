# SuperPack

This is the performance evaluation code for the dishonest-majority maliciously-secure protocol presented in the paper "SUPERPACK: Dishonest Majority MPC with Constant Online Communication".

It is built from the implementation of [TurboPack](https://github.com/deescuderoo/turbopack).

## Requirements

The requirements are

- `cmake >= 3.15`
- `lcov`
- `Catch2`
- `secure-computation-library` (included)

These are installed as follows:

### CMake

Follow the [official instructions](https://cmake.org/install/).

### LCov

Follow the [official instructions](http://ltp.sourceforge.net/coverage/lcov.php) (also available in the standard Ubuntu repository).

### Catch2

We require version 2 of Catch2 (the current version is 3).
To install version 2, proceed as follows:

```
$ git clone -b v2.x https://github.com/catchorg/Catch2.git
$ cd Catch2
$ cmake -Bbuild -H. -DBUILD_TESTING=OFF
$ sudo cmake --build build/ --target install
```

### Secure Computation Library

This self-contained library handles communication, finite field types, polynomial evaluation/interpolation, and other primitives required for our protocol.
The library is [open source](https://github.com/anderspkd/secure-computation-library) under the GNU Affero General Public License.
For `SuperPack`, an earlier version of SCL was used, which is included in this repository under `secure-computation-library/`.
This version includes ad-hoc support for packed secret-sharing, which may be included into the main SCL repository in the future for more general use.
To compile it, first enter the `secure-computation-library` directory and then run

```
$ cmake . -DCMAKE_BUILD_TYPE=Release -B build
$ cd build
$ make
```

## Installing

With all the pre-requisites in place, go to the main directory and run

```
$ cmake . -DCMAKE_BUILD_TYPE=Release -B build
$ cd build
$ make
```

This creates two executables in the main directory: `ours.x` and `ours_online.x`. `ours.x` includes our end-to-end protocol (except for the calls to VOLE and OLE) and `ours_online.x` only benchmarks for the online phase.

A script `throttle.py` from [EMP-toolkit](https://github.com/emp-toolkit/emp-readme) can be used to simulate network conditions with variable bandwidth and delay.

## Running

Each executable corresponds to one party in either protocol, and they are set to be connected locally through localhost.
To run our protocol, call

```
$ ./build/ours.x n_parties id size depth percentage
```

Where `n_parties` is the number of parties, `id` is the id of the current party (starting at zero), `size` is the number of multiplication gates, `depth` is the desired depth (this number must divide the number of multiplications, and the multiplications will be spread evenly across all layers), and `percentage` is the proportion of assumed malicious parties among `n_parties`.

There is a script that automates spawning these parties. Run

```
$ ./run_local.sh
```

or

```
$ ./run_online_local.sh
```

for the benchmark of end-to-end execution or the online phase only. To note that the results other than online phase demonstrated in `run_online_local.sh` are dummy executions.
