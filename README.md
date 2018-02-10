# MTLS

#### C++ implementation of the More-Thuente linesearch method as described in paper 'Line Search Algorithms with Guaranteed Sufficient Decrease', doi:10.1145/192115.192132, by Jorge J. More and David J. Thuente.

## Requirement
#### [Eigen](https://www.eigen.tuxfamily.org/) C++ library for linear algebra.

## Install
#### Just include the header, MTLS.h, in your project and compile as usual.

## Usage
#### Firstly, have a look at the examples in example.cpp.
#### The simplest way to initialize and run the linesearch is:
```C++
step = linesearchMoreThuente(p, x0, f0, F, G, ...Args)
```
1. p: search direction given as a vector
2. f0: function value
3. x0: initial set of parmeters
4. F: function pointer to function (needs to take x0 in as first parameter)
5. G: function pointer to get derivative (must be without arguments)

#### This will run linesearch with a default set of parameters. Additional arguments possibly used by F can be sent in as the last arguments.
#### To control the parameters used create a struct with the following syntax
```C++
struct Params {
    /* struct of default parameters */
    T maxIterations = 100; // maximum number of iterations
    T mu = 0.5; // step scaling factor, (0<mu<=1/2<eta)
    T eta = 1.0; // termination parameter (0<eta<1)
    T delta = 4.0; // delta: scaling of step updating ([1.1,4.0])
    T bisectWidth = 0.66; // extrapolation tolerance for bounds
    T bracketTol = 1e-14; // termination tolerance for brackets
    T aMin0 = 0.0; // lower bound for step (aMin>=0 and aMin<aMax)
    T aMax0 = 100.0; // upper bound for step (aMax>aMin)
} params; // end struct Params
```

### (the values given are the default ones in this case, just change them to your liking) and call as
```C++
step = linesearchMoreThuente(&params, p, x0, f0, F, G)
```

### Support for member functions is implemented, but a pointer to object is required. Syntax is as follows(with an instance of FOO set to foo)
```C++
step = linesearchMoreThuente(p, x0, f0, foo, &FOO::F, &FOO::G)
step = linesearchMoreThuente(&params, p, x0, f0, foo, &FOO::F, &FOO::G)
```
