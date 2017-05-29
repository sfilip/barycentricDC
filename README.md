barycentricDC
==============================================

## Information
This folder contains an implementation of the
differential correction algorithm for computing
best rational approximations on a discrete set
of points. The particularity of the code is the
use of an adaptive barycentric representation
for representing the rational approximations
computed during executions. For more information,
see TODO.

The code is written in MATLAB and requires that
CVX is also installed. The testing was done using
the Mosek backend for CVX. If available, it can
be set by using the command
```
cvx_solver Mosek
```

## Licensing

The provided code is MIT licensed.

## References
TODO
