Comparative benchmark of QP solvers
===================================

Solvers:
--------

- `qpmad` -> https://github.com/asherikov/qpmad
- `qpOASES` -> https://github.com/coin-or/qpOASES
- `eiquadprog` -> https://github.com/stack-of-tasks/eiquadprog


Problems
--------

- `qpOASES` problem set
  https://github.com/coin-or/qpOASES/tree/misc/testingdata/cpp, converted to
  `JSON` in https://github.com/asherikov/qp_collection/. Only positive-definite
  problems are used.


Notes
-----

### `qpOASES`

* MPC option preset is used, this gives ~50% performance boost compared to default preset.

### `eiQuadProg`

* Has some precision issues on some of the problems, used with increased tolerances.
