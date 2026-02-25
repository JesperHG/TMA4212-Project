# TMA4212-Project
Compulsory project in TMA4212 Numerical solution of differential equations by difference methods

## AI usage
We used AI to:
* Make the project file structure
* Make plotting and utils
* Find two errors in the code
    * Forgot to multiply alpha_dot by tdiff to remain dimensionless
    * Forgot to multiply c*adot by dx**2

The important parts of the project, `discretization.py` and `solvers.py`, were written by us.

## Project structure

```
README.md
answers/               Written answers (task1.md, task2.md, …)
src/                   Python package – all reusable code
    params.py            Physical / model parameters
    discretization.py    Laplacian matrix & αdot
    solvers.py           Heat and Vulcanization PDE solvers
    plotting.py          Plotting functions (animations, snapshots, analysis)
    utils.py             Save/load solutions, caching helpers
tasks/                 Runnable scripts, one per task
    task2.py             Solve & generate all Task 2 plots
output/                Generated files (git-ignored)
    solutions/           Cached .npz solution files
    plots/               Saved figures and GIFs
```

## Usage

Run a task from the project root:
```bash
python -m tasks.task2
```

Cached solutions are stored in `output/solutions/`. Delete them (or pass
`force=True` to `get_solution`) to re-solve.

TODO:
- More task scripts as new tasks arrive
