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
main.py                Main script – solve and generate all plots
src/                   Python package – all reusable code
    params.py            Physical / model parameters
    discretization.py    Laplacian matrix & αdot
    solvers.py           Heat and Vulcanization PDE solvers
    plotting.py          Plotting functions (animations, snapshots, analysis)
    utils.py             Save/load solutions, caching helpers
tasks/                 Written answers and task-specific scripts
    task1.md             Task 1 answers
    task2.md             Task 2 answers
output/                Generated files (git-ignored)
    solutions/           Cached .npz solution files
    plots/               Saved figures and GIFs
```

## Usage

Run from the project root:
```bash
python main.py
```

Cached solutions are stored in `output/solutions/`. Delete them (or pass
`force=True` to `get_solution`) to re-solve.

# Things to consider:
* Plot snapshots of heat as well

* Plot cure at edge of domain
* Plot cure and temp in same plot

* Why is curing seemingly uniform?
* Why does curing have an S shape?
    * When does alpha^n = (1-alpha)^m ?

laplacian!
