# AI Declaration

We used github copilot for AI assistance.
It was used to:
* Make the project file structure
* Write `README.md`
* Save, load and plot our results
    * See `plotting.py` and `utils.py`
* Write the function `solve_vulcanization_parametric()`. It is based on our own solvers, and used to do the parameter analysis.
* Find two errors when troubleshooting our solvers:
    * Forgot to make the sourceterm function dimensionless when implementing it in the code. So, multiply alpha_dot by tdiff.
    * Forgot to multiply our update step of b by dx**2.
* Write out the system matrix for the 3D case in markdown based on how it was implemented by us in the code (it was found without use of AI).
* Write the function `solve_vulcanization_parametric()`, it is the same as our solvers, but takes gamma, c and E as input for parameter study.

The solvers and discretization of the problem was done by us. We would argue that this is the course relevant part of the project.