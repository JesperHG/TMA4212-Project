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

The solvers and discretization of the problem was done purely by us. We would argue that this is the relevant part of the project to the course.