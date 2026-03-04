# TMA4212-Project
Compulsory project in TMA4212 Numerical solution of differential equations by difference methods.

Written by Vegard Røttingen Enden, Jesper Hassel Grotthing and Fredrik Marcelius Eriksen

## Project structure

```
README.md
ai_declaration.md      AI declaration
report.md              Theoretical groundwork
main.py                Main script – solve and generate all plots
src/                   Python package – all reusable code
    params.py               Physical / model parameters
    discretization.py       System matrix & αdot
    solvers.py              Heat and Vulcanization PDE solvers
    error_convergence.py    Numerical error convergence test of heat equation solver
    plotting.py             Plotting functions (animations, snapshots, analysis)
    utils.py                Save/load solutions, caching helpers
output/                Generated files (git-ignored)
    solutions/           Cached .npz solution files
    plots/
        gifs/               Animated GIFs (temperature, curing, overheating)
        snapshots/          Static field images (T/α grids, cross-sections, isosurface)
        analysis/           Line / scalar analysis plots (cure curves, max T, centreline)
```

## Output plots

### GIFs (`output/plots/gifs/`)
| File | Description |
|------|-------------|
| `heat_T.gif` | Pure heat equation – animated temperature |
| `vulc2d_T.gif` | 2-D vulcanization – animated temperature |
| `vulc2d_alpha.gif` | 2-D vulcanization – animated curing field |
| `vulc2d_T_above_Tm.gif` | 2-D vulcanization – only regions where T > Tm |
| `vulc3d_T_slices.gif` | 3-D vulcanization – animated T at z-slices |
| `vulc3d_alpha_slices.gif` | 3-D vulcanization – animated α at z-slices |
| `vulc3d_T_above_Tm_slices.gif` | 3-D vulcanization – only regions where T > Tm at z-slices |

### Snapshots (`output/plots/snapshots/`)
| File | Description |
|------|-------------|
| `heat_vs_vulc_T.png` | Heat vs vulcanization temperature comparison |
| `vulc2d_alpha.png` | 2-D curing snapshots over time |
| `vulc3d_T_slices.png` | 3-D temperature snapshots (z-slices × time) |
| `vulc3d_alpha_slices.png` | 3-D curing snapshots (z-slices × time) |
| `vulc3d_isosurface.png` | 3-D isosurface scatter (3×3: α levels × time) |

### Analysis (`output/plots/analysis/`)
| File | Description |
|------|-------------|
| `heat_error_convergence.png` | Heat equation error convergence |
| `vulc2d_cure_center.png` | 2-D cure at centre & edge vs time |
| `vulc2d_max_T.png` | 2-D max temperature vs time |
| `vulc2d_parameter_sensitivity.png` | 2-D vulcanization parameter sensitivity |
| `vulc3d_centerline_z.png` | 3-D centre-line profiles along z |
| `vulc3d_cure_center.png` | 3-D cure at centre, face & corner vs time |
| `vulc3d_max_T.png` | 3-D max temperature vs time |

## Usage

Run from the project root:
```bash
python main.py
```

Cached solutions are stored in `output/solutions/`. Delete them (or pass
`force=True` to `get_solution`) to re-solve.
