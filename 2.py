import numpy as np
import matplotlib.pyplot as plt

from Heat import solveHeat
from Vulcanization import solveVulcanization
from plotting import plotU, plotAlpha

def showSol(problem, M=51, N=101):
    if problem == "Heat":
        x, y, t, U = solveHeat(M=M, N=N)

        # Create animated plot
        fig, ax, anim = plotU(x, y, t, U)
        plt.show() 

    elif problem == "Vulcanization":
        x, y, t, U, alpha = solveVulcanization(M=M, N=N)

        # Create animated plot
        fig, ax, anim = plotU(x, y, t, U)
        plt.show()

        fig, ax, anim = plotAlpha(x, y, t, alpha)
        plt.show()