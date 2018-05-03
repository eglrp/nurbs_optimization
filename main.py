import os
import sys
import math
import numpy as np
from geomdl import BSpline
from geomdl import utilities
from geomdl import exchange 

# Try to load the visualization module
try:
    render_curve = True
    from geomdl.visualization import VisMPL
except ImportError:
    render_curve = False

# Fix file path
os.chdir(os.path.dirname(os.path.realpath(__file__)))

# Create a B-Spline curve instance
curve = BSpline.Curve()

# Set up curve
curve.ctrlpts = [[1.0,1.0], [2.5,0.5], [5.0,1.0]]
curve.degree = 1

# Auto-generate knot vector
curve.knotvector = utilities.generate_knot_vector(curve.degree, len(curve.ctrlpts))

# Length computation
xi=[-0.57735, 0.57735]
wi=[1, 1]
g1= lambda t: curve.derivatives(0.25*t+0.25, 1)
dcurve1 = np.array([g1(xi[i])[1] for i in range(len(xi))])
length1 = 0.25*np.dot(wi, np.linalg.norm(dcurve1, axis=1))
print length1
l1 = np.linalg.norm(np.array(curve.ctrlpts[1]) - np.array(curve.ctrlpts[0]))
l2 = np.linalg.norm(np.array(curve.ctrlpts[2]) - np.array(curve.ctrlpts[1]))
print 'l1 = ', l1
print 'l2 = ', l2
print 'ana lengh=', l1+l2

sys.exit(0)
# Set evaluation delta
curve.delta = 0.01

# Evaluate curve
curve.evaluate()

# Draw the control point polygon and the evaluated curve
if render_curve:
    vis_comp = VisMPL.VisCurve2D()
    curve.vis = vis_comp
    curve.render()