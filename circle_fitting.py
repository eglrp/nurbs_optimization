import os
import sys
import math
import numpy as np
from geomdl import BSpline
from geomdl import utilities
from geomdl import exchange 
from geomdl import utilities as utils

# Try to load the visualization module
try:
    render_curve = True
    from geomdl.visualization import VisMPL
except ImportError:
    render_curve = False

r=3.0
Cx = lambda t: r + r*math.cos(0.5*math.pi*t + math.pi)
Cy = lambda t: r + r*math.sin(0.5*math.pi*t + math.pi)

# choose curve parameter
degree = 2
nCtrlpts = 3

# gauss quadrature
xi=[-0.90618, -0.538469, 0, 0.538469, 0.90618]
wi=[0.236927, 0.478629, 0.568889, 0.478629, 0.236927]

# solve
knotvector = utils.generate_knot_vector(degree, nCtrlpts)
span = lambda xi: utils.find_span(degree, tuple(knotvector), nCtrlpts, xi*0.5+0.5)
bfuns = lambda xi: utils.basis_functions(2, tuple(knotvector), span(xi), xi*0.5 + 0.5)

t=-0.5
print knotvector
print span(t)
print bfuns(t)

A = np.zeros((nCtrlpts, nCtrlpts))
A_temp = np.zeros((nCtrlpts, nCtrlpts))

B = np.zeros(3)
for g in range(0,len(xi)):
    for i in range(0,3):
        for j in range(0,3):
            A_temp[i,j] = bfuns(xi[g])[j]*bfuns(xi[g])[i]
    A += wi[g]*A_temp
    A_temp[:,:] = 0

    for k in range(0, 3):
        B[k] += wi[g] * Cx(0.5*xi[g]+0.5) * bfuns(xi[g])[k]

print np.linalg.inv(A)
print B
xCtrpts = np.dot(np.linalg.inv(A), B) 

B = np.zeros(3)
for g in range(0,len(xi)):
    for k in range(0, 3):
        B[k] += wi[g] * Cy(0.5*xi[g]+0.5) * bfuns(xi[g])[k]

yCtrpts = np.dot(np.linalg.inv(A), B) 

print xCtrpts
print yCtrpts
ctrpts = zip(xCtrpts, yCtrpts)

print Cx(0)
print bfuns(-1)
print bfuns(1)
print xCtrpts[0] * np.array(bfuns(-1)) 

# visualize results

# Create a B-Spline curve instance
curve = BSpline.Curve()

# Set up curve
curve.ctrlpts = [ctrpts[0], ctrpts[1], ctrpts[2]]
curve.degree = 2

# Auto-generate knot vector
curve.knotvector = utilities.generate_knot_vector(curve.degree, len(curve.ctrlpts))

# Length computation
xi=[-0.90618, -0.538469, 0, 0.538469, 0.90618]
wi=[0.236927, 0.478629, 0.568889, 0.478629, 0.236927]
g1= lambda t: curve.derivatives(0.25*t+0.25, 1)
dcurve1 = np.array([g1(xi[i])[1] for i in range(len(xi))])
length1 = 0.25*np.dot(wi, np.linalg.norm(dcurve1, axis=1))
print length1
print 'ana lengh=', 2*3.14*3/4/2
print 'knot vector = ', curve.knotvector

# Set evaluation delta
curve.delta = 0.01

# Evaluate curve
curve.evaluate()

# Draw the control point polygon and the evaluated curve
if render_curve:
    vis_comp = VisMPL.VisCurve2D()
    curve.vis = vis_comp
    curve.render()