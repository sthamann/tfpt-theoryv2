import sys
sys.path.append('/Users/stefanhamann/Projekte/Q5/Pyrate3/pyrate/results/E8Cascade2Loop/E8Cascade2LoopGravity/PythonOutput')

from E8Cascade2LoopGravity import RGEsolver

##############################################
# First, create an instance of the RGEsolver #
##############################################

rge = RGEsolver('rge', tmin=0, tmax=20, initialScale=0)


##########################################################
# We fix the running scheme and initial conditions below #
##########################################################

# Running scheme :

rge.loops = {'GaugeCouplings': 2,
             'Yukawas': 2,
             'QuarticTerms': 2,
             'ScalarMasses': 2,
             'Vevs': 2}

# Gauge Couplings

rge.g1.initialValue = 0
rge.g2.initialValue = 0
rge.g3.initialValue = 0

# Yukawa Couplings

rge.Yu.initialValue = [[0., 0., 0.],
                       [0., 0., 0.],
                       [0., 0., 0.]]

rge.Yd.initialValue = [[0., 0., 0.],
                       [0., 0., 0.],
                       [0., 0., 0.]]

rge.Ye.initialValue = [[0., 0., 0.],
                       [0., 0., 0.],
                       [0., 0., 0.]]

rge.yN1.initialValue = [[0.],
                        [0.],
                        [0.]]

rge.yN2.initialValue = [[0.],
                        [0.],
                        [0.]]

rge.yN3.initialValue = [[0.],
                        [0.],
                        [0.]]


# Quartic Couplings

rge.lambda_.initialValue = 0
rge.lPhi.initialValue = 0
rge.lHphi.initialValue = 0

# Scalar Mass Couplings

rge.mu2.initialValue = 0
rge.MPhi.initialValue = 0

# Vacuum-expectation Values

rge.vSM.initialValue = 0
rge.vPQ.initialValue = 0


############################
# Solve the system of RGEs #
############################

rge.solve(step = .05)

# Another way to call rge.solve() :
# rge.solve(Npoints = 500)

####################
# Plot the results #
####################

rge.plot(subPlots=True, printLoopLevel=True)


#############################################
# Possibly save the results for a later use #
#############################################

# Save the results in some file

# rge.save('rgeResults.save')

# Later, load the rge object with :

# rge = RGEsolver.load('rgeResults.save')

