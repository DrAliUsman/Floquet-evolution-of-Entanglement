from qutip import *
from numpy import *
from numpy import mean


## Return the cosine of an operator ##
def cosine(Op):
    return 0.5*( (i*Op).expm() + (-i*Op).expm() )


## Parameters ##
i = 1.0j 				# imaginary unit
q = 2                                   # dimensionless ratio between frequencie           
alpha = 2*pi/q			        # scaled time (kick to kick period)
heff = 1.0                              # scaled Planck's constant
etta = 1/sqrt(2)                        # classicality parameter for kicked oscillators
V0 = 0.25   		                # quantum kick strength
N = 64 		                        # number of states of quantum system

## Basis State ##
for n in range (1, 50):                 # number of kicks to simulate
    psi = (coherent(N, 5)).unit()       # initial state
    psic = psi.dag()

    ## Operators ##
    a = destroy(N) 			# lowering operator
    x = position(N)                     # Position operator
    v = 0.5*(a + a.dag())
    u = 0.5*i*(a.dag() - a)

    ## Evolution Operator ##
    F0 = (-i*alpha*(a.dag()*a + 0.5)).expm()
    F1 = (-i*(V0/heff)*(cosine((etta)*(a + a.dag())))).expm()
    F = F0*F1				# Floquet operator
    Fd = F.dag()

    S = 1-((abs((((Fd)**n)*((F0)**n)).matrix_element(psic, psi)))**2)
    print S                 #linear entropy as a measure of LQU


