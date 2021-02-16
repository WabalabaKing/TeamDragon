import nlopt
from numpy import *
from run_single import *

def myfunc(x,grad):
    AR = x[0]
    A_wing = x[1]
    weight = run_single(AR,A_wing)
    if grad.size > 0:
        grad[0] = (run_single(AR+0.01,A_wing)-weight)/0.01
        grad[1] = (run_single(AR,A_wing+0.01)-weight)/0.01
    print(count)
    return weight

def myConstr(x,grad):
    AR = x[0]
    A_wing = x[1]
    weight = run_single(AR,A_wing)
    if grad.size > 0:
        grad[0] = -0.5*weight*AR**(-3/2)*A_wing**(-1/2)
        grad[1] = -0.5*weight*AR**(-1/2)*A_wing**(-3/2)
    return 388-(weight/sqrt(AR*A_wing))
            
count = 0;  
opt = nlopt.opt(nlopt.LD_MMA, 2)
opt.set_lower_bounds([0, 0])
opt.set_min_objective(myfunc)
opt.add_inequality_constraint(lambda x,grad: myConstr(x,grad), 1e-8)
opt.verbose = 1
opt.set_xtol_rel(1e-4)
x = opt.optimize([5,200])
minf = opt.last_optimum_value()
print("optimum at ", x[0], x[1])
print("minimum value = ", minf)
print("result code = ", opt.last_optimize_result())