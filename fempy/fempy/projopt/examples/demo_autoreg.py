import numpy as np
import matplotlib.pyplot as plt
from fempy.projopt import Opt, Sense, Relation, OptParams, Method

rng = np.random.default_rng(0)
n = 80

A = rng.standard_normal((n, 30))
Q = A @ A.T + 1e-2 * np.eye(n)
c = rng.standard_normal(n)

x0 = np.zeros(n)

quad = Opt.create_response(
    "quad",
    value=lambda x: 0.5 * x @ (Q @ x) + c @ x,
    grad =lambda x: Q @ x + c
)

vol = Opt.create_response(
    "mean_x",
    value=lambda x: float(np.mean(x)),
    grad =lambda x: np.ones_like(x) / n
)

a = rng.standard_normal(n)
alpha = -0.1
lin = Opt.create_response(
    "aTx",
    value=lambda x, a=a: float(a @ x),
    grad =lambda x, a=a: a
)

obj   = Opt.objective(quad, Sense.MIN)
c_vol = Opt.constraint(vol, Relation.LT, rhs=0.2)
c_lin = Opt.constraint(lin, Relation.GT, rhs=alpha)

# Pick solver + simple uniform bounds (e.g., densities in [0,1])
params = OptParams(
    method=Method.SLSQP,
    maxiter=400,
    ftol=1e-12,
    disp=False,
    lb_val=0.0,   # set all lower bounds to 0.0
    ub_val=1.0,   # set all upper bounds to 1.0
    verbose=True
)

opt = Opt.new(params)
opt.set_objective(obj)
opt.add_constraint(c_vol)
opt.add_constraint(c_lin)

info = opt.solve(x0)
print("\nDone:", info)

# Stats after solve
stats = opt.get_stats()
print("Stats:", stats)

# Plot log(loss) (shift if needed)
vals = np.asarray(opt.history_f)
shift = 0.0 if np.all(vals > 0) else -vals.min() + 1e-16
plt.figure()
plt.semilogy(vals + shift, marker='o')
plt.xlabel("Iteration"); plt.ylabel("Objective (shifted, log)")
plt.grid(True, which="both"); plt.tight_layout(); plt.show()
