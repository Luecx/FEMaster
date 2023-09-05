import numpy as np
import model
from scipy.optimize import minimize

model = model.Model(width=10, height=10, length=20)

corners, mids = model.get_special_node_ids()

model.set_support(corners[0], 0 ,0 ,0)
model.set_support(corners[1], 0 ,0 ,0)
model.set_support(corners[2], 0 ,0 ,0)
model.set_support(corners[3], 0 ,0 ,0)

model.set_force(corners[4], 0 ,1 ,1)

model.set_material(youngs=210000, nu=0.3)

model.set_exponent(3)
model.set_proximity_radius(3)
model.set_solver(use_cpu=True, use_direct=True)

target_density = 0.5
n_elements = model.width * model.height * model.length

last_iter_result = None

# Objective function
def objective(x):
    global last_iter_result
    last_iter_result = model.run(x)
    return last_iter_result['COMPLIANCE_ADJ'].sum()

# Gradient of the objective function
def gradient(x):
    return np.array(last_iter_result['DENS_GRAD'])

# Constraints
def constraints(x):
    avg_density = np.mean(x)
    return target_density - avg_density

# Jacobian of the constraints
def cons_jacobian(x):
    return np.ones_like(x) * -1/len(x)

# Initial guess
n_elements = model.width * model.height * model.length
initial_densities = np.ones(n_elements) * target_density

# Bounds
bounds = [(0.01, 1) for _ in range(n_elements)]

# Constraints (equality)
cons = {'type': 'eq', 'fun': constraints, 'jac': cons_jacobian}

print("beginning optimization")

# Call the optimizer
res = minimize(fun=objective,
               x0=initial_densities,
               jac=gradient,
               bounds=bounds,
               constraints=cons,
               method='SLSQP',
               options={'disp': True})

print("Optimized Densities:", res.x)
print("Final Objective Value:", res.fun)

model.plot_cuboid(densities=res.x, threshold=0.5)