import torch
import pygad
import pygad.torchga
from gnn import GNNNetwork
from gnn_geom import create_geom


# Define the GNN
values_per_layer = [4, 16, 16, 1]  # Adjust as necessary for your problem
model = GNNNetwork(values_per_layer)

def fitness_func(ga_instance, solution, sol_idx):
    """
    Fitness function for the genetic algorithm.
    Minimizes the maximum stress in the problem.
    """
    global model
    problem = create_geom()
    adj_matrix = torch.tensor(problem.adjacency_matrix(), dtype=torch.float32).unsqueeze(0)  # Add batch dimension

    # Load the weights into the model
    model_weights_dict = pygad.torchga.model_weights_as_dict(model=model, weights_vector=solution)
    model.load_state_dict(model_weights_dict)

    # compute initial stresses
    problem.run()
    initial_stress = problem.max_mises()

    # Get inputs from the problem
    inputs = torch.tensor(problem.get_inputs(), dtype=torch.float32).unsqueeze(0)  # Add batch dimension

    with torch.no_grad():
        movements = model(inputs, adj_matrix).squeeze(0).numpy()  # Remove batch dimension

    # adjust movements so that mean = 0
    movements -= movements.mean()
    # clamp to [-1, 1]
    movements = movements.clip(-1, 1)
    # mean again
    movements -= movements.mean()
    # scale to [-0.1, 0.1]
    movements *= 0.1

    # Update the geometry and rerun FEM analysis
    try:
        problem.update_geom({point.id: movements[point.id][0] for point in problem})
        problem.run()
    except Exception as e:
        print(f"Error updating geometry: {e}")
        return 0

    # Fitness value: Negative max Mises stress (as we want to minimize it)
    max_stress = problem.max_mises()

    # print with 3 decimal places

    # print("evaluating solution {0} with stress {1}, (reduction={2})".format(sol_idx, max_stress, initial_stress - max_stress))
    print("evaluating solution {0} with stress {1:.3f}, (reduction={2:.3f})".format(sol_idx, max_stress, initial_stress - max_stress))

    return -max_stress

def callback_generation(ga_instance):
    print("Generation = {generation}".format(generation=ga_instance.generations_completed))
    print("Fitness    = {fitness}".format(fitness=ga_instance.best_solution()[1]))

# Create an instance of the pygad.torchga.TorchGA class to build the initial population.
torch_ga = pygad.torchga.TorchGA(model=model, num_solutions=10)

# Define the genetic algorithm
ga_instance = pygad.GA(
    num_generations=50,             # Number of generations
    num_parents_mating=5,           # Number of parents selected for mating
    initial_population=torch_ga.population_weights,  # Initial population of network weights
    fitness_func=fitness_func,      # Fitness function
    on_generation=callback_generation,  # Callback after each generation
    mutation_percent_genes=10,      # Percentage of genes to mutate
    parent_selection_type="rank",  # Parent selection method
    keep_parents=2,                 # Number of parents to keep for the next generation
    mutation_type="random",        # Mutation strategy
    crossover_type="single_point"  # Crossover strategy
)

# Run the genetic algorithm
ga_instance.run()

# Extract the best solution
solution, solution_fitness, solution_idx = ga_instance.best_solution()
print("Best Solution Fitness:", -solution_fitness)

# Save the best weights to the model
best_solution_weights = pygad.torchga.model_weights_as_dict(model=model, weights_vector=solution)
model.load_state_dict(best_solution_weights)

print("GNN weights updated with the best solution.")

# Visualize the results
problem.plot()
