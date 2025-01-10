import torch.optim as optim
from gnn import *
from gnn_env import *
from gnn_geom import *

# Initialize the environment
problem = create_geom()
env = CustomGNNEnv(problem)

adj_matrix = problem.adjacency_matrix()

# Initialize the GNN
gnn = GNNNetwork([4, 32, 1])  # Adjust input/output layers based on your env
optimizer = optim.Adam(gnn.parameters(), lr=0.001)
loss_fn = nn.MSELoss()  # Example loss function

# Training loop
num_episodes = 100
for episode in range(num_episodes):
    state, _ = env.reset()  # Get initial state

    total_reward = 0
    done = False
    while not done:
        state_tensor      = torch.tensor(state, dtype=torch.float32).unsqueeze(0)  # Add batch dim
        adj_matrix_tensor = torch.tensor(adj_matrix, dtype=torch.float32).unsqueeze(0)

        # GNN predicts actions
        action = gnn(state_tensor, adj_matrix_tensor).detach().numpy()[0]

        # Take a step in the environment
        next_state, reward, terminated, truncated, _ = env.step(action)
        done = terminated or truncated

        total_reward += reward

        reward_tensor = torch.tensor(reward, dtype=torch.float32)
        loss = -reward_tensor * action.sum()  # Example: Encourage higher rewards with a scaling term

        # Backpropagation
        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        state = next_state  # Update state

    print(f"Episode {episode + 1}/{num_episodes}, Total Reward: {total_reward}")