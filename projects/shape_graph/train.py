from geom import *

import torch
import torch.nn.functional as F
from torch_geometric.nn import GCNConv
from torch_geometric.data import Data

class FEMGraphNN(torch.nn.Module):
    def __init__(self, input_dim=2, hidden_dim=64, output_dim=1):
        super(FEMGraphNN, self).__init__()
        # Define graph convolutional layers
        self.conv1 = GCNConv(input_dim, hidden_dim)
        self.conv2 = GCNConv(hidden_dim, hidden_dim)
        self.conv3 = GCNConv(hidden_dim, output_dim)

    def forward(self, x, edge_index):
        x = F.relu(self.conv1(x, edge_index))
        x = F.relu(self.conv2(x, edge_index))
        return self.conv3(x, edge_index)  # Output displacement per node


def build_graph(model):

    # fix curvature at boundary
    model.point_group.points[0 ].curvature = model.point_group.points[ 1].curvature
    model.point_group.points[-1].curvature = model.point_group.points[-2].curvature

    # normalise curvature and stress
    curvatures = np.array([p.curvature for p in model.point_group.points])
    stresses = np.array([p.stress for p in model.point_group.points])

    curvatures = curvatures / np.max(np.abs(curvatures))
    stresses   = stresses / np.max(np.abs(stresses))

    # Node features: curvature and initial stress
    x = np.array([[c, s] for c, s in zip(curvatures, stresses)], dtype=np.float32)

    # Edge index: connections between neighboring nodes (e.g., a chain structure)
    edge_index = []
    num_points = len(model.point_group.points)
    for i in range(num_points - 1):
        edge_index.append([i, i + 1])
        edge_index.append([i + 1, i])

    edge_index = torch.tensor(edge_index, dtype=torch.long).t().contiguous()

    print(x)

    data = Data(x=x, edge_index=edge_index)
    return data

# Initialize GNN model, optimizer, and loss function
gnn_model = FEMGraphNN(input_dim=2, hidden_dim=64, output_dim=1)
optimizer = torch.optim.Adam(gnn_model.parameters(), lr=0.001)

# print overview
print(gnn_model)

import pickle
models = pickle.load(open("models100.pkl", "rb"))

# Training loop
num_epochs = 100
for epoch in range(num_epochs):
    total_loss = 0
    for model in models:
        # Create graph data for each model
        data = build_graph(model)

        # Forward pass through the GNN
        displacements = gnn_model(data.x, data.edge_index).detach().numpy()

        # Apply displacements to each node in the model
        for i, p in enumerate(model.point_group.points):
            p.pos += displacements[i] * np.array([p.tangent[1], -p.tangent[0]])  # Adjust along normal

        # Recompute stress after displacement
        model.compute_stress()
        new_stress = np.array([p.stress for p in model.point_group.points])

        # Calculate reward based on stress reduction
        initial_stress = np.array([p.stress for p in model.point_group.points])
        reward = compute_reward(initial_stress, new_stress)

        # Backpropagation to minimize the loss based on the reward
        reward_tensor = torch.tensor(reward, dtype=torch.float32)
        optimizer.zero_grad()
        loss = -reward_tensor.mean()
        loss.backward()
        optimizer.step()
        total_loss += loss.item()

    print(f"Epoch {epoch+1}/{num_epochs}, Loss: {total_loss:.4f}")
