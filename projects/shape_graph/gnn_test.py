import torch
import torch.nn as nn
import torch.optim as optim

from gnn import GNNNetwork

# Define a dummy input graph
batch_size = 1
num_nodes = 4
node_feature_dim = 2

# Dummy node features (random values for this test)
node_feats = torch.rand(batch_size, num_nodes, node_feature_dim)

# Dummy adjacency matrix (fully connected graph with self-loops)
adj_matrix = torch.zeros(batch_size, num_nodes, num_nodes)

# Add identity to ensure self-loops
adj_matrix += torch.eye(num_nodes).unsqueeze(0)

# Define dummy target (e.g., one target value per node)
dummy_target = torch.rand(batch_size, num_nodes, 1)

# Initialize the GNN network
net = GNNNetwork([2, 24, 24, 1])

# Define a loss function and optimizer
loss_fn = nn.MSELoss()
optimizer = optim.Adam(net.parameters(), lr=0.01)

# Training loop
num_epochs = 100
for epoch in range(num_epochs):
    optimizer.zero_grad()

    # Forward pass
    output = net(node_feats, adj_matrix)

    # Compute loss
    loss = loss_fn(output, dummy_target)

    # Backward pass
    loss.backward()
    optimizer.step()

    if (epoch + 1) % 10 == 0:
        print(f"Epoch [{epoch + 1}/{num_epochs}], Loss: {loss.item():.4f}")

# Test the trained model
with torch.no_grad():
    test_output = net(node_feats, adj_matrix)
    print("\nDummy Target:")
    print(dummy_target)
    print("\nModel Output After Training:")
    print(test_output)
