from gnn_spatial import *

import matplotlib.pyplot as plt

def create_line(n=5, spacing=1.0):
    angle = random.uniform(0, 2*np.pi)
    p1 = np.array([0, 0])
    points = [p1]
    for i in range(n):
        p2 = p1 + np.array([np.cos(angle), np.sin(angle)]) * spacing
        points.append(p2)
        angle += random.uniform(-np.pi/3, np.pi/3)
        p1 = p2
    return np.asarray(points)

def augment_line(line, reps):
    # rotate and move the line around
    lines = []
    for rep in range(reps):
        angle = random.uniform(0, 2*np.pi)
        offset = np.random.randn(2)
        new_line = line.copy()
        for p in new_line:
            p0 = p[0] * np.cos(angle) - p[1] * np.sin(angle)
            p1 = p[0] * np.sin(angle) + p[1] * np.cos(angle)
            p[0] = p0 + offset[0]
            p[1] = p1 + offset[1]
        lines.append(new_line)
    return lines

def step(line):
    # moves all points to the average of the neighbors
    new_line = line.copy()
    for i in range(1, len(line)-1):
        new_line[i] = (line[i-1] + line[i+1]) / 2
    return new_line

def delta(line1, line2):
    return line2 - line1


lines = [create_line() for _ in range(50)]
augmented_lines = [augment_line(line, 10) for line in lines]
all_lines = [aug_line for augments in augmented_lines for aug_line in augments]

# Prepare data
all_line_diffs = []
all_line_positions = []
for line in all_lines:
    line_target = step(line)
    line_diff = delta(line, line_target)
    all_line_diffs.append(torch.tensor(line_diff, dtype=torch.float32))
    all_line_positions.append(torch.tensor(line, dtype=torch.float32))
node_connections = [[i-1, i+1] if 0 < i < len(line)-1 else ([i+1] if i == 0 else [i-1]) for i in range(len(lines[0]))]

# Stack into batched tensors
batched_positions = torch.stack(all_line_positions)
batched_diffs = torch.stack(all_line_diffs)
batched_features = torch.zeros((len(all_lines), len(all_lines[0]), 1), dtype=torch.float32)

class LineDataset(Dataset):
    def __init__(self, positions, features, diffs):
        self.positions = positions
        self.features = features
        self.diffs = diffs

    def __len__(self):
        return len(self.positions)

    def __getitem__(self, idx):
        return (
            self.positions[idx],
            self.features[idx],
            self.diffs[idx],
        )

dataset = LineDataset(batched_positions, batched_features, batched_diffs)

# Training Loop
# Define dataset and dataloader
dataset = LineDataset(batched_positions, batched_features, batched_diffs)
dataloader = DataLoader(dataset, batch_size=4, shuffle=True)

# Updated training loop
def train_network():
    # Define network, optimizer, and loss function
    model = SGCNNetwork([1, 16, 16, 2])
    optimizer = torch.optim.Adam(model.parameters(), lr=0.001)
    loss_fn = nn.MSELoss()

    epochs = 500
    for epoch in range(epochs):
        model.train()
        total_loss = 0.0

        for batch in dataloader:
            positions, features, diffs = batch
            optimizer.zero_grad()

            outputs = model(features,node_connections, positions)

            # Compute loss
            loss = loss_fn(outputs, diffs)

            # Backward pass and optimization
            loss.backward()
            optimizer.step()

            total_loss += loss.item()

        # Print epoch statistics
        print(f"Epoch {epoch}/{epochs}, Loss: {total_loss:.6f}")
    return model
# model = train_network()
# torch.save(model.state_dict(), 'line_model.pth')
model = SGCNNetwork([1, 16, 16, 2])
model.load_state_dict(torch.load('line_model.pth'))

# Perform 5 steps with the real function and the network
steps = 5

# Initialize real and predicted lines
new_line = create_line()
real_lines = [new_line]
predicted_lines = [new_line.copy()]

# Perform the steps
for _ in range(steps):
    # Real steps
    real_line = step(real_lines[-1])
    real_lines.append(real_line)

    # Predicted steps
    with torch.no_grad():
        input = torch.zeros((1, len(predicted_lines[-1]), 1), dtype=torch.float32)
        predicted_line_diff = model(
            input,
            node_connections,
            torch.tensor(predicted_lines[-1], dtype=torch.float32).unsqueeze(0)
        ).squeeze(0).numpy()

    # Convert deltas to the next predicted line
    predicted_line = predicted_lines[-1] + predicted_line_diff
    predicted_lines.append(predicted_line)

# Plot the results
plt.figure(figsize=(10, 6))

# Plot real lines
for i, real_line in enumerate(real_lines, start=1):
    plt.plot(real_line[:, 0], real_line[:, 1], label=f"Real Step {i}", color="blue", alpha=0.6)

# Plot predicted lines
for i, predicted_line in enumerate(predicted_lines, start=1):
    plt.plot(predicted_line[:, 0], predicted_line[:, 1], label=f"Predicted Step {i}", alpha=0.6)

plt.title("Real vs Predicted Lines (5 Steps)")
plt.xlabel("X Coordinate")
plt.ylabel("Y Coordinate")
plt.legend()
plt.grid(True)
plt.show()
