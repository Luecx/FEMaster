
import torch
import torch.nn as nn
import torch.nn.functional as F

from torch import Tensor

import numpy as np

class GNNLayer(nn.Module):
    def __init__(self, c_in, c_out):
        super().__init__()
        self.projection = nn.Linear(c_in, c_out)

    def forward(self, node_feats, adj_matrix):
        """Forward.

        Args:
            node_feats: Tensor with node features of shape [batch_size, num_nodes, c_in]
            adj_matrix: Batch of adjacency matrices of the graph. If there is an edge from i to j,
                         adj_matrix[b,i,j]=1 else 0. Supports directed edges by non-symmetric matrices.
                         Assumes to already have added the identity connections.
                         Shape: [batch_size, num_nodes, num_nodes]

        """

        # Num neighbours = number of incoming edges
        num_neighbours = adj_matrix.sum(dim=-1, keepdims=True)
        node_feats = self.projection(node_feats)
        node_feats = torch.bmm(adj_matrix, node_feats)
        node_feats = node_feats / num_neighbours
        return node_feats

class GNNNetwork(nn.Module):
    def __init__(self, values_per_layer):
        super().__init__()
        self.layers = nn.ModuleList()
        for i in range(len(values_per_layer) - 1):
            self.layers.append(GNNLayer(values_per_layer[i], values_per_layer[i + 1]))

    def forward(self, node_feats, adj_matrix):
        for layer in self.layers[:-1]:
            node_feats = F.relu(layer(node_feats, adj_matrix))
        node_feats = self.layers[-1](node_feats, adj_matrix)
        return node_feats

# net = GNNNetwork([2, 4, 4, 4, 1])
#
#
#
# node_feats = torch.arange(8, dtype=torch.float32).view(1, 4, 2)
# adj_matrix = Tensor([[[1, 1, 0, 0],
#                       [1, 1, 1, 1],
#                       [0, 1, 1, 1],
#                       [0, 1, 1, 1]]])
#
# print("Node features:\n", node_feats)
# print("\nAdjacency matrix:\n", adj_matrix)
# print("\nOutput features:\n", net(node_feats, adj_matrix))
#
# #
# # print("Node features:\n", node_feats)
# # print("\nAdjacency matrix:\n", adj_matrix)
# #
# # layer = GNNLayer(c_in=2, c_out=2)
# # layer.projection.weight.data = Tensor([[1.0, 0.0], [0.0, 1.0]])
# # layer.projection.bias.data = Tensor([0.0, 0.0])
# #
# # with torch.no_grad():
# #     out_feats = layer(node_feats, adj_matrix)
# #
# # print("Adjacency matrix", adj_matrix)
# # print("Input features", node_feats)
# # print("Output features", out_feats)