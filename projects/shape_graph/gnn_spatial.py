
import torch
import torch.nn as nn
import torch.nn.functional as F

from torch.utils.data import DataLoader, Dataset
from torch import Tensor

import random
import numpy as np

# implements the paper with spatial information during the aggregation step
# where instead of simply averaging the features, one uses the relative direction vector as well as a trainable
# weight matrix applied to the direction to scale the features
class SGCNLayer(nn.Module):
    def __init__(self, c_in, c_out, n_filters, activation=F.relu):
        super().__init__()

        self.c_in = c_in
        self.c_out = c_out
        self.n_filters = n_filters
        self.activation = activation

        # maps the position difference to the features
        self.affine_directions = nn.Linear(3, c_in * n_filters)
        self.affine_features   = nn.Linear(c_in * n_filters, c_out)

    def aggregate(self, node_feats, pos_this, pos_other):
        # node_feats is batch_size x 1 x c_in
        # pos_this   is batch_size x 3
        # pos_other  is batch_size x 3
        # check this
        if node_feats.size() != torch.Size([node_feats.size(0), 1, self.c_in]):
            raise ValueError(f"Node features have wrong shape: {node_feats.size()}, expected {torch.Size([node_feats.size(0), 1, self.c_in])}")
        if pos_this.size() != torch.Size([node_feats.size(0), 3]) and pos_this.size() != torch.Size([node_feats.size(0), 2]):
            raise ValueError(f"Node positions have wrong shape: {pos_this.size()}")
        if pos_other.size() != torch.Size([node_feats.size(0), 3]) and pos_other.size() != torch.Size([node_feats.size(0), 2]):
            raise ValueError(f"Node positions have wrong shape: {pos_other.size()}")

        # compute the relative direction vector
        direction = pos_other - pos_this
        # direction is batch_size x 3

        # make sure its a flat vector with 3 entries. if not, add a zero
        if direction.size(-1) < 3:
            direction = torch.cat([direction, torch.zeros(direction.size(0), 3-direction.size(-1))], dim=-1)

        # transform the directions to each stacked feature
        stacked_feats = node_feats.repeat(1, 1, self.n_filters)
        affine_direction = F.relu(self.affine_directions(direction))

        # scale the features with the direction
        node_feats = stacked_feats[:] * affine_direction[:,None, :]

        return node_feats

    def aggreate_nodes(self, node_id, node_feats, node_connections, node_positions):
        agg = self.aggregate(node_feats[:,node_id:node_id+1,:], node_positions[:,node_id,:], node_positions[:,node_id,:])
        for con in node_connections[node_id]:
            if con == node_id:
                continue
            agg += self.aggregate(node_feats[:,con:con+1,:], node_positions[:,node_id,:], node_positions[:,con,:])
        return agg

    def compute_node(self, node_id, node_feats, node_connections, node_positions):
        aggregation = self.aggreate_nodes(node_id, node_feats, node_connections, node_positions)
        affined     = self.activation(self.affine_features(aggregation))
        return affined
        # # size is batch_size x c_out
        # # transform to batch_size x 1 x c_out
        # return affined.unsqueeze(1)

    def forward(self, node_feats, node_connections, node_positions):
        """Forward.

        Args:
            node_feats: Tensor with node features of shape [batch_size, num_nodes, c_in]
            adj_matrix: Batch of adjacency matrices of the graph. If there is an edge from i to j,
                         adj_matrix[b,i,j]=1 else 0. Supports directed edges by non-symmetric matrices.
                         Assumes to already have added the identity connections.
                         Shape: [batch_size, num_nodes, num_nodes]

        """
        return torch.cat([self.compute_node(i, node_feats, node_connections, node_positions)
                          for i in range(node_feats.size(1))], dim=1)


class SGCNNetwork(nn.Module):
    def __init__(self, layers_config):
        super().__init__()

        self.layers = nn.ModuleList()
        for i in range(len(layers_config) - 1):
            self.layers.append(SGCNLayer(layers_config[i], layers_config[i + 1], n_filters=16))
        # remove activation from the last layer
        self.layers[-1].activation = lambda x: x

    def forward(self, node_feats, node_connections, node_positions):
        for layer in self.layers:
            if isinstance(layer, nn.BatchNorm1d):
                node_feats = layer(node_feats)
            else:
                node_feats = layer(node_feats, node_connections, node_positions)
        return node_feats




#
# # Example Usage
# layer = SGCNLayer(2, 2, 2)
#
# # create a batch of features, batch = 2, 3 nodes each with 4 features
# node_feats = torch.arange(12, dtype=torch.float32).view(2, 3, 2)
# positions = torch.tensor([[[0.0, 0.0], [1.0, 0.0], [0.0, 1.0]],
#                           [[0.0, 0.0], [2.0, 1.0], [1.0, 2.0]]], dtype=torch.float32)
#
# # Aggregate for the second node from the first node
# node_connections = [[1], [0, 2], [1]]
#
# print(layer.forward(node_feats, node_connections, positions))
# # print(layer.forward(node_feats, node_connections, positions))

