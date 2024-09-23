import torch
import torch.nn as nn
import torch.nn.functional as F
from torch_geometric.nn import GATConv

def full_block(in_features, out_features, p_drop):
    """Create a fully connected block with BatchNorm, ELU activation, and Dropout."""
    return nn.Sequential(
        nn.Linear(in_features, out_features),
        nn.BatchNorm1d(out_features),
        nn.ELU(),
        nn.Dropout(p=p_drop)
    )

class GNNVAEModel(nn.Module):
    """Graph Neural Network Variational Autoencoder Model."""
    def __init__(self, input_dim, params, num_classes=1):
        super().__init__()
        self.var = params.var
        self.num_classes = num_classes
        self.params = params

        # Encoder
        self.encoder = nn.Sequential(
            full_block(input_dim, params.feat_hidden1, params.p_drop),
            full_block(params.feat_hidden1, params.feat_hidden2, params.p_drop)
        )

        # GAT Layers
        self.gat1 = GATConv(
            in_channels=params.feat_hidden2,
            out_channels=params.gat_hidden1,
            heads=params.nheads,
            dropout=params.p_drop
        )
        self.gat2 = GATConv(
            in_channels=params.gat_hidden1 * params.nheads,
            out_channels=params.gat_hidden2,
            heads=1,
            concat=False,
            dropout=params.p_drop
        )
        if self.var:
            self.gat3 = GATConv(
                in_channels=params.gat_hidden1 * params.nheads,
                out_channels=params.gat_hidden2,
                heads=1,
                concat=False,
                dropout=params.p_drop
            )

        # Decoder
        self.decoder = nn.Sequential(
            full_block(params.gat_hidden2, params.feat_hidden2, params.p_drop),
            full_block(params.feat_hidden2, params.feat_hidden1, params.p_drop),
            nn.Linear(params.feat_hidden1, input_dim)
        )

        # Clustering Layer
        self.cluster = nn.Sequential(
            full_block(params.gat_hidden2, params.feat_hidden2, params.p_drop),
            nn.Linear(params.feat_hidden2, self.num_classes)
        )

    def encode(self, x, edge_index):
        """Encode the input features into latent space."""
        x = self.encoder(x)
        x = self.gat1(x, edge_index)
        x = F.relu(x)
        x = F.dropout(x, p=self.params.p_drop, training=self.training)

        mu = self.gat2(x, edge_index)
        if self.var:
            logvar = self.gat3(x, edge_index)
            return mu, logvar
        else:
            return mu, None

    def reparameterize(self, mu, logvar):
        """Reparameterization trick for sampling."""
        if self.training and logvar is not None:
            std = torch.exp(0.5 * logvar)
            eps = torch.randn_like(std)
            return eps * std + mu
        else:
            return mu

    def forward(self, x, edge_index):
        """Forward pass through the model."""
        mu, logvar = self.encode(x, edge_index)
        z = self.reparameterize(mu, logvar)
        x_reconstructed = self.decoder(z)
        pred_label = F.softmax(self.cluster(z), dim=1)
        return pred_label, x_reconstructed, z, mu, logvar
