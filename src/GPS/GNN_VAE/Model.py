#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul  3 11:42:44 2023

@author: songliyang
"""

import torch
import torch.nn.functional as F
from torch.nn.modules.module import Module
from torch.nn.parameter import Parameter
import torch.nn as nn
from torch_geometric.nn import GATConv


def full_block(in_features, out_features, p_drop):
    return nn.Sequential(nn.Linear(in_features, out_features),
                         nn.BatchNorm1d(out_features),
                         nn.ELU(),
                         nn.Dropout(p=p_drop))


class GNN(nn.Module):
    def __init__(self, in_features, out_features, dr=0, act=F.relu,heads=1):
        super().__init__()
        self.conv1 = GATConv(in_features, out_features,heads)
        self.act = act
        self.dr = dr
    #-
    def forward(self, x, edge_index):
        out = self.conv1(x, edge_index)
        out = self.act(out)
        out = F.dropout(out, self.dr, self.training)
        return out


    
class GNN_VAE_Model(nn.Module):
    def __init__(self, input_dim,params,num_classes=1):
        super(GNN_VAE_Model, self).__init__()
        self.var = params.var
        self.num_classes = num_classes
            
        # Encoder
        self.encoder = nn.Sequential()
        self.encoder.add_module('encoder_L1', full_block(input_dim, params.feat_hidden1, params.p_drop))
        self.encoder.add_module('encoder_L2', full_block(params.feat_hidden1, params.feat_hidden2, params.p_drop))
        
        # GNN (GAT)
        self.gn1 = GNN(params.feat_hidden2, params.gcn_hidden1, params.p_drop, act=F.relu,heads = params.nheads)
        self.gn2 = GNN(params.gcn_hidden1*params.nheads, params.gcn_hidden2, params.p_drop, act=lambda x: x)
        self.gn3 = GNN(params.gcn_hidden1*params.nheads, params.gcn_hidden2, params.p_drop, act=lambda x: x) 
        
        # Decoder
        self.decoder = nn.Sequential()
        self.decoder.add_module('decoder_L1', full_block(params.gcn_hidden2, params.feat_hidden2, params.p_drop))
        self.decoder.add_module('decoder_L2', full_block(params.feat_hidden2, params.feat_hidden1, params.p_drop))
        self.decoder.add_module('decoder_output', nn.Sequential(nn.Linear(params.feat_hidden1, input_dim)))
        
        # Cluster
        self.cluster = nn.Sequential()
        self.cluster.add_module('cluster_L1', full_block(params.gcn_hidden2, params.feat_hidden2, params.p_drop))
        self.cluster.add_module('cluster_output', nn.Linear(params.feat_hidden2, self.num_classes))
           
    def encode(self, x, adj):
        feat_x = self.encoder(x)
        hidden1 = self.gn1(feat_x, adj)
        mu = self.gn2(hidden1, adj)
        if self.var:
            logvar = self.gn3(hidden1, adj)
            return mu, logvar
        else:
            return mu, None
     
    def reparameterize(self, mu, logvar):
        if self.training and logvar is not None:
            std = torch.exp(logvar)
            eps = torch.randn_like(std)
            return eps.mul(std).add_(mu)
        else:
            return mu
    
    def forward(self, x, adj):
        mu, logvar = self.encode(x, adj)
        gnn_z = self.reparameterize(mu, logvar)
        x_reconstructed = self.decoder(gnn_z)
        pred_label = F.softmax(self.cluster(gnn_z),dim=1)
        return pred_label, x_reconstructed, gnn_z, mu, logvar