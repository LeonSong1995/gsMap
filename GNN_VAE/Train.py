#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul  4 19:58:58 2023

@author: songliyang
"""
import time
import torch
import torch.nn.functional as F
from progress.bar import Bar
from Model import GNN_VAE_Model


def reconstruction_loss(decoded, x):
    loss_fn = torch.nn.MSELoss()
    loss = loss_fn(decoded, x)
    return loss


def label_loss(pred_label, true_label):
    loss_fn = torch.nn.CrossEntropyLoss()
    loss = loss_fn(pred_label, true_label)
    return loss


class Model_Train:
    def __init__(self, node_X, graph_dict, params, label=None):
        device = 'cuda' if torch.cuda.is_available() else 'cpu'
        torch.cuda.empty_cache()

        self.params = params
        self.device = device
        self.epochs = params.epochs
        self.node_X = torch.FloatTensor(node_X.copy()).to(device)
        self.adj_norm = graph_dict["adj_norm"].to(device).coalesce()
        self.label = label
        self.num_classes = 1
        
        if not self.label is None:
            self.label = torch.tensor(self.label).to(self.device)
            self.num_classes = len(self.label.unique())
        
        # Set Model 
        self.model = GNN_VAE_Model(self.params.feat_cell,self.params,self.num_classes).to(device)
        self.optimizer = torch.optim.Adam(params = list(self.model.parameters()),
                                          lr = self.params.gcn_lr, weight_decay = self.params.gcn_decay)               
    
    # Train    
    def run_train(self):
        self.model.train()
        prev_loss = float('inf')
        
        bar = Bar('GAT-AE model train:', max = self.epochs)
        bar.check_tty = False 
        for epoch in range(self.epochs):
            start_time = time.time()
            self.model.train()
            self.optimizer.zero_grad()
            pred_label, de_feat, latent_z, mu, logvar = self.model(self.node_X, self.adj_norm)
            loss_rec = reconstruction_loss(de_feat, self.node_X)
            
            # Check whether annotation was provided
            if not self.label is None:
                loss_pre = label_loss(pred_label, self.label)
                loss = (self.params.rec_w * loss_rec) + (self.params.label_w * loss_pre)
            else:
                loss = loss_rec
                
            loss.backward()
            self.optimizer.step()
            
            # Update process
            end_time = time.time()
            batch_time = end_time - start_time
            
            
            bar_str = '{} / {} | Left time: {batch_time:.2f} mins| Loss: {loss:.4f}'
            bar.suffix = bar_str.format(epoch + 1,self.epochs,
                                        batch_time = batch_time * (self.epochs - epoch) / 60, loss=loss.item())
            bar.next()
            
            # Check convergence
            if abs(loss.item() - prev_loss) <= self.params.convergence_threshold and epoch >= 200:
                print('\nConvergence reached. Training stopped.')
                break

            prev_loss = loss.item()
            
        bar.finish()
    #-    
    def get_latent(self):
        self.model.eval()
        pred, de_fea, latent_z, mu, logvar = self.model(self.node_X, self.adj_norm)
        latent_z = latent_z.data.cpu().numpy()
        return latent_z
