import logging
import time

import torch
import torch.nn.functional as F
from tqdm import tqdm

from gsMap.GNN.model import GATModel

logger = logging.getLogger(__name__)


def reconstruction_loss(decoded, x):
    """Compute the mean squared error loss."""
    return F.mse_loss(decoded, x)


def label_loss(pred_label, true_label):
    """Compute the cross-entropy loss."""
    return F.cross_entropy(pred_label, true_label)


class ModelTrainer:
    def __init__(self, node_x, graph_dict, params, label=None):
        """Initialize the ModelTrainer with data and hyperparameters."""
        self.device = torch.device('cuda' if torch.cuda.is_available() else 'cpu')
        self.params = params
        self.epochs = params.epochs
        self.node_x = torch.FloatTensor(node_x).to(self.device)
        self.adj_norm = graph_dict["adj_norm"].to(self.device).coalesce()
        self.label = label
        self.num_classes = 1

        if self.label is not None:
            self.label = torch.tensor(self.label).to(self.device)
            self.num_classes = len(torch.unique(self.label))

        # Set up the model
        self.model = GATModel(self.params.feat_cell, self.params, self.num_classes).to(self.device)
        self.optimizer = torch.optim.Adam(
            self.model.parameters(),
            lr=self.params.gat_lr,
            weight_decay=self.params.gcn_decay
        )

    def run_train(self):
        """Train the model."""
        self.model.train()
        prev_loss = float('inf')
        logger.info('Start training...')
        pbar = tqdm(range(self.epochs), desc='GAT-AE model train:', total=self.epochs)
        for epoch in range(self.epochs):
            start_time = time.time()
            self.optimizer.zero_grad()
            pred_label, de_feat, latent_z, mu, logvar = self.model(self.node_x, self.adj_norm)
            loss_rec = reconstruction_loss(de_feat, self.node_x)

            if self.label is not None:
                loss_pre = label_loss(pred_label, self.label)
                loss = self.params.rec_w * loss_rec + self.params.label_w * loss_pre
            else:
                loss = loss_rec

            loss.backward()
            self.optimizer.step()

            batch_time = time.time() - start_time
            left_time = batch_time * (self.epochs - epoch - 1) / 60  # in minutes

            pbar.set_postfix({'Left time': f'{left_time:.2f} mins', 'Loss': f'{loss.item():.4f}'})
            pbar.update(1)

            if abs(loss.item() - prev_loss) <= self.params.convergence_threshold and epoch >= 200:
                pbar.close()
                logger.info('Convergence reached. Training stopped.')
                break
            prev_loss = loss.item()
        else:
            pbar.close()
            logger.info('Max epochs reached. Training stopped.')


    def get_latent(self):
        """Retrieve the latent representation from the model."""
        self.model.eval()
        with torch.no_grad():
            _, _, latent_z, _, _ = self.model(self.node_x, self.adj_norm)
        return latent_z.cpu().numpy()
