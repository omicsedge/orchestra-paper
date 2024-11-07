import logging
import os
import shutil
from copy import deepcopy
from pathlib import Path
from timeit import default_timer

import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import torch.nn.functional as F
import torch.optim as optim
from torch.utils.data import DataLoader

from .dataset import LaiDataset
from .multiheadmodel import MHAttentionModel

logging.basicConfig(level=logging.DEBUG)


class FocalLoss(nn.CrossEntropyLoss):
    """Focal loss for classification tasks on imbalanced datasets"""

    def __init__(self, gamma, alpha=None, ignore_index=-100, reduction="none"):
        super().__init__(weight=alpha, ignore_index=ignore_index, reduction="none")
        self.reduction = reduction
        self.gamma = gamma

    def forward(self, input_, target):
        cross_entropy = super().forward(input_, target)
        # Temporarily mask out ignore index to '0' for valid gather-indices input.
        # This won't contribute final loss as the cross_entropy contribution
        # for these would be zero.
        target = target * (target != self.ignore_index).long()
        input_prob = torch.gather(F.softmax(input_, 1), 1, target.unsqueeze(1))
        loss = torch.pow(1 - input_prob, self.gamma) * cross_entropy
        if self.reduction == "mean":
            return torch.mean(loss)
        elif self.reduction == "sum":
            return torch.sum(loss)
        else:
            return loss


def load_data(batch_size, training_dir, validation_dir, device):
    # datasets
    train_dataset = LaiDataset(training_dir, device=device)
    val_dataset = LaiDataset(validation_dir, device=device)

    # data loaders
    train_loader = DataLoader(
        train_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=0,
    )
    val_loader = DataLoader(
        val_dataset,
        batch_size=batch_size,
        shuffle=True,
        num_workers=0,
    )

    return train_loader, val_loader


def get_populatin_weights(train_loader):
    # get flatten label matrix
    all_labels = np.array(train_loader.dataset.y.flatten().cpu())

    # inverse count weightage
    all_labels_counts = 1 / pd.Series(all_labels).value_counts()
    all_labels_counts = all_labels_counts / all_labels_counts.sum()

    # reorder to align with torch classes
    all_labels_counts = all_labels_counts[
        np.array(train_loader.dataset.y.unique().cpu())
    ]
    all_labels_counts = torch.Tensor(all_labels_counts.values).float()

    return all_labels_counts


def train_batch(network, optimizer, criterion, batch):
    # zero the parameter gradients
    optimizer.zero_grad()

    # get data
    data, labels = batch["inputs"], batch["labels"]

    # forward
    preds = network(data)

    # backward
    loss = criterion(preds, labels)
    loss.backward()
    optimizer.step()

    return loss


def evaluate_populations(loader, network):
    test_pred_ = []
    test_label_ = []
    for batch in loader:
        test_pred_.append(
            torch.argmax(network(batch["inputs"]), 1).cpu(),
        )
        test_label_.append(batch["labels"].cpu())

    test_pred_ = torch.cat(test_pred_, 0).numpy()
    test_label_ = torch.cat(test_label_, 0).numpy()

    pop_accuracies = []
    for i in range(9):
        idx = test_label_ == i
        full_test_acc = (test_pred_[idx] == test_label_[idx]).mean()
        pop_accuracies.append(f"{full_test_acc:0.3f}")

    return pop_accuracies


def test_model(network, criterion, loader):
    """Test the model with test data."""
    network.eval()
    with torch.no_grad():
        loss, n_batches = 0, 0
        true_preds, count = 0, 0
        pop_true_preds = None
        pop_count = None

        for batch in loader:
            data, labels = batch["inputs"], batch["labels"]
            pred_proba = network(data)
            preds = pred_proba.argmax(dim=1)

            loss += criterion(pred_proba, labels)
            n_batches += 1

            true_check = preds == labels
            true_preds += true_check.sum().item()
            count += labels.numel()

            n_populations = data.shape[1]
            if pop_true_preds is None:
                pop_true_preds = [0] * n_populations
                pop_count = [0] * n_populations

            for i in range(n_populations):
                idx = labels == i
                pop_true_preds[i] += true_check[idx].sum().item()
                pop_count[i] += idx.sum().item()

        acc = true_preds / count
        loss = loss / n_batches

    return loss, acc, [tp / c for tp, c in zip(pop_true_preds, pop_count)]


def train_model(
    output_dir,
    training_dir,
    validation_dir,
    population_map_file,
    hparams,
    device,
    level,
):
    # get loaders
    train_loader, val_loader = load_data(
        hparams["batch_size"], training_dir, validation_dir, device
    )

    # get population size
    if level == 3:
        id_to_pop = (
            pd.read_csv(population_map_file, sep="\t")
            .set_index("id")["level_3_code"]
            .to_dict()
        )
    else:
        df = pd.read_csv(population_map_file, sep="\t")
        id_to_pop = dict(
            enumerate(df[f"level_{level}_name"].drop_duplicates().to_list())
        )

    n_populations = len(id_to_pop)
    logging.info(f"Level: {level}")
    logging.info(f"Number of populations: {n_populations}")

    # network
    network = MHAttentionModel(
        n_populations=n_populations,
        device=device,
    ).to(device)
    logging.info(network)

    # optimizer
    optimizer = optim.SGD(
        network.parameters(),
        lr=hparams["lr"],
        momentum=hparams["momentum"],
        weight_decay=hparams["weight_decay"],
    )
    scheduler = optim.lr_scheduler.StepLR(
        optimizer,
        step_size=50,
        gamma=0.8,
        last_epoch=-1,
        verbose=False,
    )
    # optimizer = optim.Adam(network.parameters(), lr=0.05, weight_decay=0.01)

    # loss function
    criterion = FocalLoss(
        gamma=5,
        reduction="mean",
    )

    # train parameters
    best_test_acc = -np.inf
    # best_train_acc = -np.inf
    best_model = deepcopy(network)
    start_timer = default_timer()
    all_pop_accuracies = []

    for epoch in range(hparams["n_epochs"]):  # loop over the dataset multiple times
        network.train()
        train_loss = 0.0

        # train
        n_batches = 0
        # true_preds = 0
        # count = 0
        for batch in train_loader:
            loss = train_batch(network, optimizer, criterion, batch)
            train_loss += loss.item()
            n_batches += 1

        # decay learning rate
        scheduler.step()

        # epoch metrics
        train_loss = train_loss / n_batches
        # train_acc = true_preds/count
        test_loss, test_acc, pop_accuracies = test_model(network, criterion, val_loader)
        all_pop_accuracies.append([epoch] + pop_accuracies)
        weighted_acc = np.mean(pop_accuracies)
        # train_loss, train_acc = test_model(network, criterion, train_loader)

        time_taken = int(default_timer() - start_timer)
        start_timer = default_timer()

        if best_test_acc < test_acc:
            best_test_acc = test_acc
            # best_train_acc = train_acc
            best_epoch = epoch
            best_model = deepcopy(network)

        best_msg = f"Best epoch: #{best_epoch} best_test_acc: {best_test_acc:0.3f} {'[NEW BEST]' if epoch == best_epoch else ''}"
        msg = (
            f"[{time_taken}s]  Epoch: #{epoch} loss: {train_loss:.3f} "
            f"test_loss: {test_loss:.3f}\t"
            f"test_acc: {test_acc:.3f} weighted_acc: {weighted_acc:.3f}"
            f"\t{best_msg}"
        )

        logging.info(msg)

        if (epoch - best_epoch) >= hparams["early_stop"]:
            logging.info("Training terminated by early stopping")
            break

            # track population accuracy
    df_pop_accuracy = pd.DataFrame(
        all_pop_accuracies,
        columns=["epoch"] + [id_to_pop[i] for i in range(n_populations)],
    )
    df_pop_accuracy.to_csv(
        os.path.join(output_dir, "pop_accuracies.tsv.gz"), sep="\t", index=False
    )

    if epoch == hparams["n_epochs"] - 1:
        logging.info("Training terminated - finished max epoch")

    # remove memory load
    network = network.cpu()

    # save model
    logging.info(f"Best epoch: #{best_epoch} test_acc: {best_test_acc:0.3f}")
    model_path = output_dir / "model.pt"
    logging.info(f"Model path: {model_path}")
    model_scripted = torch.jit.script(best_model)
    model_scripted.save(model_path)


def train_smooth_layer(
    training_dir: Path,
    validation_dir: Path,
    pop_map: Path,
    output_dir: Path,
    level: int = 3,
    epochs: int = 1000,
    label_smoothing: float = 0.1,
    early_stop: int = 100,
    learning_rate: float = 0.0005,
    batch_size: int = 256,
):
    output_dir.mkdir(parents=True, exist_ok=True)

    hyperparameters = {
        "n_epochs": epochs,
        "label_smoothing": label_smoothing,
        "lr": learning_rate,
        "batch_size": batch_size,
        "early_stop": early_stop,
        "momentum": 0.9,
        "weight_decay": 0.005,
    }

    # copy parameters file to model
    shutil.copyfile(training_dir / "parameters.json", output_dir / "parameters.json")

    # train model
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    train_model(
        output_dir,
        training_dir,
        validation_dir,
        pop_map,
        hyperparameters,
        device,
        level,
    )
