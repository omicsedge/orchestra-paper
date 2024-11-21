import torch
import torch.nn as nn


def PADDING_REQ(KERNEL_SIZE):
    assert KERNEL_SIZE % 2 == 1, "Use kernels of odd number length"
    return int((KERNEL_SIZE - 1) // 2)


def get_CONV1d_layers(param_tuples):
    # =============================================================================
    #     for ind,(inp,_,_) in enumerate(param_tuples[1:]):
    #         assert param_tuples[ind][1] == inp/reduction,\
    #             f"In channels and outchannel dont match in {ind+1} [{param_tuples[ind]} -> {param_tuples[ind+1]}]"
    #
    # =============================================================================
    assert all([k % 2 == 1 for _, _, k in param_tuples]), "Use only odd kernel sizes"

    layers = [
        nn.Sequential(
            nn.Conv1d(
                in_channels=inp,
                out_channels=out,
                kernel_size=kernel,
                padding=PADDING_REQ(kernel),
            ),
            nn.BatchNorm1d(out),
            nn.ReLU(),
        )
        for inp, out, kernel in param_tuples
    ]

    return layers


class MHScaledDotProductSoftmaxAttention(nn.Module):
    def __init__(self, device, nheads):
        super(MHScaledDotProductSoftmaxAttention, self).__init__()

        self.bias = torch.rand((1)).to(device)
        self.bias = self.bias.requires_grad_()
        self.activation = nn.ReLU()
        self.nheads = nheads

    def forward(self, X):

        ns, npop, nw = X.shape

        X1 = X.clone()
        X2 = X.clone()
        X3 = X.clone()

        X1 = torch.reshape(X1, (X1.shape[0], self.nheads, -1, X1.shape[2]))

        X2 = torch.reshape(X2, (X2.shape[0], self.nheads, -1, X2.shape[2]))

        X3 = torch.reshape(X3, (X3.shape[0], self.nheads, -1, X3.shape[2]))

        # calculate batch attention
        batch_attention = (
            torch.einsum("nhpl,nhpw->nhlw", X1, X2)
            / torch.sqrt(torch.tensor([npop / self.nheads])).item()
            + self.bias
        )
        batch_attention = self.activation(batch_attention)
        batch_attention = torch.softmax(batch_attention, dim=2)

        # weigthed values
        X_transform = torch.einsum("nhab,nhcb->nhca", batch_attention, X3)
        X_transform = X_transform.reshape(
            X_transform.shape[0], -1, X_transform.shape[3]
        )

        return torch.cat((X, X_transform), 1)


def force_cudnn_initialization():
    s = 32
    dev = torch.device("cuda")
    torch.nn.functional.conv2d(
        torch.zeros(s, s, s, s, device=dev), torch.zeros(s, s, s, s, device=dev)
    )


class MHAttentionModel(nn.Module):
    """
    only convolutional 1d layers, with depthwise convolutions at the head
    """

    def __init__(self, n_populations, device):
        super(MHAttentionModel, self).__init__()

        self.inp_dimension = n_populations
        self.dropout = nn.Dropout(0.35)
        self.n_populations = n_populations
        self.device = device

        # channels
        self.N1d = 256

        # nheads has to be a divisor of the N1d
        self.nheads = 16

        # declare containers
        self.layers = nn.ModuleList()

        # setup architecture
        self.add_layers()

    def add_layers(self):
        self.param_tuples = [
            (self.inp_dimension, self.N1d, 1),
            (self.N1d, self.N1d, 1),
            (self.N1d, self.N1d, 3),
            (self.N1d * 2, self.N1d, 9),
            (self.N1d * 2, self.n_populations, 19),
        ]

        self.n_tail = len(self.param_tuples) - 1

        for ind, layer in enumerate(get_CONV1d_layers(self.param_tuples)):
            self.layers.append(layer)

            if ind > 1:
                # CNN_receptive_field = self.param_tuples[ind-1][2] if ind >= 1 else 1
                if ind != self.n_tail:
                    self.layers.append(
                        MHScaledDotProductSoftmaxAttention(
                            device=self.device, nheads=self.nheads
                        )
                    )

    def forward(self, x):
        out = x
        for ind, layer in enumerate(self.layers):
            out = layer(out)

            if ind != self.n_tail:
                out = self.dropout(out)

        return out
