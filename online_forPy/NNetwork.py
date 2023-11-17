import torch
import torch.nn as nn
import numpy as np
from collections import OrderedDict

def DAML(XA,XB,dt_slow):

    class CNN(nn.Module):
        """
        Fully CNN with 4 convolutional layers
        The input 'args' should be a dictionary containing
        details of the network hyperparameters and architecture
        """

        def __init__(self,args):
            super(CNN, self).__init__()
            torch.manual_seed(args['seed'])

            self.conv_net = nn.Sequential(OrderedDict([
                ('C1', nn.Conv2d(args['n_channels'], args['h_channels'][0], kernel_size=args['kernel_size'],\
                                 padding=args['zero_padding'],stride=args['stride'], bias=args['bias'])),
                ('Relu1', nn.ReLU()),
                ('C2', nn.Conv2d(args['h_channels'][0], args['h_channels'][1], kernel_size=args['kernel_size'],\
                                 padding=args['zero_padding'],stride=args['stride'],bias=args['bias'])),
                ('Relu2', nn.ReLU()),
                ('C3', nn.Conv2d(args['h_channels'][1], args['h_channels'][2], kernel_size=args['kernel_size'],\
                                 padding=args['zero_padding'],stride=args['stride'],bias=args['bias'])),
                ('Relu3', nn.ReLU()),
                ('C4', nn.Conv2d(args['h_channels'][2], args['n_classes'], kernel_size=args['kernel_size'],\
                                 padding=args['zero_padding'],stride=args['stride'],bias=args['bias'])),
            ]))

        def forward(self, x):
            return self.conv_net(x)

    def Net(x,args,path=None):
        """
        Function to train CNN and/or perform inference.
        Can either provide a 'path' to the pre-computed network weights,
        or leave as None to perform optimization from scratch (the latter
        will be very slow). If path is None, then the optimization will be
        performed using training data 'x_train' and 'y_train'. Inference
        will be done on 'x_valid' and 'y_valid' in either case.
        The input 'args' should be a dictionary containing
        details of the network hyperparameters and architecture
        """

        torch.manual_seed(args['seed'])

        x = torch.from_numpy(x.astype(np.float32))
        args['n_channels'] = x.shape[1]

        model = CNN(args)
        model.load_state_dict(torch.load(path,map_location=torch.device('cpu')))
        model.eval()
        with torch.no_grad():
            predictions = model(x)
        return predictions.numpy()

    argsA = {
    'kernel_size':3,
    'zero_padding':0,
    'h_channels':[32,64,128],
    'n_classes':1,
    'stride':1,
    'bias':False,
    'seed':711,
    }

    argsB = {
    'kernel_size':1,
    'zero_padding':0,
    'h_channels':[32,64,128],
    'n_classes':5,
    'stride':1,
    'bias':False,
    'seed':711,
    }

    NetworkA_weights = '/lustre/f2/dev/William.Gregory/CNNForpy/NetworkA_weights_1982-2017_notend_P2.pt'
    NetworkB_weights = '/lustre/f2/dev/William.Gregory/CNNForpy/NetworkB_weights_1982-2017_notend_P2.pt'

    NetworkA_stats = np.load('/lustre/f2/dev/William.Gregory/CNNForpy/NetworkA_statistics_1982-2017_allsamples_notend.npz')
    NetworkB_stats = np.load('/lustre/f2/dev/William.Gregory/CNNForpy/NetworkB_statistics_1982-2017_allsamples_notend.npz')

    sec5day = 432000
    scale = sec5day/dt_slow
    halo = 4

    XA = np.transpose(XA[None],(0,1,3,2)) #(1,9,nj+halo,ni+halo)
    XB = np.transpose(XB[None,:,halo:-halo,halo:-halo],(0,1,3,2)) #(1,6,nj,ni)
    nj,ni = XA.shape[2],XA.shape[3]
    mj,mi = XB.shape[2],XB.shape[3]
    
    for N in range(XA.shape[1]-1):
        XA[:,N] = (XA[:,N]-NetworkA_stats['mu'][N])/NetworkA_stats['sigma'][N]
        XA[:,N][XA[:,-1]==0] = 0
        
    XA[np.isnan(XA)] = 0
    argsA['n_channels'] = XA.shape[1]
    dSIC = Net(XA,argsA,path=NetworkA_weights)
    dSIC[XB[:,-1:]==0] = 0

    XB = np.hstack((dSIC,XB))
    for N in range(XB.shape[1]-1):
        XB[:,N] = (XB[:,N]-NetworkB_stats['mu'][N])/NetworkB_stats['sigma'][N]
        XB[:,N][XB[:,-1]==0] = 0

    XB[np.isnan(XB)] = 0
    argsB['n_channels'] = XB.shape[1]
    dSICN = Net(XB,argsB,path=NetworkB_weights)
    for k in range(dSICN.shape[1]):
        dSICN[:,k][XB[:,-1]==0] = 0
    
    dSICN = np.transpose(dSICN[0]/scale,(2,1,0)) #(ni,nj,5)
    dSICN = np.vstack((np.vstack((np.zeros((halo,mj,argsB['n_classes'])),dSICN)),np.zeros((halo,mj,argsB['n_classes'])))) #(ni+halo,nj,5)
    dSICN = np.hstack((np.hstack((np.zeros((ni,halo,argsB['n_classes'])),dSICN)),np.zeros((ni,halo,argsB['n_classes'])))) #(ni+halo,nj+halo,5)
    
    return dSICN
