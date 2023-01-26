import random
import numpy as np
import matplotlib.pyplot as plt
from scipy.stats import rv_continuous
from scipy.optimize import curve_fit
import math
import time
import warnings
import sys
import os
import os.path
import multiprocessing as mp
from itertools import accumulate as _accumulate, repeat as _repeat
from bisect import bisect as _bisect
dir = os.getcwd()


def choices(population, weights=None, *, cum_weights=None, k=1):
    """Return a k sized list of population elements chosen with replacement.
    If the relative weights or cumulative weights are not specified,
    the selections are made with equal probability.
    """
    n = len(population)
    if cum_weights is None:
        if weights is None:
            _int = int
            n += 0.0    # convert to float for a small speed improvement
            return [population[_int(random.random() * n)] for i in _repeat(None, k)]
        cum_weights = list(_accumulate(weights))
    elif weights is not None:
        raise TypeError('Cannot specify both weights and cumulative weights')
    if len(cum_weights) != n:
        raise ValueError('The number of weights does not match the population')
    bisect = _bisect
    total = cum_weights[-1] + 0.0   # convert to float
    hi = n - 1
    return [population[bisect(cum_weights, random.random() * total, 0, hi)]
            for i in _repeat(None, k)]


'''
N = int(sys.argv[1])
amostras = int(sys.argv[2])
d = int(sys.argv[3])
aA = float(sys.argv[4])
aG = float(sys.argv[5])
w0 = float(sys.argv[6])
eta = float(sys.argv[7])
'''
N = 1000
amostras = 10000
d = 3
aA = 3.00
aG = 1.00
w0 = 1.00
eta = 2.00
a = [0.2, 0.3, 0.4, 0.2, 0.3, 0.4]
a_interaction = [1, 1, 1, 3, 3, 3]

class Geo_Node:

    '''
    For some algorithms, we will need networks that take into account the euclidean
    positions and distances between the nodes, which is not (to my knowledge)
    supported by networkx, so this class and the Geo_Network class will be created
    in order to create such networks. The present class will be used for each
    individual node.
    '''

    def __init__(self, position, label, weight=1.0):
        '''
        The position has to be any 2 valued (x,y) set (a tuple, a list, whatever)
        containing the x and y component of the position of the node. Maybe
        in the future I can work with z components. The weight of the node
        is set to 1 (in case your network in unweighted), otherwise you can
        set a weight yourself.
        '''

        self.label=label #every node has a label to distinguish it from the others,
        #it can be a number for example.
        self.position= np.array(position)
        self.energy=weight
        self.connections=[]

    def update_energy(self, edge_dict):
        #the edge_dict has to be the dictionary of edges available with the
        #Geo_Network class.

        new_energy=0
        for edge in edge_dict:
            if self in edge:
                new_energy+=edge_dict[edge]/2
        self.energy = new_energy

class Geo_Node_1d:
    def __init__(self, position, label, weight=1.0):
        self.label=label
        self.x = position
        self.position= (self.x)
        self.energy=weight
        self.connections=[]

    def update_energy(self, edge_dict):
        new_energy=0

        for edge in edge_dict:
            if self in edge:
                new_energy+=edge_dict[edge]/2

        self.energy = new_energy

class Geo_Node_2d:
    def __init__(self, position, label, weight=1.0):
        self.label = label
        self.x = position[0]
        self.y = position[1]
        self.position = (self.x, self.y)
        self.energy = weight
        self.connections=[]

    def update_energy(self, edge_dict):
        new_energy=0

        for edge in edge_dict:
            if self in edge:
                new_energy+=edge_dict[edge]/2
        self.energy = new_energy

class Geo_Node_3d:
    def __init__(self, position, label, weight=1.0):
        self.label=label
        self.x = position[0]
        self.y = position[1]
        self.z = position[2]
        self.position = (self.x, self.y, self.z)
        self.energy = weight
        self.connections=[]
        self.energy_with_correction = weight
        self.total_energy = weight

    def update_energy(self, edge_dict):
        new_energy=0

        for edge in edge_dict:
            if self in edge:
                new_energy+=edge_dict[edge]/2
        self.energy = new_energy

class Geo_Network:
    '''
    Read the Geo_Nodes class description first. This class will store the network
    itself of that kind of network, and will have functions that concern the whole
    network.
    '''

    def __init__(self, edges_list=None):
        '''
        The network starts completely empty if you don't initiate it with any
        arguments, and then you can build the network yourself. If, on the other
        hand, you have a pre-made network and want to use it, you can do it, the
        edges_list argument receives a list on the format [(A,B,w1), (A,C,w2), (B,D,w3), etc]
        containing all the edges of the network, where A, B, C, D, etc, are nodes
        (Geo_Node format) of your network and (A,B,w) represents an edge between
        nodes A and B with weight w.
        '''
        if edges_list!=None:

            self.edges=[(edge[0], edge[1]) for edge in edges_list]
            self.nodes=[]
            self.edges_weights={} #dict storing the weights of the edges

            for edge in edges_list:
                if len(edge)==3:
                    self.edges_weights[edge]=edge[3]

            for tup in self.edges:

                if tup[0] not in self.nodes:
                    self.nodes.append(tup[0])

                if tup[1] not in self.nodes:
                    self.nodes.append(tup[1])

            positionsum = np.zeros(len(node[0].position))
            weight_sum=0

            for node in self.nodes:
                positionsum += node.energy * node.position
                weight_sum+=node.energy

            self.center_mass= positionsum/weight_sum

        else:
            self.nodes=[]
            self.edges=[]
            self.center_mass=None
            self.edges_weights={}

    def add_node(self, node): #node has to be an object of the type Geo_Node

        self.nodes.append(node)

    def add_edge(self, node1, node2, weight):
        '''
        This function can create an edge between two preexisting nodes, or you
        can create new nodes from it. If you do self.add_edge(A,B), both A and
        B still not in the network, they will automatically be put in the network,
        along with the edge between them. node1 and node2 obviously have to be
        of the class Geo_Node.
        '''
        self.edges_weights[(node1, node2)]=weight

        if (node1, node2) not in self.edges:
            self.edges.append((node1, node2))

        if node1 not in self.nodes:
            self.nodes.append(node1)

        if node2 not in self.nodes:
            self.nodes.append(node2)

    def euclid_distance(self, node1, node2):
        #calculate the euclidean distance between node1 and node2
        subtraction = node1.position - node2.position
        distance= math.sqrt(np.sum(subtraction*subtraction))

        return distance

    def update_center(self):
        #update the center of mass of the network
        positionsum = np.zeros(d)
        weight_sum=0

        for node in self.nodes:
            positionsum += node.energy * node.position
            weight_sum+= node.energy

        self.center_mass= positionsum/weight_sum

class Geo_Network_1d:
    def __init__(self, edges_list=None):
        if edges_list!=None:

            self.edges=[(edge[0], edge[1]) for edge in edges_list]
            self.nodes=[]
            self.edges_weights={} #dict storing the weights of the edges

            for edge in edges_list:
                if len(edge)==3:
                    self.edges_weights[edge]=edge[3]

            for tup in self.edges:

                if tup[0] not in self.nodes:
                    self.nodes.append(tup[0])

                if tup[1] not in self.nodes:
                    self.nodes.append(tup[1])

            xsum=0
            weight_sum=0

            for node in self.nodes:
                xsum+=node.energy*node.x
                weight_sum+=node.energy

            self.center_mass=(xsum/weight_sum)

        else:
            self.nodes=[]
            self.edges=[]
            self.center_mass=None
            self.edges_weights={}

    def add_node(self, node): #node has to be an object of the type Geo_Node
        self.nodes.append(node)

    def add_edge(self, node1, node2, weight):

        self.edges_weights[(node1, node2)]=weight

        if (node1, node2) not in self.edges:
            self.edges.append((node1, node2))

        if node1 not in self.nodes:
            self.nodes.append(node1)

        if node2 not in self.nodes:
            self.nodes.append(node2)

    def euclid_distance(self, node1, node2):
        distance= ((node1.x-node2.x)**2)**0.5
        return distance

    def update_center(self):
        xsum=0
        weight_sum=0

        for node in self.nodes:
            xsum+=node.energy*node.x
            weight_sum+=node.energy

        self.center_mass= (xsum/weight_sum)

class Geo_Network_2d:
    def __init__(self, edges_list=None):

        if edges_list!=None:

            self.edges=[(edge[0], edge[1]) for edge in edges_list]
            self.nodes=[]
            self.edges_weights={} #dict storing the weights of the edges

            for edge in edges_list:
                if len(edge)==3:
                    self.edges_weights[edge]=edge[3]

            for tup in self.edges:

                if tup[0] not in self.nodes:
                    self.nodes.append(tup[0])

                if tup[1] not in self.nodes:
                    self.nodes.append(tup[1])

            xsum=0
            ysum=0
            weight_sum=0

            for node in self.nodes:
                xsum+=node.energy*node.x
                ysum+=node.energy*node.y
                weight_sum+=node.energy

            self.center_mass=(xsum/weight_sum, ysum/weight_sum)

        else:
            self.nodes=[]
            self.edges=[]
            self.center_mass=None
            self.edges_weights={}

    def add_node(self, node): #node has to be an object of the type Geo_Node
        self.nodes.append(node)

    def add_edge(self, node1, node2, weight):

        self.edges_weights[(node1, node2)]=weight

        if (node1, node2) not in self.edges:
            self.edges.append((node1, node2))

        if node1 not in self.nodes:
            self.nodes.append(node1)

        if node2 not in self.nodes:
            self.nodes.append(node2)

    def euclid_distance(self, node1, node2):
        distance= ((node1.x-node2.x)**2+(node1.y-node2.y)**2)**0.5
        return distance

    def update_center(self):
        xsum=0
        ysum=0
        weight_sum=0

        for node in self.nodes:
            xsum+=node.energy*node.x
            ysum+=node.energy*node.y
            weight_sum+=node.energy

        self.center_mass=(xsum/weight_sum, ysum/weight_sum)

class Geo_Network_3d:
    def __init__(self, edges_list=None):

        if edges_list!=None:

            self.edges=[(edge[0], edge[1]) for edge in edges_list]
            self.nodes=[]
            self.edges_weights={} #dict storing the weights of the edges

            for edge in edges_list:
                if len(edge)==3:
                    self.edges_weights[edge]=edge[3]

            for tup in self.edges:

                if tup[0] not in self.nodes:
                    self.nodes.append(tup[0])

                if tup[1] not in self.nodes:
                    self.nodes.append(tup[1])

            xsum=0
            ysum=0
            zsum=0
            weight_sum=0

            for node in self.nodes:
                xsum+=node.energy*node.x
                ysum+=node.energy*node.y
                zsum+=node.energy*node.z
                weight_sum+=node.energy

            self.center_mass=(xsum/weight_sum, ysum/weight_sum, zsum/weight_sum)

        else:
            self.nodes=[]
            self.edges=[]
            self.center_mass=None
            self.edges_weights={}

    def add_node(self, node): #node has to be an object of the type Geo_Node
        self.nodes.append(node)

    def add_edge(self, node1, node2, weight):

        self.edges_weights[(node1, node2)]=weight

        if (node1, node2) not in self.edges:
            self.edges.append((node1, node2))

        if node1 not in self.nodes:
            self.nodes.append(node1)

        if node2 not in self.nodes:
            self.nodes.append(node2)

    def euclid_distance(self, node1, node2):
        distance= ((node1.x-node2.x)**2+(node1.y-node2.y)**2+(node1.z-node2.z)**2)**0.5
        return distance

    def update_center(self):
        xsum=0
        ysum=0
        zsum=0
        weight_sum=0

        for node in self.nodes:
            xsum+=node.energy*node.x
            ysum+=node.energy*node.y
            zsum+=node.energy*node.z
            weight_sum+=node.energy

        self.center_mass=(xsum/weight_sum, ysum/weight_sum, zsum/weight_sum)

class exponential(rv_continuous):
    '''
    This is a VERY confusing part of the code. This is supposed to create a custom
    probability distribution function, equation 3 in the article. This is what is
    (apparently) called an abstract base class (ABC), its kind of a base class where
    you can build the "rest" of the class by yourself. In this case, we give this class
    a distribution function, and the class takes care of making it work with every
    other function and class in the module and in python. Here is the link for the class:
    https://docs.scipy.org/doc/scipy/reference/generated/scipy.stats.rv_continuous.html
    '''
    def _pdf(self, x):

        return eta*math.exp(-(x/w0)**eta)/(w0*math.gamma(1/eta))

stretched_exponential=exponential(name='stretched_exponential', a=0)

class alpha_G_prob(rv_continuous):

    def _pdf(self, x, alpha_G, d):

        return (d+alpha_G-1)*(1/x**(d+alpha_G))

alpha_dist=alpha_G_prob(name='alpha_dist', a=1)

def TN_model_generate(alpha_A, alpha_G, N, d, a, a_interaction, index):

    '''
    This function will create a network based on a network model presented by
    Constantino Tsallis in the paper available in the references text. The model
    is still unamed, so I named it for the purposes of this code as Tsallis Network
    model (explaining why TN_model). The parameters alpha are, obviously, the same
    alphas from the article used as reference and N is the number of iterations.
    The dimension, for now, will be always 2 for simplicity.
    '''
    if d == 1:
        nk=Geo_Network_1d()
        node1=Geo_Node_1d((0), 1)
        nk.add_node(node1)
        nk.update_center()

        r=alpha_dist.rvs(alpha_G, 1)
        o=choices((-1,1))

        node2_x=nk.center_mass+r*o[0]

        node2=Geo_Node_1d((node2_x), 2)
        nk.add_node(node2)
        edge12_weight=stretched_exponential.rvs()
        nk.add_edge(node1, node2, edge12_weight)
        node1.connections.append(node2.label)
        node2.connections.append(node1.label)
        node1.update_energy(nk.edges_weights)
        node2.update_energy(nk.edges_weights)
        nk.update_center()

        for i in range(2,N):
            r=alpha_dist.rvs(alpha_G, 1)
            o=choices((-1,1))
            nodei_x=nk.center_mass+r*o[0]
            nodei=Geo_Node_1d((nodei_x),i+1)

            connections_list = choices(nk.nodes, weights=[node.energy/nk.euclid_distance(nodei, node) for node in nk.nodes], k=1)
            nk.add_node(nodei)
            for connection in connections_list:
                edgeij_weight=stretched_exponential.rvs()
                nk.add_edge(nodei, connection, edgeij_weight)
                nodei.connections.append(connection.label)
                connection.connections.append(nodei.label)
                nodei.update_energy(nk.edges_weights)
                connection.update_energy(nk.edges_weights)


            nk.update_center()
            #print('Iteration', i+1, 'is complete!')


    elif d == 2:
        nk=Geo_Network_2d()
        node1=Geo_Node_2d((0,0), 1)
        nk.add_node(node1)
        nk.update_center()

        r=alpha_dist.rvs(alpha_G, 2)
        o=random.uniform(0,1)
        phi = 2*math.pi*o
        node2_x=nk.center_mass[0]+r*math.cos(phi)
        node2_y=nk.center_mass[1]+r*math.sin(phi)

        node2=Geo_Node_2d((node2_x, node2_y), 2)
        nk.add_node(node2)
        edge12_weight=stretched_exponential.rvs()

        nk.add_edge(node1, node2, edge12_weight)
        node1.connections.append(node2.label)
        node2.connections.append(node1.label)
        node1.update_energy(nk.edges_weights)
        node2.update_energy(nk.edges_weights)

        nk.update_center()

        for i in range(2,N):
            r=alpha_dist.rvs(alpha_G, 2)
            u=random.uniform(0,1)
            phi=2*math.pi*u
            nodei_x=nk.center_mass[0]+r*math.cos(phi)
            nodei_y=nk.center_mass[1]+r*math.sin(phi)
            nodei=Geo_Node_2d((nodei_x, nodei_y),i+1)

            connections_list = choices(nk.nodes, weights=[node.energy/nk.euclid_distance(nodei, node) for node in nk.nodes], k=1)
            nk.add_node(nodei)
            for connection in connections_list:
                edgeij_weight=stretched_exponential.rvs()
                nk.add_edge(nodei, connection, edgeij_weight)
                nodei.connections.append(connection.label)
                connection.connections.append(nodei.label)
                nodei.update_energy(nk.edges_weights)
                connection.update_energy(nk.edges_weights)

            nk.update_center()
            #print('Iteration', i+1, 'is complete!')

        # Correção da curva

    elif d == 3:
        nk=Geo_Network_3d()
        node1=Geo_Node_3d((0,0,0), 1)
        nk.add_node(node1)
        nk.update_center()

        r=alpha_dist.rvs(alpha_G, 3)
        u=random.uniform(0,1)
        v=random.uniform(0,1)
        theta= 2*math.pi*u
        phi=np.arccos(1-2*v)
        node2_x=nk.center_mass[0]+r*math.cos(theta)*math.sin(phi)
        node2_y=nk.center_mass[1]+r*math.sin(theta)*math.sin(phi)
        node2_z=nk.center_mass[2]+r*math.cos(phi)

        node2=Geo_Node_3d((node2_x, node2_y, node2_z), 2)
        nk.add_node(node2)
        edge12_weight=stretched_exponential.rvs()

        nk.add_edge(node1, node2, edge12_weight)
        node1.connections.append(node2.label)
        node2.connections.append(node1.label)
        node1.update_energy(nk.edges_weights)
        node2.update_energy(nk.edges_weights)

        nk.update_center()

        for i in range(2,N):
            r=alpha_dist.rvs(alpha_G, 3)
            u=random.uniform(0,1)
            v=random.uniform(0,1)
            theta= 2*math.pi*u
            phi=np.arccos(1-2*v)
            nodei_x=nk.center_mass[0]+r*math.cos(theta)*math.sin(phi)
            nodei_y=nk.center_mass[1]+r*math.sin(theta)*math.sin(phi)
            nodei_z=nk.center_mass[2]+r*math.cos(phi)
            nodei=Geo_Node_3d((nodei_x, nodei_y, nodei_z),i+1)

            connections_list = choices(nk.nodes, weights=[node.energy/nk.euclid_distance(nodei, node) for node in nk.nodes], k=1)
            nk.add_node(nodei)
            for connection in connections_list:
                edgeij_weight=stretched_exponential.rvs()
                nk.add_edge(nodei, connection, edgeij_weight)
                nodei.connections.append(connection.label)
                connection.connections.append(nodei.label)
                nodei.update_energy(nk.edges_weights)
                connection.update_energy(nk.edges_weights)

            nk.update_center()
            #print('Iteration', i+1, 'is complete!')

    else:
        nk=Geo_Network()
        position_node1 = np.zeros(d)
        node1=Geo_Node(position_node1, 1)
        nk.add_node(node1)
        nk.update_center()
        r=alpha_dist.rvs(alpha_G, d) #take a look at the pareto variate prob
        #distribution to understand why I'm using it, I'm also assuming the module
        #uses xm=1 to be compatible with the r>=1 of the article.

        v = np.random.randn(d)
        v_norm = np.sqrt(np.sum(v*v))
        v_hat = v / v_norm
        position_node2 = (v_hat * r) + nk.center_mass

        node2=Geo_Node(position_node2, 2)
        nk.add_node(node2)
        edge12_weight=stretched_exponential.rvs()

        nk.add_edge(node1, node2, edge12_weight)
        node1.connections.append(node2.label)
        node2.connections.append(node1.label)
        node1.update_energy(nk.edges_weights)
        node2.update_energy(nk.edges_weights)

        nk.update_center()

        for i in range(2,N):
            r=alpha_dist.rvs(alpha_G, d)
            v = np.random.rand(d)
            v_norm = np.sqrt(np.sum(v*v))
            v_hat = v / v_norm
            position_nodei = nk.center_mass + (r*v_hat)
            nodei=Geo_Node(position_nodei,i+1)
            connections_list=choices(nk.nodes, weights=[(node.energy/nk.euclid_distance(nodei, node)) for node in nk.nodes], k=1)
            nk.add_node(nodei)

            for connection in connections_list:
                edgeij_weight=stretched_exponential.rvs()
                nk.add_edge(nodei, connection, edgeij_weight)
                nodei.connections.append(connection.label)
                connection.connections.append(nodei.label)
                nodei.update_energy(nk.edges_weights)
                connection.update_energy(nk.edges_weights)
            nk.update_center()
            #print('Iteration', i+1, 'is complete!')

    
    #Energy correction (not working)
    correction_array = []
    for i in range (N):
        a
        N_star = 0
        energy_product_sum = 0
        for connection in nk.nodes[i].connections:
            N_star += nk.euclid_distance(nk.nodes[i], nk.nodes[connection-1])**(-a_interaction)
            energy_product_sum += nk.nodes[i].energy *  nk.nodes[connection-1].energy
        correction = (a/N_star)*energy_product_sum
        correction_array.append(correction)

    for i in range (N):
        nk.nodes[i].energy += correction_array[i]

    file = 'N_%d_d_%d_aA_%.2f_aG_%.2f_w0_%.2f_eta_%.2f_a_%.3f_a_interaction_%.3f_amostra_%d.txt' %(N, d, aA, aG, w0, eta, a, a_interaction, index+1)
    path = 'N_%d_d_%d_aA_%.2f_aG_%.2f_w0_%.2f_eta_%.2f_a_%.6f_a_interaction_%.3f' %(N, d, aA, aG, w0, eta, a, a_interaction)
    name = os.path.join(path, file)

    f = open(name, "w")
    energy_list = [node.energy for node in nk.nodes]
    label_list = [node.label for node in nk.nodes]
    for i in range (N):
        f.write(str(energy_list[i]) + '\t' + str(label_list[i]) + '\n')
    f.close()
    return nk

def nk_generation(n, alpha_A, alpha_G, N, d, a, a_interaction):
    '''
    This will create n networks with the TN_model, using the parameters shown,
    and print the graph of p(e) x e where p(e) is the probability to find a node
    with energy e in the final network.
    '''
    path = 'N_%d_d_%d_aA_%.2f_aG_%.2f_w0_%.2f_eta_%.2f_a_%.6f_a_interaction_%.3f' %(N, d, aA, aG, w0, eta, a, a_interaction)
    if not os.path.exists(path):
        os.mkdir(path)

    energies=[]

    for i in range(n):
        nk=TN_model_generate(alpha_A, alpha_G, N, d, a, i)
        energies+=[node.energy for node in nk.nodes]
        print('Simulation', i+1, 'of', n, 'is complete!')

    return energies

def nk_generation_multi(n, alpha_A, alpha_G, N, d, a, a_interaction):
    '''
    This will create n networks with the TN_model, using the parameters shown,
    and print the graph of p(e) x e where p(e) is the probability to find a node
    with energy e in the final network.
    '''

    path = 'N_%d_d_%d_aA_%.2f_aG_%.2f_w0_%.2f_eta_%.2f_a_%.6f_a_interaction_%.3f' %(N, d, aA, aG, w0, eta, a, a_interaction)
    if not os.path.exists(path):
        os.mkdir(path)

    cores=mp.cpu_count()
    pool=mp.Pool(processes=cores)
    energies=[]
    pool_list=[]

    for i in range(n):
        p=pool.apply_async(TN_model_generate, (alpha_A, alpha_G, N, d, a, a_interaction, i))
        pool_list.append(p)
        print('Simulation', 100*len(pool_list)/n, '% completed!')

    networks=[p.get() for p in pool_list]

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#

#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#---#

if __name__=='__main__':
    t0=time.time()
    for i in range(len(a)):
        #nk_generation(amostras, aA, aG, N, d, a[i])
        nk_generation_multi(amostras, aA, aG, N, d, a[i], a_interaction[i])
    #analise(N, d, aA, aG, w0, eta, a, 50, q_fit=True)
    #print('Execution time:', time.time()-t0)
