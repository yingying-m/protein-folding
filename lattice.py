import networkx as nx
import numpy as np
from matplotlib import pyplot as plt

# Rescale the coordinates and return the new coordinates and max x,y-dimensions
def rescale(coordinate_list):
    max_x = coordinate_list[0].x
    max_y = coordinate_list[0].y
    min_x = coordinate_list[0].x
    min_y = coordinate_list[0].y

    for c in coordinate_list:
        if c.x > max_x: max_x = c.x
        if c.x < min_x: min_x = c.x
        if c.y > max_y: max_y = c.y
        if c.y < min_y: min_y = c.y

    for i in range(len(coordinate_list)):
        coordinate_list[i] = (coordinate_list[i].x-min_x, coordinate_list[i].y-min_y)
    
    return coordinate_list, max_x, max_y

# Create the graph using NetworkX
def create_graph(coordinate_list, dictionary, max_x, max_y):
    M = max_x + 1
    N = max_y + 1
    size = max(N,M)
    
    G = nx.grid_2d_graph(N, M)

    pos = {(i,j): np.array([i,j]) for i in range(N) for j in range(M)}
    elist = []

    for i in range(len(coordinate_list)):
        if i < len(coordinate_list) - 1:
            elist.append([coordinate_list[i], coordinate_list[i + 1]])

    for i in range(N):
        for j in range(M):
            if (i,j) not in coordinate_list: G.remove_node((i,j))

    node_colour = []
    for node in pos:
        for key in dictionary:
            if node == key:
                if dictionary[key] == 'H':
                    node_colour.append('k')
                else:
                    node_colour.append('white')

    nodes = nx.draw_networkx_nodes(G, pos, node_color=node_colour, node_size=300/(0.05*size))
    nodes.set_edgecolor('k')

    nx.draw_networkx_edges(G, pos, elist, edge_color='k', width=1.0)
    plt.gca().invert_yaxis()
    plt.box(False)
    plt.show()

# Visualise a conformation on a lattice
def show_conformation(seq, coordinate_list):
    coordinate_list, max_x, max_y = rescale(coordinate_list)

    seq_array = [char for char in seq]
    dictionary = {coordinate_list[i]: seq_array[i] for i in range(len(coordinate_list))}

    create_graph(coordinate_list, dictionary, max_x, max_y)