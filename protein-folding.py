import time
import pandas as pd
import numpy as np
from copy import deepcopy
from queue import PriorityQueue
import lattice

class Point:
    def __init__(self, x, y):
        self.x = x
        self.y = y
    
    # Print the coordinates of some point
    def __repr__(self):
        return "".join(["(", str(self.x), ",", str(self.y), ")"])

    # Find the (un)occupied positions adjacent to some point
    def find_positions(self, grid, relative_position_list, is_occupied):
        position_list = []
        for relative_position in relative_position_list:
            adjacent_position = Point(self.x + relative_position[0], self.y + relative_position[1])
            if adjacent_position.x > (len(grid) - 1) or adjacent_position.x < 0 or adjacent_position.y > (len(grid[len(grid)-1]) - 1) or adjacent_position.y < 0:
                continue
            if is_occupied:
                if grid[adjacent_position.y][adjacent_position.x] == ' ': continue
            else:
                if grid[adjacent_position.y][adjacent_position.x] != ' ': continue
            position_list.append(adjacent_position)
        return position_list

class Node:
    def __init__(self, coordinates, grid, g):
        self.coordinates = coordinates
        self.grid = grid
        self.g = g
        self.f = g

    # Define the criteria on which nodes must be sorted
    def __lt__(self, other):
        if self.f < other.f: return True
        elif self.f == other.f:
            if len(self.coordinates) < len(other.coordinates): return True
        return False

    # Returns a list of children for a node
    def find_children(self, seq_array, index):
        current_position = self.coordinates[-1]
        relative_position_list = [(-2, 0), (2, 0), (0, -2), (0, 2)]
        is_occupied = False
        child_positions = current_position.find_positions(self.grid, relative_position_list, is_occupied)
        children = []
        amino_acid = seq_array[index]

        for position in child_positions:
            coordinates_list = deepcopy(self.coordinates)
            coordinates_list.append(position)
            new_grid = deepcopy(self.grid)
            new_grid[position.y][position.x] = amino_acid
            draw_bond(new_grid, current_position, position)
            if len(coordinates_list) < len(seq_array):
                free_positions = position.find_positions(new_grid, relative_position_list, is_occupied)
                if len(free_positions) == 0: continue
            child = Node(coordinates_list, new_grid, self.f)
            children.append(child)
        return children
        
    # Calculate the change in stability after placing a new amino acid on the grid
    def calculate_stability(self):
        stability_change = 0
        current_position = self.coordinates[-1]
        current_val = self.grid[current_position.y][current_position.x]
        if (current_val == 'H' and current_position.y > 0 and self.grid[current_position.y - 2][current_position.x] == 'H' and self.grid[current_position.y - 1][current_position.x] == ' '):
            self.grid[current_position.y - 1][current_position.x] = ':'
            stability_change += 1
        if (current_val == 'H' and current_position.y < rows - 2 and self.grid[current_position.y + 2][current_position.x] == 'H' and self.grid[current_position.y + 1][current_position.x] == ' '):
            self.grid[current_position.y + 1][current_position.x] = ':'
            stability_change += 1
        if (current_val == 'H' and current_position.x > 0 and self.grid[current_position.y][current_position.x - 2] == 'H' and self.grid[current_position.y][current_position.x - 1] == ' '):
            self.grid[current_position.y][current_position.x - 1] = '"'
            stability_change += 1
        if (current_val == 'H' and current_position.x < cols - 2 and self.grid[current_position.y][current_position.x + 2] == 'H' and self.grid[current_position.y][current_position.x + 1] == ' '):
            self.grid[current_position.y][current_position.x + 1] = '"'
            stability_change += 1
        return -stability_change

class Output:
    def __init__(self, node, runtime):
        self.structure = node.grid
        self.stability = node.f
        self.runtime = runtime

# Read amino acid sequence from .csv file
def read_sequence(input_file, id):
    df = pd.read_csv(input_file)
    seq = str(df.iloc[id].sequence)
    return seq

# Print the grid
def display(grid):
    print("\n")
    for rows in range(len(grid)):
        empty_row = grid[rows].count(grid[rows][0]) == len(grid[rows])
        if not empty_row:
            print("%s" % ' '.join(map(str, grid[rows][:])))
    print("\n")

# Draw a bond between a position and an adjacent position
def draw_bond(grid, current_position, adjacent_position):
    distance = Point((adjacent_position.x - current_position.x)/2, (adjacent_position.y - current_position.y)/2)
    bond_position = Point(int(adjacent_position.x - distance.x), int(adjacent_position.y - distance.y))
    if distance.x == 0:
        grid[bond_position.y][bond_position.x] = "|"
    else:
        grid[bond_position.y][bond_position.x] = "-"

# Print output for an amino acid sequence
def display_output(seq_id, seq, best_node, num_nodes, computation_time):
    print("---" * 30)
    print("Amino acid sequence", seq_id, seq)
    print("\nOptimal folding:")
    display(best_node.grid)
    print("Stability:", best_node.g)
    print("Number of iterations:", num_nodes)
    print("Amino acids: {}/{}".format(len(best_node.coordinates), len(seq)))
    print("Runtime: %0.3f seconds" % computation_time)
    print("---" * 30)

# Fold the sequence
def fold(beam_width):
    counter = 0
    amino_acid_index = 0
    start_coordinate = Point(int(cols/2), int(rows/2))
    grid[start_coordinate.y][start_coordinate.x] = seq_array[amino_acid_index]
    coordinate_list = [start_coordinate]
    start_node = Node(coordinate_list, grid, 0)
    best_node = deepcopy(start_node)
    q1 = PriorityQueue()
    q1.put(start_node)
    q2 = PriorityQueue()

    while not (q1.empty() and q2.empty()):
        iter = 0
        while not q1.empty():
            current_node = q1.get()
            if iter < beam_width:
                if iter == 0: best_node = deepcopy(current_node)
                amino_acid_index = len(current_node.coordinates)
                if amino_acid_index < len(seq_array):
                    
                    children_list = current_node.find_children(seq_array, amino_acid_index)
                    for child in children_list:
                        child.g = current_node.g + child.calculate_stability()
                        child.f = child.g
                        q2.put(child)
                counter += 1
            iter += 1
        
        while not q2.empty():
            item = q2.get()
            q1.put(item)
    
    computation_time = (time.time() - start_time)
    output = Output(best_node, computation_time)
    display_output(seq_id, seq, best_node, counter, computation_time)
    
    # # Uncomment to visualise the conformation using NetworkX
    # lattice.show_conformation(seq, best_node.coordinates)

    exit_code = 0
    if len(best_node.coordinates) < len(seq_array): exit_code = -1
    return exit_code, output

start_time = time.time()
input_file = 'HP50.csv'
num_sequences = 1
beam_width = 3

# To execute Dijkstra's algorithm, set to True
is_exhaustive = False

# # Uncomment to create an output file
# output_file = os.path.splitext(input_file)[0] + '_BS' + str(beam_width) + '.txt'
# sys.stdout = open(output_file, "w")

invalid_ids = []
structure_list = []
stability_list = []
runtime_list = []

# Loop through the input data and run the folding algorithm
for i in range(num_sequences):
    start_time = time.time()
    seq_id = i
    seq = read_sequence(input_file, seq_id)
    seq_array = [char for char in seq]
    cols, rows =  2 * len(seq_array), 2 * len(seq_array)
    grid = [[' ' for i in range(cols)] for j in range(rows)]

    # Set beam width to a value larger than the number of possible conformations
    if is_exhaustive: beam_width = 4**len(seq)

    exit_code, output = fold(beam_width)
    if exit_code == -1: invalid_ids.append(seq_id)

    structure_list.append(output.structure)
    stability_list.append(output.stability)
    runtime_list.append(output.runtime)

# Display the IDs of sequences which could not be folded due to collisions
if len(invalid_ids) > 0:
    print("Invalid IDs:")
    print(invalid_ids)

# Print a list of stability values and a list of runtime values of all input sequences
print("\nStability list:")
print(stability_list)
print("\nRuntime list:")
print(runtime_list)