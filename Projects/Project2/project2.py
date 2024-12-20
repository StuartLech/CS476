import numpy as np

def read_distance_matrix(file_path):
    with open(file_path, 'r') as file:
        lines = file.readlines()
    num_taxa = int(lines[0].strip())
    matrix = []
    for line in lines[2:]:
        matrix.append(list(map(float, line.strip().split())))
    return np.array(matrix)

def upgma(distance_matrix): # UPGMA Implementation
    num_taxa = distance_matrix.shape[0]
    clusters = {i: [i] for i in range(num_taxa)}
    newick = {i: f"{i}" for i in range(num_taxa)}
    active_clusters = list(range(num_taxa))
    
    while len(active_clusters) > 1:
        min_dist = float('inf')
        for i in active_clusters:
            for j in active_clusters:
                if i < j and distance_matrix[i, j] < min_dist:
                    min_dist = distance_matrix[i, j]
                    min_pair = (i, j)
        
        i, j = min_pair
        print(f"Merging clusters {newick[i]} and {newick[j]} with average distance {min_dist}")
        
        new_cluster = clusters[i] + clusters[j]
        new_cluster_index = max(clusters.keys()) + 1
        clusters[new_cluster_index] = new_cluster
        
        newick[new_cluster_index] = f"({newick[i]},{newick[j]})"
        
        new_dist = np.zeros(len(distance_matrix) + 1)
        for k in active_clusters:
            if k != i and k != j:
                new_dist[k] = (distance_matrix[i, k] + distance_matrix[j, k]) / 2
        
        new_distance_matrix = np.zeros((len(distance_matrix) + 1, len(distance_matrix) + 1))
        new_distance_matrix[:-1, :-1] = distance_matrix
        new_distance_matrix[-1, :-1] = new_dist[:-1]
        new_distance_matrix[:-1, -1] = new_dist[:-1]
        
        distance_matrix = new_distance_matrix
        
        active_clusters.remove(i)
        active_clusters.remove(j)
        active_clusters.append(new_cluster_index)
    
    return newick[active_clusters[0]]

def neighbor_joining(distance_matrix): # Neighbor Joining Implementation
    num_taxa = distance_matrix.shape[0]
    clusters = {i: [i] for i in range(num_taxa)}
    newick = {i: f"{i}" for i in range(num_taxa)}
    active_clusters = list(range(num_taxa))
    
    while len(active_clusters) > 2:
        total_distances = np.sum(distance_matrix, axis=1)
        Q_matrix = np.zeros_like(distance_matrix)
        
        for i in active_clusters:
            for j in active_clusters:
                if i != j:
                    Q_matrix[i, j] = (len(active_clusters) - 2) * distance_matrix[i, j] - total_distances[i] - total_distances[j]
        
        min_Q = float('inf')
        for i in active_clusters:
            for j in active_clusters:
                if i < j and Q_matrix[i, j] < min_Q:
                    min_Q = Q_matrix[i, j]
                    min_pair = (i, j)
        
        i, j = min_pair
        print(f"Merging clusters {newick[i]} and {newick[j]} with Q value {min_Q}")
        
        new_cluster = clusters[i] + clusters[j]
        new_cluster_index = max(clusters.keys()) + 1
        clusters[new_cluster_index] = new_cluster
        
        total_i = total_distances[i]
        total_j = total_distances[j]
        dist_ij = distance_matrix[i, j]
        branch_i = (dist_ij + (total_i - total_j) / (len(active_clusters) - 2)) / 2
        branch_j = dist_ij - branch_i
        newick[new_cluster_index] = f"({newick[i]}:{branch_i},{newick[j]}:{branch_j})"
        
        new_dist = np.zeros(len(distance_matrix) + 1)
        for k in active_clusters:
            if k != i and k != j:
                new_dist[k] = (distance_matrix[i, k] + distance_matrix[j, k] - distance_matrix[i, j]) / 2
        
        new_distance_matrix = np.zeros((len(distance_matrix) + 1, len(distance_matrix) + 1))
        new_distance_matrix[:-1, :-1] = distance_matrix
        new_distance_matrix[-1, :-1] = new_dist[:-1]
        new_distance_matrix[:-1, -1] = new_dist[:-1]
        
        distance_matrix = new_distance_matrix
        
        print(f"Updated distance matrix:\n{distance_matrix}") # Handle Matrix Updates
        
        active_clusters.remove(i)
        active_clusters.remove(j)
        active_clusters.append(new_cluster_index)
    
    i, j = active_clusters
    branch_i = distance_matrix[i, j] / 2
    branch_j = distance_matrix[i, j] - branch_i
    final_newick = f"({newick[i]}:{branch_i},{newick[j]}:{branch_j})"
    
    return final_newick
# Print Console Outputs
def main():
    files = ['DM-p127.txt', 'DM-p139.txt'] # Read Sample Files
    
    for file in files:
        print(f"\nTesting with {file}:\n")
        
        distance_matrix = read_distance_matrix(file)
        
        print("Select algorithm:")
        print("1. UPGMA (pysip)")
        print("2. Neighbor Joining (oyop)")
        choice = input("Enter 1 or 2: ").strip()
        
        if choice == "1":
            print("UPGMA:")
            upgma_tree = upgma(distance_matrix)
            print("Final UPGMA Newick format tree:", upgma_tree)
        elif choice == "2":
            print("Neighbor Joining:")
            neighbor_joining_tree = neighbor_joining(distance_matrix)
            print("Final Neighbor Joining Newick format tree:", neighbor_joining_tree)
        else:
            print("Invalid choice. Skipping this file.")

if __name__ == "__main__":
    main()


