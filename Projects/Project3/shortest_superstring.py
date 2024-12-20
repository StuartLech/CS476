def read_fragments(filename):
    with open(filename, 'r') as f:
        fragments = [line.strip() for line in f if line.strip()]
    return fragments

def remove_substrings(fragments):
    unique_fragments = set(fragments)
    fragments_to_remove = set()
    for frag1 in unique_fragments:
        for frag2 in unique_fragments:
            if frag1 != frag2 and frag1 in frag2:
                fragments_to_remove.add(frag1)
                break
    return list(unique_fragments - fragments_to_remove)

def compute_overlap(s1, s2):
    max_overlap = min(len(s1), len(s2))
    for i in range(max_overlap, 0, -1):
        if s1[-i:] == s2[:i]:
            return i
    return 0

def build_overlap_graph(fragments):
    edges = []
    for i, s1 in enumerate(fragments):
        for j, s2 in enumerate(fragments):
            if i != j:
                overlap = compute_overlap(s1, s2)
                if overlap > 0:
                    edges.append((overlap, i, j))
    return edges

def find(parent, i):
    if parent[i] != i:
        parent[i] = find(parent, parent[i])
    return parent[i]

def union(parent, rank, x, y):
    xroot = find(parent, x)
    yroot = find(parent, y)
    if xroot == yroot:
        return False
    if rank[xroot] < rank[yroot]:
        parent[xroot] = yroot
    else:
        parent[yroot] = xroot
        if rank[xroot] == rank[yroot]:
            rank[xroot] += 1
    return True

def greedy_hamiltonian_path(fragments, edges):
    n = len(fragments)
    parent = [i for i in range(n)]
    rank = [0] * n
    in_degree = [0] * n
    out_degree = [0] * n
    selected_edges = {}

    edges.sort(reverse=True)

    for overlap, u, v in edges:
        if out_degree[u] == 0 and in_degree[v] == 0:
            if find(parent, u) != find(parent, v):
                union(parent, rank, u, v)
                out_degree[u] = 1
                in_degree[v] = 1
                selected_edges[u] = (v, overlap)
                if len(selected_edges) == n - 1:
                    break

    # Find the starting node (node with in_degree 0)
    start_nodes = [i for i in range(n) if in_degree[i] == 0]
    if len(start_nodes) != 1:
        raise Exception("Multiple starting nodes found, invalid Hamiltonian path")

    path = []
    current = start_nodes[0]
    while True:
        path.append(current)
        if current in selected_edges:
            current = selected_edges[current][0]
        else:
            break

    return path, selected_edges

def reconstruct_superstring(fragments, path, selected_edges):
    superstring = fragments[path[0]]
    for i in range(1, len(path)):
        prev = path[i - 1]
        curr = path[i]
        overlap = selected_edges[prev][1]
        superstring += fragments[curr][overlap:]
    return superstring

def main():
    filename = input("Enter the filename: ")
    fragments = read_fragments(filename)
    fragments = remove_substrings(fragments)
    edges = build_overlap_graph(fragments)
    path, selected_edges = greedy_hamiltonian_path(fragments, edges)
    superstring = reconstruct_superstring(fragments, path, selected_edges)
    print("Superstring:", superstring)

if __name__ == "__main__":
    main()
