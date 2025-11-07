def greedy_pd(V, E, w, L, n, root):

    #Function for getting maximized evolutionary diversity from subset of species of phylo tree

    '''Returns (S, totalPD) where:
      S = set of selected leaves (size min(n, len(L)))
      totalPD = total phylogenetic diversity gained'''

    #Build adjacency list graph (used graph not matrix bec should be sparse as it's a tree)
    adjacency = {v: [] for v in V}
    for u, v in E:
        adjacency[u].append(v)
        adjacency[v].append(u)

    #Helper func made to read an undirected edge weight
    def edge_weight(u, v):
        if (u, v) in w:
            return w[(u, v)]
        if (v, u) in w:
            return w[(v, u)]

    #"Root" the tree at abstract "root" of phylo tree, set curr weight = 0 -> parent[] and weight_to_parent[]
    parent = {root: None}
    weight_to_parent = {root: 0.0}

    #Stack of visited vertices
    stack = [root]

    while stack:
        u = stack.pop()
        for v in adjacency[u]:
            #If it hasn't been visited yet
            if v not in parent:
                parent[v] = u
                weight_to_parent[v] = edge_weight(u, v)
                stack.append(v)

    #Precompute the reverse path edges i.e. child->parent, for each candidate leaf
    leaves = list(L)
    path_edges = {}
    for x in leaves:
        path = []
        c = x
        while parent[c] is not None:
            p = parent[c]
            path.append((c, p))  #Child -> Parent
            c = p
        path_edges[x] = path  #order = leaf->root

    # Greedy loop
    covered = set()  #Directed edges = already covered
    selected = set()
    total_pd = 0.0
    target = min(n, len(leaves))

    while len(selected) < target:
        best_leaf = None
        best_gain = float("-inf")

        #Evaluate gain in PD for each leaf.
        for x in leaves:
            if x in selected:
                continue
            gain = 0.0
            for child, par in path_edges[x]:
                if (child, par) in covered:
                    break
                gain += weight_to_parent[child]
            if gain > best_gain:
                best_gain = gain
                best_leaf = x

        #Make the choice and update coverage/gain
        if best_leaf is None:
            break
        selected.add(best_leaf)
        total_pd += best_gain
        covered.update(path_edges[best_leaf])

    return selected, total_pd  #Return the subset & total phylo divv


