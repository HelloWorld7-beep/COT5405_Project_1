from Bio import Phylo
import Greedy_PD
import time, statistics
import time
import statistics
import csv
import matplotlib.pyplot as plt


#Load the shark data file of 500 trees (lets just take the first one)
tree_path = "Chond.610sp.10Cal.500TreeSet.tre"
trees = list(Phylo.parse(tree_path, "newick"))
tree = trees[0]

#Some error checking
#print(f"Loaded {len(trees)} trees; using the first one.")
#print(f"Number of tips (species): {len(tree.get_terminals())}")

#Convert to my (V, E, w) representation
V = set()
E = []
w = {}

#Iterate through all clades (i.e. nodes in tree i.e. graph)
for clade in tree.find_clades(order="level"):
    V.add(clade)
    if clade.clades:
        for child in clade.clades:
            E.append((clade, child))
            w[(clade, child)] = child.branch_length or 0.0
            V.add(child)

print(f"Total nodes: {len(V)}, edges: {len(E)}")

#Identify the leaves and root
L = [term for term in tree.get_terminals()]
root = tree.root

print(f"Ex leaves: {[t.name for t in L[:5]]}")

#Map clade objects to simple string IDs for labellig purposes
def simplify_node(x):
    return x.name if x.name else f"id_{id(x)}"

#So I can use strings instead of complex objects
V_named = {simplify_node(v) for v in V}  #Set of string node names
E_named = [(simplify_node(u), simplify_node(v)) for u, v in E]  #List of string edge tuples
w_named = {(simplify_node(u), simplify_node(v)): wt for (u, v), wt in w.items()}  #Weight dict keyed by string edge tuples
L_named = [simplify_node(l) for l in L]  #List of string leaves
root_named = simplify_node(root)  #string root id

#Now calle the greedy pd func
S, total_PD = Greedy_PD.greedy_pd(V_named, E_named, w_named, L_named, n=50, root=root_named)

#print(S)


# -------Runtime Analysis experiment------- #

#Eveything already built above, so
def time_once(n: int) -> float:
    start = time.perf_counter()
    S, total_PD = Greedy_PD.greedy_pd(V_named, E_named, w_named, L_named, n=n, root=root_named)
    t = time.perf_counter() - start
    return t

#Values of n to test
n_values = [10, 25, 50, 100, 200]
results = []

#Run trials and collect data, do 3 trials cause science
for n in n_values:
    trials = 3
    times = [time_once(n) for _ in range(trials)]
    mean = statistics.mean(times)
    tmin, tmax = min(times), max(times)
    results.append({"n": n, "mean": mean, "min": tmin, "max": tmax})

#Save results to CSV
csv_path = "shark_runtime_vs_n.csv"
with open(csv_path, "w", newline="") as f:
    writer = csv.writer(f)
    writer.writerow(["n", "mean_seconds", "min_seconds", "max_seconds"])
    for r in results:
        writer.writerow([r["n"], r["mean"], r["min"], r["max"]])

print(f'Saved CSV: {csv_path}')

#-------Minimalist Runtime Plot from CSV------- #

csv_path = "shark_runtime_vs_n.csv"

#Load CSV data
ns, means = [], []
with open(csv_path, "r") as f:
    reader = csv.DictReader(f)
    for row in reader:
        ns.append(int(row["n"]))
        means.append(float(row["mean_seconds"]))

#Plot
plt.figure()
plt.plot(ns, means, "s--", label="Mean")

plt.xlabel("n (number of selected species)")
plt.ylabel("Runtime (seconds)")
plt.title("Greedy_PD runtime vs n (shark dataset)")
plt.legend()
plt.grid(True, alpha=0.3)

#Save plot
plt.savefig("shark_runtime_plot.png", dpi=300, bbox_inches="tight")
print("Saved plot: shark_runtime_plot.png")
