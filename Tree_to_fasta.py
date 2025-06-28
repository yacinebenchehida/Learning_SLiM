import tskit  # Library for tree sequence data
import random  # For random sampling

def sample_individuals_per_population(ts, sample_size=3):
    pop_counts = {pop.id: 0 for pop in ts.populations()}  # Count individuals per population
    pop_individuals = {pop.id: [] for pop in ts.populations()}  # Store individuals by population
    for ind in ts.individuals():  # Iterate all individuals
        first_node = ind.nodes[0]  # Get first node ID for individual
        pop_id = ts.node(first_node).population  # Population of this node
        pop_counts[pop_id] += 1  # Increment population count
        pop_individuals[pop_id].append(ind.id)  # Append individual ID to population list
    sampled_inds_per_pop = {}  # Dict for sampled individuals
    for pop_id, count in pop_counts.items():  # For each population and count
        if count == 0:  # Skip empty populations
            continue
        pop_name = ts.population(pop_id).metadata.get("name", f"pop{pop_id}")  # Get pop name
        print(f"Population {pop_id} ({pop_name}): {count} individuals")  # Print pop info
        n = min(sample_size, count)  # Number to sample
        sampled_inds = random.sample(pop_individuals[pop_id], n)  # Sample individuals
        sampled_inds_per_pop[pop_id] = sampled_inds  # Store sampled inds
        print(f"Sampled individuals from population {pop_id}: {sampled_inds}")  # Print sampled inds
    return sampled_inds_per_pop  # Return dict

def subset_ts_by_individuals(ts, sampled_inds):
    all_individuals = []  # Flat list of all sampled individual IDs
    for inds in sampled_inds.values():  # Iterate lists per population
        all_individuals.extend(inds)  # Add to flat list
    nodes_to_keep = []  # List of nodes to keep
    for ind_id in all_individuals:  # For each sampled individual
        ind = ts.individual(ind_id)  # Get individual object
        nodes_to_keep.extend(ind.nodes)  # Add individual's nodes
    subsetted_ts = ts.simplify(samples=nodes_to_keep)  # Simplify TS to these nodes
    return subsetted_ts  # Return subsetted TS

def fasta_all_haplosomes_with_names(ts, pop_names=None):  # Define function with optional population names
    nucleotides = ["A", "T", "G", "C"]  # List of possible bases
    seq_len = int(ts.sequence_length)  # Total length of the simulated sequence

    pop_to_inds = {pop.id: [] for pop in ts.populations()}  # Map each population to its individuals
    for ind in ts.individuals():  # Loop over all individuals
        first_node = ind.nodes[0]  # Take the first haplotype node of the individual
        pop_id = ts.node(first_node).population  # Get the population ID of that node
        pop_to_inds[pop_id].append(ind.id)  # Add the individual to the corresponding population list

    non_empty_pops = [pop_id for pop_id, inds in pop_to_inds.items() if inds]  # Keep only populations that have individuals

    pop_id_to_name = {}  # Dictionary mapping population ID to name
    for i, pop_id in enumerate(non_empty_pops):  # For each non-empty population
        if pop_names and i < len(pop_names):  # If user provided a name for this population
            pop_id_to_name[pop_id] = pop_names[i]  # Use user-provided name
        else:
            pop_id_to_name[pop_id] = f"pop{pop_id}"  # Otherwise, assign default name

    variant_dict = {int(var.site.position): var for var in ts.variants()}  # Map position to variant object
    samples = ts.samples()  # Get all haploid nodes (haplosomes)

    nt_map = {}  # Map each position to allele index → base (e.g. {0: 'A', 1: 'G'})
    for pos in range(seq_len):  # Loop over all positions
        var = variant_dict.get(pos)  # Check if a variant exists at this position
        if var is None:  # No mutation recorded at this position
            base = random.choice(nucleotides)  # Random base for invariant position
            nt_map[pos] = {0: base}  # All haplosomes will get this base
        else:
            alleles_present = set(var.genotypes[samples])  # Get allele values in current haplosomes
            if len(alleles_present) == 1:  # Monomorphic variant (same allele everywhere)
                base = random.choice(nucleotides)  # Random base
                nt_map[pos] = {0: base}  # Map that allele to the base
            else:  # True variant site
                ref = random.choice(nucleotides)  # Random reference base
                alt = random.choice([n for n in nucleotides if n != ref])  # Alternate base ≠ reference
                nt_map[pos] = {0: ref, 1: alt}  # Map alleles to ref/alt bases

    # Choose one haplotype per individual and build sequences for them
    sequences = {}  # Dictionary key: (pop_id, ind_id), value: list of bases
    for pop_id, inds in pop_to_inds.items():  # For each population and its individuals
        if pop_id not in pop_id_to_name:  # Skip empty populations
            continue
        for ind_id in inds:  # For each individual in population
            ind = ts.individual(ind_id)  # Get individual object
            chosen_node = random.choice(ind.nodes)  # Pick one haplotype node randomly
            sequences[(pop_id, ind_id)] = []  # Initialize sequence list for this individual

            for pos in range(seq_len):  # Loop over all positions
                var = variant_dict.get(pos)  # Get variant object at position (if any)
                if var is None:  # No variant here
                    base = nt_map[pos][0]  # Invariant base
                else:
                    gt = var.genotypes[chosen_node]  # Allele index for chosen haplotype
                    base = nt_map[pos].get(gt, nt_map[pos][0])  # Corresponding base or fallback to ref
                sequences[(pop_id, ind_id)].append(base)  # Append base to sequence

    for (pop_id, ind_id), seq_list in sequences.items():  # For each individual sequence
        pop_name = pop_id_to_name[pop_id]  # Get population name
        header = f">{pop_name}_ind{ind_id}"  # FASTA header with pop name and individual id
        seq_str = "".join(seq_list)  # Join bases into string
        print(header)  # Print header
        print(seq_str)  # Print sequence
           

if __name__ == "__main__":
    ts = tskit.load("output.trees")  # Load tree sequence from file
    sampled_inds = sample_individuals_per_population(ts)  # Sample 3 individuals per population
    subsetted_ts = subset_ts_by_individuals(ts, sampled_inds)  # Subset TS to sampled individuals only
    pop_names = ["NAT", "BS", "IB", "MA", "NBB"]
    fasta_str = fasta_all_haplosomes_with_names(subsetted_ts,pop_names)  # Generate FASTA for sampled data only
