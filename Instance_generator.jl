using Random
using Printf

function generate_pooling_problem(num_sources::Int, num_pools::Int, num_terminals::Int, num_k::Int, 
                                       density_prob::Float64, 
                                       output_file::String; 
                                       acyclic::Bool=true,
                                       no_direct_st::Bool=false,
                                       min_inter_pool_arcs::Int=1,
                                       seed::Union{Int, Nothing}=nothing)
    """
    Generates a GAMS file for a randomized pooling problem instance 
    
    Parameters:
    - num_sources: Number of sources (|S|)
    - num_pools: Number of intermediate pools (|I|)
    - num_terminals: Number of terminals (|T|)
    - num_k: Number of quality attributes (|K|)
    - density_prob: Expected network density (probability for each possible arc)
    - output_file: Path to the output .gms file
    - acyclic: If true, enforce acyclicity by avoiding backward arcs in pool ordering (groups A-D); if false, allow cycles (group E)
    - min_inter_pool_arcs: Minimum number of arcs between pools (I → I) to guarantee in the instance (default: 1).
                         Ensures a "general" pooling problem with at least this many inter-pool connections.
                         Added deterministically after probabilistic generation, respecting acyclicity if enabled.
                         If num_pools < 2, this is ignored (no possible I → I arcs).
    - seed: Optional random seed for reproducibility
    
    Generation Rules :
    - Total nodes |i| = |S| + |I| + |T|, ordered: S (1 to |S|), I (|S|+1 to |S|+|I|), T (last |T|).
    - Possible arcs: from i ∈ S ∪ I to j ∈ I ∪ T, i ≠ j.
    - For acyclic=true: Pools ordered by index; from pool ik (k≥2) avoid arcs to pools i1..i_{k-1} (no backward arcs).
    - Each possible arc exists independently with probability density_prob.
    - After probabilistic generation, if fewer than min_inter_pool_arcs exist between pools (I → I),
      additional arcs are added deterministically: random selection from remaining possible I → I pairs
      (forward-only if acyclic=true; any direction if false). This guarantees at least min_inter_pool_arcs
      inter-pool arcs, ensuring a general pooling structure with multi-stage blending.
    - Arc costs: c_ij = d_i - d_j if arc exists (0 otherwise).
      - d_i = 0 for i ∈ I.
      - d_i ~ Uniform{Integer} {0,...,10} for i ∈ S.
      - d_j ~ Uniform{Integer} {10,...,30} for j ∈ T.
    - Node capacities: bl_i = 0 for all i (listed only for S ∪ T in output).
      - bu_i ~ Uniform{Integer} {20,...,80} for all i ∈ i.
    - Qualities q(i,k):
      - For i ∈ S: q_{i,k} ~ Uniform{Integer} {0,...,10} for each k.
      - For i ∈ T: q_{i,k} ~ Uniform{Integer} {2,...,6} for each k (quality upper bounds).
      - Undefined (0) for i ∈ I.
    - All values are integers; formatted as floats (e.g., 5.00) in tables for GAMS compatibility.
    - Tables structured as in the example: a/c rows for S∪I, columns for I∪T; q/bl for S∪T; bu for all.
    
    Returns: Nothing (writes to file); prints confirmation, including the final number of inter-pool arcs.
    """
    
    if seed !== nothing
        Random.seed!(seed)
    end
    
    total_i = num_sources + num_pools + num_terminals
    if total_i < 1 || num_sources < 1 || num_pools < 1 || num_terminals < 1
        error("All counts must be positive integers")
    end
    
    if num_pools < 2 && min_inter_pool_arcs > 0
        @warn "num_pools < 2: Cannot add inter-pool arcs; ignoring min_inter_pool_arcs."
        min_inter_pool_arcs = 0
    end
    
    # Node ranges
    sources = 1:num_sources
    pools = (num_sources + 1):(num_sources + num_pools)
    terminals = (total_i - num_terminals + 1):total_i
    from_nodes = vcat(sources, pools)  # S ∪ I
    to_nodes = vcat(pools, terminals)  # I ∪ T
    
    # Generate d_i
    d = zeros(Int, total_i)
    for i in sources
        d[i] = rand(0:10)
    end
    # pools: d=0 already
    for i in terminals
        d[i] = rand(10:30)
    end
    
    # Generate qualities q(i,k)
    q = zeros(Int, total_i, num_k)
    for i in sources
        for kk in 1:num_k
            q[i, kk] = rand(0:10)
        end
    end
    for i in terminals
        for kk in 1:num_k
            q[i, kk] = rand(2:6)
        end
    end
    
    # Generate capacities
    bl = zeros(Int, total_i)  
    bu = [rand(20:80) for _ in 1:total_i]
    
    # Generate adjacency a(i,j) and costs c(i,j)
    a = zeros(Int, total_i, total_i)
    c = zeros(Float64, total_i, total_i)  
    for i in from_nodes
        # Determine possible j based on no_direct_st
        if i in sources
            possible_j = no_direct_st ? pools : to_nodes  
        else  # i in pools
            possible_j = to_nodes  # I ∪ T
        end
        
        for j in possible_j
            if i == j
                continue  # No self-loops
            end
            
            # Acyclicity check: if acyclic and both in pools and j < i (backward)
            if acyclic && i in pools && j in pools && j < i
                continue
            end
            
            # Add arc with probability density_prob
            if rand() < density_prob
                a[i, j] = 1
                c[i, j] = Float64(d[i] - d[j])
            end
        end
    end
    
    # Ensure minimum inter-pool arcs (I → I) - unchanged
    if min_inter_pool_arcs > 0
        # Count existing inter-pool arcs
        existing_inter_pool = 0
        for i in pools
            for j in pools
                if i != j && a[i, j] == 1
                    existing_inter_pool += 1
                end
            end
        end
        
        needed = max(0, min_inter_pool_arcs - existing_inter_pool)
        if needed > 0
            # Collect possible I → I pairs that don't have an arc yet
            possible_pairs = Vector{Tuple{Int, Int}}()
            for i in pools
                for j in pools
                    if i != j && a[i, j] == 0
                        # Respect acyclicity: only forward if acyclic
                        if acyclic && j < i
                            continue
                        end
                        push!(possible_pairs, (i, j))
                    end
                end
            end
            
            if length(possible_pairs) < needed
                @warn "Only $(length(possible_pairs)) possible I → I pairs available; adding all $(length(possible_pairs)) (less than requested $needed)."
                needed = length(possible_pairs)
            end
            
            # Shuffle for randomness and add the first 'needed' pairs
            shuffle!(possible_pairs)
            for (idx, (i, j)) in enumerate(possible_pairs[1:needed])
                a[i, j] = 1
                c[i, j] = Float64(d[i] - d[j])
                println("Added deterministic inter-pool arc: $i → $j (to meet minimum)")
            end
        end
        
        # Recount final inter-pool arcs for reporting
        final_inter_pool = 0
        for i in pools
            for j in pools
                if i != j && a[i, j] == 1
                    final_inter_pool += 1
                end
            end
        end
        println("Inter-pool arcs: probabilistic generated $existing_inter_pool, final total $final_inter_pool (min requested: $min_inter_pool_arcs)")
    end
    
    # Build the GAMS code as a vector of strings
    lines = String[]
    
    push!(lines, "# Declare sets")
    push!(lines, "set i    / 1*$(total_i)  /;")
    push!(lines, "set s / 1*$(num_sources)  /;")
    term_start = total_i - num_terminals + 1
    push!(lines, "set t / $(term_start)*$(total_i)  /;")
    push!(lines, "set k    / 1*$(num_k)  /;")
    
    push!(lines, "")
    push!(lines, "alias (i,j);")
    
    push!(lines, "")
    push!(lines, "# The arc unit cost c_{ij}")
    push!(lines, "table c(i,j)")
    
    # Headers for c: columns for to_nodes
    header_width = 7
    header = join([lpad(string(j), header_width, " ") for j in to_nodes], "")
    push!(lines, "          $(header)")
    
    # Rows for from_nodes
    for i in from_nodes
        row_vals = [isapprox(c[i,j], 0.0) ? " 0.00" : @sprintf("%7.2f", c[i,j]) for j in to_nodes]
        row_str = join(row_vals, "")
        push!(lines, "  $(i)   $(row_str)")
    end
    push!(lines, "  ;")
    
    push!(lines, "")
    push!(lines, "# The adjacency matrix (the arcs set A)")
    push!(lines, "table a(i,j)")
    
    # Headers for a: shorter spacing
    header_width_a = 3
    header = join([lpad(string(j), header_width_a, " ") for j in to_nodes], "")
    push!(lines, "      $(header)")
    
    for i in from_nodes
        row_vals = [" $(a[i,j])" for j in to_nodes]
        row_str = join(row_vals, "")
        push!(lines, "  $(i)  $(row_str)")
    end
    push!(lines, "  ;")
    
    push!(lines, "")
    push!(lines, "# Source qualities/terminal quality upper bounds")
    push!(lines, "table q(i,k)")
    
    # Headers for k
    header_width_q = 7
    header = join([lpad(string(kk), header_width_q, " ") for kk in 1:num_k], "")
    push!(lines, "          $(header)")
    
    # Rows for sources and terminals
    st_nodes = vcat(sources, terminals)
    for i in sort(st_nodes)  # Sort to match example order (sources then terminals)
        row_vals = [@sprintf("%7.2f", Float64(q[i,kk])) for kk in 1:num_k]  # Cast to float for formatting
        row_str = join(row_vals, "")
        push!(lines, "  $(i)   $(row_str)")
    end
    push!(lines, "  ;")
    
    push!(lines, "")
    push!(lines, "# Node capacity lower bound")
    push!(lines, "parameter bl(i) /")
    
    # List only for S ∪ T as 0.00
    bl_strs = []
    for i in sort(vcat(sources, terminals))
        push!(bl_strs, "                   $(i)     0.00")
    end
    push!(lines, join(bl_strs, "\n"))
    push!(lines, " / ;")
    
    push!(lines, "")
    push!(lines, "# Node capacity upper bound")
    push!(lines, "parameter bu(i) /")
    
    bu_strs = []
    for i in 1:total_i
        push!(bu_strs, "                   $(i)    $(@sprintf("%6.2f", Float64(bu[i])))")
    end
    push!(lines, join(bu_strs, "\n"))
    push!(lines, " / ;")
    
    # Write to file
    gams_code = join(lines, "\n")
    write(output_file, gams_code)
    println("Generated GAMS file: $output_file (acyclic=$(acyclic), no_direct_st=$(no_direct_st), density=$(density_prob), min_inter_pool_arcs=$min_inter_pool_arcs, sizes: S=$num_sources, I=$num_pools, T=$num_terminals, K=$num_k)")
    
    return nothing
end

#Sources , Pools, Terminals, K, Density, Output file, Acyclic, Min inter-pool arcs, Seed

for i in 1:5
    generate_pooling_problem(30, 30, 30, 5, 0.2,"P$i.gms", acyclic=true, seed=i*100)

end

#I : 50 sources, 25 pools, 10 terminals, 5 quality attributes, density 0.6, acyclic=true, min_inter_pool_arcs=1, seed=100-500 (5 instances with different seeds)
#S : 20 sources, 50 pools, 10 terminals, 5 quality attributes, density 0.7, acyclic=true, min_inter_pool_arcs=1, seed=100-500 (5 instances with different seeds)
#A : 30 sources, 25 pools, 18 terminals, 8 quality attributes, density 0.6, acyclic=true, min_inter_pool_arcs=1, seed=100-500 (5 instances with different seeds)
#K : 30 sources, 40 pools, 30 terminals, 4 quality attributes, density 0.4, acyclic=true, min_inter_pool_arcs=1, seed=100-500 (5 instances with different seeds)