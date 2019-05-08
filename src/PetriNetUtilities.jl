module PetriNetUtilities

using LightGraphs, GraphPlot, DiffEqBiological

function make_petri_net(rs)
    # extract nodes
    species_nodes  = species(rs)
    reaction_nodes = [Symbol("R$j") for j in eachindex(rs.reactions)]

    number_species = length(species_nodes)

    # initialize edge set
    edge_set = Tuple{Int, Int}[]

    # labels
    nodelabel = vcat(species_nodes, reaction_nodes)
    edgelabel = []

    dep_graph = rxtospecies_depgraph(rs)

    for j in eachindex(dep_graph)
        # reaction -> species
        for (k, stoich) in productstoich(rs, j)
            push!(edge_set, (number_species + j, k))
            if stoich > 1
                push!(edgelabel, stoich)
            else
                push!(edgelabel, "")
            end
        end

        # species -> reaction
        for (k, stoich) in substratestoich(rs, j)
            push!(edge_set, (k, number_species + j))
            if stoich > 1
                push!(edgelabel, stoich)
            else
                push!(edgelabel, "")
            end
        end
    end

    number_nodes = length(species_nodes) + length(reaction_nodes)

    g = SimpleDiGraph(number_nodes)

    for edge in edge_set
        add_edge!(g, edge)
    end

    return gplot(g,
        nodelabel = nodelabel,
        edgelabel = edgelabel,
        layout    = spectral_layout)
end

export make_petri_net

end # module
