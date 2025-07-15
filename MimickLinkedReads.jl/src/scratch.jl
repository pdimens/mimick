mimick(
    ["test.hap1.fa", "test.hap2.fa"],
    format = "stlfr",
    prefix = "simulated/SIM",
    coverage = 50,
    n_molecules = 2,
    mol_cov = 0.2,
    mol_len = 80000,
    insert_length = 500,
    insert_stdev = 50,
    read_len = [150,150],
    singletons = 0.35,
    circular = false,
    attempts = 25,
    seed = 0
)

mimick(
    "resources/fasta.fa",
    "resources/test.vcf",
    format = "stlfr",
    prefix = "simulated/",
    coverage = 30,
    n_molecules = 2,
    mol_cov = 0.2,
    mol_len = 80000,
    insert_length = 500,
    insert_stdev = 50,
    read_len = [150,150],
    singletons = 0.35,
    circular = false,
    attempts = 25,
    seed = 0
)


#=
julia --project=MimickLinkedReads.jl ./MimickLinkedReads.jl/src/Main.jl test.hap1.fa test.hap2.fa
=#
