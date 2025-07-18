mimick(
    ["test/hap1.fa", "test/hap2.fa"],
    "stlfr",
    prefix = "simulated/SIM",
    coverage = 100,
    n_molecules = 2,
    mol_coverage = 0.2,
    mol_length = 80000,
    insert_length = 500,
    insert_stdev = 50,
    read_length = [150,150],
    singletons = 0.35,
    circular = false,
    attempts = 25,
    seed = 0,
    quiet = true
)

mimick(
    "test/fasta.fa",
    "test/test.vcf.gz",
    "standard:haplotagging",
    prefix = "simulated/",
    coverage = 10,
    n_molecules = 2,
    mol_coverage = 0.2,
    mol_length = 80000,
    insert_length = 500,
    insert_stdev = 50,
    read_length = [150,150],
    singletons = 0.35,
    circular = false,
    attempts = 25,
    seed = 0
)