schema = setup_schema(["test.hap1.fa", "test.hap2.fa"], 10.,[150,120])

params = SimParams("this/that", 0.001, 500, 100, 150, 131, 80000, 3.0, 0.05; circular = false, attempts = 50)

molsize = get_molecule_size(params, length(schema["Contig1_1"].sequence))

frags = calculate_fragments(params, molsize)

molecule = get_sequences(schema["Contig1_1"], params, "ATGCAGA", molsize, frags)

convert2fastq(molecule, "I")

#mimick("this/that", ["test.hap1.fa", "test.hap2.fa"], 10.0, 4.0, 80000, 600, 150, [120,150], 0.35, true, 50)