rule generate_group_by_position:
    input:
        "Input/{file}/reference.fasta",
        "Input/{file}/quasispecies.fasta"
    output:
        "Output/{file}/group_by_position.csv"
    script:
        "Pipeline/generate_group_by_position.py"

rule generate_variability_map:
    input:
        "Output/{file}/group_by_position.csv"
    output:
        "Output/{file}/variability_map.png",
        "Output/{file}/variability_map.csv"
    script:
        "Pipeline/generate_variability_map.py"
    