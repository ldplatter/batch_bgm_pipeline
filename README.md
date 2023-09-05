*format_bgm_json_result.py, parse_fasta.py, and prune_tree.py were integrated into the run_bgm.py thus **not needed** if 
running the whole test.*

## Batch_BGM_pipeline
The goal of this pipeline is to carry out co-evolutionary analysis 
of proteins of interest utilizing HyPhy Bayesian Graphical Model package, assisted with the previous ERC pipeline.

### Input files
1. Two protein multiple sequences alignment (MSA) files, **could be 
either raw or aligned files**. Some common taxa are needed for the BGM 
to function and provide valid results.
2. A complete version of phylogenetic master tree that contains all species
(or more) in MSAs. The tree will be automatically pruned so that it could 
match with two proteins' **shared taxa**.

A basic knowledge about BGM method (see hyphy webpage) is favorable, including
the parameter settings. The BGM pipeline uses General Time Reversible (GTR)
model as the baseline model when inferring substitutions and evolution rates,
but could be modified in the script (_"call_bgm func"_) along with other 
parameters, if necessary.

### Output files
1. A json file encoded with the raw data of conditional dependencies between 
every tested nodes. Each amino acid residue that met the minimum substitution
requirements were regarded as one node.
2. Processed summary CSV files filtering interactions
based on posterior probabilities cutoffs (primarily 50% and 90%).
4. Bayesian graphs visualizing interactions between selected nodes using ```networkx```.

