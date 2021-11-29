"""
####################
##                ##
##    ANNOTATE    ##
##                ##
####################

Snakemake rule file to preprocess Illumina sequence data for viral contigg analysis.

What is accomplished with this script?
    - Taxonomic Annotation

Additional Reading:
- Decatomb GitHub: https://github.com/shandley/decatomb
- Hecatomb GitHub: https://github.com/shandley/hecatomb
- Official Snakemake documentation: https://snakemake.readthedocs.io/en/stable/

Written by: Scott Handley (handley.scott@gmail.com), March 2021
"""

import os
import sys
    
# NOTE: bbtools uses "threads=auto" by default that typically uses all threads, so no need to specify. 
# -Xmx is used to specify the memory allocation for bbtools operations
# Set your -Xmx specifications in your configuration file 

rule calculate_contig_dictionary_properties:
"""

Step 22: Calculate contig sequence properties properties (ie. GC-content, tetramer frequencies) per sequence

"""
input:
    os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "assembly.fasta")
output:
    gc = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary_properties.gc"),
    tetramer = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary_properties.tetramer"),
    seq_properties = os.path.join(ASSEMBLY, "CONTIG_DICTIONARY", "FLYE", "contig_dictionary_properties.tsv")
benchmark:
    os.path.join(BENCH, "PREPROCESSING", "s22.calculate_contig_dictionary_properties.txt")
log:
    log1 = os.path.join(LOGS, "step_20", "s20.gc.log"),
    log2 = os.path.join(LOGS, "step_20", "s20.tetramer.log")
resources:
    mem_mb=60000,
    cpus=40
conda:
    "../envs/bbmap.yaml"
shell:
    """
    # Calcualate per sequence GC content
    countgc.sh in={input} format=2 ow=t | awk 'NF' > {output.gc};
    sed -i '1i id\tGC' {output.gc} 2> {log.log1};
    
    # Calculate per sequence tetramer frequency
    tetramerfreq.sh in={input} w=0 ow=t | \
    tail -n+2 | \
    cut --complement -f2 > {output.tetramer} 2> {log.log2};
    
    sed -i 's/scaffold/id/' {output.tetramer};
    
    # Combine
    csvtk join -f 1 {output.gc} {output.tetramer} -t -T > {output.seq_properties};
    
    """