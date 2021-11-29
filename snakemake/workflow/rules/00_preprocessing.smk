"""
##########################
##                      ##
##    PREPROCESSING     ##
##                      ##
##########################

Snakemake rule file to preprocess Illumina sequence data for viral contigg analysis.

What is accomplished with this script?
    - Non-biological sequence removal (primers, adapters)
    - Common laboratory contaminant removal (NCBI UniVec: https://www.ncbi.nlm.nih.gov/tools/vecscreen/univec/)
    - Host sequence removal (user selected)
    - Removal of redundant sequences (clustering)
        - Creation of sequence count table
        - Calculation of sequence properties (e.g. GC content, tetramer frequencies)

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

rule remove_5prime_primer:
    """
    
    Step 01: Remove 5' primer.
    
    Default RdA/B Primer sequences are provided in  the file primerB.fa. If your lab uses other primers you will need to place them in CONPATH (defined in the Snakefile) and change the file name from primerB.fa to your file name below.
    
    """
    input:
        r1 = os.path.join(READDIR, PATTERN_R1 + file_extension),
        r2 = os.path.join(READDIR, PATTERN_R2 + file_extension),
        primers = os.path.join(CONPATH, "primerB.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_01", PATTERN_R1 + ".s1.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_01", PATTERN_R2 + ".s1.out.fastq")),
        stats = os.path.join(STATS, "step_01", "{sample}.s1.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s01.removeprimerB_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_01", "s01_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=l restrictleft=20 \
            removeifeitherbad=f trimpolya=10 ordered=t rcomp=f ow=t \
            -Xmx{config[System][Memory]}g 2> {log};
        """

rule remove_3prime_contaminant:
    """
    
    Step 02: Remove 3' read through contaminant. This is sequence that occurs if the library fragment is shorter than 250 bases and the sequencer reads through the the 3' end. We use the full length of primerB plus 6 bases of the adapter to detect this event and remove everything to the right of that molecule when detected.
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_01", PATTERN_R1 + ".s1.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_01", PATTERN_R2 + ".s1.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_02", PATTERN_R1 + ".s2.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_02", PATTERN_R2 + ".s2.out.fastq")),
        stats = os.path.join(STATS, "step_02", "{sample}.s2.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s02.remove_3prime_contaminant_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_02", "s02_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=11 ktrim=r removeifeitherbad=f ordered=t rcomp=f ow=t \
            -Xmx{config[System][Memory]}g 2> {log};
        """

rule remove_primer_free_adapter:
    """
    
    Step 03: Remove primer free adapter (both orientations). Rarely the adapter will be seen in the molecule indpendent of the primer. This removes those instances as well as everything to the right of the detected primer-free adapter.
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_02", PATTERN_R1 + ".s2.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_02", PATTERN_R2 + ".s2.out.fastq"),
        primers = os.path.join(CONPATH, "nebnext_adapters.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_03", PATTERN_R1 + ".s3.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_03", PATTERN_R2 + ".s3.out.fastq")),
        stats = os.path.join(STATS, "step_03", "{sample}.s3.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s03.remove_primer_free_adapter_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_03", "s3_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=1 mink=10 ktrim=r removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{config[System][Memory]}g 2> {log};
        """

rule remove_adapter_free_primer:
    """
    
    Step 04: Remove adapter free primer (both orientations). Rarely the primer is detected without the primer. This removes those instances as well as everything to the right of the detected adapter-free primer. 
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_03", PATTERN_R1 + ".s3.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_03", PATTERN_R2 + ".s3.out.fastq"),
        primers = os.path.join(CONPATH, "rc_primerB_ad6.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_04", PATTERN_R1 + ".s4.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_04", PATTERN_R2 + ".s4.out.fastq")),
        stats = os.path.join(STATS, "step_04", "{sample}.s4.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s04.remove_adapter_free_primer_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_04", "s4_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=16 hdist=0 removeifeitherbad=f ordered=t rcomp=t ow=t \
            -Xmx{config[System][Memory]}g 2> {log};
        """

rule remove_vector_contamination:
    """
    
    Step 05: Vector contamination removal (PhiX + NCBI UniVecDB)
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_04", PATTERN_R1 + ".s4.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_04", PATTERN_R2 + ".s4.out.fastq"),
        primers = os.path.join(CONPATH, "vector_contaminants.fa")
    output:
        r1 = temp(os.path.join(TMPDIR, "step_05", PATTERN_R1 + ".s5.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, "step_05", PATTERN_R2 + ".s5.out.fastq")),
        stats = os.path.join(STATS, "step_05", "{sample}.s5.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s05.remove_vector_contamination_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_05", "s5_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            ref={input.primers} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            k=31 hammingdistance=1 ordered=t ow=t \
            -Xmx{config[System][Memory]}g 2> {log};
        """
        
rule remove_low_quality:
    """
    
    Step 06: Remove remaining low-quality bases and short reads. Quality score can be modified in config.yaml (QSCORE).
    
    """
    input:
        r1 = os.path.join(TMPDIR, "step_05", PATTERN_R1 + ".s5.out.fastq"),
        r2 = os.path.join(TMPDIR, "step_05", PATTERN_R2 + ".s5.out.fastq")
    output:
        r1 = temp(os.path.join(TMPDIR, PATTERN_R1 + ".clean.out.fastq")),
        r2 = temp(os.path.join(TMPDIR, PATTERN_R2 + ".clean.out.fastq")),
        stats = os.path.join(STATS, "step_06", "{sample}.s6.stats.tsv")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s06.remove_low_quality_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_06", "s6_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbduk.sh in={input.r1} in2={input.r2} \
            out={output.r1} out2={output.r2} \
            stats={output.stats} \
            ordered=t \
            qtrim=r maxns=2 \
            entropy={config[ENTROPY]} \
            entropywindow=25 \
            trimq={config[QSCORE]} \
            minlength={config[READ_MINLENGTH]} \
            -Xmx{config[System][Memory]}g 2> {log};
        """

rule host_removal_mapping:
    """
    
    Step 07a: Host removal. Must define host in config file (see Paths: Host: in config.yaml). Host should be masked of viral sequence.
    
    Step 7 has several substeps (7a, 7b, 7c and 7d) to subset and fix read pairing issues
    
    If your reference is not available post an issue on GitHub requesting it to be added (https://github.com/shandley/hecatomb)
    
    """
    input:
        r1 = os.path.join(TMPDIR, PATTERN_R1 + ".clean.out.fastq"),
        r2 = os.path.join(TMPDIR, PATTERN_R2 + ".clean.out.fastq"),
        hostpath = HOSTPATH
    output:
        sam = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".sam")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s07.host_removal_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_07a", "s7a_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/minimap2.yaml"
    shell:
        """
        # Map with minimap
        # Remove supplementary alignments (-F 2048)
        # Isolate unmapped (non-host) reads (-f 4)
        minimap2 -ax sr -t {config[System][Threads]} --secondary=no {input.hostpath} {input.r1} {input.r2} 2> {log} \
            | samtools view -f 4 -h --threads {resources.cpus} > {output.sam} 2> {log};
        """

rule extract_host_unmapped:
    """
    
    Step 07b: Extract unmapped (non-host) sequences from sam files
    
    """
    input:
        sam = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".sam")
    output:
        r1 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq")),
        r2 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq")),
        singletons = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s07b.extract_host_unmapped_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_07b", "s07b_{sample}.log")

    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/samtools.yaml"
    shell:
        """
        samtools fastq --threads {resources.cpus} -NO -1 {output.r1} -2 {output.r2} \
        -0 /dev/null \
        -s {output.singletons} \
        {input.sam} 2> {log};
        """

rule nonhost_read_repair:
    """
    
    Step 07c: Parse R1/R2 singletons (if singletons at all)
    
    """
    input:
        singletons = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.singletons.fastq")
    output:
        r1 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq")),
        r2 = temp(os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq"))
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s07c.nonhost_read_repair_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_07c", "s07c_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        reformat.sh in={input.singletons} out={output.r1} out2={output.r2} \
        -Xmx{config[System][Memory]}g 2> {log};
        """

rule nonhost_read_combine:
    """
    
    Step 07d: Combine R1+R1_singletons and R2+R2_singletons to create single R1 and R2 sequence files.
    
    """
    input:
        r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq"),
        r1s = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
        r2s = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
    output:
        r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".all.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".all.fastq")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s07d.nonhost_read_combine_{sample}.txt")
    resources:
        mem_mb=60000,
        cpus=40
    shell:
        """
        cat {input.r1} {input.r1s} > {output.r1};
        cat {input.r2} {input.r2s} > {output.r2};
        """

rule remove_exact_dups:
    """
    
    Step 08: Remove exact duplicates
    
    - Exact duplicates are considered PCR generated and not accounted for in the count table (seqtable_all.tsv)
    
    """
    input:
        os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".all.fastq")
    output:
        os.path.join(QC, "CLUSTERED", PATTERN_R1 + ".deduped.out.fastq")
    priority: 5
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s08.remove_exact_dupes_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_08", "s08_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        dedupe.sh in={input} out={output} \
        ac=f ow=t \
        -Xmx{config[System][Memory]}g 2> {log};
        """
          
rule assembly_kmer_normalization:
    """
    
    Step 14: k-mer normalization. Data reduction for assembly improvement
    
    """
    input:
        r1 = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".unmapped.fastq"),
        r2 = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".unmapped.fastq"),
        r1s = os.path.join(QC, "HOST_REMOVED", PATTERN_R1 + ".singletons.fastq"),
        r2s = os.path.join(QC, "HOST_REMOVED", PATTERN_R2 + ".singletons.fastq")
    output:
        r1_norm = os.path.join(ASSEMBLY, PATTERN_R1 + ".norm.fastq"),
        r2_norm = os.path.join(ASSEMBLY, PATTERN_R2 + ".norm.fastq")
    benchmark:
        os.path.join(BENCH, "PREPROCESSING", "s14.kmer_normalization_{sample}.txt")
    log:
        log = os.path.join(LOGS, "step_14", "s14_{sample}.log")
    resources:
        mem_mb=60000,
        cpus=40
    conda:
        "../envs/bbmap.yaml"
    shell:
        """
        bbnorm.sh in={input.r1} in2={input.r2} \
        extra={input.r1s},{input.r2s} \
        out={output.r1_norm} out2={output.r2_norm} \
        target=100 \
        ow=t \
        -Xmx{config[System][Memory]}g 2> {log}
        """