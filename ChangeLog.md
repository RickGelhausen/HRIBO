### version 2.0.0-dev [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de)
 * Moved to the standard Snakemake layout (workflow/ + config/), so the workflow is
   relocatable and no longer has to be cloned into a directory named HRIBO
 * Replaced the 73 hardcoded HRIBO/scripts/... shell paths with a SCRIPTS global
 * Config is now read from config/config.yaml in the project directory rather than
   from inside the clone
 * Replaced the hand-written config validation with JSON schemas
   (workflow/schemas/config.schema.yaml, samples.schema.yaml)
 * Added preflight input validation that runs before the DAG is built, reporting
   every problem at once with a suggested fix rather than failing later inside an
   unrelated tool. Notably it detects genome/annotation sequence identifier
   mismatches and explains the likely intended mapping
 * Added checks for embedded ##FASTA sections, out-of-bounds features, missing CDS
   or rRNA/tRNA features, unreadable or truncated fastq files, unusable
   differential expression setups and impossible metagene windows
 * Fixed parse_read_lengths rejecting the documented "22,23,27,34-35" syntax
 * Fixed equalize_dictionary_keys dropping the stop profile for chromosomes seen
   only on the start side, and aliasing all filled-in read lengths onto one array
 * Restored the missing report/ directory; all 13 caption files were absent, so the
   report: directive and --report were broken
 * Added envs/bed.yaml, which was referenced but missing
 * Added a pytest suite (41 tests) covering the validation and metagene helpers
 * The TIS advisor now evaluates both the 5' and the 3' read end and recommends
   whichever the protocol defines more precisely, rather than analysing only the
   end named in the config. Which end carries the cleaner signal is organism and
   nuclease dependent, and ORFBounder accepts either
 * The two ends are compared on the consistency of the estimated offset across
   read lengths, not on peak height. For a read of fixed length the 5' and 3'
   profiles are the same profile shifted, so they are equally sharp by
   construction and peak height cannot separate them; a difference in pooled
   sharpness mostly records how many read lengths each end's search pooled
 * Both ends are reported side by side, so a large difference between them is
   visible rather than hidden behind the winner
 * A read length now joins the recommended set only if it improves the pooled peak
   by at least 2%, which stops a length with no real signal being swept in on a
   rounding difference
 * tisAdvisorSettings.mappingMethod became mappingMethods, a list

 * Fixed the DeepRibo A-site occupancy track. The reverse-strand A-site was
   computed from `read.pos - read_length`, placing it roughly a full read length
   outside the alignment: for a 30 nt read at position 1000 the A-site landed at
   982, 18 nt before the read even starts. Every reverse-strand gene therefore fed
   DeepRibo a misplaced signal. The forward strand was 2 nt off, from mixing a
   1-based position into a 0-based bedgraph. Both strands now use DeepRibo's
   documented convention, a 12 nt offset from the 3' end of the read
 * The A-site offset is now configurable as predictionSettings.deepriboASiteOffset.
   DeepRibo's published value of 12 was derived from E. coli and does not
   necessarily transfer to other organisms or digestion protocols. Note that it is
   not interchangeable with the offset the TIS advisor reports: that one is
   measured from the 5' end to the P-site, this one from the 3' end to the A-site,
   and converting between them depends on read length
 * A-site strand detection now uses the reverse flag rather than testing
   `flag == 0` / `flag == 16`. SAM FLAG is a bit field, so the equality test only
   matches reads with no other bit set. In HRIBO's own pipeline the alignments
   reaching this step carry only flags 0, 16 and 4, so this was fragility rather
   than data loss; it matters for externally produced BAMs
 * A-site bedgraph output is now sorted, and an empty track fails immediately
   rather than letting DeepRibo fail later and obscurely
 * Fixed parameter_estimation.R writing every library's S-curve diagnostic to the
   same relative path, so concurrent jobs overwrote each other. The destination is
   now per library
 * read_parameters now raises on a truncated parameters file instead of passing
   "failed" to DeepRibo as an RPKM cutoff
 * Added log directives and resource declarations to the DeepRibo rules, and
   documented that create_deepribo_gff.py writes `dist` into the GFF phase column
   deliberately, since merge_duplicates_deepribo.py reads it back from there

 * Rewrote the metagene figures. The previous version drew every read length as a
   line on one pair of shared axes, cycling 10 colours and 6 dash patterns, which
   with the default 10 read lengths is 20 overlapping traces. Replaced by a
   read-length-against-position heatmap, small multiples per read length, a
   reading frame composition chart and a read length distribution
 * The heatmap encodes enrichment over each read length's own background, so a
   sparse read length stays legible beside a deep one. Normalising each row by its
   maximum, the obvious choice, inverts the figure: noise saturates and real peaks
   wash out
 * Added lib/theme.py, a single colour-vision-deficiency safe plotly theme; no
   plotting code contains a literal colour any more
 * Interactive reports now embed plotly.js once per page rather than once per
   figure, which was adding several megabytes per plot
 * Added a TIS caller advisor (workflow/scripts/tis_advisor.py) that estimates a
   P-site offset per read length, scores each for usability, greedily searches
   read length combinations, and emits an ORFBounder-ready config block alongside
   an HTML report, a JSON document and a TSV of the evidence
 * The advisor reports "no recommendation" with reasons rather than inventing a
   setup when no read length carries a usable initiation signal, and grades
   confidence by which evidence it rests on. Peak detection is a z-test against the
   upstream background rather than a bare ratio: a ratio test accepts pure Poisson
   noise, which routinely reaches 4x its own median
 * New tisAdvisorSettings section in the config

 * Collapsed the duplicated rule layer: 159 rules in 2739 lines became 105 rules
   in 2102 lines, with no change to any output path
   - 38 near-identical coverage track rules became 2 generic ones
   - 14 read counting rules became 3
   - the 4 createOverviewTable* variants became 1
 * Moved helper functions and lookup tables into workflow/rules/common.smk
 * Fixed CONTRASTS being auto-populated after the diffex rule files were already
   included, which left poolxtail, poolriborex, pooldeltate and contrastInput with
   empty input lists whenever contrasts were not set explicitly in the config.
   xtail_all.csv, riborex_all.csv and deltate_all.csv were pooled from no files at
   all, so the differential expression columns of overview.xlsx were wrong in the
   default configuration
 * Replaced the FTP/HTTP storage plugin downloads with plain curl, removing the
   snakemake-storage-plugin-ftp and -http requirements. Both rules also moved their
   downloads out of the storage cache with mv, which broke re-runs
 * Removed 78 of 85 redundant mkdir calls; the remaining 7 create directories that
   are not output parents and are genuinely needed
 * Removed dead code: the samstrandswap NOTSET branch that no sample sheet could
   select, the segemehl params.fastq branch referencing undefined inputs, a
   duplicated get_inputs_paired definition, and a commented-out copy of merge_fastq
 * Input helpers now raise on a layout mismatch rather than returning None, which
   previously produced a rule with no inputs
 * Declared mem_mb and runtime on 39 rules, so cluster requirements live with the
   rules instead of in the profile
 * Rewrote the SLURM profile for Snakemake 9 (executor: slurm). It referenced
   "slurm-jobscript.sh" while the file was named slurm_jobscript.sh, and 11 of its
   set-threads/set-resources entries named rules that do not exist. The custom job
   script is unsupported by the SLURM executor plugin, so its module loads moved
   into slurm_run.sh
 * Replaced the deprecated singularity: directive with container:
 * Removed 7 scripts that no rule referenced
 * Switched the JSON schemas to draft 2020-12 to match the installed validator
 * Bumped the minimum Snakemake version to 9.0.0 and updated the CLI usage

### version 1.8.1 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 22.05.25
 * fixed off-by-one error in rna filtering rule. 

### version 1.8.0 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 19.04.24
 * Added full support for paired-end data and mixed data (paired-end samples mixed with single-end samples)
 * paired-end data is now correctly trimmed, but will still be converted to single-end data due to lack of paired-end bam support of many downstream tools
 * Updated to version 8.10.7 of snakemake
 * Changed FileProviders to Storage
 * Fixed bug with blastDB generation in reparation
 * Added labels to the PCA plots
 * updated ReadTheDocs

### version 1.7.1 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 23.06.2023
 * exchanged ftp server for swiss prot fasta download

### version 1.7.0 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 24.02.2023
 * improved Snakefile and added validation for config file
 * updated config file structure
 * added 3D PCA plots for quality control
 * added 2D PCA plots + hierarchical clustering of read counts for quality control
 * fixed error that caused crashes when annotation had extra plasmids not present in genome file
 * updated manual
 * updated ReadTheDocs

### version 1.6.2 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 04.10.2022
 * fixed bug where the new input tables were not correctly created depending on the input contrast.
 * fixed bug with wrong condition vector labeling causing incorrect behavior in riborex log2FC calculation
 * ensured that all three tools, deltaTE, riborex and xtail use the same orientation: treated vs untreated

### version 1.6.1 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 29.09.2022
 * moved differential expression input creation to python to avoid compatability errors with new R versions.
 * fixed bug that caused overview_excel to fail when using differential expression

### version 1.6.0 [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) 23.09.2022
 * fixed bug in metagene-profiling that caused negative strand not to be used correctly * fixed bug in metagene-profiling that caused negative strand not to be used correctly.
 * complete restructering of differential expression analysis
 * added deltaTE and re-added riborex
 * fixed issue with sorting between python and R
 * changed .xlsx format for differential expression, using sheets rather than multiple files.
 * cleaned differential expression scripts
 * added option to customize the input contrasts
 * updated tool versions to newest available

### version 1.5.2 (added to 1.6.0) [Rick Gelhausen](mailto:gelhausr@informatik.uni-freiburg.de) -
 * removed riborex support and improved xtail support
 * fixed issue with missing xtail results in overview tables
 * added all annotated features to overview.xlsx
 * fixed bug where the min normalization was not correctly using the minimal value.
 * changed excel_writer to create xlsx files with freezed header
 * ensured correct sorting of xtail sorted files
 * fixed issue where xtail result tables had wrong naming
 * Changed naming of makereport.sh output files to ensure correct sorting

### version 1.5.1 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 28.07.2021
 * fixed issue which sometimes caused crashes in coverage file generation.
 * fixed versioning issues with snakemake and multiQC.
 * updated documentation + Manual

### version 1.5.0 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 09.06.2021
 * added improved meta gene profiling figures
 * added automatic detection of peak read lengths
 * added detection of best offset for individual readlengths
 * added pip packaging
 * functions organized in library
 * added old_locus_tag information to annotation excel files.
 * fixed a bug in the metagene-profiling that would cause a strange peak at the start/end for 5 or 3prime mappings.
 * improved the plot output for the metagene-profiling when checking many read lengths.
 * slightly changed the metagene-profiling read-counting to make it more comparable to the build in coverage files.
 * added metagene profiling for stop codons, this is automatically used on TTS and RNATTS libraries.

### version 1.4.4 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 18.09.2020
 * added overlap column to the overview table, to show whether an entry overlaps with an annotated gene.
 * improved the handling of TE wildcards in excel scripts, allowing input combinations that were previously not possible.
 * minor fix for mapping script that ensures the header is set correctly for the "mil" mapping.
 * minor fix that ensures the gene_name is set correctly and not sometimes replaced by the locus_tag.
 * fixed a bug where tranlation of nt_seqs were done only up to the first stop codon.
 * improved gtf2gff3 script to use gene_biotype to find RNA features.
 * improved start/stop gff files (including motifs and frames)
 * added 15nt upstream of the start codon to each .xlsx table.
 * updated reparation environment
 * some code clean-up

### version 1.4.3 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 15.05.2020
 * added genome-browser identifier to overview table
 * completely reworked gff2 support by transforming gff2 to gff3. This is unavoidable as some tools require gff3 format
 * added script to convert gff2 to gff3
 * improved overview table by adding and reordering columns
 * improved TE calculation to be less confusing, using NaN if read-counts are 0 and would lead to division by 0
 * improved excel file generation
 * resolved redundancies in excel files
 * updated manual + online documentation

### version 1.4.2 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 25.04.2020
 * fixed bug causing generate_overview_excel.py to crash for certain deepribo files
 * updated manual

### version 1.4.1 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 20.04.2020
 * added support for multi-condition results for differential expression in the overview table
 * updated manual

### version 1.4.0 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 17.04.2020
 * added overview table that aggregates all important information about ORFs
 * updated readcounting to consitently use the same method
 * updated makereport script
 * updated manual

### version 1.3.2 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 24.03.2020
 * updated sge.yaml
 * updated torque.yaml
 * bugfixed retrieveAnnotation
 * added xlsx output for xtail and riborex
 * added updated annotation containing reparation predictions and the orginal annotation
 * updated manual

### version 1.3.1 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 13.02.2020

 * added makereport script

### version 1.3.0 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 30.01.2020

 * integration of deepribo
 * added metagene profiling
 * added pseudogene read counting
 * refactoring of scripts folder
 * renamed summary.xlsx -> predictions_reparation.xlsx
 * removal of unneeded dependencies
 * Refactored Snakefiles

### version 1.2.0 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 22.01.2020

 * updated reparation to v1.0.9
 * added differential expression analysis with Riborex and Xtail
 * added TIS support to generate_excel.py
 * standardized output annotation to GFF3 format
 * updated HRIBO Manual

### version 1.1.0 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 22.08.2019

  * added paired-end support
  * paired-end mapping (currently no tools allowing this exist)
  * combined paired-end mapping
  * improved samples.xlsx (for paired end)
  * improved count_tables (rounding numbers to 2 floating points)
  * fix for annotation generation

### version 1.0.0 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 08.08.2019

 + initial commit
