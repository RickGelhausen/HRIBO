# Changelog

## 2.0.0

HRIBO 2.0 is a workflow-only release. It modernizes the Snakemake
layout and environments, adds stage-based execution and strict preflight
validation, repairs scientific and failure-signalling defects, and brings the
maintained documentation into this repository.

Notable changes since 1.8.1 include:

- migrated to the standard `workflow/`, `config/`, `tests/`, and `docs/`
  layout, with portable local and SLURM launchers for Snakemake 9;
- added schema plus semantic input validation and explicit output stages,
  including canonical sample identifiers, full-stream gzip/FASTQ checks, and
  a stage-aware DeepRibo reference-alphabet check;
- corrected start/stop metagene geometry, exact read-length output, raw-count
  P-site markers, CIGAR-block-aware global metagenes and global/centered
  browser tracks, empty-profile handling, RPKM/CPM and track normalization,
  including fractional `1/NH` multi-mapper totals;
- repaired differential contrast handling, current workbook pooling schemas,
  PCA rank handling, and hard failure validation for external engines;
- hardened Reparation, DeepRibo, and deltaTE container boundaries with pinned
  artifacts, atomic output contracts, recovery receipts, and production CI
  smokes;
- made prediction/annotation aggregation deterministic and strict GFF3, exposed
  the combined updated annotation, and declared every overview side product;
- pinned the complete Linux launcher, development, and per-rule environments;
- added a semantic real-data result comparator; consolidated and updated the
  former `HRIBO_ReadTheDocs` material in the maintained Sphinx documentation;
  fully locked the documentation build; and preserved that repository's
  complete history and release tags on the namespaced
  `archive/hribo-readthedocs` ref.

The automated suite covers workflow DAGs, focused executed rules, golden
scientific outputs, crash recovery, and production container boundaries. The
exact current test count is enforced by CI.

## Development history leading to 2.0.0

### Work leading to 2.0.0 — Rick Gelhausen
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
 * Updated the conda environments. Versions were determined from the package
   indexes rather than by solving environments locally, and every Python pin is a
   version the test suite was actually run against
   - pandas 0.23.4/1.4.1/1.5.2 -> 2.3.3, numpy -> 2.5.2, pysam 0.19.1 -> 0.24.0,
     biopython 1.79 -> 1.88, plotly 5.11.0 -> 6.9.0, xlsxwriter 3.0.3 -> 3.2.9,
     openpyxl 3.0.9 -> 3.1.5. The suite passes on pandas 2.2.1/numpy 1.26,
     pandas 2.3.3/numpy 2.5.2 and pandas 3.0.5/numpy 2.5.2
   - samtools 1.9/1.18 -> 1.24, bedtools 2.27.1/2.30.0 -> 2.31.1,
     cutadapt 4.4 -> 5.2, multiqc 1.18 -> 1.35, deeptools 3.2.0 -> 3.5.6,
     subread 2.0.1 -> 2.1.1, the UCSC tools 377 -> 482, bedops 2.4.41 -> 2.4.42,
     gawk 5.0.1 -> 5.4.1, curl 8.5.0 -> 8.21.0, blast 2.15.0 -> 2.17.0,
     pear pinned at 0.9.11, bioconductor-deseq2 1.42.0 -> 1.50.2
   - python-kaleido is deliberately held at 0.2.1: version 1 removed the bundled
     Chrome and requires a browser on the machine, which would break SVG and PDF
     export on a compute node. Verified that plotly 6.9.0 still exports with it
   - xtail moves from R 3.5.1 to its R 4.0 build, the newest that exists
   - riborex cannot be updated at all: both bioconda builds require
     r-base >=3.4.1,<3.4.2, so R 3.4.1 is the only version it installs against
   - reparation stays on Python 3.7 because reparation_blast pins biopython 1.77
     and pysam 0.16; only its blast dependency could be moved
   - segemehl 0.3.4 and fastqc 0.12.1 were already current
   - removed the unused imagemagick, normalization and xtailcounts environments
   - replaced the 184-line frozen environment.yaml export with a readable spec
 * Made every DataFrame sort deterministic. pandas sorts with an unstable
   quicksort by default, so rows with an equal sort key came out in a different
   order depending on the pandas and numpy build; merge_duplicates_deepribo.py
   assigns prediction ranks from that order. Ties are now broken explicitly and
   all sorts use kind="stable"

 * Audited the remaining scripts. Fixes:
   - samples_to_xlsx.py shortened the fastq paths with str[1], the second path
     component, so an absolute path became "data" and a bare file name became
     NaN; and it tested for a column named "Fastqfile2" while the sample sheet
     calls it "fastqFile2", so the second read file was never shortened at all
   - call_featurecounts.py ignored the featureCounts exit code and reused the
     same temporary file for every feature, so a failed run silently re-read the
     previous feature's counts. It also chose the annotation identifier from the
     first row alone, which picks the wrong attribute when the file opens with a
     region or source feature. Temporary files are now cleaned up
   - motif_to_gff.py emitted a stray colon in the forward-strand identifiers
     (ID=chr:1-3:+:) while the reverse strand had none
   - create_reparation_gff.py compared the strand with "is" against a string
     literal, which Python warns about and which is not guaranteed to work
   - mapping.py parsed --clip_length and then passed a hardcoded 11 instead
   - enrich_annotation.py had two bare except clauses that would have swallowed
     unrelated errors
   - preparePCAinput.py printed the whole sample sheet twice as leftover debug
     output
 * read_length_statistics.py now uses the shared plot theme and page shell rather
   than its own inline HTML and size-24 fonts, so it matches the metagene figures

 * Added golden-output tests for the seven annotation-transforming scripts, which
   had none, and consolidated their GFF attribute parsing into gff_utils.py.
   Twelve hand-written copies became three; the remaining three parse differently
   on purpose
 * Fixed reannotate_orfs.py crashing on "ORF_type=;". Reparation emits an empty
   ORF type routinely, and splitting the attributes on both ";" and "=" while
   dropping empty fields left an odd number of items, so rebuilding them pairwise
   ran off the end of the list. That took down the whole workflow at
   reannotatedORFs
 * Fixed enrich_annotation.py raising KeyError on a Parent attribute naming a
   feature absent from the file. It looked the parent up before checking that it
   existed. This one is on the critical path too: every read counting rule depends
   on auxiliary/enriched_annotation.gff
 * Fixed the merged prediction tracks being non-deterministic. The Evidence field
   was joined from a set, and Python randomises string hashing per process, so the
   same input produced a different file on every run
 * Fixed merge_duplicates_deepribo.py filling old_locus_tag from the "locus_tag"
   attribute, so the merged track and everything downstream carried the current
   locus tag where the old one belonged. Also fixed old_locus_tag being read
   without ever being assigned when a feature had no resolvable parent, which
   silently carried the previous row's value

 * Standardised the spreadsheet column headers across every workbook. Same
   information, one naming convention:
   - fixed the misspelled "identifer" header in the DeepRibo table
   - "15nt upstream" -> "Upstream_15nt" and "Feature count" -> "Feature_count",
     so no header contains a space
   - the differential expression tables now lead with Identifier, matching the
     column order of the prediction and annotation tables
   - log2FoldChange -> log2FC, lfcSE -> log2FC_SE and padj -> pvalue_adjusted in
     the riborex and deltaTE tables, which is what xtail and the overview table
     already called them
   - Pred_probability -> Reparation_probability, and Pred_value / Pred_rank ->
     Deepribo_score / Deepribo_rank, matching the overview table
   Verified that every data value is unchanged; only headers were renamed and the
   Identifier column moved

 * Deduplicated the excel generating scripts. generate_excel.py,
   generate_excel_reparation.py and generate_excel_deepribo.py shared roughly 90%
   of their bodies, differing only in which columns they emit; they now declare a
   column list against one shared table builder and shrank from 363 to 185 lines
   combined. The riborex, xtail and deltaTE tables likewise share one builder
 * Collapsed the duplication inside excel_utils: the GFF attribute parsing block
   appeared seven times and is now one function, the three differential
   expression readers became one parameterised by column names, the two
   prediction readers share their row parsing, and generate_annotation_dict and
   generate_non_cds_dict share theirs
 * All of this is verified by golden-output tests: every workbook is byte
   identical to what the previous implementation produced

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

 * Made the pipeline runnable in parts. workflowSettings.workflow, which offered
   the three fixed choices full/preprocessing/trimming, is replaced by
   workflowSettings.stages: a list of the thirteen output groups the run should
   produce (trimming, mapping, qc, tracks, genome_tracks, readcounts, metagene,
   tis_advisor, correlation, pca, predictions, differential_expression,
   overview), or one of the presets "full" and "preprocessing". A stage only says
   what is requested, so everything it depends on is still built; commenting one
   out is what stops a run at, for instance, the BAM files. A single run can
   override the config file with --config stages=mapping,tracks
   - differentialExpressionSettings.differentialExpression and
     tisAdvisorSettings.tisAdvisor are gone; they are now the
     "differential_expression" and "tis_advisor" stages, so there is one place
     that decides what runs rather than two that can contradict each other
   - the preflight only checks the settings a run actually depends on, and warns
     about stages that need Ribo-seq libraries when the sample sheet has none
     instead of failing
   - the shipped config lists the stages explicitly with differential_expression
     commented out, which reproduces the previous default exactly: the target set
     of a default run, of an RNA-only run and of a run with differential
     expression enabled are each unchanged

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
