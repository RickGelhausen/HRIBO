### version 1.4.1 [Rick Gelhausen & Florian Eggenhofer](mailto:gelhausr@informatik.uni-freiburg.de) 20.04.2020
 * added support for multi-condition results for differential expression in the overview table
 
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
