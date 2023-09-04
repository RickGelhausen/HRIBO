"""
Contains validation scripts necessary for the pipeline.
Author: Rick Gelhausen
"""
import os
import re
import pandas as pd

def validate_config(conf, unique_conditions):
    """
    Validate the config file.
    """

    # Biological Settings
    if "biologySettings" not in conf:
        raise ValueError("Missing 'biologySettings' section in conf file.")
    else:
        if "adapterS3" not in conf["biologySettings"]:
            raise ValueError("Missing 'adapterS3' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["adapterS3"], str):
                raise ValueError("'adapterS3' in 'biologySettings' section in conf file must be a string.")

        if "adapterS5" not in conf["biologySettings"]:
            raise ValueError("Missing 'adapterS5' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["adapterS5"], str):
                raise ValueError("'adapterS5' in 'biologySettings' section in conf file must be a string.")

        if "adapterP3R1" not in conf["biologySettings"]:
            raise ValueError("Missing 'adapterP3R1' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["adapterP3R1"], str):
                raise ValueError("'adapterP3R1' in 'biologySettings' section in conf file must be a string.")

        if "adapterP3R2" not in conf["biologySettings"]:
            raise ValueError("Missing 'adapterP3R2' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["adapterP3R2"], str):
                raise ValueError("'adapterP3R2' in 'biologySettings' section in conf file must be a string.")

        if "adapterP5R1" not in conf["biologySettings"]:
            raise ValueError("Missing 'adapterP5R1' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["adapterP5R1"], str):
                raise ValueError("'adapterP5R1' in 'biologySettings' section in conf file must be a string.")

        if "adapterP5R2" not in conf["biologySettings"]:
            raise ValueError("Missing 'adapterP5R2' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["adapterP5R2"], str):
                raise ValueError("'adapterP5R2' in 'biologySettings' section in conf file must be a string.")

        if "genome" not in conf["biologySettings"]:
            raise ValueError("Missing 'genome' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["genome"], str):
                raise ValueError("'genome' in 'biologySettings' section in conf file must be a string.")
            else:
                if not os.path.isfile(conf["biologySettings"]["genome"]):
                    raise ValueError("'genome' in 'biologySettings' section in conf file must be a valid file path.")

        if "annotation" not in conf["biologySettings"]:
            raise ValueError("Missing 'annotation' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["annotation"], str):
                raise ValueError("'annotation' in 'biologySettings' section in conf file must be a string.")
            else:
                if not os.path.isfile(conf["biologySettings"]["annotation"]):
                    raise ValueError("'annotation' in 'biologySettings' section in conf file must be a valid file path.")

        if "samples" not in conf["biologySettings"]:
            raise ValueError("Missing 'samples' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["samples"], str):
                raise ValueError("'samples' section in conf file must be a string path.")
            else:
                if not os.path.isfile(conf["biologySettings"]["samples"]):
                    raise ValueError("'samples' section in conf file must be a valid file path.")

        if "alternativeStartCodons" not in conf["biologySettings"]:
            raise ValueError("Missing 'alternativeStartCodons' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["alternativeStartCodons"], list):
                raise ValueError("'alternativeStartCodons' in 'biologySettings' section in conf file must be a list.")
            else:
                for start_codon in conf["biologySettings"]["alternativeStartCodons"]:
                    if not isinstance(start_codon, str):
                        raise ValueError("All elements in 'alternativeStartCodons' in 'biologySettings' section in conf file must be strings.")
                    else:
                        if len(start_codon) != 3:
                            raise ValueError("All elements in 'alternativeStartCodons' in 'biologySettings' section in conf file must be 3 characters long.")
                        else:
                            if not re.match("^[ATGC]+$", start_codon):
                                raise ValueError("All elements in 'alternativeStartCodons' in 'biologySettings' section in conf file must be valid DNA sequences.")


    # Differential Expression Settings
    if "differentialExpressionSettings" not in conf:
        raise ValueError("Missing 'differentialExpressionSettings' section in conf file.")
    else:
        if "differentialExpression" not in conf["differentialExpressionSettings"]:
            raise ValueError("Missing 'differentialExpression' in 'differentialExpressionSettings' section in conf file.")
        else:
            if not isinstance(conf["differentialExpressionSettings"]["differentialExpression"], str):
                raise ValueError("'differentialExpression' in 'differentialExpressionSettings' section in conf file must be a string.")
            else:
                if conf["differentialExpressionSettings"]["differentialExpression"].lower() not in ["on", "off"]:
                    raise ValueError("'differentialExpression' in 'differentialExpressionSettings' section in conf file must be either 'on' or 'off'.")

        if "contrasts" not in conf["differentialExpressionSettings"]:
            raise ValueError("Missing 'contrasts' in 'differentialExpressionSettings' section in conf file.")
        else:
            if not isinstance(conf["differentialExpressionSettings"]["contrasts"], list):
                raise ValueError("'contrasts' in 'differentialExpressionSettings' section in conf file must be a list.")
            else:
                if len(conf["differentialExpressionSettings"]["contrasts"]) != 0:
                    for contrast in conf["differentialExpressionSettings"]["contrasts"]:
                        if not isinstance(contrast, str):
                            raise ValueError("All elements in 'contrasts' in 'differentialExpressionSettings' section in conf file must be a string.")
                        else:
                            if not "-" in contrast:
                                raise ValueError("All elements in 'contrasts' in 'differentialExpressionSettings' section in conf file must be in the format 'treated1-untreated1'.")
                            else:
                                if len(contrast.split("-")) != 2:
                                    raise ValueError("All elements in 'contrasts' in 'differentialExpressionSettings' section in conf file must be in the format 'treated1-untreated1'.")
                                else:
                                    if not set(contrast.split("-")).issubset(unique_conditions):
                                        raise ValueError("All elements in 'contrasts' in 'differentialExpressionSettings' section in conf file must be in the format 'treated1-untreated1' and must be a valid condition.")

        if "features" not in conf["differentialExpressionSettings"]:
            raise ValueError("Missing 'features' in 'differentialExpressionSettings' section in conf file.")
        else:
            if not isinstance(conf["differentialExpressionSettings"]["features"], list):
                raise ValueError("'features' in 'differentialExpressionSettings' section in conf file must be a list.")
            else:
                if len(conf["differentialExpressionSettings"]["features"]) != 0:
                    for feature in conf["differentialExpressionSettings"]["features"]:
                        if not isinstance(feature, str):
                            raise ValueError("All elements in 'features' in 'differentialExpressionSettings' section in conf file must be a string.")
                        else:
                            if feature.lower() not in ["cds", "srna"]:
                                print("WARNING: All elements in 'features' in 'differentialExpressionSettings' section in conf file are accepted but make sure they are valid and suitable for differential expression analysis. (suggested features: CDS, sRNA)")
                else:
                    raise ValueError("'features' in 'differentialExpressionSettings' section in conf file must contain atleast one feature (suggested: CDS)")

        if "padjCutoff" not in conf["differentialExpressionSettings"]:
            raise ValueError("Missing 'padjCutoff' in 'differentialExpressionSettings' section in conf file.")
        else:
            if not isinstance(conf["differentialExpressionSettings"]["padjCutoff"], float):
                raise ValueError("'padjCutoff' in 'differentialExpressionSettings' section in conf file must be a float.")
            else:
                if conf["differentialExpressionSettings"]["padjCutoff"] <= 0 or conf["differentialExpressionSettings"]["padjCutoff"] >= 1:
                    raise ValueError("'padjCutoff' in 'differentialExpressionSettings' section in conf file must be between 0 and 1.")

        if "log2fcCutoff" not in conf["differentialExpressionSettings"]:
            raise ValueError("Missing 'log2fcCutoff' in 'differentialExpressionSettings' section in conf file.")
        else:
            if not isinstance(conf["differentialExpressionSettings"]["log2fcCutoff"], float):
                raise ValueError("'log2fcCutoff' in 'differentialExpressionSettings' section in conf file must be a float.")
            else:
                if conf["differentialExpressionSettings"]["log2fcCutoff"] < 0:
                    raise ValueError("'log2fcCutoff' in 'differentialExpressionSettings' section in conf file must be greater than 0. (the cutoff will be adapted for downregulated genes!)")

    if "predictionSettings" not in conf:
        raise ValueError("Missing 'predictionSettings' section in conf file.")
    else:
        if "deepribo" not in conf["predictionSettings"]:
            raise ValueError("Missing 'deepribo' in 'predictionSettings' section in conf file.")
        else:
            if not isinstance(conf["predictionSettings"]["deepribo"], str):
                raise ValueError("'deepribo' in 'predictionSettings' section in conf file must be a string.")
            else:
                if conf["predictionSettings"]["deepribo"].lower() not in ["on", "off"]:
                    raise ValueError("'deepribo' in 'predictionSettings' section in conf file must be either 'on' or 'off'.")

    if "readstatSettings" not in conf:
        raise ValueError("Missing 'readstatSettings' section in conf file.")
    else:
        if "readLengths" not in conf["readstatSettings"]:
            raise ValueError("Missing 'readLengths' in 'readstatSettings' section in conf file.")
        else:
            if not isinstance(conf["readstatSettings"]["readLengths"], str):
                raise ValueError("'readLengths' in 'readstatSettings' section in conf file must be a string.")
            else:
                if conf["readstatSettings"]["readLengths"] == "":
                    raise ValueError("'readLengths' in 'readstatSettings' section in conf file must be a non-empty string.")

    if "metageneSettings" not in conf:
        raise ValueError("Missing 'metageneSettings' section in conf file.")
    else:
        if "positionsOutsideORF" not in conf["metageneSettings"]:
            raise ValueError("Missing 'positionsOutsideORF' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["positionsOutsideORF"], int):
                raise ValueError("'positionsOutsideORF' in 'metageneSettings' section in conf file must be an integer.")
            else:
                if conf["metageneSettings"]["positionsOutsideORF"] < 0:
                    raise ValueError("'positionsOutsideORF' in 'metageneSettings' section in conf file must be a positive integer.")

        if "positionsInORF" not in conf["metageneSettings"]:
            raise ValueError("Missing 'positionsInORF' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["positionsInORF"], int):
                raise ValueError("'positionsInORF' in 'metageneSettings' section in conf file must be an integer.")
            else:
                if conf["metageneSettings"]["positionsInORF"] < 0:
                    raise ValueError("'positionsInORF' in 'metageneSettings' section in conf file must be a positive integer.")

        if "normalizationMethods" not in conf["metageneSettings"]:
            raise ValueError("Missing 'normalizationMethods' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["normalizationMethods"], list):
                raise ValueError("'normalizationMethods' in 'metageneSettings' section in conf file must be a list.")
            else:
                for method in conf["metageneSettings"]["normalizationMethods"]:
                    if not isinstance(method, str):
                        raise ValueError("All elements in 'normalizationMethods' in 'metageneSettings' section in conf file must be a string.")
                    else:
                        if method.lower() not in ["raw", "cpm", "window"]:
                            raise ValueError("All elements in 'normalizationMethod' in 'metageneSettings' section in conf file must be either 'raw', 'cpm', and/or 'window'.")

        if "filteringMethods" not in conf["metageneSettings"]:
            raise ValueError("Missing 'filteringMethods' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["filteringMethods"], list):
                raise ValueError("'filteringMethods' in 'metageneSettings' section in conf file must be a list.")
            else:
                for method in conf["metageneSettings"]["filteringMethods"]:
                    if not isinstance(method, str):
                        raise ValueError("All elements in 'filteringMethods' in 'metageneSettings' section in conf file must be a string.")
                    else:
                        if method.lower() not in ["overlap", "length", "rpkm"]:
                            raise ValueError("All elements in 'filteringMethods' in 'metageneSettings' section in conf file must be either 'overlap', 'length', and/or 'rpkm'.")

        if "neighboringGenesDistance" not in conf["metageneSettings"]:
            raise ValueError("Missing 'neighboringGenesDistance' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["neighboringGenesDistance"], int):
                raise ValueError("'neighboringGenesDistance' in 'metageneSettings' section in conf file must be an integer.")
            else:
                if conf["metageneSettings"]["neighboringGenesDistance"] < 0:
                    raise ValueError("'neighboringGenesDistance' in 'metageneSettings' section in conf file must be a positive integer.")

        if "rpkmThreshold" not in conf["metageneSettings"]:
            raise ValueError("Missing 'rpkmThreshold' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["rpkmThreshold"], float):
                raise ValueError("'rpkmThreshold' in 'metageneSettings' section in conf file must be a float.")
            else:
                if conf["metageneSettings"]["rpkmThreshold"] < 0:
                    raise ValueError("'rpkmThreshold' in 'metageneSettings' section in conf file must be a positive float.")

        if "lengthCutoff" not in conf["metageneSettings"]:
            raise ValueError("Missing 'lengthCutoff' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["lengthCutoff"], int):
                raise ValueError("'lengthCutoff' in 'metageneSettings' section in conf file must be an integer.")
            else:
                if conf["metageneSettings"]["lengthCutoff"] < 0:
                    raise ValueError("'lengthCutoff' in 'metageneSettings' section in conf file must be a positive integer.")

        if "mappingMethods" not in conf["metageneSettings"]:
            raise ValueError("Missing 'mappingMethods' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["mappingMethods"], list):
                raise ValueError("'mappingMethods' in 'metageneSettings' section in conf file must be a list.")
            else:
                for method in conf["metageneSettings"]["mappingMethods"]:
                    if not isinstance(method, str):
                        raise ValueError("All elements in 'mappingMethods' in 'metageneSettings' section in conf file must be a string.")
                    else:
                        if method.lower() not in ["threeprime", "fiveprime", "centered", "global"]:
                            raise ValueError("All elements in 'mappingMethods' in 'metageneSettings' section in conf file must be either 'threeprime', 'fiveprime', 'centered', and/or 'global'.")

        if "readLengths" not in conf["metageneSettings"]:
            raise ValueError("Missing 'readLengths' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["readLengths"], str):
                raise ValueError("'readLengths' in 'metageneSettings' section in conf file must be a string.")
            else:
                if conf["metageneSettings"]["readLengths"] == "":
                    raise ValueError("'readLengths' in 'metageneSettings' section in conf file must be a non-empty string.")

        if "normalizationMethods" not in conf["metageneSettings"]:
            raise ValueError("Missing 'normalizationMethods' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["normalizationMethods"], list):
                raise ValueError("'normalizationMethods' in 'metageneSettings' section in conf file must be a list.")
            else:
                for method in conf["metageneSettings"]["normalizationMethods"]:
                    if not isinstance(method, str):
                        raise ValueError("All elements in 'normalizationMethods' in 'metageneSettings' section in conf file must be a string.")
                    else:
                        if method.lower() not in ["raw", "cpm", "window"]:
                            raise ValueError("All elements in 'normalizationMethods' in 'metageneSettings' section in conf file must be either 'raw', 'cpm', and/or 'window'.")

        if "outputFormats" not in conf["metageneSettings"]:
            raise ValueError("Missing 'outputFormats' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["outputFormats"], list):
                raise ValueError("'outputFormats' in 'metageneSettings' section in conf file must be a list.")
            else:
                for method in conf["metageneSettings"]["outputFormats"]:
                    if not isinstance(method, str):
                        raise ValueError("All elements in 'outputFormats' in 'metageneSettings' section in conf file must be a string.")
                    else:
                        if method.lower() not in ["interactive", "svg", "png", "pdf", "jpg"]:
                            raise ValueError("All elements in 'outputFormats' in 'metageneSettings' section in conf file must be either 'interactive', 'svg', 'png', 'pdf', and/or 'jpg'.")

        if "includePlotlyJS" not in conf["metageneSettings"]:
            raise ValueError("Missing 'includePlotlyJS' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["includePlotlyJS"], str):
                raise ValueError("'includePlotlyJS' in 'metageneSettings' section in conf file must be a string.")
            else:
                if conf["metageneSettings"]["includePlotlyJS"] not in ["integrated", "online", "local"]:
                    raise ValueError("'includePlotlyJS' in 'metageneSettings' section in conf file must be either 'integrated', 'online', or 'local'.")

        if "colorList" not in conf["metageneSettings"]:
            raise ValueError("Missing 'colorList' in 'metageneSettings' section in conf file.")
        else:
            if not isinstance(conf["metageneSettings"]["colorList"], list):
                raise ValueError("'colorList' in 'metageneSettings' section in conf file must be a list.")

    print("Config file validated! No errors found!")


def validate_sample_sheet(samples):
    """
    Validate the contents of the sample sheet
    """

    # Check if the sample sheet is empty
    if len(samples) == 0:
        raise ValueError("Sample sheet is empty!")

    # Check if the sample sheet contains the correct columns
    required_columns = ["method", "condition", "replicate", "fastqFile", "fastqFile2"]
    for column in required_columns:
        if column not in samples.columns:
            raise ValueError(f"Missing column {column} in sample sheet!")

    # Check method
    try:
        if not samples["method"].astype(str).isin(["RIBO", "TIS", "TTS", "RNA"]).all():
            raise ValueError("Column 'method' in sample sheet must contain either 'RIBO', 'TIS', 'TTS', or 'RNA'!")
    except Exception as e:
        raise ValueError(f"Error while validating 'method': {e}")

    # Check condition
    try:
        if not samples["condition"].astype(str).str.isalnum().all():
            raise ValueError("Column 'condition' in sample sheet must contain only alphanumeric characters!")
    except Exception as e:
        raise ValueError(f"Error while validating 'condition': {e}")

    # Check replicate
    try:
        if samples["replicate"].astype(int).lt(1).any():
            raise ValueError("Column 'replicate' in sample sheet must contain only positive integers!")
    except Exception as e:
        raise ValueError(f"Error while validating 'replicate': {e}")

    # Check fastqFile
    try:
        if not samples["fastqFile"].astype(str).apply(os.path.isfile).all():
            raise ValueError("Column 'fastqFile' in sample sheet must contain valid file paths!")
    except Exception as e:
        raise ValueError(f"Error while validating 'fastqFile': {e}")

    # Check fastqFile2
    try:
        non_empty_values = samples["fastqFile2"].astype(str).str.strip().ne("").fillna(True)
        non_na_values = samples["fastqFile2"].astype(str).apply(lambda x: not pd.isna(x))
        valid_values = non_empty_values & non_na_values
        if not valid_values.all():
            if not samples.loc[valid_values, "fastqFile2"].apply(os.path.isfile).all():
                raise ValueError(f"Column 'fastqFile2' in sample sheet must contain valid file paths or be NaN/empty! Current value: {samples['fastqFile2']}")
    except Exception as e:
        raise ValueError(f"Error while validating 'fastqFile2': {e}")

    print("Sample sheet validated! No errors found!")