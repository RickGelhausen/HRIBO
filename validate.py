"""
Contains validation scripts necessary for the pipeline.
Author: Rick Gelhausen
"""
import os
import re

def validate_config(conf, unique_conditions):
    """
    Validate the config file.
    """

    # Biological Settings
    if "biologySettings" not in conf:
        raise ValueError("Missing 'biologySettings' section in conf file.")
    else:
        if "adapter" not in conf["biologySettings"]:
            raise ValueError("Missing 'adapter' in 'biologySettings' section in conf file.")
        else:
            if not isinstance(conf["biologySettings"]["adapter"], str):
                raise ValueError("'adapter' in 'biologySettings' section in conf file must be a string.")
            else:
                if conf["biologySettings"]["adapter"].strip() == "":
                    print("WARNING: 'adapter' in 'biologySettings' section in conf file is empty. Skipping adapter trimming.")

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


    print("Config file validated! No errors found.")