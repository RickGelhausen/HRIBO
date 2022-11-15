"""
Contains validation scripts necessary for the pipeline.
Author: Rick Gelhausen
"""
import os
import re

def validate_config(conf):
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


    print("Config file validated! No errors found.")