"""
A Python3 wrapper for the STAR RNASeq aligner.
"""

import os
import sys
import yaml


def align_rnaseq(config, star="STAR"):
    """
    Produce the command for aligning FASTQ reads to a reference using the
    STAR aligner.

    Usage:
        align_rnaseq()

    Input:
          * config: dictionary containing STAR options to use in alignment (required)
          * star: full path to the STAR aligner (default: STAR)

    Output:
        Returns a dictionary containing:
          * cmd: command to execute for alignment
          * output: name of output file alignment is written to

    """
    # flatten the config dictionary into an options string that can be
    # combined with the STAR program to create an executable command line
    # string
    options = ' '.join([' '.join([''.join(['--', key]), str(value)]) for key, value in config.items()])
    cmd = " ".join([star, options])
    return cmd


def get_default_config(file):
    """
    Get the default configuration for the STAR aligner as a dictionary.  Note
    that the values can be overwritten for custom configurations.

    Usage:
        get_default_config(file)

    Input:
        * file: default configuration file in YAML format containing all STAR
                options (required)

    Output:
        Returns a dictionary containing standard STAR options

    """
    with open(file, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return config
