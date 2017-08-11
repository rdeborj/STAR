"""
A Python3 wrapper for the STAR RNASeq aligner.
"""
import yaml

def align_rnaseq(config,
                 star="STAR"):
    """
    Produce the command for aligning FASTQ reads to a reference using the
    STAR aligner.

    USAGE:
        align_rnaseq()

    INPUT:
          * config: dictionary containing STAR options to use in alignment (required)
          * star: full path to the STAR aligner (default: STAR)

    OUTPUT:
        Returns a dictionary containing:
          * cmd: command to execute for alignment
          * output: name of output file alignment is written to

    """
    # flatten the config dictionary into an options string that can be
    # combined with the STAR program to create an executable command line
    # string
    options = ' '.join([
        ' '.join([
            ''.join(['--', key]), str(value)
            ]) for key, value in config.items()
        ])
    cmd = " ".join([star, options])

    return cmd


def get_default_config(file):
    """
    Get the default configuration for the STAR aligner as a dictionary.  Note
    that the values can be overwritten for custom configurations.  The config
    file will be limited to rna->align->star as several other parameters for
    additional tools will be present.

    Usage:
        get_default_config(file)

    Input:
        * file: default configuration file in YAML format containing all STAR
                options (required)

    Output:
        Returns a dictionary containing standard STAR alignment options.

    """
    with open(file, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)

    return config["rna"]["align"]["star"]


def flatten_config_arrays(config):
    """
    Flatten list configuration parameters as a space separated string.

    USAGE:
        object.flatten_config_arrays(
            config
            )

    INPUT:
        * config: YAML configuration option

    OUTPUT:
        Returns the config data structure but with lists/arrays flattened as
        strings.
    """
    for key, value in config.items():
        if isinstance(value, list):
            flattened_value = ' '.join(value)
            config[key] = flattened_value

    return config
