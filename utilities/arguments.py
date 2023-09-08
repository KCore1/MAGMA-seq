import argparse
import ast
import configparser

from utilities.tile import Tile


def bounded_number(converter, low=None, high=None):
    def f(s, converter=converter):
        try:
            x = converter(s)
        except:
            raise argparse.ArgumentTypeError(f"Failed to convert '{s}'")
        if (low is not None and x < low) or (high is not None and x > high):
            raise argparse.ArgumentTypeError(f"Invalid value: {x}")
        return x

    return f


non_negative_int = bounded_number(int, low=0)


def yes_or_no(s):
    s = s.lower()
    if s in ["true", "t", "yes", "y"]:
        return True
    elif s in ["false", "f", "no", "n"]:
        return False
    else:
        raise argparse.ArgumentTypeError(
            f"Invalid value: {s}. Must be one of: True, T, Yes, Y,"
            " False, F, No, N. (Not case sensitive.)"
        )


def maybe_quoted_string(s):
    if len(s) > 0 and s[0] in ["'", '"']:
        return ast.literal_eval(s)
    else:
        return s


CLI_ONLY_ARGUMENTS = {
    "config": {"help": "Configuration file to use.", "type": str, "required": True},
}

ARGUMENTS = {
    "use_multiprocessing": {
        "help": (
            "Should the multiprocessing module be used to process FASTQ"
            " files using multiple cores?"
        ),
        "type": yes_or_no,
        "default": False,
    },
    "fastq_file_dir": {
        "help": "Directory where FASTQ files are located.",
        "type": maybe_quoted_string,
        "default": "",
    },
    "output_dir": {
        "help": "Directory to place output files.",
        "type": maybe_quoted_string,
        "default": "",
    },
    "max_mismatches": {
        "help": (
            "Maximum number of mismatches a read can have" " without being discarded."
        ),
        "type": non_negative_int,
        "default": None,
    },
    "min_quality": {
        "help": "Reads with any quality score lower than this are discarded.",
        "type": non_negative_int,
        "default": None,
    },
    "gene_min_quality": {
        "help": (
            "Reads with any quality score lower than this in the gene"
            " region are discarded."
        ),
        "type": non_negative_int,
        "default": None,
    },
    "min_ref_counts": {
        "help": (
            "Variants with fewer than this number of reference counts"
            " will be discarded."
        ),
        "type": non_negative_int,
        "default": 5,
    },
    "pseudocount": {
        "help": (
            "Pseudocount to use for mutations that were observed in the"
            " reference but not the selected population."
        ),
        "type": non_negative_int,
        "default": 1,
    },
}

MATCH_ARGUMENTS = {
    "use_multiprocessing": {
        "help": (
            "Should the multiprocessing module be used to process FASTQ"
            " files using multiple cores?"
        ),
        "type": yes_or_no,
        "default": False,
    },
    "fastq_file_dir": {
        "help": "FASTQ (can be gzipped) file directory with samples from all label concentrations and bins",
        "type": maybe_quoted_string,
        "default": "",
    },
    "map_file_path": {
        "help": "Path to barcode map file",
        "type": maybe_quoted_string,
        "default": "",
    },
    "output_dir": {
        "help": "Directory to place output files.",
        "type": maybe_quoted_string,
        "default": "",
    },
    "barcode_start": {
        "help": "First index of barcode",
        "type": int,
        "default": 0,
    },
    "template": {
        "help": "Template for the barcode",
        "type": maybe_quoted_string,
        "default": "",
    },
    "limit_csv_path": {
        "help": "Path to csv file containing high and low limits for each concentration and bin",
        "type": maybe_quoted_string,
        "default": "",
    },
    "B": {
        "help": "Variant independent autofluorescence",
        "type": float,
        "default": 350,
    },
    "sigma": {
        "help": "Standard deviation of the clonal population",
        "type": float,
        "default": 1.00,
    },
    "gene": {
        "help": "Name of the gene to analyze. This should be the same as what is before ",
        "type": maybe_quoted_string,
        "default": "",
    },
}

BARCODE_ARGUMENTS = {
    "amplen": {
        "help": "Length of the amplicon.",
        "type": non_negative_int,
        "default": 0,
    },
    "barcode": {
        "help": "Barcode sequence encoded in mixed base notation.",
        "type": maybe_quoted_string,
    },
    "barcode_start": {
        "help": "Index of first base of barcode.",
        "type": non_negative_int,
        "default": 0,
    },
    "barcode_end": {
        "help": "Index of last base of barcode.",
        "type": non_negative_int,
        "default": 1,
    },
    "barcode_min_quality": {
        "help": "Minimum quality score for barcode bases.",
        "type": non_negative_int,
        "default": 15,
    },
    "merge_only": {
        "help": "Should only the merging step be performed?",
        "type": yes_or_no,
        "action": "store_true",
    },
    "vh_or_vl": {
        "help": "VH or VL gene",
        "type": maybe_quoted_string,
        "default": "VH",
    },
    "gene_start": {
        "help": "Index of first base of gene.",
        "type": non_negative_int,
    },
    "gene_end": {
        "help": "Index of last base of gene.",
        "type": non_negative_int,
    },
    "wt_csv": {
        "help": "File that contains WT V_H and V_L genes.",
        "type": maybe_quoted_string,
    },
    "count_threshold": {
        "help": "Minimum number of reads for a variant to be considered.",
        "type": non_negative_int,
        "default": 0,
    },
    "frequency_threshold": {
        "help": "Minimum frequency for a variant to be considered.",
        "type": float,
        "default": 0.0,
    },
    "mutations": {
        "help": 'Consider all mutations or only non-silent mutations. This should be "all" or "non-silent".',
        "type": maybe_quoted_string,
        "default": "all",
    },
    "encoded_positions": {
        "help": "Remove positions of mutations NOT encoded in the library based on a list of mutations encoded in the WT CSV.",
        "type": yes_or_no,
        "default": False,
    },
}

MAPPING_ARGUMENTS = {
    "use_multiprocessing": {
        "help": (
            "Should the multiprocessing module be used to process FASTQ"
            " files using multiple cores?"
        ),
        "type": yes_or_no,
        "default": False,
    },
    "config_file_dir": {
        "help": "Directory where config files are located.",
        "type": maybe_quoted_string,
        "default": "",
    },
    "output_dir": {
        "help": "Directory to place output files.",
        "type": maybe_quoted_string,
        "default": "",
    },
}

PARAMS_NAME = "Parameters"
TILE_PREFIX = "Tile:"
SAMPLES_NAME = "Samples"
EXPERIMENTS_NAME = "Experiments"
PROTEINS_NAME = "Proteins"
BARCODE_NAME = "Barcode"
CONFIGS_NAME = "Configs"


def parse_params(arguments, config):
    params = argparse.Namespace()
    if not config.has_section(PARAMS_NAME):
        return params
    for param in config.options(PARAMS_NAME):
        if param not in arguments:
            raise ValueError(f"unknown option: {param}")
        try:
            raw = config.get(PARAMS_NAME, param, fallback=arguments[param]["default"])
            value = arguments[param]["type"](raw)
        except argparse.ArgumentTypeError as e:
            print(f"Invalid value for option {param}: {raw}")
            raise e
        setattr(params, param, value)
    for argument in arguments:
        if argument not in config.options(PARAMS_NAME):
            setattr(params, argument, arguments[argument]["default"])
    return params


def parse_tile(config, section):
    kwargs = dict(
        wt_seq=ast.literal_eval(config.get(section, "wt_seq")),
        cds_start=config.getint(section, "cds_start"),
        cds_end=config.getint(section, "cds_end"),
        first_aa=config.getint(section, "first_aa"),
    )
    if config.has_option(section, "positions"):
        kwargs["positions"] = ast.literal_eval(config.get(section, "positions"))
    return Tile(**kwargs)


def parse_tiles(config):
    tiles = {}
    for section in config.sections():
        if section.startswith(TILE_PREFIX):
            tiles[section[len(TILE_PREFIX) :]] = parse_tile(config, section)
    if len(tiles) == 0:
        if not config.has_section(BARCODE_NAME):
            raise argparse.ArgumentTypeError("config does not define any tiles.")
    return tiles


def parse_samples(config, tiles=None):  # Tiles stuff is not fully tested yet
    if not config.has_section(SAMPLES_NAME):
        raise argparse.ArgumentTypeError("config does not have a [Samples]" " section.")
    samples = {}
    for name, value in config.items(SAMPLES_NAME):
        elements = ast.literal_eval(value)
        if len(elements) != 3:
            if not config.has_section(BARCODE_NAME) or tiles:
                raise ValueError(
                    f"sample {name} does not specify a tile and two" f" files."
                )
            elif len(elements) != 2:
                raise ValueError(f"sample {name} does not specify two files.")
        if not config.has_section(BARCODE_NAME) or tiles:
            tile, file1, file2 = elements
        else:
            file1, file2 = elements
        if not config.has_section(BARCODE_NAME) and tile not in tiles:
            raise ValueError(f"sample {name} specifies undefined tile {tile}.")
        if not config.has_section(BARCODE_NAME) or tiles:
            samples[name] = (tile, (file1, file2))
        else:
            samples[name] = (file1, file2)
    return samples


def parse_samples_match(config):  # matching not tested FIXME
    if not config.has_section(SAMPLES_NAME):
        raise argparse.ArgumentTypeError("config does not have a [Samples]" " section.")
    samples = {}
    for name, value in config.items(SAMPLES_NAME):
        elements = ast.literal_eval(value)
        if isinstance(elements, str):
            file1 = elements
            samples[name] = file1
        else:
            if len(elements) != 2:
                raise ValueError(f"sample {name} does not specify a file.")
            file1 = elements[0]
            samples[name] = file1
    return samples


def parse_configs(config):
    if not config.has_section(CONFIGS_NAME):
        raise argparse.ArgumentTypeError("config does not have a [Configs]" " section.")
    configs = {}
    for name, value in config.items(CONFIGS_NAME):
        elements = ast.literal_eval(value)
        if len(elements) != 2:
            raise ValueError(
                f"[Configs] section {name} does not specify both a VH and VL config."
            )
        file1, file2 = elements
        configs[name] = (file1, file2)
    return configs


def parse_experiments(config, samples):
    if not config.has_section(EXPERIMENTS_NAME):
        raise argparse.ArgumentTypeError(
            "config does not have an [Experiments]" " section."
        )
    experiments = {}
    dms_flag = False
    for section in config.sections():
        if section.startswith(TILE_PREFIX):
            dms_flag = True
    for name, value in config.items(EXPERIMENTS_NAME):
        pair = ast.literal_eval(value)
        pair = (pair,) if type(pair) != tuple else pair
        if dms_flag:
            if len(pair) != 2:
                raise ValueError(f"experiment {name} does not specify two" f" samples.")
            if samples[pair[0]][0] != samples[pair[1]][0]:
                raise ValueError(
                    f"experiment {name} specifies samples with" " different tiles."
                )
        else:
            for sample in pair:
                if sample not in samples:
                    raise ValueError(f"sample {pair} is not defined.")
        experiments[name] = tuple(pair)
    return experiments


def parse_proteins(config, tiles, samples, experiments):
    if not config.has_section(PROTEINS_NAME):
        raise argparse.ArgumentTypeError(
            "config does not have a [Proteins]" " section."
        )
    proteins = {}
    for name, value in config.items(PROTEINS_NAME):
        exps = ast.literal_eval(value)
        if len(exps) == 0:
            raise ValueError(f"protein {name} does not specify any" " experiments.")
        proteins[name] = tuple(exps)
    if len(proteins) == 0:
        raise ValueError("config does not specify any proteins.")
    return proteins


def parse_barcode(arguments, config):
    barcode_params = argparse.Namespace()
    if not config.has_section(BARCODE_NAME):
        return None
    for param in config.options(BARCODE_NAME):
        if param not in arguments:
            raise ValueError(f"unknown option: {param}")
        try:
            raw = config.get(BARCODE_NAME, param)
            value = arguments[param]["type"](raw)
        except argparse.ArgumentTypeError as e:
            print(f"Invalid value for option {param}: {raw}")
            raise e
        setattr(barcode_params, param, value)
    # Setting defaults, does not appear to be possible using both configparser and argparse as of 10/27/2022
    if not config.has_option(BARCODE_NAME, "merge_only"):
        setattr(barcode_params, "merge_only", False)
    if not config.has_option(BARCODE_NAME, "barcode"):
        setattr(barcode_params, "barcode", "")
    if not config.has_option(BARCODE_NAME, "count_threshold"):
        setattr(barcode_params, "count_threshold", 0)
    if not config.has_option(BARCODE_NAME, "frequency_threshold"):
        setattr(barcode_params, "frequency_threshold", 0.0)
    if not config.has_option(BARCODE_NAME, "mutations"):
        setattr(barcode_params, "mutations", "all")
    if not config.has_option(BARCODE_NAME, "encoded_positions"):
        setattr(barcode_params, "encoded_positions", False)
    return barcode_params


def parse_config(path, args, barcode_args=None):
    """Parse a configuration file.

    Returns (params, barcode, tiles, samples, experiments). params is a dict
    mapping param_name -> value. barcode is a dict mapping barcode_param_name -> value. tiles is a dict mapping tile_name ->
    Tile. samples is a dict mapping sample_name -> (tile_name, (path1,
    path2)), where path1 and path2 are the forward and reverse
    paired-end reads for a single sample. experiments is a dict
    mapping experiment_name -> (ref_sample_name, sel_sample_name).
    """
    config = configparser.ConfigParser(strict=True)
    config.optionxform = str  # make the parser case-sensitive
    config.read(path)
    if len(config.sections()) == 0:
        raise ValueError(
            "config file does not contain any sections. Please verify that config exists and is not empty."
        )
    params = parse_params(args, config)
    barcode = parse_barcode(barcode_args, config)
    tiles = parse_tiles(config)
    samples = parse_samples(config, tiles)
    experiments = parse_experiments(config, samples)
    proteins = parse_proteins(config, tiles, samples, experiments)
    return params, barcode, tiles, samples, experiments, proteins


def parse_config_match(path, args=MATCH_ARGUMENTS):
    """Parse a configuration file, overloaded for matching config file.

    Returns (params, samples, experiments). params is a dict
    mapping param_name -> value. samples is a dict mapping sample_name -> (tile_name, (path1,
    path2)) (NOTE: remove tile from here if it works later), where path1 and path2 are the forward and reverse
    paired-end reads for a single sample. experiments is a dict
    mapping experiment_name -> (ref_sample_name, sel_sample_name).
    """
    config = configparser.ConfigParser(strict=True)
    config.optionxform = str  # make the parser case-sensitive
    config.read(path)
    if len(config.sections()) == 0:
        raise ValueError(
            "config file does not contain any sections. Please verify that config exists and is not empty."
        )
    params = parse_params(args, config)
    samples = parse_samples_match(config)
    return params, samples  # , experiments


def parse_config_mapping(path, args=MAPPING_ARGUMENTS):
    config = configparser.ConfigParser(strict=True)
    config.optionxform = str  # make the parser case-sensitive
    config.read(path)
    if len(config.sections()) == 0:
        raise ValueError(
            "config file does not contain any sections. Please verify that config exists and is not empty."
        )
    params = parse_params(args, config)
    configs = parse_configs(config)
    return params, configs


def make_arg_parser(cli_only_args, args):
    parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in cli_only_args.items():
        parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    for name, kwargs in args.items():
        parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    return parser


def parse_args_and_read_config(
    argv,
    cli_only_arguments=CLI_ONLY_ARGUMENTS,
    arguments=ARGUMENTS,
    barcode_arguments=BARCODE_ARGUMENTS,
):
    # First, parse CLI-only arguments.
    cli_only_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in cli_only_arguments.items():
        cli_only_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    namespace, remaining_argv = cli_only_parser.parse_known_args(argv)
    # Next, get the configuration file and parse that.
    params, barcode, tiles, samples, experiments, proteins = parse_config(
        namespace.config, arguments, barcode_arguments
    )
    # Finally, read everything else from the commandline, overriding
    # anything also specified in the config file Parameters block.
    arg_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in arguments.items():
        arg_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    arg_parser.parse_args(remaining_argv, params)
    return params, barcode, tiles, samples, experiments, proteins


def parse_args_and_read_config_match(
    argv, cli_only_arguments=CLI_ONLY_ARGUMENTS, match_arguments=MATCH_ARGUMENTS
):
    """Parse command-line arguments and read config file, overloaded for matching config file."""
    # First, parse CLI-only arguments.
    cli_only_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in cli_only_arguments.items():
        cli_only_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    namespace, remaining_argv = cli_only_parser.parse_known_args(argv)
    # Next, get the configuration file and parse that.
    params, samples = parse_config_match(
        namespace.config, match_arguments
    )  # , experiments
    # Finally, read everything else from the commandline, overriding
    # anything also specified in the config file Parameters block.
    arg_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in match_arguments.items():
        arg_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    arg_parser.parse_args(remaining_argv, params)
    return params, samples  # , experiments


def parse_args_and_read_config_mapping(
    argv, cli_only_arguments=CLI_ONLY_ARGUMENTS, arguments=MAPPING_ARGUMENTS
):
    """Parse command-line arguments and read config file, overloaded for mapping config file."""
    # First, parse CLI-only arguments.
    cli_only_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in cli_only_arguments.items():
        cli_only_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    namespace, remaining_argv = cli_only_parser.parse_known_args(argv)
    # Next, get the configuration file and parse that.
    params, configs = parse_config_mapping(namespace.config, arguments)  # , experiments
    # Finally, read everything else from the commandline, overriding
    # anything also specified in the config file Parameters block.
    arg_parser = argparse.ArgumentParser(allow_abbrev=False)
    for name, kwargs in arguments.items():
        arg_parser.add_argument(f'--{name.replace("_", "-")}', **kwargs)
    arg_parser.parse_args(remaining_argv, params)
    return params, configs
