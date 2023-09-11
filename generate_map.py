import sys
import os
import numpy as np
import pandas as pd
import time
from utilities.logo import print_logo

from main import do_barcoding
from merge_vhvl import merge_vhvl
from utilities.arguments import (
    parse_args_and_read_config_mapping,
    parse_config,
    BARCODE_ARGUMENTS,
    ARGUMENTS,
)

if __name__ == "__main__":
    print_logo()
    map_params, configs = parse_args_and_read_config_mapping(sys.argv[1:])

    # do haplotyping for each of the paired configs
    start = time.time()
    for config in configs:
        outputs_1 = []  # VH
        outputs_2 = []  # VL
        print(
            "Haplotyping config: ",
            os.path.join(map_params.config_file_dir, configs[config][0]),
        )
        print(
            "Output directory: ",
            os.path.join(map_params.output_dir),
        )
        params, barcode, tiles, samples_VH, experiments, proteins = parse_config(
            os.path.join(map_params.config_file_dir, configs[config][0]),
            ARGUMENTS,
            BARCODE_ARGUMENTS,
        )

        if not os.path.exists(map_params.output_dir):
            os.makedirs(map_params.output_dir)

        outputs_1.extend(
            [
                os.path.join(map_params.output_dir, sample, "map_success.csv")
                for sample in samples_VH
            ]
        )
        do_barcoding(params, samples_VH, barcode)

        params, barcode, tiles, samples_VL, experiments, proteins = parse_config(
            os.path.join(map_params.config_file_dir, configs[config][1]),
            ARGUMENTS,
            BARCODE_ARGUMENTS,
        )
        print("Finished VH in", round((time.time() - start) / 60, 2), "min")
        start2 = time.time()
        if len(samples_VH) != len(samples_VL):
            raise ValueError("Number of samples in VH and VL configs must be the same")

        outputs_2.extend(
            [
                os.path.join(map_params.output_dir, sample, "map_success.csv")
                for sample in samples_VL
            ]
        )
        do_barcoding(params, samples_VL, barcode)
        print("Finished VL in", round((time.time() - start2) / 60, 2), "min")
        # merge together each config pair from the output files

        for output_1, output_2, sample_VH, sample_VL in zip(
            outputs_1, outputs_2, samples_VH, samples_VL
        ):
            input_df1 = pd.read_csv(output_1)
            input_df2 = pd.read_csv(output_2)
            df, failure_df = merge_vhvl(input_df1, input_df2)
            print("Merging VH and VL for sample: ", sample_VH, sample_VL)
            df.to_csv(
                f"{map_params.output_dir}/{sample_VH}_{sample_VL}_vhvl.csv", index=False
            )
            # failure_df.to_csv(f"{map_params.output_dir}/{sample_VH}_{sample_VL}_vhvl_failures.csv", index=False)
        print("Total time:", round((time.time() - start) / 60, 2), "min")
