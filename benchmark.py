import pyperf
import main
import utilities.haplotype as haplotype
import merge_vhvl
import scan_barcode
import utilities.arguments
import os
import pandas as pd
import argparse
import matplotlib.pyplot as plt
import pylab
import scipy.stats as stats


def display_histogram_scipy(bench, mean, bins):
    values = bench.get_values()
    values = sorted(values)
    name = bench.get_name()

    if mean:
        fit = stats.norm.pdf(values, bench.mean(), bench.stdev())
        pylab.plot(values, fit, "-o", label="mean-stdev")
    else:
        fit = stats.norm.pdf(values, bench.mean(), bench.stdev())
        pylab.plot(values, fit, "-o", label="mean-stdev")

    plt.legend(loc="upper right", shadow=True, fontsize="x-large")
    pylab.hist(values, bins=bins)
    pylab.title(name)
    # pylab.show()
    print(name)
    pylab.savefig(f"./figures/{name}.png")


def plot_benchmark(args: dict = None):
    if args == None:
        parser = argparse.ArgumentParser()
        parser.add_argument(
            "-n",
            "--bins",
            type=int,
            default=25,
            help="Number of histogram bars (default: 25)",
        )
        parser.add_argument(
            "--mean", action="store_true", help="Use mean-stdev, instead of median-mad"
        )
        parser.add_argument("-b", "--benchmark")
        parser.add_argument("filename")
        args = parser.parse_args()

        if args.benchmark:
            suite = pyperf.BenchmarkSuite.load(args.filename)
            bench = suite.get_benchmark(args.benchmark)
        else:
            bench = pyperf.Benchmark.load(args.filename)

        display_histogram_scipy(bench, args.mean, args.bins)
    else:
        bench = pyperf.Benchmark.load(args["filename"])
        display_histogram_scipy(bench, args["mean"], args["bins"])


def bench_write_merged_reads(runner: pyperf.Runner, argv):
    (
        params,
        barcode,
        tiles,
        samples,
        experiments,
        proteins,
    ) = arguments.parse_args_and_read_config(argv)
    # get arguments
    inputs = []
    for sample, filenames in samples.items():
        path1, path2 = [os.path.join(params.fastq_file_dir, f) for f in filenames]
        inputs.append(
            (
                f"{params.output_dir}/{sample}/merged.csv",
                path1,
                path2,
                barcode.wt_csv,
                barcode.vh_or_vl,
                params.max_mismatches,
                params.min_quality,
                barcode.barcode_min_quality,
            )
        )  # benchmarking VH only in benchmark_haplotype.config
    inputs = inputs[0]
    inputs = [*inputs]
    # args = runner.parse_args(['--track-memory']) # remember to include --track-memory # FIXME: Broken for now
    runner.metadata["description"] = "Benchmark write_merged_reads"
    bench = runner.bench_func("merge_reads", haplotype.write_merged_reads, *inputs)
    return bench


def bench_barcode_to_variant_map(runner: pyperf.Runner, argv):
    params, barcode, _, samples, _, _ = arguments.parse_args_and_read_config(argv)
    # get arguments
    inputs = []
    sample = list(samples.keys())[0]
    df = pd.read_csv(f"{params.output_dir}/{sample}/merged.csv")
    reads = df.loc[:, "Seq"]
    amplens = df.loc[:, "Amplen"]
    inputs = [reads, amplens, barcode.wt_csv, barcode]
    # args = runner.parse_args(['--track-memory']) # remember to include --track-memory # FIXME: Broken for now
    runner.metadata["description"] = "Benchmark barcode_to_variant_map"
    bench = runner.bench_func(
        "barcode_to_variant_map", haplotype.barcode_to_variant_map, *inputs
    )
    return bench


def bench_choose_variant(argv):  # FIXME not sure how to do this yet
    params, barcode, _, samples, _, _ = arguments.parse_args_and_read_config(argv)
    # get arguments
    inputs = []
    for sample, filenames in samples.items():
        path1, path2 = [os.path.join(params.fastq_file_dir, f) for f in filenames]
        inputs.append(
            (
                f"{params.output_dir}/{sample}/merged.csv",
                path1,
                path2,
                barcode.wt_csv,
                barcode.vh_or_vl,
                params.max_mismatches,
                params.min_quality,
                barcode.barcode_min_quality,
            )
        )  # benchmarking VH only in benchmark_haplotype.config
    inputs = inputs[0]
    inputs = [*inputs]
    runner = pyperf.Runner()
    runner.args.track_memory = True
    runner.metadata["description"] = "Benchmark write_merged_reads"

    return runner.bench_func("merge_reads", haplotype.write_merged_reads, *inputs)


def bench_merge_vhvl(argv):  # FIXME
    params, barcode, _, samples, _, _ = arguments.parse_args_and_read_config(argv)
    # get arguments
    inputs = []
    for sample, filenames in samples.items():
        path1, path2 = [os.path.join(params.fastq_file_dir, f) for f in filenames]
        inputs.append(
            (
                f"{params.output_dir}/{sample}/merged.csv",
                path1,
                path2,
                barcode.wt_csv,
                barcode.vh_or_vl,
                params.max_mismatches,
                params.min_quality,
                barcode.barcode_min_quality,
            )
        )  # benchmarking VH only in benchmark_haplotype.config
    inputs = inputs[0]
    inputs = [*inputs]
    runner = pyperf.Runner(processes=3)
    runner.args.track_memory = True
    runner.metadata["description"] = "Benchmark write_merged_reads"

    return runner.bench_func("merge_reads", haplotype.write_merged_reads, *inputs)


def bench_scan_barcode(runner: pyperf.Runner, argv):  # FIXME
    runner.metadata["description"] = "Benchmark scan_barcode main function"

    bench = runner.bench_func("scan_barcode", scan_barcode.main, argv)
    return bench


if __name__ == "__main__":
    #     runner = pyperf.Runner()
    #     bench = bench_write_merged_reads(runner, ['--config', './configs/benchmark_haplotype.config'])
    #     bench.dump('bench_write_merged_reads.json', replace=True)

    #     bench = bench_barcode_to_variant_map(runner, ['--config', './configs/benchmark_haplotype.config'])
    #     bench.dump('bench_barcode_to_variant_map.json', replace=True)
    #
    #     # bench = bench_choose_variant(['--config', './configs/benchmark_haplotype.config'])
    #     # bench.dump('bench_choose_variant.json', replace=True)
    #     # args = {'filename': 'bench_choose_variant.json', 'mean': True, 'bins': 25}
    #     # plot_benchmark(args)
    #     # bench = bench_merge_vhvl(['--config', './configs/benchmark_haplotype.config'])
    #     # bench.dump('bench_merge_vhvl.json', replace=True)
    #     # args = {'filename': 'bench_merge_vhvl.json', 'mean': True, 'bins': 25}
    #     # plot_benchmark(args)
    #     bench = bench_scan_barcode(runner, ['--config', './configs/benchmark_matching.config'])
    #     bench.dump('bench_scan_barcode.json', replace=True)
    args = {"filename": "bench_write_merged_reads.json", "mean": True, "bins": 25}
    plot_benchmark(args)
    args = {"filename": "bench_barcode_to_variant_map.json", "mean": True, "bins": 25}
    plot_benchmark(args)
#     args = {'filename': 'bench_scan_barcode.json', 'mean': True, 'bins': 25}
#     plot_benchmark(args)
