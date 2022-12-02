# Run results: use python os module
# Date: July 7, 2021

# Script to run this executable: nohup python KEP_RO_Carvalho.py 2>&1 &

import argparse
import multiprocessing
import subprocess
import typing
import os
import os.path

import pydantic
import yaml


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the test for 2StageROThisWork",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument(
        "config_file",
        type=str,
        help="Path to the config file",
    )
    parser.add_argument(
        "--num-processes",
        type=int,
        default=max(1, multiprocessing.cpu_count() - 2),
        help="Number of processes to use",
    )
    parser.add_argument(
        "--dry-run",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    print(f"Loading configuration file {args.config_file}")
    with open(args.config_file, "rt") as f:
        config = RunConfig.parse_obj(yaml.safe_load(f))

    instance_path = {}
    name = None
    for name in config.run_setup.InstanceFolders:  # Folder
        instances_folders = []
        for s in range(0, 30):  # Seeds
            for n in range(0, 10):  # NDDs
                instance_file = os.path.join(
                    config.run_setup.InstancesDir,
                    name,
                    f"CarvalhoRO2021_{s}_{name}_{n}.txt",
                )
                if os.path.isfile(instance_file):
                    instances_folders.append(
                        (instance_file, f"CarvalhoRO2021_{s}_{name}_{n}")
                    )
        assert len(instances_folders) > 0
        instance_path[name] = instances_folders
    assert name is not None, "No instance folders found"
    print(
        f"setup {len(instance_path[name])} instances for {config.run_setup.InstanceFolders} folders"
    )

    commands = []
    sweep_options = [[]]
    for option_list in [
        config.sweep_options.Cycle_Length,
        config.sweep_options.Chain_Length,
        config.sweep_options.VertexBudget,
        config.sweep_options.ArcBudget,
        config.sweep_options.Policy,
        config.sweep_options.Formulation,
        config.sweep_options.LowerBound,
    ]:
        new_options = []
        assert len(option_list) > 0, "Empty sweep options not allowed"
        for option in option_list:
            for prev_lst in sweep_options:
                new_options.append(prev_lst + [str(option)])
        sweep_options = new_options
    print(f"setup {len(sweep_options)} commands:")
    print("\n".join([str(s) for s in sweep_options]))
    os.makedirs(config.run_setup.LogPrintFolder + "/stderr", exist_ok=True)
    os.makedirs(config.run_setup.LogPrintFolder + "/stdout", exist_ok=True)
    os.makedirs(os.path.dirname(config.run_setup.OutputPath), exist_ok=True)

    f_num = 0
    for options_list in sweep_options:
        for n in config.run_setup.InstanceFolders:
            for full_path, file_name in instance_path[n]:
                commands.append(
                    (
                        [
                            config.run_setup.exe,
                            full_path,
                            *options_list,
                            config.run_setup.OutputPath,
                            str(config.run_setup.TimeLimit),
                        ],
                        os.path.join(
                            config.run_setup.LogPrintFolder,
                            "stdout",
                            f"{n}_{file_name}_{'_'.join(options_list)}.txt",
                        ),
                        os.path.join(
                            config.run_setup.LogPrintFolder,
                            "stderr",
                            f"{n}_{file_name}_{'_'.join(options_list)}.txt",
                        ),
                        args.dry_run,
                    )
                )
    print(f"setup {len(commands)} commands")

    with multiprocessing.Pool(args.num_processes) as pool:
        pool.starmap(run_command, commands, chunksize=1)
    print("Done")


def run_command(
    command: typing.List[str],
    output_path_stdout: str,
    output_path_stderr: str,
    dry_run: bool,
) -> None:
    print(
        f"Writing to: {command[-2]}, stdout: {output_path_stdout}, Running command: {' '.join(command)}"
    )
    if not dry_run:
        with open(output_path_stdout, "wb") as fstd, open(
            output_path_stderr, "wb"
        ) as ferr:
            subprocess.run(command, stderr=ferr, stdout=fstd, shell=False)


class SweepOptions(pydantic.BaseModel):
    Cycle_Length: typing.List[int]
    Chain_Length: typing.List[int]
    VertexBudget: typing.List[int]
    ArcBudget: typing.List[int]
    Policy: typing.List[str]
    Formulation: typing.List[str]
    LowerBound: typing.List[str]


class RunSetup(pydantic.BaseModel):
    OutputPath: str
    TimeLimit: int
    InstanceFolders: typing.List[str]
    InstancesDir: str
    exe: str
    LogPrintFolder: str


class RunConfig(pydantic.BaseModel):
    sweep_options: SweepOptions
    run_setup: RunSetup


if __name__ == "__main__":
    main()
