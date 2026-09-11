"""Exercise embedded and relocated mockturtle CMake consumers without dependencies."""

import argparse
import os
from pathlib import Path
import shutil
import subprocess
import tempfile


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--generator", default="Ninja")
    parser.add_argument("--build-dir", type=Path)
    parser.add_argument("--parallel", default="2")
    parser.add_argument("--define", action="append", default=[])
    args = parser.parse_args()
    source = Path(__file__).resolve().parents[2]
    if args.build_dir:
        args.build_dir.mkdir(parents=True, exist_ok=True)
    work = Path(tempfile.mkdtemp(prefix="mockturtle-cmake-", dir=args.build_dir))
    print(f"Contract artifacts: {work}", flush=True)
    definitions = [f"-D{value}" for value in args.define]
    ctest = str(Path(shutil.which(args.cmake) or args.cmake).with_name("ctest.exe" if os.name == "nt" else "ctest"))

    def run(*command):
        print("+ " + " ".join(map(str, command)), flush=True)
        subprocess.run(list(map(str, command)), check=True)

    def configure(src, build, *options):
        run(args.cmake, "-S", src, "-B", build, "-G", args.generator,
            "-DCMAKE_BUILD_TYPE=Release", *definitions, *options)

    def build_and_test(build):
        run(args.cmake, "--build", build, "--config", "Release", "--parallel", args.parallel)
        run(ctest, "--test-dir", build, "-C", "Release", "--output-on-failure")

    embedded = work / "embedded"
    configure(source / "test/cmake", embedded, f"-DMOCKTURTLE_SOURCE={source.as_posix()}", "-DBASE_ONLY=ON")
    # IPO checks may compile probe archives during configuration.
    configured_archives = {path for path in embedded.rglob("*") if path.suffix.lower() in {".a", ".lib"}}
    build_and_test(embedded)
    archives = [path for path in embedded.rglob("*")
                if path.suffix.lower() in {".a", ".lib"} and path not in configured_archives]
    if archives:
        raise RuntimeError(f"Base-only default build produced archives: {archives}")
    configure(source / "test/cmake", embedded, f"-DMOCKTURTLE_SOURCE={source.as_posix()}", "-DBASE_ONLY=OFF")
    build_and_test(embedded)

    # Build a copy so removing the original package paths cannot affect a checkout.
    producer = work / "producer-source"
    producer.mkdir()
    for directory in ("include", "lib", "cmake"):
        shutil.copytree(source / directory, producer / directory)
    for name in ("CMakeLists.txt", "LICENSE"):
        shutil.copy2(source / name, producer / name)
    install_build = work / "producer-build"
    prefix = work / "prefix"
    configure(producer, install_build, "-DMOCKTURTLE_INSTALL=ON",
              "-DMOCKTURTLE_BUILD_EXAMPLES=OFF", "-DMOCKTURTLE_BUILD_TESTS=OFF",
              "-DMOCKTURTLE_BUILD_EXPERIMENTS=OFF", f"-DCMAKE_INSTALL_PREFIX={prefix.as_posix()}")
    run(args.cmake, "--build", install_build, "--config", "Release", "--parallel", args.parallel)
    run(args.cmake, "--install", install_build, "--config", "Release")
    relocated = work / "relocated"
    prefix.rename(relocated)
    producer.rename(work / "unavailable-source")
    install_build.rename(work / "unavailable-build")
    installed = work / "installed"
    configure(source / "test/cmake", installed, f"-DCMAKE_PREFIX_PATH={relocated.as_posix()}")
    build_and_test(installed)
    if any(value.upper() == "BILL_Z3=ON" for value in args.define):
        configure(source / "test/cmake", work / "missing-z3",
                  f"-DCMAKE_PREFIX_PATH={relocated.as_posix()}", "-DMISSING_Z3=ON")
    print("All embedded and relocated package contracts passed.", flush=True)


if __name__ == "__main__":
    main()
