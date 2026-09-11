"""Build the contract consumer against embedded and installed mockturtle."""

import argparse
import os
import shutil
import subprocess
import tempfile
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--cmake", default="cmake")
    parser.add_argument("--generator", default="Ninja")
    args = parser.parse_args()

    source = Path(__file__).resolve().parents[2]
    consumer = source / "test/package"
    work = Path(tempfile.mkdtemp(prefix="mockturtle-package-"))
    ctest = Path(shutil.which(args.cmake) or args.cmake).with_name("ctest.exe" if os.name == "nt" else "ctest")

    def run(*command):
        print("+ " + " ".join(map(str, command)), flush=True)
        subprocess.run(list(map(str, command)), check=True)

    def configure(src, build, *options):
        run(args.cmake, "-S", src, "-B", build, "-G", args.generator, "-DCMAKE_BUILD_TYPE=Release", *options)

    def build_and_test(build):
        run(args.cmake, "--build", build, "--config", "Release", "--parallel", "2")
        run(ctest, "--test-dir", build, "-C", "Release", "--output-on-failure")

    # Embedded: mockturtle must not disturb the parent project.
    configure(consumer, work / "embedded", f"-DMOCKTURTLE_SOURCE={source.as_posix()}")
    build_and_test(work / "embedded")

    # Installed: build from a copy, then move the prefix and delete the sources
    # so that nothing can resolve through a build-tree or source-tree path.
    producer = work / "producer-source"
    for directory in ("include", "lib", "cmake"):
        shutil.copytree(source / directory, producer / directory)
    for name in ("CMakeLists.txt", "LICENSE"):
        shutil.copy2(source / name, producer / name)
    prefix = work / "prefix"
    configure(producer, work / "producer-build", "-DMOCKTURTLE_INSTALL=ON",
              "-DMOCKTURTLE_BUILD_EXAMPLES=OFF", f"-DCMAKE_INSTALL_PREFIX={prefix.as_posix()}")
    run(args.cmake, "--build", work / "producer-build", "--config", "Release", "--parallel", "2")
    run(args.cmake, "--install", work / "producer-build", "--config", "Release")
    relocated = work / "relocated"
    prefix.rename(relocated)
    shutil.rmtree(producer)
    shutil.rmtree(work / "producer-build")

    configure(consumer, work / "installed", f"-DCMAKE_PREFIX_PATH={relocated.as_posix()}")
    build_and_test(work / "installed")
    print("Embedded and installed package contracts passed.", flush=True)


if __name__ == "__main__":
    main()
