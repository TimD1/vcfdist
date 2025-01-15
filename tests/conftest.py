import os
from pathlib import Path
import pytest
import subprocess


# get .../vcfdist/ from current path: .../vcfdist/tests/conftest.py
VCFDIST_REPO_PATH = Path(__file__).parent.parent.resolve()
os.environ["VCFDIST_REPO_PATH"] = str(VCFDIST_REPO_PATH)


def print_result(result: subprocess.CompletedProcess, process_name: str, padding: int = 30) -> None:
    """Helper function to pretty-print stdout and stderr from a process."""
    print("="*padding + f" {process_name} " + '='*padding)
    print("-"*padding + " stdout " + '-'*padding)
    print(result.stdout.decode())
    print("-"*padding + " stderr " + '-'*padding)
    print(result.stderr.decode())

def build_cpp_project(directory: str) -> None:
    """Helper function to rebuild a C++ project before running PyTest tests."""
    result = subprocess.run(
        f"""
            cd {directory} &&
            make clean &&
            make &&
            cd {VCFDIST_REPO_PATH}/tests
        """,
        capture_output=True,
        shell=True
    )
    print_result(result, f"rebuilding vcfdist C++ project '{directory}'")


def pytest_sessionstart(session):
    """Automatically rebuild vcfdist and Google Test framework before running tests."""
    build_cpp_project(f"{VCFDIST_REPO_PATH}/src")
    build_cpp_project(f"{VCFDIST_REPO_PATH}/tests/unit/build")
