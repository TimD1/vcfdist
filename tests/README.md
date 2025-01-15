## Testing
Tests are run using `pytest`, `pytest-workflow`, and Google's C++ testing framework.
From the `vcfdist/test` directory, `pytest -vv` will run all tests (with verbose output).
The GRCh38 reference `GCA_000001405.15_GRCh38_no_alt_analysis_set.fasta` is required to run several tests, and should be located in the `vcfdist/data/refs/` directory.
