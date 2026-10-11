# Contributing

Thanks for your interest. This is an exploratory re-analysis, so corrections and well-supported improvements are welcome.

1. Fork the repository and create a branch (`git checkout -b fix/short-description`).
2. Install dependencies: `pip install -r requirements.txt`.
3. Make your change, then run `make test` (rebuilds the processed data and results, then runs `pytest`). All tests must pass.
4. If you change a result, figure or table, update the numbers quoted in `README.md` and `docs/`.
5. Open a pull request describing what changed and why.

Please keep the statement that results are **exploratory and TPM-based** (no raw counts) wherever results are reported.
