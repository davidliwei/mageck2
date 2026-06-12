# Contributing to MAGeCK2

Thank you for your interest in improving MAGeCK2. Contributions of all kinds are welcome — bug
reports, feature requests, documentation, tests, and code.

## Ways to contribute

- **Report a bug or request a feature** by opening an
  [issue](https://github.com/davidliwei/mageck2/issues). For bugs, please include your MAGeCK2
  version, operating system, Python version, the exact command you ran, and the full error output.
  A small example that reproduces the problem is the single most helpful thing you can provide.
- **Ask or answer questions** on the [MAGeCK user forum](https://groups.google.com/g/mageck).
- **Improve documentation** in this repository or in
  [mageck2-doc](https://github.com/davidliwei/mageck2-doc).
- **Contribute code** via a pull request, as described below.

Note: issues about the example datasets or the documentation should be filed in the
[mageck2-demo](https://github.com/davidliwei/mageck2-demo) and
[mageck2-doc](https://github.com/davidliwei/mageck2-doc) repositories, respectively.

## Development setup

```bash
git clone https://github.com/davidliwei/mageck2.git
cd mageck2
pip install -e .
```

MAGeCK2 includes a C component (RRA) that is compiled during installation, so a working C compiler
(e.g. `gcc` or `clang`) and `make` are required.

## Pull request process

1. **Open an issue first** for anything beyond a small fix, so the approach can be discussed before
   you invest significant effort.
2. **Fork the repository** and create a topic branch from `main` (e.g. `fix/count-umi-edgecase`).
3. **Make focused changes.** Keep each pull request to a single logical change where possible.
4. **Add or update tests** that cover your change, and make sure the existing tests pass.
5. **Update documentation** (docstrings, README, or the docs repo) when behavior changes.
6. **Open a pull request** against `main` with a clear description of what changed and why, linking
   the related issue.

A maintainer will review your pull request. Please be responsive to review feedback; small,
well-described pull requests are reviewed fastest.

## Coding guidelines

- Target **Python 3** and keep the public command-line behavior backward compatible unless a change
  is explicitly agreed in an issue.
- Follow [PEP 8](https://peps.python.org/pep-0008/) style. Prefer clear, readable code and add
  comments only where the intent is not obvious from the code itself.
- Match the style and structure of the surrounding code.
- Do not commit generated artifacts, large data files, or editor/build byproducts.

## License

By contributing, you agree that your contributions will be licensed under the same
[BSD 3-Clause License](LICENSE) that covers this project.

## Code of conduct

Participation in this project is governed by the [Code of Conduct](CODE_OF_CONDUCT.md). By
participating, you are expected to uphold it.
