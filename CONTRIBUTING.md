# Contributing to Angelica

Thank you for your interest in contributing. This document explains how to
report bugs, suggest improvements, and submit code changes.

## Code of conduct

Be respectful and constructive in all interactions. Harassment or personal
attacks of any kind will not be tolerated.

## Reporting bugs

Open an issue on [GitHub](https://github.com/ScimonCFD/Angelica/issues) and
include:

- A minimal description of the network that triggers the problem
- The Angelica version (shown in the title bar)
- What you expected and what actually happened
- Any error message or traceback

If possible, attach the `.gui.json` file for the case that fails.

## Suggesting features

Open an issue with the `enhancement` label and describe:

- What you are trying to do that is currently not possible
- Why it would be useful to others
- Any references (papers, standards) that describe the desired behaviour

## Setting up a development environment

Requires Python 3.8 or later.

```bash
git clone https://github.com/ScimonCFD/Angelica.git
cd Angelica
pip install -e ".[dev]"
```

## Running the tests

```bash
pytest tests/
```

All tests must pass before submitting a pull request.

## Submitting a pull request

1. Fork the repository and create a branch from `main`.
2. Make your changes. Keep each PR focused on a single fix or feature.
3. Add or update tests if the change affects solver behaviour.
4. Confirm the full test suite passes.
5. Open the pull request with a clear description of what changes and why.

PRs that break existing tests or lack a description will not be merged.

## Code style

- Follow the conventions already present in the file you are editing.
- No unnecessary comments. Code should be self-explanatory through clear
  naming.
- Keep functions short and focused.
