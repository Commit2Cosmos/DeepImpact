# Contributing

## Getting Started
- Install project dependencies following `README.md`.

## Branching & Workflow
- Never commit directly to `main`; treat it as protected.
- Create a new branch (e.g., `aiburst/<branch_name>` or `damage/<branch_name>`).
- Fetch from `main` frequently.

## Pull Requests
- Run the full test suite (`pytest`, black, ruff) before pushing.
- Update documentation and add tests when behavior change.
- Assign `ada-ab2125` as the primary reviewer on every PR.

## Code Style & Quality
- Use type hints and docstrings for classes and functions.

## Commits
Write meaningful commit messages and use consistent prefixes:
- `feature` for new implementations
- `fix` for bug fixes
- `docs` for documentation
- `style` for code style
- `refactor` for code changes
- `perf` for performance
- `test` for tests
- `build` for builds
- `ci` for CI
- `chore` for other changes


## Testing
- Add or update unit/integration tests alongside code changes.
- Ensure tests are deterministic and clean up any external resources they create.
- Include reproduction steps in the PR if fixing a bug.

## Communication
- Use GitHub Issues/Discussions for feature ideas or large changes before coding.
- Flag blockers or open design questions early so reviewers can help.
- Be respectful and clear in all collaboration; we're all here to build together.

Thank you for contributing!

