# Contributing to MuRAT

Thank you for your interest in MuRAT. Contributions of all kinds are welcome:
bug reports, fixes, documentation, new examples, and new features. The people
who have given substantial contributions so far are listed in
[AUTHORS.md](AUTHORS.md).

Please note that this project is governed by a [Code of Conduct](CODE_OF_CONDUCT.md).
All contributors are expected to adhere to it.

## Getting help

- **Questions about using MuRAT:** first check the [README](README.md),
  `Documentation_MuRAT.pdf` and the [wiki](https://github.com/LucaDeSiena/MuRAT/wiki).
  If your question is not answered, open an issue.
- **Anything private or sensitive:** contact the principal developer at
  `lucadesiena80@gmail.com`.

## Reporting a bug

Open an issue at <https://github.com/LucaDeSiena/MuRAT/issues> and include:

1. What you expected to happen and what happened instead, with the complete
   MATLAB error message.
2. Your MATLAB release (`version`), operating system, and the output of `ver`
   (installed toolboxes).
3. The input file you used (`Murat_input*.m`), or the changes you made to one
   of the examples.
4. If possible, a minimal dataset that reproduces the problem. Please do not
   attach large data files; link to them instead.

## Requesting a feature

Open an issue describing the scientific use case, the method (with references
if it is published), and what you would expect as input and output. Discussing
the idea first avoids duplicated work.

## Contributing code or documentation

1. Fork the repository and create a branch from `master` (for example
   `fix-sac-header-check`).
2. Make your change. Please:
   - follow the existing naming convention (`Murat_*.m` functions in `bin/`)
     and give every new function a help header (purpose, inputs, outputs,
     example);
   - keep units consistent with the README (event depth in km, station
     elevation in m);
   - avoid adding large data files to the repository.
3. Add or update tests in `Tests/`. Tests are class-based (`matlab.unittest`).
   To run them locally from the repository root:

   ```matlab
   results = runtests('Tests');
   assertSuccess(results);
   ```

   The same tests run automatically on GitHub Actions for every push and pull
   request.
4. Update `CHANGELOG.md` and the documentation if your change affects
   behaviour or inputs.
5. Open a pull request and describe what changed and why. A maintainer will
   review it; please be prepared to revise.

## Third-party code

`Utilities_Matlab/` bundles several third-party components under their own
licences. Before adding new third-party code:

1. Check that the licence permits redistribution alongside MuRAT's MIT licence.
   **GNU GPL v2/v3 code is incompatible** with MIT redistribution; do not add
   GPL-licensed code without prior discussion.
2. Add the licence file to the component's subfolder.
3. Add a row to the table in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md)
   with the folder, file list, authors/copyright, and licence.
4. If the component has no stated licence, add a clarifying notice file (see
   the existing `LICENSE_inferno.txt`, `LICENSE_redblue.txt`, and
   `LICENSE_inpaintn.txt` for the pattern to follow) and document the
   outstanding uncertainty in the Notes section of `THIRD_PARTY_NOTICES.md`.
5. Prefer MIT- or BSD-licensed alternatives where functionally equivalent
   options exist.

## Licence

By contributing your own code, you agree that your contribution is released
under the repository's [MIT licence](LICENSE.md). Contributions that consist
of third-party code remain under their original licence as documented in
[THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md).
