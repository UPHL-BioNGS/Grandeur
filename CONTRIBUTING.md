# Contributing to Grandeur

First off, thank you for considering contributing to Grandeur! It's people like you that make Grandeur a great tool for the bioinformatics and public health communities.

## How to Contribute

### Reporting Bugs

If you find a bug or run into a pipeline error, please open an issue on our [GitHub Issues page](https://github.com/UPHL-BioNGS/Grandeur/issues). 

When reporting a bug, please include:
* A detailed description of the issue.
* The exact Nextflow command you ran.
* The Nextflow `.nextflow.log` file.
* Information about your compute environment (e.g., local, Slurm, AWS) and container system (Docker or Singularity/Apptainer).

### Suggesting Enhancements

We welcome suggestions for new features, organism-specific subtyping tools, or general pipeline improvements. Please submit an issue on the [GitHub Issues page](https://github.com/UPHL-BioNGS/Grandeur/issues) detailing your request.

### Pull Requests

We gladly accept pull requests (PRs) for bug fixes and new features. 

1. **Fork the repository** and create your branch from `main`.
2. **Write clean code** and document your changes.
3. **Test your changes:** Ensure that your changes pass the existing Continuous Integration (CI) tests in GitHub Actions. Grandeur relies on `nf-test` and several dataset-specific CI workflows (e.g., `ecoli.yml`, `salmonella.yml`, etc.) to verify functionality.
4. **Lint your code:** Grandeur utilizes `nf-core` guidelines. Please ensure your code passes `nf-core pipelines lint`.
5. **Update documentation:** If you are adding a new parameter or process, please update the `nextflow_schema.json` and relevant Wiki/Markdown documentation.
6. **Submit the PR** with a clear description of why the changes were made and how they were tested.

## Code of Conduct

Please note that this project is released with a Contributor Code of Conduct. By participating in this project you agree to abide by its terms.
