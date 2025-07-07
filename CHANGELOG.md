# nf-core/grandeur: Changelog

All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.0.0dev - 2025-06-17

Initial release of nf-core/grandeur, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- update project with the latest commits.
- added gitlab templates for merge requests, feature requests, and reporting bugs.
- instantiated nf-core pipeline template.

### `Fixed`

### `Dependencies`

### `Deprecated`

Moved old pipeline structure to temporary location.

## 2025-06-18

Updated project config files.

### `Added`

- updated nextflow.config with the legacy nextflow.config configurations
- updated modules.config with legacy process configs
- added legacy configs to configs
- added legacy assets to assets

## 2025-06-20

Updated pipleine schema.

### `Added`

- Added Grandeur params to schema. Organized parameters in schema.

## 2025-06-23

### `Added`

- Added pipeline scripts to bin directory. 
- Added pipeline-specific input validation to the nf-core initialization tasks.

## 2025-06-25

### `Added`

- Added preprocessing subworkflow.
- Added de novo assembly subworkflow.
- Added relevant modules for subworkflows.
- Updated grandeur workflow with the new subworkflows.

### `Fixed`

- Fixed a couple files for end-to-end integration testing.

## 2025-06-30

### `Added`

- Added quality assessment subworkflow
- Added relevant modules for running the quality assessment subworkflow.
- Updated grandeur workflow with the new subworkflow.

## 2025-07-06

### `Added`

- Added the species identification subworkflow
- Added relevant modules for the new subworkflow.
- Updated grandeur.sh with the new subworkflow.

### `Fixed`

- fixed the MASH_DIST process in mash.nf to work with nf-core

## 2025-07-07

- Added the read identification subworkflow.
- Added relevant modules for the new subworkflow.
- Updated grandeur.sh with the new subworkflow
- Added blobtools.sif file to the project directory
