# SCWorkflow AI Coding Instructions

## Project Overview
SCWorkflow is an R package for single-cell RNA-seq analysis built on the Seurat framework. It's designed for analyzing multimodal 10x Genomics data, with support for CITE-Seq, cell hashing, and TCR-seq data. The package is deployed as both an R package and Docker container for use in NIDAP (Palantir Foundry) and Biowulf HPC environments.

## Core Architecture

### Sequential Workflow Pattern
Functions follow a numbered workflow sequence:
1. **processRawData()** - Process H5 files into Seurat objects
2. **filterQC()** - Quality control and filtering  
3. **combineNormalize()** - Merge samples, normalize, dimension reduction
4. **Harmony integration** (optional) - Batch correction
5. **annotateCellTypes()** - Automatic cell type annotation via SingleR
6. **Analysis functions** - DEG, visualization, clustering

### Seurat Object as Central Data Structure
- All functions expect/return Seurat objects as primary data containers
- Functions often modify objects in-place and return modified versions
- Metadata is heavily used for sample tracking and analysis parameters
- Multiple assays supported: RNA, SCT (SCTransform), ADT (CITE-seq), etc.

### Function Design Patterns
- **Comprehensive parameter lists**: Most functions have 15-30+ parameters with sensible defaults
- **Conditional workflows**: Functions check input types (H5 vs Seurat) and adapt behavior
- **Multi-sample handling**: Functions can process lists of Seurat objects or file paths
- **Plotting integration**: Most analysis functions return both data and visualization outputs

## Key Conventions

### File Organization
- `R/` contains one function per file, named descriptively (e.g., `Process_Raw_Data.R`)
- Function names use camelCase: `processRawData()`, `annotateCellTypes()`
- File names use Snake_Case with capitalization: `Process_Raw_Data.R`

### Parameter Naming Patterns
- `object` - Primary Seurat object input
- `samples.to.include` - Character vector for sample subsetting
- `reduction.type` - Visualization method ("umap", "tsne", "pca")
- `organism` - Species specification ("Human" or "Mouse")
- Boolean parameters use `.` separator: `do.normalize.data`, `draw.umap`

### Documentation Standards
- Extensive roxygen2 documentation with `@details` sections explaining workflow step numbers
- `@importFrom` statements for specific function imports
- `@export` for all user-facing functions
- Parameter descriptions include defaults and valid options

## Development Workflows

## Output standards
- Generate R CMD check friendly code.
- Use roxygen2 for exported functions.
- Add utils::globalVariables() for NSE variables.
- Prefer explicit namespaces (dplyr::, ggplot2::).
- Avoid side effects in tests.

### Template conventions
- JSON templates live under inst/extdata/NIDAPjson/.
- JSON is the source-of-truth for arguments, defaults, and behavior.
- orderedMustacheKeys defines argument order.

### Testing Structure
- Comprehensive `tests/testthat/` with both unit tests and integration tests
- `fixtures/` directory contains real Seurat objects for testing
- Helper functions in `helper-*.R` files for test setup
- Tests follow naming convention: `test-Function_Name.R`

### Deliverables
- Update CHANGELOG.md, decision_log.md, docs/session_notes.md.
- Add testthat tests for each generated function.

### Docker Development
- Base image: `nciccbr/ccbr_ubuntu_22.04:v4`
- Conda-based R environment (R 4.3.2)
- Multi-stage build supporting both development and production
- Container designed for HPC environments (Biowulf) and cloud platforms (NIDAP)

### CI/CD Integration  
- GitFlow-based R package workflow via `gitflow-R-action.yml`
- Automatic NIDAP deployment on successful builds
- pkgdown documentation generation
- Docker image building and deployment

## Critical Dependencies and Integration Points

### Bioconductor/Seurat Ecosystem
- **Seurat 4.1.1+**: Core framework for single-cell analysis
- **SingleR**: Automated cell type annotation
- **MAST/limma/edgeR**: Differential expression analysis
- **Harmony**: Batch effect correction and integration

### External API Integration
- `palantir_api_call.R`: Custom integration with Palantir Foundry APIs
- Designed for NIDAP cloud environment deployment
- Handles authentication and data transfer patterns specific to NIH infrastructure

### Memory Management Patterns
- `cell.count.limit` parameters (default: 35000) trigger memory conservation
- `only.var.genes` option reduces memory footprint for large datasets
- SCTransform normalization levels: sample-wise vs merged strategies

## Common Debugging Areas

### Data Type Handling
- Functions check for H5 vs Seurat object inputs: `is(input, "Seurat")`
- Cell filtering can drastically reduce object sizes - verify cell counts
- Assay availability: functions may require specific assays (RNA, SCT, ADT)

### Visualization Rendering
- Plot functions return both Seurat objects AND plot objects
- `reduction.type` must match available reductions in object
- Color palettes automatically generated but can be overridden

### Environment-Specific Issues
- Package designed for both local R and containerized environments
- File path handling differs between NIDAP, Biowulf, and local development
- Some functions have environment-specific conditional logic

## Quick Start for AI Agents
When working with this codebase:
1. Always check if input is a Seurat object: `class(object)`
2. Verify required reductions exist: `object@reductions`
3. Check available assays: `names(object@assays)`
4. Use `str(object@meta.data)` to understand metadata structure
5. Test functions start with small datasets from `tests/testthat/fixtures/`