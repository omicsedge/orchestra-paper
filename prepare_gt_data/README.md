# Ancestry Analysis in Admixed Populations: A Toy Example

## Overview

This repository contains scripts for analyzing the ancestry of an admixed population, focusing on genetic data from Mexicans in Los Angeles (1000 Genomes Project). The analysis uses a reference panel comprising three populations:

* **Spanish**: Iberians from 1000 Genomes Project
* **Sub-Saharan African**: Yoruba from 1000 Genomes Project
* **Native American**: Maya, Mixe, Zapotec, and Pima populations from:
  - Human Genome Diversity Project
  - Simons Genome Diversity Project

Our primary objectives are to:
1. Perform local ancestry inference (LAI) using Orchestra (Lerga-Jaso et al. 2023)
2. Identify signals of natural selection in the admixed population

## Workflow

The analysis follows a structured pipeline:

### 1. Data Retrieval
- Download genotype data for target populations
- Prepare reference and target panels

### 2. Local Ancestry Inference
- Implement Orchestra tool for LAI analysis
- Follow methodology from Lerga-Jaso et al. (2023)

### 3. Selection Analysis
- Calculate admixture-based statistics:
  - Fadm scores
  - LAD scores
- Integrate scores to identify beneficial mutations
- Explore adaptive admixture in Mexican genomes

## Purpose

This repository serves as an educational example demonstrating:
- Orchestra's capabilities as a LAI tool
- Reproduction of selection signals from Cuadros-Espinoza et al. (2022)
- A complete, reproducible workflow using Orchestra instead of RFmix

## References

1. Lerga-Jaso, J., et al. (2023). "Retracing Human Genetic Histories and Natural Selection Using Precise Local Ancestry Inference." *bioRxiv* 2023.09.11.557177. 
   DOI: [10.1101/2023.09.11.557177](https://doi.org/10.1101/2023.09.11.557177)

2. Cuadros-Espinoza, S., et al. (2022). "The genomic signatures of natural selection in admixed human populations." *Am. J. Hum. Genet.* 109, 710-726.
   DOI: [10.1016/j.ajhg.2022.02.011](https://doi.org/10.1016/j.ajhg.2022.02.011)
   PMID: 35259336
