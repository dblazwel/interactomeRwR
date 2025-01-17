# The Human Interactome: Mapping Disease Co-Occurrences and Drug Repurposing

## Project Overview
This repository contains the scripts and resources developed during my Master's Thesis (TFM) as part of the Master's in Bioinformatics and Data Science in Personalized Precision Medicine and Health. This project was conducted in collaboration with the Barcelona Supercomputing Center (BSC) from June 2024 to January 2025, under the supervision of Jon Sánchez and Iker Núñez.

The central aim of the project was to leverage network biology approaches to explore disease co-occurrences and identify opportunities for drug repurposing. Specifically, the project focused on constructing and analyzing a multiscale human interactome, integrating protein-protein interactions, functional similarities, and co-citations from scientific literature.

## Objectives
- Develop a multiscale human interactome that integrates biological networks to map disease-gene and drug-gene associations.
- Use diffusion-based algorithms (e.g., Random Walk with Restart, RWR) to calculate disease and drug profiles across the interactome.
- Identify latent similarities between diseases and drugs, providing a foundation for precision drug repurposing and insights into disease comorbidities.

## Repository Structure
The repository is organized as follows:

- data/: Contains datasets used or generated during the project (e.g., interactome layers, gene-disease associations, diffusion profiles).
- scripts/: Includes Python and R scripts used for network construction, diffusion analysis, and result visualization.
- results/: Output files, including similarity matrices, visualizations, and performance metrics.

## Highlights of the Project
- Multiscale Human Interactome: A multilayer network integrating:
  - Protein-protein interactions (PPIs).
  - Functional gene relationships.
  - Co-citation data from biomedical literature.
- Diffusion Profiles: Calculated using RWR to simulate how diseases and drugs propagate across the interactome.
- Drug Repurposing: Identified potential therapeutic uses for existing drugs by analyzing profile similarities.
- Comorbidity Insights: Explored molecular mechanisms underlying disease co-occurrences.

## Technologies Used
- Programming Language: R.
- Libraries and Tools:
  - ggplot2, Cairo: Plotting
  - data.table, Dplyr, Reshape2, lsa, pROC, parallel: Data processing, analysis, and addtional tools
  - Multinet, RandomWalkRestartMH, igraph: Multilayer construction
  - AnnotationDbi, STRINGdb: For annotation and database management.
- High-Performance Computing: Analysis conducted on the MareNostrum supercomputer at the BSC.

## How to Use
1. Clone the repository:
```bash
git clone https://github.com/your-username/human-interactome-project.git
cd human-interactome-project
```
2. Run individual scripts from the scripts/ folder for specific tasks (e.g. network construction, diffusion analysis)

## Acknowledgments
This work was carried out as part of my Master's Thesis at the Barcelona Supercomputing Center (BSC). 
I am deeply grateful to my supervisors, Jon Sánchez and Iker Núñez, for their guidance and expertise throughout this project. 
Special thanks to the BSC team for facilitating access to computational resources.

## Contact
For questions, feedback, or collaboration opportunities, feel free to reach out:

- Email: [dblazquezg99@gmail.com](dblazquezg99@gmail.com)
- LinkedIn: [David Blázquez García](https://www.linkedin.com/in/dblazquezg/)




