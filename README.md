# Decoding the naked mole-rat genome, telomere to telomere

Code for the analyses in the manuscript (Sokolowski, Nater, Miclăuș et al.), kept as run.

| Directory | Analysis | Figure | Author |
|---|---|---|---|
| `annotations` | gene ([GenAnT](https://github.com/BaderLab/GenAnT)), repeat, satellite, pseudogene and epigenetic annotation | 1, 2 | Dustin Sokolowski |
| `assembly_QC` | assembly evaluation and overview | 1 | Dustin Sokolowski |
| `SatelliteFreeCentromeres` | gene and L1 content of centromeres, centromere models | 2 | Dustin Sokolowski |
| `chr15` | chromosome 15 structural polymorphism | 2 | Dustin Sokolowski |
| `Sperm_loss` | sperm gene loss, RNA-seq and proteomics | 3 | Dustin Sokolowski |
| `Sperm_relaxed_selection` | RELAX selection pipeline (see `RUNBOOK.md`) | 3 | Wing Yui Ngan |
| `gene_duplications` | gene copy-number events | 4 | Dustin Sokolowski |

Code in separate repositories:

| Repository | Analysis | Figure | Author |
|---|---|---|---|
| [KMCG](https://github.com/Oujiang-Laboratory-Bioinformatics/KMCG) | k-mer multiplicity consensus graph, used to evaluate and curate the assembly | 1 | Oujiang Laboratory |

Data:

- Sequencing reads and assemblies: ENA umbrella project [PRJEB108765](https://www.ebi.ac.uk/ena/browser/view/PRJEB108765)
- Annotation files: [Zenodo](https://zenodo.org/records/23217020)
