# Structure

`1HK9.pdb` is the *Escherichia coli* Hfq crystal structure (PDB 1HK9; Sauter et al., 2003). It is a homologous reference structure, not a Salmonella-specific model. The download script validates the PDB header before saving, and the analysis independently validates headers before generating Fig10.

To download or refresh it (internet required), run `python scripts/fetch_hfq_structure.py`. The script prints the PDB header and refuses to save a file unless the header identifies Hfq / host factor-I.

Reference: Sauter C, Basquin J, Suck D. SM-like proteins in Eubacteria: the crystal structure of the Hfq protein from *Escherichia coli*. Nucleic Acids Res. 2003;31(14):4091-4098. [doi:10.1093/nar/gkg480](https://doi.org/10.1093/nar/gkg480).
