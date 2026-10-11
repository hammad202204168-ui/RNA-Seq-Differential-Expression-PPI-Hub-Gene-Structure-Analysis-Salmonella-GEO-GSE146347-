# Publishing this folder on GitHub (step by step)

1. Create an **empty** repository on https://github.com/new
   * Name: `RNA-Seq-Differential-Expression-PPI-Hub-Gene-Structure-Analysis-Salmonella-GEO-GSE146347-`, visibility Public, do **not** add README/.gitignore/license (they already exist here).
2. In a terminal, inside this folder:
```bash
git init -b main            # skip if already initialised
git add .
git commit -m "Initial release v2.2.0"
git remote add origin https://github.com/hammad202204168-ui/RNA-Seq-Differential-Expression-PPI-Hub-Gene-Structure-Analysis-Salmonella-GEO-GSE146347-.git
git push -u origin main
```
3. The GitHub username and repository link are already filled in README.md and CITATION.cff.
4. Check the **Actions** tab: the `pipeline` workflow should turn green.
5. Optional: Settings > Pages is not needed; add topics (`rna-seq`, `salmonella`, `hfq`, `bioinformatics`); create a release `v2.2.0`.
   For a DOI, link the repository to Zenodo (zenodo.org > GitHub) before creating the release.

Delete this file after publishing if you prefer.
