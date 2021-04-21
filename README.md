# Philip_and_Newton_2021_PsychiatricResearch
Analytical code used in the publication: Transcriptional markers of excitation-inhibition balance in germ-free mice show region-specific dysregulation and rescue after bacterial colonization.
Available at: https://doi.org/10.1016/j.jpsychires.2021.01.021

Analyses shown for all brain regions largely, some scripts only for hippocampus.
Raw data provided in /Data directory.

Guide to scripts for analysis and visualizations
1. `<ANOVAs and heatmaps.R>`
   * One-way ANOVAs (significant at 10% FDR), with Bonferroni-corrected *post-hoc* t-tests
   * Heatmap and boxplot matrices visualizing changes

2. `<PCA and 3D plotting>`
   * Principle component analysis and 3D plotting
 
3. `<WGCNA module generation and export.R>`
   * Data QC, pre-analytical steps, and generation of modules with weighted gene co-expression analysis (WGCNA)
   * Data export for permutation testing and visualization in Cytoscape

4. `<WGCNA permutation testing.R>`
   * Permutation testing of observed differences in module composition usign permutation testing
   * Performed between all groups
 
5. `<preservation p-value meta-analysis.R>`
   * Resolving module-level preservation statistics to group-level statistics
