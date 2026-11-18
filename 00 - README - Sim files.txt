Simulation files for the Equal-vs-Added paper (numbering as in the text)
=====================================================================================
Set the working directory to the folder of the script you run (RStudio: Session > Set Working
Directory > To Source File Location). All paths are relative, so the folders can be moved together.

Sec 6 - Sim evSyn vs GORICA on meta-an est
    Sec 6 - Tab 2 & App D - Tab 6 -- Sim - evSyn vs MetaAn_GORICA.R
        the simulation (S = 20 and S = 4; FE, RE, Equal, Added); output sim_results.rds and
        sim_summary.csv -> Table 2 (Section 6) and Table 6 (Appendix D).
    Sec 6 - Fig 6 & 7 -- Plots sim (K&C style).R
        uses sim_results.rds -> Figures 6 and 7 (images/sim_weight_true_compl.png and
        images/sim_weight_false_compl.png; PDF versions in the folder 'Figures (pdf)' in the main
        folder of the R files: Fig6_sim_weight_true_compl.pdf and Fig7_sim_weight_false_compl.pdf), plus additional plots (images/*.pdf) and
        summary_sim_section6_per_condition.csv (THR, mean/median/5th/95th percentile weight of the
        true hypothesis, and proportion of weights < .05, < .25, < .45) and
        summary_sim_section6_selected_hypothesis.csv (weight of the selected hypothesis for all
        replications and for the true and false selection cases; numbers in the text on Figure 7).
    functions_analysis_Equal.R: plot functions used by the plot script.

App C - Sim K&C (adj & ext)
    Orig code K&C: the (needed) original simulation and analysis code and data
        of the Added paper (main simulation 'Sim - Added - Final' and 'Sim - extra (review)').
        This code is kept as it was released and is not changed here; only the names of the folders
        were shortened ('Original code Kuiper & Clapper' -> 'Orig code K&C' and 'simulation - extra
        based on review' -> 'Sim - extra (review)'). Two notes:
        - In estimate_betas() (sim/functions.R), the latent error variances of the probit (should be 1)
          and logit (should be pi^2/3) models are swapped. This is corrected in the re-runs in folders
          1 and 2 (rerun_functions_hetero_corrected.R replaces estimate_betas(); see also
          'README - Equal results.txt').
        - Some comments in this code speak of 'tested' hypotheses (sets). Since hypotheses are not
          tested but evaluated with the GORIC(A) and GORIC(A) evidence synthesis, 'evaluated' is the
          better word there ('test' only applies to null hypothesis significance testing).
    1 Main sim (S = 4): re-run of the main simulation (N = 25-150) keeping both the
        Added and Equal weights (the original run kept Added only), the re-run of the hetero
        condition with the correct probit/logit error variances, and the processing/analyses
        (Equal versions of the Added-paper plots; Added vs Equal comparisons).
    2 Extra sim (S = 20) [plots zipped: images.zip, images_corrected.zip]: processing/analyses of
        the Equal results of the extra simulation, and its re-run with the correct probit/logit error variances
        (set corrected <- TRUE in the processing/analysis scripts for those results).
    3 App C - Fig 8: the script for Figure 8 (images/sim_KC_THR.png; PDF version
        Fig8_sim_KC_THR.pdf in the folder 'Figures (pdf)' in the main folder of the R files).
    functions_analysis_Equal.R: plot functions used by the scripts.
    See also 'README - Equal results.txt' and the 'Findings ...' files in folders 1 and 2.

Other (not in paper)
    Extra sim - Equal vs Added vs FE: a small additional exploration (S, n, effect size).
