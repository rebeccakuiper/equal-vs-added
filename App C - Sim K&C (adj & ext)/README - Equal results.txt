Equal results
=============
Summaries and plots of the equal-evidence results of the simulations of the
Added-evidence paper (Kuiper & Clapper), in the same way as was done for the
added-evidence results there (Appendix C of the Equal-vs-Added paper).
The plots for the simulation of the Equal-vs-Added paper (Section 6) are in the
folder 'Sec 6 - Sim evSyn vs GORICA on meta-an est'.

The scripts use the original code in 'Orig code K&C' (folders
'Sim - Added - Final' and 'Sim - extra (review)'); set the
working directory to the folder of the script you run (RStudio: Session > Set
Working Directory > To Source File Location).

functions_analysis_Equal.R
    Extra plot/summary functions (used next to the original analysis/functions_analysis.R):
    effect_on_hf() (the original code for this plot is not available), THR_sets_plot(),
    HSR_plot() (as in post_analysis.R), compare_marginal() (Added vs Equal), and
    marginal_effect_Eq() = marginal_effect() with a y-axis for the support panel that also
    shows the 'all' and 'false' cases (with the original y-axis these fall off the plot for Equal).

1 Main sim (S = 4)
    NB: The equal-evidence weights were computed in the original run, but not saved
    (only results_added). Therefore, the simulation is re-run (same design and code;
    row-specific seeds), saving Added and Equal (based on the same data sets).
    1 - rerun_simulation_main.R      re-run (with rerun_functions.R: goric() call for the current restriktor)
    2 - process_results_main_Equal.R  = process_results.R for Equal (and Added re-run), incl. check
                                        of the Added re-run against the original Added results
    3 - analyse_results_main_Equal.R  = analyse_results.R + post_analysis.R for Equal, plus comparisons
    4 - rerun_simulation_main_hetero_corrected.R  hetero condition with the CORRECT probit/logit error
                                        variances (rerun_functions_hetero_corrected.R); N <= 150
    5 - compare_hetero_corrected.R     original vs corrected variances, Added and Equal
    results/  results_{added,equal}_rerun.RData (output of 1), results_{added,equal}_hetero_corrected.RData (output of 4)
    data/     data_to_analyse_main_*.RData, avg_per_condition_main_*.RData (output of 2 and 5)
    images/   *_equal.pdf: Equal versions of the plots of the Added paper;
              *_added.pdf: same for the Added re-run; compare_*.pdf / THR_sets_added_vs_equal.pdf: Added vs Equal
    *_output*.txt: printed output (tables) of the scripts
    Findings main simulation - Equal vs Added.txt: summary
    NB: the re-run covers N = 25-150 (not 200-1000).

2 Extra sim (S = 20)
    process_results_extra_Equal.R   = process_results.R for Equal (and Added; reproduces data_to_analyse_2024_09_22)
    analyse_results_extra_Equal.R   = analyse_results.R for Equal (N_*.pdf, R2_*.pdf), plus Added vs Equal
                                      comparisons, and also the case where the informative hypothesis is
                                      NOT true (inverse_ratio = 1; *_Hfalse.pdf), which was not inspected
                                      in Kuiper & Clapper.
    rerun_simulation_extra_corrected.R  re-run with the CORRECT probit/logit error variances;
                                      process/analyse with corrected <- TRUE (output: *_corrected)
    images.zip, images_corrected.zip: plots (original and corrected variances)
    THR_overview_extra_Added_vs_Equal(_corrected).csv: THR and median weight true hypothesis per set, N, R2
    Findings extra simulation - Equal vs Added.txt: summary

3 App C - Fig 8
    App C - Fig 8 -- Equal vs Added in sim K&C.R
        Figure 8 (Appendix C) of the Equal-vs-Added paper (images/sim_KC_THR.png), based on the
        results of folders 1 and 2 with the corrected probit/logit error variances.

NB (probit/logit error variances): in estimate_betas() of both Added-paper simulations, the latent
error variances of the probit (should be 1) and logit (should be pi^2/3) models were swapped. This
changes the effect sizes of these studies (not which hypothesis is true); correcting it changes the
THRs by at most .03.
