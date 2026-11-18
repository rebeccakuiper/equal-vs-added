Scripts for the Equal-vs-Added paper, named using section/table/figure numbering:

Showcases = Section 5, 
Simulation = Section 6, 
Appendix A = Equal vs pooling,
Appendix B = GORICA on meta-analyzed estimates (Tables 3-4; Appendix B-4: the example of Berkey et al. (1998)
             in subsection 'When Do FE and RE Differ?'), 
Appendix C = Kuiper & Clapper simulations (Figure 8),
Appendix D = unconstrained as failsafe (Tables 5-6)). 

The simulation files are in the folders 'Sec 6 - Sim evSyn vs GORICA on meta-an est' (Section 6) and
'App C - Sim K&C (adj & ext)' (Appendix C); see '00 - README - Sim files.txt'. The folder
'Other (not in paper)' contains a small additional exploration that is not reported in the paper.

Abbreviations in the file and folder names:
  Sec = Section, App = Appendix, Fig = Figure, Tab = Table, Sim = simulation,
  K&C = Kuiper & Clapper (2025), S4/S20 = 4/20 studies, N150/N1000 = sample size of 150/1000 per
  study (otherwise 85), opp = opposite ordering (H_123 is false), indep = independent samples per study
  (otherwise the same data set for each study), evSyn = equal- and added-evidence synthesis
  (evSyn function of restriktor), meta-an est = meta-analyzed estimates, adj & ext = adjusted &
  extended.

The folder 'Figures (pdf)' contains the PDF versions (vector graphics; for the journal) of all
figures in the paper (Figures 1-8), named by figure number (e.g., 'Fig1a_...pdf'), together with
'00 - Figure descriptions.txt' (short description of each figure and which script creates it).
The scripts that create the figures save these PDF versions automatically (besides the png versions).

The script 'Sec 5-2-3 -- evSyn - S20 - indep - N1000.R' (20 independent datasets with 1000 observations
per study; same population and seed as for 85 and 150 observations) gives the example with n_s = 1000 in
the text of Section 5.2.3; it has no figure in the paper. It also saves the PDF plots used on the slide
'Showcase: 20 independent datasets, large sample size' in the presentations GORIC_evSyn and
GORIC_evSyn_MultipleStudies.
