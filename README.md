# Temporal Credit Assignment Project

Codes for analysis and drawing figures in temporal credit assignment project.

## Authors
Dong-Hyun Lim 1,2, Young-Ju Yoon 1,2, Eunsil Her 3, Suehee Huh 2 & Min Whan Jung 1,2,4*

1 Department of Biological Sciences, Korea Advanced Institute of Science and Technology, Daejeon, 34141, Korea

2 Center for Synaptic Brain Dysfunctions, Institute for Basic Science, Daejeon, 34141, Korea

3 Department of Psychiatry, Ajou University School of Medicine, Suwon, 16499, Korea

4 Lead Contact

*Correspondence: mwjung@kaist.ac.kr

## Task description
For detailed task description please note our previous study in the reference (Her et al., 2016).

Rats were trained to perform a dynamic two-armed bandit (TAB) task.

## Code description
The analysis was done with MATLAB R2013a, and drawing figures with R2013a and R2017a for some cases.
*****
### Raw data sets
#### Data description
The raw data sets are described below:
1. txt files of behavior results including trial number, choice, reward and trial information (error, duration, block) for each session
2. nev files of beam-break times
3. nvt files of recorded trajectory of the animal
4. t files of isolated single units
#### Data path description
Each sessions's .nev, .nvt, .t files are assumed to be stored in <code>E:/Data/Neural analysis/"BRAIN_REGION"/"YYYY-MM-DD_HH-MM-SS"</code>

txt files are assumed to be stored as above and also in <code>E:/Data/Behavioral analysis/"ANIMAL_NAME"</code>
*****
### MATLAB Path setting
    FindFiles.m
    get_Spiketime.m
    figq_Q_Y2.m
    QLfun_Q.m
    AxesPoint_1raster_new.m
    
should be on the MATLAB path. (Or download and adding path <code>/basic</code> folder.)

*****
### Making event .mat files
The .nev and .txt files are first converted to mat files.
1. Run <code>Behavior code/1_Event.m</code> with proper file directories. This converts nev files to <code>Event.mat</code> files.
2. Run <code>fig1C_1220_conversing_mean_yj.m</code> with proper file directories. This makes <code>Event_New_YJ_1220.mat</code> files.
Here the memory stage onset (convergence point) is calculated as described on the reference (Her et al., 2016) with 100ms (3 points) of criteria.
3. Run <code>diversing_mean_1220_yj.m</code> with proper file directories. This makes <code>Event_AC_YJ_1220.mat</code> files, the mat file with event timings. In this analysis, only event #4 and #5 (Memory and Reward stage) was used. Here the divergence point is calculated as described on the reference (Her et al., 2016).

*****
### Single unit files
From all sorted signle units, putative pyramidal cells were used in the analysis.

Run <code>clustering_3_1220_yj.m</code> varying the variable <code>region = 'IL';</code> among 'PRL','ACC' and 'IL' to classify and make a (.mat file) list of putative pyamidal cells.

Their properties distribution plot (fig 2B) is drawn by <code>fig2B_1220_clustering_plot_mpfc_yj.m</code> file.

*****
### Behavioral information
1. Run <code>Behavior code/1_Beh.m</code> with proper file directories. This reads the raw txt files and converts to <code>beh_.mat</code> file.
<code>beh_.mat</code> file consists of dummy variables regarding on choice (C[t]), reward (R[t]) and their interaction (X[t]).
2. Example trajectory figure is made by <code>fig1B_1220_exampletrace_yj.m</code>.
3. Figures (1C) with the x-position data are made by the codes with prefix 'fig1C'.
4. Example recording session's performance figure (1D) is made by <code>fig1D_1220_P_left_recording_sess_yj.m</code>.

*****
### Calculating Q-values
For each animal, from the behavioral result of all sessions the Q-values were calculated (Her et al., 2016).
1. Make the input data matrix by running <code>Behavior code/total_cond.m</code> for each animal.
2. Fit the behavior data on our model, with <code>Behavior code/total_animal_Q.m</code> for each animal.
3. Run <code>Q_value_per_rat_RL_ws_yj.m</code> to calculate Q_Left, Q_Right and Q_Chosen for each session as described in the reference.
*****
### Regression analysis
For regression analysis in figure 3 and 4, follow the codes with the prefix 'fig3' or 'fig4'.

First run 'regression' and then 'Fig' codes.

*****
### Example neurons
For choice signal example, neuron number 60 of ILC was selected.

For reward signal example, neuron number 53 of ACC was selected.

*****
### Other calculations
Behavioral performance calculations and mean discharging rate calculations mentioned in the results are done with <code>means_1220_yj.m</code> and <code>P_high_P_rwd_1220_yj.m</code>.

*****
## References
* Her ES, Huh N, Kim J, Jung MW (2016) Neuronal activity in dorsomedial and dorsolateral striatum under the requirement for temporal credit assignment. Sci Rep 6:27056.


https://github.com/in0565/TemporalCreditAssignment
