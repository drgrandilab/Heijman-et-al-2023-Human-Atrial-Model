Matlab code for Heijman et al. (Circ Res. 2023) human atrial model.

This model is based on the framework of the Grandi et al. human atrial model (Circ Res. 2011;
109:1055-66, available for download on this website), previously modified to reproduce Na-dependent
regulation of IK1 and IK,ACh, and here extended to integrate a previously published model of IK2P
and a novel formulation of ISK.

_____________________________________________________________________________________________________
Contents:

human_atrial_main.m		loads initial conditions and runs the simulation
human_atrial_model.m		excitation-contraction coupling model
human_atrial_calcCurrents.m	supporting function for simulation output analysis

.mat files contain steady-state initial conditions in nSR or cAF with normal or blocked SK current 
obtained at different pacing frequencies

Executing the following codes allows reproducing the published figures:
Fig_2_Linear_regression.m
Fig_S6_Freq_analysis.m
Fig_S8_Caffeine_protocol.m
Fig_S9_ICaL_block_protocol.m
Fig_S10_Freq_analysis_SK_block.m
_____________________________________________________________________________________________________

References:

Heijman J, Zhou X, Morotti S, Molina CE, Abu-Taha IH, Tekook M, Jespersen T, Zhang Y, Dobrev S,
Milting H, Gummert J, Karck M, Kamler M, El-Armouche A, Saljic A, Grandi E, Nattel S, Dobrev D.
Enhanced Ca2+-Dependent SK-Channel Gating and Membrane Trafficking in Human Atrial Fibrillation.
Circ Res. 2023. https://doi.org/10.1161/CIRCRESAHA.122.321858

Grandi E, Pandit SV, Voigt N, Workman AJ, Dobrev D, Jalife J, Bers DM.
Human Atrial Action Potential and Ca2+ Model: Sinus Rhythm and Chronic Atrial Fibrillation.
Circ Res. 2011 Oct 14;109(9):1055-66. doi: https://doi.org/10.1161/CIRCRESAHA.111.253955

Please cite the above papers when using this model.
