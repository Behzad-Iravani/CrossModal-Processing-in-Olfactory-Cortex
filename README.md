# CrossModal-Processing-in-Olfactory-Cortex
 We explored the dynamic functional connectivity between the posterior piriform cortex (PPC) and visual/auditory object processing areas.  We investigate the relevance of object-related information in the olfactory system. Our code utilizes BOLD activity from 272 ROIs, applying t-tests to generate group-level t-maps and identifying the strongest connections. Explore the shortest path between auditory/visual areas and PPC

## Methods 
To implement this analysis, BOLD activities from 272 spherical ROIs with a radius of 18 mm were extracted across the whole brain. The BOLD signals were denoised through a standard denoising pipeline, which involved the regression of confounding variables. Specifically, these variables encompassed motion parameters along with their first-order derivatives (12 factors) as well as adjustments for task-related effects (2 factors: pictures and sounds). The denoised BOLD were windowed to a five-time-point segments using a Hanning sliding window with 50% overlap. Subsequently, we computed the Pearson correlation between brain 272 brain ROIs within each window yielding time-series data representative of the dFC. 
On the individual level, we constructed a linear regression analysis to assess the task modulation of the dFC. A design matrix with intervals labeled either auditory or visual was created. These intervals represent the consecutive time points where each corresponding stimulus was presented to individuals. To account for the hemodynamic response, the design matrix was convolved with HRF. The linear regression determined the extent that the task (i.e., picture and sounds) affected the brain-wide dFCs, yielding so-called individual connectivity beta maps. At the group level, the beta maps underwent a t-transformation and were subjected to averaging across cerebral hemispheres. To eliminate false positive connections, the t-transformed maps were thresholded based on two criteria: 1) any average gray matter probability falling below 40% was discarded; 2) only the top 20% of the strongest connections were retained.
We opted to focus on LOC and hAC based on the assumption that olfactory system preferentially integrates object-related information rather than the low-level features processed by the primary visual and auditory areas (V1 and A1). In general, the PPC receives more input from other cortical areas than does APC (Cohen et al., 2017), and we therefore consider only PPC in this analysis. To assess the functional connectivity between PPC and object processing areas within visual (LOC) and auditory (hAC) cortices, we converted the group t-transformed and thresholded dFC connectivity matrix and convert it to distances matrices, by inverse transform. The inverse transformation was required to be able to use  Floyd–Warshall algorithms (Cf. Fransson et al., 2011) to find the shortest path between each pair of nodes, i.e., between auditory and visual object processing areas and the PPC.

## Results 
We demonstrated above that PPC processes unisensory visual and auditory objects. We next sought to understand which cerebral node might mediate the information to the PPC. We extracted the BOLD signal from 272 ROIs. Figure 5A. The denoised and task-related adjusted bold signals were windowed using a 50% overlapping Hanning window and correlated to created dFCs, Figure 5B. A linear regression was used to estimate the task modulation of dFC at individual level, Figure 5C. T-transformed and thresholded group level dFCs were calculated for the two conditions, Figure 5D. Using a Floyd–Warshall algorithm we assessed the task relevant shortest path between the auditory and visual object-oriented cortices and PPC.  . We found that both auditory and visual object areas during associated tasks (i.e., hAC cortex during sounds and LOC during pictures) are connected to PPC via amygdala.

![Figure1](https://github.com/Behzad-Iravani/CrossModal-Processing-in-Olfactory-Cortex/assets/7909726/155521f1-a1c7-4c6d-995e-005cd66e3e7b)
Figure1. 
