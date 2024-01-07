# CrossModal-Processing-in-Olfactory-Cortex
 We explored the dynamic functional connectivity between the posterior piriform cortex (PPC) and visual/auditory object processing areas.  We investigate the relevance of object-related information in the olfactory system. Our code utilizes BOLD activity from 272 ROIs, applying t-tests to generate group-level t-maps and identifying the strongest connections. Explore the shortest path between auditory/visual areas and PPC

## Methods 
To implement this analysis, BOLD activities from 272 spherical ROIs with a radius of 18 mm were extracted across the whole brain. The BOLD signals were denoised through a standard denoising pipeline, which involved the regression of confounding variables. Specifically, these variables encompassed motion parameters along with their first-order derivatives (12 factors) as well as adjustments for task-related effects (2 factors: pictures and sounds). The denoised BOLD were windowed to a five-time-point segments using a Hanning sliding window with 50% overlap. Subsequently, we computed the Pearson correlation between 272 brain ROIs within each window, yielding time-series data representative of the dFC. 
On the individual level, we constructed a linear regression analysis to assess the task modulation of the dFC. A design matrix with intervals labeled either sound or picture was created. These intervals represent the consecutive time points where each corresponding stimulus was presented to individuals. To account for the hemodynamic response, the design matrix was convolved with HRF. The linear regression determined the extent that the task (i.e., sounds and pictures) affected the brain-wide dFCs, yielding so-called individual connectivity beta maps. At the group level, the beta maps underwent a t-transformation and were subjected to averaging across cerebral hemispheres. To eliminate false positive connections, the t-transformed maps were thresholded based on two criteria: 1) any ROI with averaged gray matter probability falling below 40% was discarded; 2) only the top 20% of the strongest connections were retained.
We opted to focus on LOC and hAC based on the assumption that olfactory system preferentially integrates object-related information rather than the low-level features processed by the primary visual and auditory areas (V1 and A1). In general, the PPC receives more input from other cortical areas than does APC (Cohen et al., 2017), and we therefore consider only PPC in this analysis. To assess the functional connectivity between PPC and object processing areas within visual (LOC) and auditory (hAC) cortices, we converted the group t-transformed and thresholded dFC connectivity matrix and convert it to distances matrices, by inverse transform. The inverse transformation was required to be able to use  Floyd–Warshall algorithms (Cf. Fransson et al., 2011) to find the shortest path between each pair of nodes, i.e., between auditory and visual object processing areas and the PPC.


## Results 
We demonstrated above that PPC processes unisensory visual and auditory objects. We next sought to understand which cerebral node might mediate the information to the PPC. We extracted the BOLD signal from 272 ROIs. Figure 5A. The denoised and task-related adjusted bold signals were windowed using a 50% overlapping Hanning window and correlated to create dFCs, Figure 5B. A linear regression was used to estimate the task modulation of dFC at individual level, Figure 5C. T-transformed and thresholded group level dFCs were calculated for the two conditions, Figure 5D. Using a Floyd–Warshall algorithm we assessed the task relevant shortest path between the auditory and visual object-oriented cortices and PPC. We found that both auditory and visual object areas during their associated tasks (i.e., hAC cortex during sounds and LOC during pictures) are connected to PPC via amygdala [x ± 26.25, y -8.25, z -26.25],  (Pictures: LOC--AMY t(46) = 17, p = 0, AMY-- PPC t(46) = 17, p = 0; Sound: hAC--AMY t(46) = 13, p = 0, AMY--PPC, t(46) = 17, p = 0), Figure 5E.

![Figure1](https://github.com/Behzad-Iravani/CrossModal-Processing-in-Olfactory-Cortex/assets/7909726/13021af4-626e-4afb-a106-69391138921a)
Figure 5. The dFC analysis determined the task relevant shortest functional path between object-oriented cortices and PPC. A) An atlas with 272 ROIs were selected to cover the whole cerebral. B) BOLD signals were extracted from the 272 ROIs and were denoise prior to being windowed by a 5-sample-long Hanning sliding window with 50% overlapping (alternating between 2 and 3 samples overlaps). Pearson correlation was used to construct the functional connectivity for each window, cumulatively creating brain-wide dynamic functional connectivity (dFC). C) Linear regression was set up to assess the task modulation of the dFC at individual level. The design matrix indicates samples relevant to each stimuli (i.e., sounds and pictures). By fitting the regression the effect tasks were estimated as individual Beta values. D) On the group level the Beta values were t-transformed and thresholded based on the average gray matter of ROIs in the cohort. The small brain map in the middle indicates the averaged gray matter probability that was used to threshold (discarded ROIs with GM < 40%) the group level t-transformed dFCs. E) The 20% strongest dFC were subjected to inverse transformation tprior to fed to Floyd–Warshall algorithm to find the task relevant shortest path between the audio (hAC) and visual (LOC) object oriented cortices and PPC . 
