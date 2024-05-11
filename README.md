This is the step-to-step note by running my own data under manual of brainGraph.

The brainGraph is more recommanded to use in your work than BCT or GRETNA, from myself experience.

Anyway, I have learned much from the brainGraph, like GLM anlysis, mediation analysis, causal mediation analysis, and some basical information about graph theory. Hope you (someone) can also enjoy it like me.

This tutorial included my data and code. I have run it successfully in my Windows computer, and I believed it will also work in Linux.

**Notices:**
1. Here are [some notes](https://editor.csdn.net/md/?articleId=128781466) written by me, but with Chinese lauguage.
2. Here are the [result validations](https://editor.csdn.net/md/?articleId=129985875) between GRETNA (another commonly used software) and brainGraph. The conclusion is that results of both two are **totally same** (much great!), in spite of their different names of multiple network building. Specifically, the "Network Sparsity" in GRETNA responses to "density" in brainGraph, while the "Value of matrix element" in GRETNA responses to "threshold" in brainGraph.

**update May 11 2024**

I have uploaded a script 'lateralization_analysis.R'.
Previously, I have succesfully run fmriprep + qsiprep in my clinic data (Chronic Fatigue Syndrome) and the big population data (ISYB https://doi.org/10.11922/sciencedb.00740), getting the files of MNI sapce preprocessed Bold file and whole-brain-track tck file.

Then I used brainnectome atlas (246 ROIs) to compute the functional network and fiber network. Additionally, I have split those two
networks to three types of absolute network, positive network, and negative network, respectively.

After all preparing, I used this lateralization_analysis.R to compute network properties with brainGraph.

We know that the brainGraph reports errors frequently and is hard for new researchers who are not familiar with R software, especially 
the author of brainGraph seems not focusing on this package (from my guess). So I hope this tutorial finds you well. This tutorial works
for my data and I have already preparing my article. Hope you can also get helps from my script if you are also interested to brainGraph.
Despite many shortages of brainGraph, I think this package is greater than GREATNA and BrainConnectome Toolbox and others similar tools.

Sincerely yours,
Kang.
Department of Neurology, University of Iowa.
kangwu@uiowa.edu
clancy_wu@126.com
