The flagme package, is designed to process, visualise and analyze sets of GC-MS samples. The ideas discussed here were originally designed with GC-MS-based metabolomics in mind, but indeed some of the methods and visualizations could be useful for LC-MS data sets. The fragment-level analysis though, takes advantage of the rich fragmentation patterns observed from electron
interaction (EI) ionization. There are many aspects of data processing for GC-MS data. Generally, algorithms are run separately on each sample to detect features, or peaks (e.g. AMDIS). Due to retention time shifts from run-to-run, an alignment
algorithm is employed to allow the matching of the same feature across multiple samples. Alternatively, if known
standards are introduced to the samples, retention indices can be computed for each peak and used for alignment.
After peaks are matched across all samples, further processing steps are employed to create a matrix of abundances,
leading into detecting differences in abundance.
Many of these data processing steps are prone to errors and they often tend to be black boxes. But, with effective
exploratory data analysis, many of the pitfalls can be avoided and any problems can be fixed before proceeding to
the downstream statistical analysis. The package provides various visualizations to ensure the methods applied are
not black boxes.
The flagme package gives a complete suite of methods to go through all common stages of data processing. In
addition, R is especially well suited to the downstream data analysis tasks since it is very rich in analysis tools and
has excellent visualization capabilities.
