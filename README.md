# big_malaria_mapping_methods_comparison

This analysis is in the archive drawer.
I'm not really doing anything with it.
If you find it interesting feel free to
take code or email me to discuss elements of the analysis.

I started working on this as a big comparison of different malaria mapping methods.
I was aiming to examine the benefits of 1) adding new, better covariates and 2)
the benefits of various modelling techniques.

One angle I was looking at was the different way space is used by different models.
Model based geostatistics and spatial machine learning models use space as a special 
covariate. 
Spatially weighted regression "controls" for space by fitting models to local data only.
Varying coefficient models use space only as a modifier for other covariate effects.

There was one particular model that I thought was novel but never coded up or fitted.
I wanted to fit ML models and then ensemble them in a spatially varying way.
So in one area, the elastic net will be the dominant model while in another area the
RandomForest would be. 
This was to be done with a number of Guassian Processes (1 for each ML model) that
summed to one.


Overall, this repo is a mess and I'm not really taking this analysis forward.
I still sometimes use it as a test for a new method. 
I've got some nice benchmarks etc. so I just check a method at the data.
If it underperforms other methods there's little reason to consider that
model further.



