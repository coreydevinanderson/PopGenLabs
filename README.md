# PopGenLabs

This is a general (public repository) for an undergraduate laboratory in population genetics with R. I developed this lab for my undergraduate Ecology and Evolution course at Valdosta State University. It is a student favorite (when the real candy is used). The R part makes the point, but could be improved in many ways.

The lab is a demonstration of genetic drift, focusing on how finite sampling results in a loss of genetic variation each generation through time; eventually only one allele (M&M color) will remain due to sampling error alone.

The students start with 20 M&Ms, four for each of the five colors chosen (let them eat one color). They sample 10 M&Ms from the brown bag of 20 M&Ms. They spill out what was not sampled into a side pool. They then repopulate the bag back up to a count of 20 M&Ms, based on their relative frequencies for the sample of 10. For example, if 2 of 10 were Blue, then 4 get returned to the brown bag for the next generation of sampling; if there was one Red, two Red would be returned to the brown bag. Just tell them to double the counts of each color when re-populating from 10 back up to 20. It is a good idea to keep on an eye on the students as they get going. A common mistake is to not empty the bag after each generation.

Procedure:

* Sample 10 from bag of 20; example sample: 
 
 2-Blue (0.2), 1-Green (0.1), 3-Orange (0.3), 2-Yellow (0.2), 2-Red (0.2) -> 

* Pour unsampled 10 into "slush pile." 

* Repopulate brown bag at count of 20, based on relative frequency out of 10 (above). Based on the sample from step 1 (above):
  
  4-Blue (0.2), 2-Green (0.1), 6-Orange (0.3), 4-Yellow (0.2), 4-Red (0.2)

* Repeat up to 20 generations, or until only one color is left. 


In the most basic form of the lab, just have them keep track of which generation they are on. If you want to get fancy, you could have them record the frequency of a certain color and then make histograms at different time points to demonstrate how drift increases in the variance in the frequency of alleles in demes through time, and how this relates to the measure of variance effective size (which measures the rate of loss of genetic variation in a deme due to drift).

Some groups will lose colors quickly, others will not reach fixation by generation 20. While not all groups will have bags that go to fixation in 20 generations, all groups seem to get the basic idea: they will lose colors through time until only one color is left.

The R code (below) makes plots to simulate the process and show the outcome after each simulation through time, including how many generations it took until only one color was left. The exercise can be twisted in many cool ways (bottle-necks, gene flow and drift, selection, etc.).

Adding a nice interface, allowing students to change the number of alleles (or other basic parameters) could improve learning and make it more accesible to students that are not learning about R programming.

---
Repository manager:

Corey Devin Anderson

Professor of Biology

Valdosta State University

coreanderson@valdosta.edu
