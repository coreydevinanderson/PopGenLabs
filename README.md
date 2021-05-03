# PopGenLabs

* This is a general (public repository) for an undergraduate laboratory in population genetics with R. I developed this lab for my undergraduate Ecology and Evolution course at Valdosta State University. It is a student favorite (when the real candy is used). The R part makes the point, but could be improved in many ways.

The lab is a demonstration of genetic drift, focusing on how drift results in a loss of genetic variation each generation through time; eventually only one allele (M&M color) will remain due to sampling error alone.

The students start with 20 M&Ms, four for each of the five colors chosen (let them eat one color). They sample 10 M&Ms from the brown bag of 20 M&Ms. They spill out what was not sampled into a side pool. They then repopulate the bag back up to a count of 20 M&Ms, based on their relative frequencies for the sample of 10. For example, if 2 of 10 were Blue, then 4 get returned to the brown bag for the next generation of sampling; if their was one Red, two "Red" would be returned to the brown bag. The easiest solution for them is to just double the counts when re-populating.

Procedure:

# Sample of 10 from bag of 20
2-Blue (0.2), 1-Green (0.1), 3-Orange (0.3), 2-Yellow (0.2), 2-Red (0.2) ->

# Pour unsampled 10 into dead pool.

# Repopulate brown bag at count of 20, based on relative frequency out of 10 (above):
 4-Blue (0.2), 2-Green (0.1), 6-Orange (0.3), 4-Yellow (0.2), 4-Red (0.2)

# Repeat up to 20 generations, or until only one color is left.

In the most basic form of the lab, just have them keep track of which generation they are on.

Some groups will lose colors quickly, others will not reach fixation by generation 20. While not all students will go to fixation in 20 generations, it is a cool magic trick that makes the point.

The plots simulate the process and show the outcome after each simulation through time, including how many generations it took until only one color was left. The exercise can be twisted in many cool ways (bottle-necks, gene flow and drift, selection, etc.).

Adding a nice interface, allowing students to change the number of alleles (or other basic parameters) could improve learning and make it more accesible to students that are not learning about R programming.

---
Repository manager:

Corey Devin Anderson

Associate Professor of Biology

Valdosta State University

coreanderson@valdosta.edu
