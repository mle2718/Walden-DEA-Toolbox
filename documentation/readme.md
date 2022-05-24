# Overview

This project is a repository for most of the DEA models I have worked with over the past 20+ years. Most of the models are written in R, although some older models I developed in GAMS will also be included. I am also including documents I have published over the years which used some of the DEA models found here. If you just want to use a well developed package without doing all the detailed programming, I highly recommend the Benchmarking package. If you have a problem that isn't easily accomplished using the Benchmarking package, then the programs found here should allow you to build whatever DEA model you like. For many of the programs, I use the Charnes1981 data supplied with the Benchmarking package because it is a publicly available data set and allows replication of results.

# A digression about solvers.

When I first moved from GAMS to R for solving DEA models, I used the LpSolveAPI package. Subsequentle, I switched to the Rglpk package because I found the Syntax easier. However, for really large data files, I found that the LpSolveAPI package ran significantly faster. 

Recently, a package OMPR was developed and released which I believe makes programming DEA models much easier. I provide an example of this in the Output oriented Technical Efficiency section. The package allows the user to program the DEA problem using equations. The user can export the structure of the model and then use other solvers to solve the DEA model. The reason I recommend exporting the structure of the model and solving outside of the OMPR package is so that the model structure doesn't need to be revised at each iteration of the DEA model.  


# Description of the folders

Output Oriented Technical Efficiency - Holds programs which solve the output oriented technical efficiency problem.

Input Oriented Technical Efficiency - Holds programs which solve the input oriented technical efficiency problem.
