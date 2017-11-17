# Kernel Equating Without Pre-Smoothing

Maintainer: Cambridge Assessment, Assessment, Research and Development (ARD) <benton.t@cambridgeassessment.org.uk>

Authors: Tom Benton

## About

An R package to carry out kernel equating without needing to choose a loglinear model for pre-smoothing.
This is desirable since, although pre-smoothing may possibly improve accuracy, skipping this step allows kernel equating to become 
more straightforward and easier to automate as there are no decisions to make about the order of loglinear models.
In addition most of the functions in this package are designed to work in cases where non-integer scores occur.
In other words, for most functions, it is not necessary for test scores to be integers.

As an alternative to pre-smoothing the package relies on the kernel within kernel equating itself to supply sufficient smoothing 
to provide accurate results.
This may be done either using the plug-in formula of Andersson and von Davier (2014) or by using the method of cross-validation suggested
by Liang and von Davier (2014).

For the puposes of post-stratification equating (PSE) pre-smoothing would require not only fitting a model to each marginal score
distribution but also to the joint relationship between test and anchor forms. In this package this is avoided.
Instead, data from different populations is weighted so that the first n moments of scores on the anchor test match.
This is achieved using the Adjustment by Minimum Discriminant Information method also known as MDIA (Haberman, 1984).
Functions to perform MDIA weighting are included within this package and may potentially be used in contexts aside from equating.

Aside from post-stratification and chained equating, this package also implements the hybrid equating method decribed 
in von Davier and Chen (2013). 
 
Finally this package implements two experimental equating methods that were devised by training neural networks
to replicate known "true" equating functions within the nonequivalent groups with anchor test (NEAT) design.
See Benton (2017) for more details.

**References**

Andersson, B., & von Davier, A. A. (2014). Improving the bandwidth selection in kernel equating. 
*Journal of Educational Measurement, 51*(3), 223-238.

Benton, T. (2017). Can AI learn to equate?, 
*presented at the International Meeting of the Psychometric Society, Zurich, 2017*. Cambridge, UK: Cambridge Assessment.

Haberman, S. J. (1984). Adjustment by minimum discriminant information. 
*The Annals of Statistics, 12*(3), 971-988.

Liang, T., & von Davier, A. A. (2014). Cross-validation: An alternative bandwidth-selection method in kernel equating. 
*Applied Psychological Measurement, 38*(4), 281-295.

von Davier, A. A., & Chen, H. (2013). The Kernel Levine Equipercentile Observed-Score Equating Function. ETS Research Report Series.

## Installation
This package can be installed using the devtools package using the commands below.



library(devtools)



install_github("CambridgeAssessmentResearch/KernEqWPS")



In order to run all of the examples provided with this package you will also need to install the ggplot2 package.
To see examples of how the functions in the package work type.

library(KernEqWPS)

help("KernEqWPS-package")

## License

The MIT License (MIT)

Copyright (c) 2017 Cambridge Assessment

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
