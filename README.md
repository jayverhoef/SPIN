
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

[![minimal R version](https://img.shields.io/badge/R%3E%3D-4.0.0-6666ff.svg)](https://cran.r-project.org/) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/kotzeb0912)](https://cran.r-project.org/package=kotzeb0912) [![packageversion](https://img.shields.io/badge/Package%20version-1.0.0-orange.svg?style=flat-square)](commits/master)

[![Last-changedate](https://img.shields.io/badge/last%20change-2023--06--06-yellowgreen.svg)](/commits/master)

# SPIN
## An R package in support of publication, "Indexing and Partitioning the Spatial Linear Model for Large Data Sets." 

#### Jay M. Ver Hoef<sup>a</sup>, Michael Dumelle, Matt Higham, Erin Peterson, and Daniel Isaak

#### <sup>a</sup>NOAA Fisheries (NMFS) Alaska Fisheries Science Center, Marine Mammal Laboratory, Seattle, WA

As a scientific work, and in keeping with common scientific practicies, I kindly request that you cite my research project and applicable publications if you use my work(s) or data in your publications or presentations. Additionally, I strongly encourage and welcome collaboration to promote use of these data in the proper context and scope.  The publication is currently submitted for publication in <em>PLOS ONE</em>:

#### Ver Hoef, J.M., Dumelle, M., Higham, M., Peterson, E.E., and Isaak, D.J. 2023. Species Density Models from Opportunistic Citizen Science Data. Submitted, <em>PLOS ONE</em>.


Abstract
-----------------

We consider four main goals when fitting spatial linear models: 1) estimating covariance parameters, 2) estimating fixed effects, 3) kriging (making point predictions), and 4) block-kriging (predicting the average value over a region).  Each of these goals can present different challenges when analyzing large spatial data sets.  Current research uses a variety of methods, including spatial basis functions (reduced rank), covariance tapering, etc, to achieve these goals.  However, spatial indexing, which is very similar to composite likelihood, offers some advantages.  We develop a simple framework for all four goals listed above by using indexing to create a block covariance structure and nearest-neighbor predictions while maintaining a coherent linear model. We show exact inference for fixed effects under this block covariance construction. Spatial indexing is very fast, and simulations are used to validate methods and compare to another popular method. We study various sample designs for indexing and our simulations showed that indexing leading to spatially compact partitions are best over a range of sample sizes, autocorrelation values, and generating processes.  Partitions can be kept small, on the order of 50 samples per partition.  We use nearest-neighbors for kriging and block kriging, finding that 50 nearest-neighbors is sufficient.  In all cases, confidence intervals for fixed effects, and prediction intervals for (block) kriging, have appropriate coverage. Some advantages of spatial indexing are that it is available for any valid covariance matrix, can take advantage of parallel computing, and easily extends to non-Euclidean topologies, such as stream networks.  We use stream networks to show how spatial indexing can achieve all four goals, listed above, for very large data sets, in a matter of minutes, rather than days, for an example data set.
Installation
------------

Installation of this R code and data package is done through the `devtools::install_github()` function or by downloading the [source package from the latest release](https://github.com/jayverhoef/POP). 


```
library("devtools")
install_github("jayverhoef/SPIN")
```


-------------
##### Disclaimer

<sub>This repository is a scientific product and is not official communication of the Alaska Fisheries Science Center, the National Oceanic and Atmospheric Administration, or the United States Department of Commerce. All AFSC Marine Mammal Laboratory (AFSC-MML) GitHub project code is provided on an ‘as is’ basis and the user assumes responsibility for its use. AFSC-MML has relinquished control of the information and no longer has responsibility to protect the integrity, confidentiality, or availability of the information. Any claims against the Department of Commerce or Department of Commerce bureaus stemming from the use of this GitHub project will be governed by all applicable Federal law. Any reference to specific commercial products, processes, or services by service mark, trademark, manufacturer, or otherwise, does not constitute or imply their endorsement, recommendation or favoring by the Department of Commerce. The Department of Commerce seal and logo, or the seal and logo of a DOC bureau, shall not be used in any manner to imply endorsement of any commercial product or activity by DOC or the United States Government.</sub>

