# NK_hikes
This repository contains supporting data and codes for the paper
### Simulating short- and long-term evolutionary dynamics on rugged landscapes
by L. Trujillo, P. Banse and G. Beslon (2021)
##### Abstract
We propose a minimal model to simulate long waiting times followed by evolutionary bursts on rugged landscapes. It combines point and inversions-like mutations as sources of genetic variation. The inversions are intended to simulate one of the main chromosomal rearrangements. Using the well-known family of NK fitness landscapes, we simulate random adaptive walks, i.e. successive mutational events con- strained to incremental fitness selection. We report the emergence of different time scales: a short-term dynamics mainly driven by point mutations, followed by a long-term (stasis- like) waiting period until a new mutation arises. This new mutation is an inversion which can trigger a burst of successive point mutations, and then drives the system to new short-term increasing-fitness period. We analyse the effect of genes epistatic interactions on the evolutionary time scales. We suggest that the present model mimics the process of evolutionary innovation and punctuated equilibrium.

The preprint is available on [arXiv](https://arxiv.org/).

#### Repository content
1 Codes:
- nk_walk: Main routine to simulate adative walks on the Stuart Kauffman's NK-fitness landscape model, with point mutations (as usual) and inversions (as new). 
-  Python codes to analyse the data and plot the figures reported in the paper.

2 Data:
-  Bursts candidates (figure 6)
- Hikes1 (figures 3 and 4)
- Hikes2 (figure 5)
