# ABEL INVERSION

Is the inverse of the Abel Transform, which is an integral transform created by Niels Henrik Able. It is used in the analysis of spherically symmetric functions. It is a commonly used in interferometry and polarimetry as plasmas often have the property that they are cylindrically summetric.

This code, draws from Principles of Plasma Diagnostics by I.H. Hutchinson and code written by Pat Carle. This along with other open resources such as matlab and wikipedia.

From my research, the two largest issues that this deals with are that the derivative tends to amplify any noise in the data. In our case, as this is experimental data, smoothing is necessary. Secondly, a singularity occurs at /( y = r)/ which can introduce errors when making adjusments. Third, it is assumed that the data and the reconstruction tends to zero as y gets very large.

## Inputs

Data of some sort (I am guessing from a PI machine)

## Outputs

2D reconstruction of the plasmas
