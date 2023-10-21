# FetalSSFSESimulation - A modified FaBiAN phantom

## Summary
A silmulation of single-shot Fast spin echo sequences on a simulated moving fetal brain for a Philip's system. This simulation is based on Hélène Lajous' FaBiAN simulation:
*Code:* 
https://github.com/Medical-Image-Analysis-Laboratory/FaBiAN https://doi.org/10.5281/zenodo.5471094 
*Paper:*
https://doi.org/10.1038/s41598-022-10335-4

We use the T1 and T2 values suggested by Lajous and colleagues and the standard atlas by Gholipour and colleagues ( Atlas: http://crl.med.harvard.edu/research/fetal_brain_atlas/  Paper: https://doi.org/10.1038/s41598-017-00525-w ) 

We deal with smaller matrices (thus lower RAM usage) by going slice by slice instead of calculating on the entire matrix. In addition, we use linear interpolation on individual tissue classes followed by max voting for transformations and we model the slice profiles using Bloch simulations and then integrate to get an overall signal of the slice rather than weighting the slice with a Gaussian function.

## Current version
We currently assume perfect recovery of all slices, with no slice cross-talk, and no in-slice motion or B1 inhomogeneites as this was not within our focus in our paper. In addition, we only use Philip's parameters. The code requires using RF pulses that will give the appropriate slice profiles. The example shows how to deal with the RF pulses. 

## Acknwoledgements 
This work was supported by funding from the EPSRC Centre for Doctoral Training in Smart Medical Imaging (EP/S022104/1), the core funding from the Wellcome/EPSRC Centre for Medical Engineering [WT203148/Z/16/Z] and by the National Institute for Health Research Clinical Research Facility. 

## Dependencies
We use S.J Malik's EPG-X instead of Weigel's (used in the original FaBiAN phantom). This library includes the code from that library, but find the full library here: https://github.com/mriphysics/EPG-X and the paper is here: https://doi.org/10.1002/mrm.27040
We also use the edit-image function from MIRTK (https://mirtk.github.io/commands/edit-image.html) to make it usable with SVRTK (https://github.com/SVRTK/SVRTK) so for best results I recommend installing this.

## Publication
This code was made to help with my personal publication:
https://doi.org/10.1007/978-3-031-45544-5_4 
If this code was of use for you, please cite the above publication as well as Hélène's original paper: https://doi.org/10.1038/s41598-022-10335-4

## Contact
Please feel free to get in contact with me with any queries regarding this code:
Email: suryava.bhattacharya@kcl.ac.uk
