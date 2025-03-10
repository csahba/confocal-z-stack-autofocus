# confocal-z-stack-autofocus
A matlab data processing pipeline to convert a stack of z slices from a confocal microscope into a single focused image, only for purposes of tracking cells that may drift in the z axis over time. As it was adequate for my needs, more accurately as it stands this a form of cell outlining rather than actual focussing.

This is an old script I wrote when I ran experiments that required long time lapses between imaging, over which cells were prone to moving slightly in the z axis. As the microscope I was using did not have any autotracking, I set it to automate a series of z slices to span the space the cells would be in, initially to allow me to manually select focused slices.

I then noticed a quirk of confocal microscopy that allowed me to apply simple but effective maths to effectively autofocus. As focus changes, the cell acts as a lens, casting a shadow or focusing light, causing distinct changes in phase contrast in reverse directions at the cell border and in the cell body. This script fits a simple linear gradient for each pixel stack in the z axis, and uses the array of these gradients to remove noise and provide a clear image.

![z stack animation](https://github.com/csahba/confocal-z-stack-autofocus/raw/master/z%20focussing.gif)
Here is an animation of a z stack, as we descend through focus you can see the change in contrast around the cells, this was the feature I exploited.

My initial intention was to use the intersection of the fitted linear regression between the cell body and the cell wall/'halo' to select z slices in focus, but this turned out to not be necessary for my needs, though there is an imperfect version in the script and a pipeline involving creating a mask and selecting clusters of cells.

Having written this in matlab, I hope to come back and port it to python as an exercise.

imreadBF and imreadBFmeta Copyright (c) 2011, Christoph Moehl copyright and license information available here https://uk.mathworks.com/matlabcentral/fileexchange/32920-imread-for-multiple-life-science-image-file-formats#license_modal

loci_tools.jar from openmicroscopy.org license information here
https://www.openmicroscopy.org/bio-formats/downloads/
