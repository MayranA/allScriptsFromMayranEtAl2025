# Quantification of CDH1 and CDH2 on membrane

First, lif images containing multiple images were splitted using [this macro](./0-SplitLIF.ijm).

Then, [this macro](./1-Macro-Analyse-Coloc.groovy) is run for each image. It uses ilastik to classify pixels as being part of the membrane or not. The ilastik project is [here](./MembraneClassifier.ilp). The outputs have been copied to [this directory](../../output.files/imaging/membrane/groovy_outputs/).

The plots were generated with [this quarto](./Pixel.Membrane.measurements.to.Plots.qmd) and the outputs can be found [here](../../output.files/imaging/membrane/).
