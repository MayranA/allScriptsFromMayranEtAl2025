# toolshed.g2.bx.psu.edu/repos/bgruening/deeptools_bigwig_average/deeptools_bigwig_average/3.5.4+galaxy0
# command_version:bigwigAverage 3.5.4
bigwigAverage --numberOfProcessors "${GALAXY_SLOTS:-4}" --bigwigs bigwig normalized per million reads in peaks.bigwig bigwig normalized per million reads in peaks.bigwig --outFileName 'bigwigAverage on data 2391 and data 2390.bigwig' --outFileFormat 'bigwig'     --scaleFactors '1' --binSize 10
