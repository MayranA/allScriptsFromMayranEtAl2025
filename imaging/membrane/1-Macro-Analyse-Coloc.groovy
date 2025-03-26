/*
* * = AUTHOR INFORMATION =
 * Code written by Nicolas Chiaruttini, EPFL - SV - PTECH - BIOP
 */
#@File imagePath
#@int channel1
#@int channel2
#@File ilastik_membrane_classifier_file
#@File(label = "", style = "directory") outputFolder

#@RoiManager rm
#@Double(label = "Pre filter images - median filter") median_radius

#@Integer fluo_bins = 255
#@Integer fluo_min = 0
#@Integer fluo_max = 255

#@CommandService cs

def image = IJ.openImage(imagePath.getAbsolutePath())

image.setOverlay(new Overlay())

// We build a synthetic membrane channel by adding both channels
def ori_image_name = image.getTitle();

image.setRoi(null)

// Makes the sum image

def image_channel1 = new Duplicator().run(image, channel1, channel1, 1, image.getNSlices(), 1, 1) 
def image_channel2 = new Duplicator().run(image, channel2, channel2, 1, image.getNSlices(), 1, 1) 

def sum_image = ImageCalculator.run(image_channel1, image_channel2, "Add create 32-bit stack");


// Which we classify : we only keep the pixels classified as belonging to the membrane
def predictions_Img = cs.run(IlastikPixelClassificationCommand.class, false,
		"projectFileName", ilastik_membrane_classifier_file,
		"inputImage", sum_image,
		"pixelClassificationType","Segmentation"
	).get().getOutput("predictions")
	
def predictions = new Duplicator().run(ImageJFunctions.wrap(predictions_Img, "predictions")) // wraps the ImgPlus as an ImagePlus, and duplicate because otherwise we can't operate on it (-1 *255)


// The mask returns value 1 = bg, 2 = membrane, 
// We want to convert it ot an IJ1 mask : 0 = BG, 255 = Object -> subtract 1 and multiply by 255
IJ.run(predictions, "Subtract...", "value=1 stack"); 
IJ.run(predictions, "Multiply...", "value=255 stack"); 

def slices_img = predictions.getNSlices()

Prefs.blackBackground = true // IJ specific thing

rm.reset() // clears ROI manager

// Median filter
if (median_radius>0) {
	rf = new RankFilters();
	for (iSlice = 1; iSlice<=slices_img; iSlice++) {
		rf.rank(image_channel1.getStack().getProcessor(iSlice), median_radius, RankFilters.MEDIAN)
		rf.rank(image_channel2.getStack().getProcessor(iSlice), median_radius, RankFilters.MEDIAN)
	}
}

def imgs = [] // To store image per slice

for (iSlice = 1; iSlice<=slices_img; iSlice++) {
	IJ.log("iSlice = "+iSlice)
	def currentProcessor = predictions.getStack().getProcessor(iSlice)
	def max_img = currentProcessor.getStatistics().max

	
	if (max_img>0) {
		
		currentProcessor.setThreshold(255, 255, ImageProcessor.NO_LUT_UPDATE);
		def roiAnalysis = new ThresholdToSelection().convert( currentProcessor )
		def current_roi = new ThresholdToSelection().convert( currentProcessor )
		current_roi.setPosition(1, iSlice, 1)
		image.getOverlay().add(current_roi)
		
		def slice_channel1 = new Duplicator().run(image_channel1, 1, 1, iSlice, iSlice, 1, 1) // run("Duplicate...", "duplicate");
		def slice_channel2 = new Duplicator().run(image_channel2, 1, 1, iSlice, iSlice, 1, 1) // run("Duplicate...", "duplicate");
		slice_channel1.setRoi(roiAnalysis)
		slice_channel2.setRoi(roiAnalysis)
		slice_channel1.setCalibration(image.getCalibration())
		slice_channel2.setCalibration(image.getCalibration())
		
		def reportPath = outputFolder.getAbsolutePath()+File.separator+"Median_"+median_radius+"_"+ori_image_name+"_slice_"+iSlice+"_Fluorogram.csv"
		
		def sliceReport = performColoc(slice_channel1, slice_channel2, roiAnalysis, iSlice, ori_image_name, reportPath)
		imgs.add(sliceReport)	
	}

}

predictions.close()

def result = Concatenator.run(imgs as ImagePlus[]);

new FileSaver(result).saveAsTiff(outputFolder.getAbsolutePath()+File.separator+"Median_"+median_radius+"_"+ori_image_name+"_Report.tiff")
new FileSaver(image).saveAsTiff(imagePath.getAbsolutePath())

public ImagePlus performColoc(slice_channel1, slice_channel2, roiAnalysis, iSlice, imageName, pathToSaveRawFluorogram) {

	slice_channel1.setRoi(roiAnalysis)
	slice_channel2.setRoi(roiAnalysis)
	
	def avg_channel1 = slice_channel1.getProcessor().getStats().mean
	def avg_channel2 = slice_channel2.getProcessor().getStats().mean
	
	def ic = new ImageColocalizer(slice_channel1, slice_channel2, slice_channel1.getCalibration())
	ic.thrA = 0
	ic.thrB = 0
	ic.SpearmanRank();
	ic.Pearson();
	ic.Areas();
	
	// Make Montage
	ArrayList<ImagePlus> imgs = new ArrayList<ImagePlus>();
	def columns = 3;
	def rows = 1;
	def stroke_width = 0f;
	
	imgs.add(ic.getRGBColocImage(stroke_width))
	imgs.add(ic.getRGBImageA(false, stroke_width))
	imgs.add(ic.getRGBImageB(false, stroke_width))
	
    ImagePlus montage;
    
    // Need to do a different thing if we are using stack or a single slice.        
    if(imgs.get(0).getNSlices() == 1) {
    	def montagestk = imgs.get(0).createEmptyStack();
        // Make a normal montage
	 	for(ImagePlus i : imgs) {
	 		montagestk.addSlice(i.getProcessor());
	 	}
	 	// Use Montage Maker
	 	MontageMaker mm = new MontageMaker();
    	montage = mm.makeMontage2(new ImagePlus("for montage",montagestk), columns, rows, 1.0, 1, imgs.size(), 1, 0, false);

    } 
	

	def is_auto_fluo = false
	
    //Add the fluorogram
    ic.addResult("Slice", iSlice);
    ic.addResult("Channel 1", channel1);
    ic.addResult("Channel 2", channel2);
    ic.addResult("Mean Channel 1", avg_channel1);
    ic.addResult("Mean Channel 2", avg_channel2);
    ic.addResult("Median_Radius", median_radius);
    ic.addResult("Ilastik Classifier", ilastik_membrane_classifier_file.getAbsolutePath())
   
    ic.CytoFluo();
   
    def fp = ic.getFluorogram() 
    
    if ((pathToSaveRawFluorogram!=null) && (pathToSaveRawFluorogram!="")) {
	    float[] valA = fp.getXValues();
		float[] valB = fp.getYValues();	    
		File f = new File(pathToSaveRawFluorogram)
		StringBuilder sb = new StringBuilder()
		for (int i = 0; i<ic.A.length; i++) {
            if (ic.M[i]>0) {
               sb.append(ic.A[i]+"\t"+ic.B[i]+"\n")
            }
        }
		
		f.write(sb.toString())		
    }
    

	ImagePlus fluo = (is_auto_fluo) ? ic.getFluorogramImage() : ic.getFluorogramImage(fluo_bins, fluo_min,fluo_max);
	
	ImagePlus scaledFluo;
	
	// Scale the fluorogram to the width or height of the image
   
    scaledFluo = Utils.scale(fluo, montage.getHeight());
	
    def flst = scaledFluo.getStack();
	// Make same number of dimensions as the montage (Z slices eventually)
	for(int i=1; i<montage.getStackSize(); i++) {
		flst.addSlice(scaledFluo.getProcessor().duplicate());
	}
	scaledFluo.setStack(flst);
	
	// Finally assemble them
	StackCombiner sc = new StackCombiner();
	
	ImageStack montagefluo;
    
	montagefluo = sc.combineHorizontally(montage.getStack(), scaledFluo.getStack());
	montage = new ImagePlus(imageName+" Report", montagefluo);
	
	
	// Make sure it's the title we want.
	montage.setTitle(imageName+" Report");
	montage.setCalibration(slice_channel1.getCalibration())
	
	ic.showResults();
	
	// And VOILA
	return montage
	
}

import ij.plugin.Concatenator
import ij.io.FileSaver
import ij.plugin.RGBStackMerge
import ij.ImageStack
import ij.plugin.MontageMaker
import ij.plugin.StackCombiner
import ch.epfl.biop.coloc.utils.Utils
import ij.gui.Overlay

import ch.epfl.biop.coloc.JACoP_B
import ch.epfl.biop.coloc.utils.ImageColocalizer
import ij.ImagePlus
import ij.plugin.Duplicator
import ij.plugin.ImageCalculator
import org.ilastik.ilastik4ij.ui.IlastikPixelClassificationCommand
import ij.IJ
import net.imglib2.img.display.imagej.ImageJFunctions
import ij.Prefs
import ij.process.ImageStatistics
import ij.process.ImageProcessor
import ij.plugin.filter.ThresholdToSelection 

import sc.fiji.coloc.algorithms.PearsonsCorrelation
import sc.fiji.coloc.algorithms.SpearmanRankCorrelation

import sc.fiji.coloc.Coloc_2
import ij.plugin.filter.RankFilters
import static ij.plugin.filter.RankFilters.MEDIAN