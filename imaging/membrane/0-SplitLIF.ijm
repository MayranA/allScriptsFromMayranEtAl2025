/*
 * = AUTHOR INFORMATION =
 * Code written by Nicolas Chiaruttini, EPFL - SV - PTECH - BIOP
 */
 
#@File(label = "Lif file to extract") file
#@File(style="directory") outputPath


// Clears images
run("Close All");

// CLears results window
run("Clear Results");

// Set standard measurements
run("Set Measurements...", "area mean standard min display redirect=None decimal=3");

run("Bio-Formats Importer", 
	"open="+file+" "
	+"color_mode=Composite "
	+"open_all_series "
	+"rois_import=[ROI manager] "
	+"view=Hyperstack stack_order=XYCZT");

while (nImages>0) {
	previousImageNumber = nImages;
	// Clears ROI Manager
	roiManager("reset");	
	currentImage3D = getTitle();
	currentImage3D = replace(currentImage3D,"/","_");
	rename(currentImage3D);
	selectImage(currentImage3D);
	titleProj = getTitle();
	Stack.setChannel(1); run("Grays");
	Stack.setChannel(2); run("Cyan");
	Stack.setChannel(3); run("Magenta");
	Stack.setChannel(4); run("Green");
	saveAs("Tiff", outputPath+File.separator+titleProj);

	close();

	while (nImages>previousImageNumber-1){ 
	  close();
	}
}




