//YggData is the ThompsonLab ImageJ macro for single cell analysis of IFA images
//	This macro first asks you to pick a nucleus picture to generate ROIs from
//	It then uses these ROIs to determine the raw NUCLEAR data from each image within the directory, and stores these CSVs in a newly created 'Nuclear' folder
//	Ygg will then back out and generate ROIs for each image independent of nuclear localization, and store these CSVs in a newly created 'WholeCell' folder
//	Follwing this macro, the imaGen() macro should be run in sirmixaplot() (version 7.2+) to combine all the data into a new dataset in which each row is a single cell

// First, lets make sure the measurements are set for appropriate analyses. If changes are made here, they will not be included in the default imaGen() compilation
run("Set Measurements...", "area mean standard modal min centroid center perimeter feret's integrated median area_fraction display redirect=None decimal=9");

// First the script lets you know which file is the nucleus, 
// and designate directory you're working in:
roi_image = File.openDialog("Choose a File");
input=getDirectory("current");
list=getFileList(input);
// The DNA image is opened, coverted to 8-bit and thresholded.
//	If the default threshold needs to be overwriten, comment out the setAutoThreshold() and uncomment the run("Threshold") and setThreshold() lines accordingly
open(roi_image);
run("8-bit");
setAutoThreshold("Default dark");
//run("Threshold...");
//setThreshold(10, 255);
setOption("BlackBackground", false);
run("Convert to Mask");

//The ROIs are rounded and split
// If you are running a high magnification (>10x) DNA image, it is recommended that you comment this out to avoid nuclear image fragementation
//run("Watershed");

//The ROIs are generated
run("Analyze Particles...", "size=.08-Infinity exclude add");
//The image is closed
close();
// Then a Nuclear and WholeCell directory is made, if not already present
if(!File.isDirectory(input+"/Nuclear")){
	File.makeDirectory(input+"/Nuclear");
};
if(!File.isDirectory(input+"/WholeCell")){
	File.makeDirectory(input+"/WholeCell");
};
inputn = input+"/Nuclear/";
inputc = input+"/WholeCell/";
// Then you generate a list of the images in that directory:
for (i=0; i<list.length; i++){
	// Each tif is opened in turn, analyzed based on the nuclear ROIs and saved in the Nuclear folder
	if(endsWith(list[i], ".tif")){
		open(input+list[i]);
		roiManager("measure");
		saveAs("Results", inputn+replace(list[i], ".tif", ".csv"));
		selectWindow("Results");
		run("Close");
		close();
	};
};
//Then everything is closed
selectWindow("ROI Manager");
run("Close");

//Now that the nuclear data is collected, Ygg will collected the nuclear indepedent data
//	Since the list of images will be the same, it simply iterates through each image and collects the total ROI pixel data
for (i=0; i<list.length; i++){
	// Each tif is opened in turn, analyzed for general ROIs and saved in the WholeCell folder
	if(endsWith(list[i], ".tif")){
		open(input+list[i]);
		run("8-bit");
		// Change this for deviations from default
		setAutoThreshold("Default dark");
		//run("Threshold...");
		setThreshold(15, 255);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		// If you are running a high magnification (>10x) DNA image, it is recommended that you comment this out to avoid nuclear image fragementation
		run("Watershed");
		// Alter this for different ROI sizes
		run("Analyze Particles...", "size=.008-Infinity add");
		//saveAs("Tiff", input+replace(list[i], ".tif", "_labeled.tif"));
		//selectWindow("Results");
		//run("Close");
		close();
		open(input+list[i]);
		roiManager("measure");
		saveAs("Results", inputc+replace(list[i], ".tif", ".csv"));
		//selectWindow("Results");
		run("Close");
		roiManager("reset")
		if (nImages>0) {
			close();
		};
	};
};
//Everything is closed
selectWindow("ROI Manager");
run("Close");