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
parent=File.getParent(input);
dlist=getFileList(parent);
for (j=0; j < dlist.length; j++){
	roi_image = parent+"/"+dlist[j]+"dna.tif";
	inputp = parent+"/"+dlist[j];
	list=getFileList(parent+"/"+dlist[j]);
// The DNA image is opened, coverted to 8-bit and thresholded.
//	If the default threshold needs to be overwriten, comment out the setAutoThreshold() and uncomment the run("Threshold") and setThreshold() lines accordingly
	open(roi_image);
	run("RGB Color");
	selectWindow("dna.tif");
	run("Close");
	run("8-bit");
	setAutoThreshold("Default dark");
	//run("Threshold...");
	setThreshold(40, 255);
	setOption("BlackBackground", false);
	run("Convert to Mask");

//The ROIs are rounded and split
// If you are running a high magnification (>10x) DNA image, it is recommended that you comment this out to avoid nuclear image fragementation
	run("Watershed");

//The ROIs are generated
	run("Analyze Particles...", "size=10-Infinity pixel exclude add");
//The image is closed
	close();
// Then a Nuclear and WholeCell directory is made, if not already present
	if(!File.isDirectory(inputp+"/Nuclear")){
		File.makeDirectory(inputp+"/Nuclear");
	};
	if(!File.isDirectory(inputp+"/WholeCell")){
		File.makeDirectory(inputp+"/WholeCell");
	};
	inputn = inputp+"/Nuclear/";
	inputc = inputp+"/WholeCell/";
// Then you generate a list of the images in that directory:
	for (i=0; i<list.length; i++){
		// Each tif is opened in turn, analyzed based on the nuclear ROIs and saved in the Nuclear folder
		if(endsWith(list[i], ".tif")){
			open(parent+"/"+dlist[j]+list[i]);
			run("RGB Color");
			selectWindow(list[i]);
			run("Close");
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
			open(parent+"/"+dlist[j]+list[i]);
			run("RGB Color");
			selectWindow(list[i]);
			run("Close");
			run("8-bit");
			// Change this for deviations from default
			setAutoThreshold("Default dark");
			//run("Threshold...");
			if (startsWith(list[i], "dna")){
				setThreshold(30, 255);
				setOption("BlackBackground", false);
				run("Convert to Mask");
				run("Watershed");
				run("Analyze Particles...", "size=10-Infinity pixel add");
			} else if (startsWith(list[i], "n")){
				setThreshold(25, 255);
				setOption("BlackBackground", false);
				run("Convert to Mask");
				run("Watershed");
				run("Analyze Particles...", "size=1-Infinity pixel add");
			} else if (startsWith(list[i], "edu")){
				setThreshold(20, 255);
				setOption("BlackBackground", false);
				run("Convert to Mask");
				run("Watershed");
				run("Analyze Particles...", "size=1-Infinity pixel add");
			} else {
				setAutoThreshold("Default dark");
				setOption("BlackBackground", false);
				run("Convert to Mask");
				run("Watershed");
				run("Analyze Particles...", "size=1-Infinity pixel add");
			}
			close();
			open(inputp+list[i]);
			run("RGB Color");
			selectWindow(list[i]);
			run("Close");
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
};
