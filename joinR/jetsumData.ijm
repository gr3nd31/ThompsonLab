// First the script lets you know which directory you're working in:

roi_image = File.openDialog("Choose a File");
input=getDirectory("current");
list=getFileList(input);

//
open(roi_image);
close();

for (i=0; i<list.length; i++){
	if(endsWith(list[i], ".tif")){
		open(input+list[i]);
		run("8-bit");
		setAutoThreshold("Default dark");
		//run("Threshold...");
		//setThreshold(41, 255);
		setOption("BlackBackground", false);
		run("Convert to Mask");
		//run("Watershed");
		run("Analyze Particles...", "size=0.05-Infinity display exclude add");
		//saveAs("Tiff", input+replace(list[i], ".tif", "_labeled.tif"));
		close();
		open(input+list[i]);
		roiManager("measure");
		saveAs("Results", input+replace(list[i], ".tif", ".csv"));
		selectWindow("Results");
		run("Close");
		close();
	};
};
selectWindow("ROI Manager");
run("Close");