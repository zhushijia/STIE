filelist = getArgument();
file = split(filelist,'#');
image_path = file[0];
feature = image_path + ".feature";
image_window = File.getName(image_path)

print(image_path);
print(feature);
print(image_window);

open(image_path);
run("DeepImageJ Run", "model=[Multi-Organ Nucleus Segmentation (StarDist 2D)] format=Tensorflow preprocessing=[per_sample_scale_range.ijm] postprocessing=[StarDist_Post-processing.ijm] axes=Y,X,C tile=3600,3600,3 logging=normal");
run("Set Measurements...", "area mean standard modal min centroid center perimeter bounding fit shape feret's integrated median skewness kurtosis area_fraction stack redirect=None decimal=3");
selectWindow(image_window);
roiManager("Measure");
saveAs("Results",feature);


newImage("mask", "16-bit black", getWidth(), getHeight(), 1);
for (index = 0; index < roiManager("count"); index++) {
	roiManager("select", index);
	setColor(index+1);
	fill();
}
resetMinAndMax();
run("glasbey");

maskJpg = image_path + ".mask.jpg";
maskTxt = image_path + ".mask.txt";

selectWindow("mask");
saveAs("Jpeg",maskJpg);
saveAs("Text Image",maskTxt);

run("Quit");
