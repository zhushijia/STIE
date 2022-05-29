https://github.com/almaan/stereoscope#reprodmb



#conda env list
#conda env remove --name myVE.python3.8

module load python/3.8.x-anaconda
#conda create -n myVE.python3.8 anaconda
source activate myVE.python3.8
#source deactivate

module load gcc/6.3.0
pip install -U loompy
cd /archive/SCCC/Hoshida_lab/s184554/Code/github
git clone https://github.com/almaan/stereoscope 
cd stereoscope
python setup.py install # require python >= 3.8
stereoscope test

cd data
mkdir raw curated
cd raw
curl -O https://storage.googleapis.com/linnarsson-lab-loom/l1_hippocampus.loom  


../../preprocess/hippocampus/create-mod-hippo-loom.py l1_hippocampus.loom .
../../preprocess/hippocampus/subsample-data.py -lf mod_l1_hippocampus.loom -o ../curated -lb 25 -ub 250 -cn "Celltype_1"


7z e mouse-st-data.zip
mv st-hippo*tsv curated/
ls -1 curated/*


# this is the result folder
cd ../../res 

# help
#stereoscope run -h

# run deconvolution
stereoscope run \
--sc_cnt ../data/curated/*cnt*.tsv \
--sc_labels ../data/curated/*mta*.tsv \
-sce 75000  -o hippo_1 -n 5000 \
--st_cnt ../data/curated/st-hippo*tsv \
-ste 75000 -stb 100 -scb 100


stereoscope run \
--sc_cnt ../data/curated/*cnt*.tsv \
--sc_labels ../data/curated/*mta*.tsv \
--st_cnt ../data/curated/st-hippo*tsv \
-o hippo_1



#########################################################################################################






























