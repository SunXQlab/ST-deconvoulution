### stereoscope usage guide
### systerm requirment: Linux/Unix, python_requires >= 3.5.0
### the package installation refer to https://github.com/almaan/stereoscope

### Step-1: enter the following in the terminal
git clone https://github.com/almaan/stereoscope 

### Step-2: enter the stereoscope folder
cd stereoscope

### Step-3: install dependencies of stereoscope
./setup.py install

###Step-4: enter the following in the terminal
cd data

cd curated

cd ../../res

stereoscope run --sc_cnt ../data/curated/*cnt*.tsv --sc_labels ../data/curated/*mta*.tsv -sce 75000  -o hippo_1 -n 5000 --st_cnt ../data/curated/st-hippo*tsv -ste 75000
