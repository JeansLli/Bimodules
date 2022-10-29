# The code is for rank decomposition and signed barcode.

There are 2 different versions for bimodule decomposition. The first one is one-to-one version and the other one is quantization version. One-to-one version is implemented by python, and quantization version is implemented by both python and C++.

## Python environment
First, install [anaconda](https://docs.anaconda.com/anaconda/user-guide/getting-started/) and create virtual environments. 

## Create data of scc2020 format
To create data, let's run
```
cd data
python rips_filtration_from_GMM.py 
```
Some of the code in this folder is from [https://bitbucket.org/mkerber/cgta_paper_2021/src/master/](https://bitbucket.org/mkerber/cgta_paper_2021/src/master/).

For more details about scc2020 format, please read `format_scc2020.pdf` in this folder.


## For cpp version code:  
First enter the folder:
```
cd cpp_version
```
Then compile it and generate an executable file `bimodules_quantization`:
```
mkdir build
cd build
cmake ..
make
```

Then the run the executable file and the command is like:  
```
./bimodules_quantization x y rank_dim input_file  
```
`(x,y)` is the grid size, `rank_dim` is the homology degree.

For example: 
```
./bimodules_quantization 20 20 1 ../../data/function_rips_GMM.txt
```

After running the code, we can visualize the result using
```
python visualize_barcodes.py --input  ./result/dim_1_barcodes_200020.txt
```
In the file `visualize_barcodes.py`, different modes of visualization could be set.

You may also want to try the dynamic visualization:
```
python visualize_barcodes_mouse_move.py --input ./result/dim_1_barcodes_200020.txt
```


