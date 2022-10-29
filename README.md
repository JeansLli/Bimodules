# It's for "On rectangle-decomposable 2-parameter persistence modules"

There are 2 different versions for bimodule decomposition. The first one is one-to-one version and the other one is quantization version. One-to-one version is implemented by python, and quantization version is implemented by both python and C++.

## Python environment
First, install [anaconda](https://docs.anaconda.com/anaconda/user-guide/getting-started/) and create virtual environments for this code. 

## Create data and scc2020 format
To create data, let's run
```
cd data
python rips_filtration_from_GMM.py 
```
Some of the code in this folder is from [https://bitbucket.org/mkerber/cgta_paper_2021/src/master/](https://bitbucket.org/mkerber/cgta_paper_2021/src/master/).

For more details about scc2020 format, please read `format_scc2020.pdf` in this folder.


## For cpp version code:  
First enter the folder using 
```
cd cpp_version
```
Then compile it:
```
mkdir build
cmake ..
make
```

Then the command is like:  
```
./bimodules_quantization x y rank_dim input_file  
```
`(x,y)` is the grid size, `rank_dim` is the homology degree.

For example: 
```
./bimodules_quantization 20 20 1 ../../data/function_rips_with_threshold_1000_1.scc
```

After running the code, we can visualize the result using
```
python visualize_barcodes.py --input  ./result/dim_0_barcodes_1025.txt
```
In the file `visualize_barcodes.py`, you can set different modes of visualization.

If you want to visualize the result dynamically, you can try
```
python visualize_barcodes_mouse_move.py --input ./result/dim_0_barcodes_121575.txtt
```


