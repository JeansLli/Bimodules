# It's for "On rectangle-decomposable 2-parameter persistence modules"

There are 2 different versions for bimodule decomposition. The first one is one-to-one version and the other one is quantization version. One-to-one version is implemented by python, and quantization version is implemented by both python and C++.


For cpp version code:  
the command is like  ./bimodules_quantization x y rank_dim input_file  
(x,y) is the grid size, rank_dim is the homology degree.

For example: ./bimodules_quantization 20 20 1 ../../data/function_rips_with_threshold_1000_1.scc

After running the code, we can visualize the result using

python visualize_barcodes.py --input  ./result/dim_0_barcodes_1025.txt


rips_filtration_from_GMM.py is to create the input data.