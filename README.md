# It's for "On rectangle-decomposable 2-parameter persistence modules"

For cpp version code:

the command is like  ./bimodules_quantization x y rank_dim input_file
(x,y) is the grid size, rank_dim is the homology degree.

For example: ./bimodules_quantization 20 20 1 ../../data/function_rips_with_threshold_1000_1.scc


For barcodes visualization:

python visualize_barcodes.py --input  ./result/dim_0_barcodes_1025.txt