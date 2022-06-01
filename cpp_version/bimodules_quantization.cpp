#include <gudhi/graph_simplicial_complex.h>
#include <gudhi/Simplex_tree.h>
#include <gudhi/Persistent_cohomology.h>
 
#include <iostream>
#include <utility>  // for pair
#include <vector>
#include <fstream>
#include <string>
#include <sstream>
#include <ctime>

 
using Simplex_tree = Gudhi::Simplex_tree<>;
using Vertex_handle = Simplex_tree::Vertex_handle;
using Simplex_handle =  Simplex_tree::Simplex_handle;
using Filtration_value = Simplex_tree::Filtration_value;
using typeVectorVertex = std::vector<Vertex_handle>;
using typePairSimplexBool = std::pair<Simplex_tree::Simplex_handle, bool>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;

typedef std::vector<std::vector<std::vector<std::vector<int> > > > rank; 
typedef std::pair<int,int> entry;
typedef std::pair<entry,entry> bar;
typedef std::map<bar, int> barcode;
typedef std::map<Simplex_handle,entry> simplex_grid_map;

using std::cout; 
using std::cerr;
using std::endl; 
using std::string;
using std::ifstream; 
using std::vector;

struct simplex_node{
	float grid_x;
	float grid_y;
	typeVectorVertex face;
};


vector<entry> create_stair(entry elbow, entry end_point){
	vector<entry> stair;
	int i,j,k,l;
	
	if(elbow.first==0 || elbow.second==end_point.second){

		for(i=0; i<=end_point.second; ++i){
            stair.push_back(entry(0,i));
		}
		--i;
		for(j=1; j<=end_point.first; ++j){
            stair.push_back(entry(j,i));
		}
	}
    else{
    	for(i=0; i<=elbow.second; ++i){
            stair.push_back(entry(0,i));
		}
		--i;
		for(j=0; j<elbow.first; ++j){
            stair.push_back(entry(j+1,i));
		}
		--j;
		for(k=i+1; k<=end_point.second; ++k){
            stair.push_back(entry(j+1,k));
		}
		--k;
		for(l=j+2; l<=end_point.first; ++l){
            stair.push_back(entry(l,k));
		}
    }
    return stair;
}


entry proj_grid_to_path(entry simplex_node, entry elbow, int x_range,int y_range){
	int x = simplex_node.first;
	int y = simplex_node.second;
	int x0 = elbow.first;
    int y0 = elbow.second;
    entry proj_p;

    if((x==0 && y<=y0) || (y==y0 && x<=x0) || (y>=y0 && x==x0) || (y==y_range && x>=x0)){
        proj_p.first = x;
        proj_p.second = y;
    }
    else if(x<=x0 and y<=y0){
        proj_p.first = x;
        proj_p.second = y0;
    }
    else if(x<=x0 and y>=y0){
        proj_p.first = x0;
        proj_p.second = y;
    }
    else{
        proj_p.first = x;
        proj_p.second = y_range;
    }
	return proj_p;
}


int multp (rank& r, int i, int j, int k, int l) {
    // check boundaries 
    if(i<0 || j<0 || k>=r.size() || l>=r[0].size())
        return 0;

    // compute multiplicity
    int res =  r[i][j][k][l];
    if (k+1 < r.size())
        res -= r[i][j][k+1][l];
    if (l+1 < r[0].size())
        res -= r[i][j][k][l+1];
    if (k+1 < r.size() && l+1 < r[0].size())
        res += r[i][j][k+1][l+1];

    return res;
}


void compute_R_S_incl_excl(rank& r, barcode& b) {
    for (int i = 0; i<r.size(); i++){
        for (int j = 0; j<r[0].size(); j++){
            for (int k = i; k<r.size(); k++){
                for (int l = j; l<r[0].size(); l++) {
                    int m = multp(r, i, j, k, l) - multp(r, i-1, j, k, l)
                            - multp(r, i, j-1, k, l) + multp(r, i-1, j-1, k, l);
                    if (m != 0) {
                        bar bb(entry(i,j), entry(k,l));
                        b[bb] = m;
                    }
                }
            }
        }
    }
}


void zero_rank (rank& r, int n, int m) {  // n columns, m rows in the grid (column-dominant)
    r.resize(n);
    for (std::vector<std::vector<std::vector<int> > >& u: r){
        u.resize(m);
        for (std::vector<std::vector<int> >& v: u){
            v.resize(n);
            for (std::vector<int>& w: v) {
                w.resize(m);
                w.assign(m, 0);
            }
        }
    }
}

int main(int argc, char* const argv[]) {
    int x_range;
    int y_range;
    int rank_dim;
    string filename;
    //string filename("../../../cgta_paper_2021/function_rips/function_rips_with_threshold_100_1.scc");
    //string filename("../../data/fig3-2.txt");

    if(argc!=5){
        cout <<"ERROR: Please input the grid size(x_range,y_range), rank_dim and input file"<<endl;
        return 0;
    }
    else{
        x_range=std::stoi(argv[1])-1;
        y_range=std::stoi(argv[2])-1;
        rank_dim=std::stoi(argv[3]);
        filename=argv[4];
    }

    //for running time
    double time_update_simplex_tree=0;
    double time_barcodes=0;
    double time_rank_inv=0;
    double time_multiplicity=0;
    double time_total=0;
    clock_t start_time;
    clock_t end_time;

	
    int line_i = 1;
    int vertice_id = 0;
    int dF2, dF1, dF0, dummy; 
    vector<int> rank_num;
    vector<float> x_values;
    vector<float> y_values;
    vector<simplex_node> simplices;
    
    Simplex_tree st;

    entry end_point(x_range,y_range);
    simplex_grid_map sg_map;

    //Read file
    string line;
    vector<string> block_lines;

    ifstream input_file(filename);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '"
             << filename << "'" << endl;
        return EXIT_FAILURE;
    }


    getline(input_file, line);
    assert(line == "scc2020");
    
    //read the file
    while (getline(input_file, line)){
        if(line[0]!='#'){
        	if(line_i<2){
        		line_i+=1;
        		continue;
        	}
        	else if(line_i==2){
                std::stringstream ss(line);
                ss >> dF2 >> dF1 >> dF0 >> dummy; 
                rank_num.push_back(dF2);
                rank_num.push_back(dF1);
                rank_num.push_back(dF0);
                rank_num.push_back(dummy);
        		
    			line_i += 1;
    		}else{
    			block_lines.push_back(line);
        		line_i += 1;
        	}
        }
    }
    input_file.close();
    
    int start_id = 0;
    int simplices_len=block_lines.size();
    
    /////////////////// Store the input data as simplices
    for(int i=rank_num.size()-1; i>=0; --i){
    	//string simplices_i;
    	//cout << "rank_num:" << i << endl;
    	for(int j=simplices_len-start_id-rank_num[i]; j<simplices_len-start_id; ++j){
            //cout<<"block_lines[j]="<<block_lines[j]<<endl;
            std::stringstream ss(block_lines[j]);
            float grid_x, grid_y;
            char drop;
            int face1,face2,face3,face4;
            simplex_node simplex;
            
            if(i==3){ // 0-dim
                ss>>grid_x>>grid_y;
                simplex.face={vertice_id};
                vertice_id += 1;
            }else if(i==2){// 1-dim
                ss>>grid_x>>grid_y>>drop>>face1>>face2;
                simplex.face = {face1,face2};
            }else if(i==1){// 2-dim
                ss>>grid_x>>grid_y>>drop>>face1>>face2>>face3;
                simplex.face = {face1,face2,face3};
            }else if(i==0){// 3-dim
                ss>>grid_x>>grid_y>>drop>>face1>>face2>>face3>>face4;
                simplex.face = {face1,face2,face3,face4};
            }

            simplex.grid_x = grid_x;
            simplex.grid_y = grid_y;
            x_values.push_back(grid_x);
            y_values.push_back(grid_y);
    		simplices.push_back(simplex);
    	}
    	start_id += rank_num[i];
    }
    ///////////////////


    
    //////////////////Convert the face number to vertex number
    start_id = rank_num[2] + rank_num[3];
    int start_id_last;
    for(int i=rank_num.size()-3; i>=0; --i){
    	start_id_last = start_id - rank_num[i+1];
    	for(int j=start_id; j<(start_id+rank_num[i]); ++j){
    		simplex_node simplex = simplices[j];
    		typeVectorVertex face = simplex.face;
    		typeVectorVertex new_face;
    		for(int k=0; k<face.size(); ++k){
    			typeVectorVertex vertices = simplices[start_id_last+face[k]].face;
    			for(int l=0; l<vertices.size(); ++l){
    				int vertex = vertices[l];
    				if (std::find(new_face.begin(), new_face.end(), vertex) == new_face.end()){
    					new_face.push_back(vertex);
    				}
    			}
    		}
    		simplices[j].face=new_face;
    	}
    	start_id += rank_num[i];
    }
    //////////////////


    float x_min = *min_element(x_values.begin(),x_values.end()); 
    float y_min = *min_element(y_values.begin(),y_values.end()); 
    float x_max = *max_element(x_values.begin(),x_values.end()); 
    float y_max = *max_element(y_values.begin(),y_values.end()); 
    float x_interval = (x_max-x_min)/x_range;
    float y_interval = (y_max-y_min)/y_range;


    ///////////////// Convert the float grid coord to integer grid coord
    vector<simplex_node> simplices_integer=simplices;
    for(int i=0; i<simplices.size();++i){
    	float grid_x = simplices[i].grid_x;
    	float grid_y = simplices[i].grid_y;

    	if(x_interval==0){
    		simplices_integer[i].grid_x=0;
    	}else{
    		simplices_integer[i].grid_x = round((grid_x-x_min)/x_interval);
    	}

    	if(y_interval==0){
    		simplices_integer[i].grid_y=0;
    	}else{
    		simplices_integer[i].grid_y = round((grid_y-y_min)/y_interval);
    	}

        //initialize the simplex tree and map
        std::pair<Simplex_handle, bool> insert_return = st.insert_simplex(simplices_integer[i].face, Filtration_value(i));//initialize the simplex tree
        if(insert_return.second==false){
            cout << "ERROR:insert fail!!!" << endl;
            return 0;
        }

    }
    //////////////////
    //std::cout << "The complex contains " << st.num_simplices() << " simplices - " << st.num_vertices() << " vertices "<< std::endl;
    

    // Build the map from Simplex handle to grid node coordinates(x,y)
    int i=0;
    for (auto f_simplex : st.filtration_simplex_range()) {
        sg_map[f_simplex] = entry(simplices_integer[i].grid_x,simplices_integer[i].grid_y);
        ++i;
    }


    /*
	//print the simplices_integer
    for(int i=0; i<simplices_integer.size(); ++i){
    	simplex_node simplex = simplices_integer[i];
    	cout << "grid "<<simplex.grid_x << " " << simplex.grid_y ;
    	typeVectorVertex face = simplex.face;
    	cout << " face ";
    	for(int j=0;j<face.size();++j){
    		cout << face[j] << " ";
    	}
    	cout <<""<<endl;
    }
    */


    rank rank_inv_tmp,rank_inv;
    zero_rank(rank_inv,x_range+1,y_range+1);

    clock_t start_total_time = clock();

    /////////////////// Loop for the elbow
    for(int i=0; i<=x_range; ++i){
    	for(int j=0; j<=y_range; ++j){
    		cout << "elbow: " << i << " " << j << endl;
    		entry elbow(i,j);
    		vector<entry> stair = create_stair(elbow, end_point);

    		start_time=clock();
            
            for(Simplex_handle sh: st.complex_simplex_range()){
                entry node_grid = sg_map[sh];
                entry proj_p = proj_grid_to_path(node_grid, elbow, x_range, y_range);
                st.assign_filtration(sh,proj_p.first+proj_p.second);
            }

            //TODO : not sure whether we should use this function
            st.initialize_filtration();//will sort the simplex by filtration order
            
            /*
            //print the filtration
            for (auto f_simplex : st.filtration_simplex_range()) {
                std::cout << "   " << "[" << st.filtration(f_simplex) << "] ";
                for (auto vertex : st.simplex_vertex_range(f_simplex)) {
                    std::cout << static_cast<int>(vertex) << " -- ";
                }
                std::cout << ";" << std::endl;
            
            }*/
            
    		end_time=clock();
            time_update_simplex_tree += ((double)(end_time-start_time)/CLOCKS_PER_SEC);
    		//cout << "time of building a simplex tree "<<(double)(end_time-start_time)/CLOCKS_PER_SEC <<"s"<<endl;

  			
            zero_rank(rank_inv_tmp,x_range+1,y_range+1);

            start_time=clock();
  			int coeff_field_characteristic=11;
  			int min_persistence=0;
  			Persistent_cohomology pcoh(st);
  			pcoh.init_coefficients(coeff_field_characteristic);
  			pcoh.compute_persistent_cohomology(min_persistence);
            auto persistent_pairs = pcoh.intervals_in_dimension(rank_dim);

            end_time=clock();
            time_barcodes+=(double)(end_time-start_time)/CLOCKS_PER_SEC ;
            //cout <<"time of calculating barcodes "<<(double)(end_time-start_time)/CLOCKS_PER_SEC <<"s"<<endl;
 
            
            start_time=clock();
            for (auto pair : persistent_pairs) {                
                int barcode_sid = pair.first;
                int barcode_eid;
                if(pair.second>x_range+y_range+2){
                    barcode_eid = x_range+y_range+1;
                }else{
                    barcode_eid = pair.second;
                }

                for(int r_s=barcode_sid; r_s<barcode_eid; ++r_s){
                    for(int r_e=r_s; r_e<barcode_eid; ++r_e){
                        int ii = stair[r_s].first;
                        int jj = stair[r_s].second;
                        int kk = stair[r_e].first;
                        int ll = stair[r_e].second;
                        rank_inv_tmp[ii][jj][kk][ll] += 1;
                    }
                }   
            }

    		for(int s_x=0;s_x<=x_range;++s_x){
    			for(int s_y=0;s_y<=y_range;++s_y){
    				for(int t_x=s_x; t_x<=x_range; ++t_x){
    					for(int t_y=s_y; t_y<=y_range; ++t_y){
    						if(rank_inv[s_x][s_y][t_x][t_y]==0){
    							rank_inv[s_x][s_y][t_x][t_y]=rank_inv_tmp[s_x][s_y][t_x][t_y];
    						}
    					}
    				}
    			}
    		}

    		end_time=clock();
            time_rank_inv += (double)(end_time-start_time)/CLOCKS_PER_SEC;
    		//cout << "time of calculating rank invariant "<<(double)(end_time-start_time)/CLOCKS_PER_SEC <<"s"<<endl;
 		}
    }
    ///////////////////


    //////////////////// Calculate the signed barcode decomposition
    barcode b;
    start_time=clock();
    compute_R_S_incl_excl(rank_inv,b);
    end_time=clock();
    ////////////////////

    clock_t end_total_time = clock();
    time_multiplicity += (double)(end_time-start_time)/CLOCKS_PER_SEC;
    time_total = (double)(end_total_time-start_total_time)/CLOCKS_PER_SEC;
    
    /*
    //print the rank invariance
    for(int i=0;i<=x_range;++i){
        for(int j=0;j<=y_range;++j){
            for(int k=0;k<=x_range;++k){
                for(int l=0;l<=y_range;++l)
                    if(rank_inv[i][j][k][l]!=0)
                        cout<<"rank[("<<i<<","<<j<<"),("<<k<<","<<l<<")]="<<rank_inv[i][j][k][l]<<endl;

            }
        }
    }*/

    //print barcodes
    //std::cout << "Computed barcode: " << std::endl << b << std::endl;
    cout<<"with st.initialize_filtration()"<<endl;
    std::cout << "Barcode size: " << b.size() << std::endl;

    std::map<bar, int>::iterator it;
    //std::ofstream f("../result/barcode.txt");
    for(it=b.begin(); it != b.end(); it++)
    {
        //f << (it->first).first.first  << " " << (it->first).first.second << " " <<  (it->first).second.first << " " << (it->first).second.second << " " <<   it->second <<  std::endl;
        //cout << (it->first).first.first  << " " << (it->first).first.second << " " <<  (it->first).second.first << " " << (it->first).second.second << " " <<   it->second <<  std::endl;
        cout<<"m[("<<(it->first).first.first<<","<<(it->first).first.second <<"),("<<(it->first).second.first<<","<<(it->first).second.second<<")]="<<it->second<<endl;
    }


    cout <<"time of updating a simplex tree per path is "<< time_update_simplex_tree/((x_range+1)*(y_range+1)) <<"s"<<endl;
    cout <<"time of calculating barcodes per path is "<< time_barcodes/((x_range+1)*(y_range+1)) <<"s"<<endl;
    cout <<"time of calculating rank invariance is "<< time_rank_inv/((x_range+1)*(y_range+1)) <<"s"<<endl;
    cout <<"time of calculating the signed barcode decomposition is "<< time_multiplicity <<"s"<<endl;
    cout << "total time is "<<time_total<<"s"<<endl;
    cout << "average running time per path is "<<time_total/((x_range+1)*(y_range+1)) <<"s"<<endl;


}


