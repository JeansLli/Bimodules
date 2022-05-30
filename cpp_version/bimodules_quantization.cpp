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
using Filtration_value = Simplex_tree::Filtration_value;
using typeVectorVertex = std::vector<Vertex_handle>;
using typePairSimplexBool = std::pair<Simplex_tree::Simplex_handle, bool>;
using Field_Zp = Gudhi::persistent_cohomology::Field_Zp;
using Persistent_cohomology = Gudhi::persistent_cohomology::Persistent_cohomology<Simplex_tree, Field_Zp >;

typedef std::vector<std::vector<std::vector<std::vector<int> > > > rank; 


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


vector<string> split(string &str, const string &pattern)
{
    char * strc = new char[strlen(str.c_str())+1];
    strcpy(strc, str.c_str());   //string转换成C-string
    vector<string> res;
    char* temp = strtok(strc, pattern.c_str());
    while(temp != NULL)
    {
        res.push_back(string(temp));
        temp = strtok(NULL, pattern.c_str());
    }
    delete[] strc;
    return res;
}

vector<int> create_stair(vector<int> elbow, vector<int> end_point){
	vector<int> stair;
	int i,j,k,l;
	
	if(elbow[0]==0 || elbow[1]==end_point[1]){

		for(i=0; i<=end_point[1]; ++i){
			stair.push_back(0);
			stair.push_back(i);
		}
		--i;
		for(j=1; j<=end_point[0]; ++j){
			stair.push_back(j);
			stair.push_back(i);
		}
	}
    else{
    	for(i=0; i<=elbow[1]; ++i){
			stair.push_back(0);
			stair.push_back(i);
		}
		--i;
		for(j=0; j<elbow[0]; ++j){
			stair.push_back(j+1);
			stair.push_back(i);
		}
		--j;
		for(k=i+1; k<=end_point[1]; ++k){
			stair.push_back(j+1);
			stair.push_back(k);
		}
		--k;
		for(l=j+2; l<=end_point[0]; ++l){
			stair.push_back(l);
			stair.push_back(k);
		}
    }
    return stair;
}


vector<int> proj_grid_to_path(simplex_node simplex, vector<int> elbow, int x_range,int y_range){
	int x = simplex.grid_x;
	int y = simplex.grid_y;
	int x0 = elbow[0];
    int y0 = elbow[1];
    vector<int> proj_p;

    if((x==0 && y<=y0) || (y==y0 && x<=x0) || (y>=y0 && x==x0) || (y==y_range && x>=x0)){
        proj_p.push_back(x);
    	proj_p.push_back(y);
    }
    else if(x<=x0 and y<=y0){
    	proj_p.push_back(x);
    	proj_p.push_back(y0);
    }
    else if(x<=x0 and y>=y0){
    	proj_p.push_back(x0);
    	proj_p.push_back(y);
    }
    else{
    	proj_p.push_back(x);
    	proj_p.push_back(y_range);
    }

	return proj_p;
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
	int x_range=20;
	int y_range=20;
	int line_i = 0;
    int vertice_id = 0;
    int rank_dim = 0;
    vector<int> rank_num;
    vector<float> x_values;
    vector<float> y_values;
    vector<simplex_node> simplices;
    vector<string> block_lines;

    vector<int> end_point;
    end_point.push_back(x_range);
    end_point.push_back(y_range);



    //Read file
	//string filename("../datasets/fig3-2.txt");
	string filename("../../cgta_paper_2021/function_rips/function_rips_with_threshold_1000000_1.scc");
	//vector<string> lines;
    string line;

    ifstream input_file(filename);
    if (!input_file.is_open()) {
        cerr << "Could not open the file - '"
             << filename << "'" << endl;
        return EXIT_FAILURE;
    }
    
    while (getline(input_file, line)){
        if(line[0]!='#'){
        	if(line_i<=1){
        		line_i+=1;
        		continue;
        	}
        	else if(line_i==2){
        		vector<string> rank_string = split(line, " ");
        		rank_string.pop_back();
        		//cout<<"rank_num"<<endl;
        		for(int i=0;i < rank_string.size();++i){
        			rank_num.push_back(std::stoi(rank_string[i]));
        		}
    			line_i += 1;
    		}else{
    			block_lines.push_back(line);
        		line_i += 1;
        	}
        }
    }
    input_file.close();
    
    int start_id = 0;

    for(int i=rank_num.size()-1; i>=0; --i){
    	//string simplices_i;
    	//cout << "rank_num:" << i << endl;
    	int simplices_len=block_lines.size();
    	for(int j=simplices_len-start_id-rank_num[i]; j<simplices_len-start_id; ++j){
    		vector<string> block = split(block_lines[j],";");
    		vector<string> grid = split(block[0]," ");
    		float grid_x = std::stof(grid[0]);
    		float grid_y = std::stof(grid[1]);
    		x_values.push_back(grid_x);
    		y_values.push_back(grid_y);
    		simplex_node simplex;
    		simplex.grid_x = grid_x;
    		simplex.grid_y = grid_y;

    		if(block.size()>1){
    			vector<string> face = split(block[1]," ");
    			for(int k=0; k<face.size();++k){
    				simplex.face.push_back(std::stof(face[k]));
    			}
    		}else{
    			simplex.face={vertice_id};
    			vertice_id += 1;
    		}
    		simplices.push_back(simplex);
    	}
    	start_id += rank_num[i];
    }


    start_id = rank_num[rank_num.size()-2] + rank_num[rank_num.size()-1];
    int start_id_last;
    for(int i=rank_num.size()-3; i>=0; --i){
    	start_id_last = start_id - rank_num[i+1];
    	for(int j=start_id; j<(start_id+rank_num[i]); ++j){
    		//cout<<"j="<<j<<endl;
    		simplex_node simplex = simplices[j];
    		typeVectorVertex face = simplex.face;
    		typeVectorVertex new_face;
    		//cout << "start_id_last" << start_id_last <<endl;
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
/*
  	//print the simplices  
    for(int i=0; i<simplices.size(); ++i){
    	simplex_node simplex = simplices[i];
    	cout << "grid "<<simplex.grid_x << " " << simplex.grid_y ;
    	typeVectorVertex face = simplex.face;
    	cout << " face ";
    	for(int j=0;j<face.size();++j){
    		cout << face[j] << " ";
    	}
    	cout <<""<<endl;
    }
*/


    float x_min = *min_element(x_values.begin(),x_values.end()); 
    float y_min = *min_element(y_values.begin(),y_values.end()); 
    float x_max = *max_element(x_values.begin(),x_values.end()); 
    float y_max = *max_element(y_values.begin(),y_values.end()); 

    float x_interval = (x_max-x_min)/x_range;
    float y_interval = (y_max-y_min)/y_range;

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

    Simplex_tree st;
    rank rank_inv_tmp,rank_inv;
    zero_rank(rank_inv,x_range+1,y_range+1);
    

    for(int i=0; i<=x_range; ++i){
    	for(int j=0; j<=y_range; ++j){
    		cout << "elbow: " << i << " " << j << endl;
    		vector<int> elbow;
    		elbow.push_back(i);
    		elbow.push_back(j);

    		vector<int> stair = create_stair(elbow, end_point);
    		//for(int stair_id=0;stair_id<stair.size();++stair_id)
    		//	cout << stair[stair_id] << " ";
    		//cout<<""<<endl;

    		clock_t start_time=clock();
    		for(int k=0; k<simplices_integer.size();++k){
    			simplex_node simplex = simplices_integer[k];
    			vector<int> proj_p = proj_grid_to_path(simplex, elbow, x_range, y_range);
    			float filtration_value=proj_p[0]+proj_p[1];
    			if(st.find(simplex.face)==st.null_simplex()){
    				st.insert_simplex(simplex.face, Filtration_value(filtration_value));
    			}
    		}
    		clock_t end_time=clock();
    		cout << "time of building a simplex tree "<<(double)(end_time-start_time)/CLOCKS_PER_SEC <<"s"<<endl;

  			start_time=clock();
  			int coeff_field_characteristic=2;
  			int min_persistence=0;
  			Persistent_cohomology pcoh(st);

  			pcoh.init_coefficients(coeff_field_characteristic);
 
  			pcoh.compute_persistent_cohomology(min_persistence);
  			end_time=clock();

  			cout <<"time of calculating barcodes "<<(double)(end_time-start_time)/CLOCKS_PER_SEC <<"s"<<endl;
 
  			// Output the diagram in filediag
  			//pcoh.output_diagram();
  			string output_file = "../datasets/barcodes_tmp";
  			std::ofstream out(output_file);
    		pcoh.output_diagram(out);
    		out.close();

    		std::ifstream in_file(output_file);
    		string line;

    		zero_rank(rank_inv_tmp,x_range+1,y_range+1);

    		start_time=clock();
    		while (getline(in_file, line)){
    			//cout << "barcode " << line <<endl;
    			vector<string> barcode = split(line," ");
    			if(std::stoi(barcode[1])==rank_dim){
    				int barcode_sid = std::stoi(barcode[2]);
    				int barcode_eid;
    				if(barcode[3]=="inf"){
    					barcode_eid = x_range+y_range+1;
    				}else{
    					barcode_eid = std::stoi(barcode[3]);;
    				}

    				for(int r_s=barcode_sid; r_s<barcode_eid; ++r_s){
    					for(int r_e=r_s; r_e<barcode_eid; ++r_e){
    						//cout <<"start "<<r_s<<endl;
    						//cout <<"end " <<r_e<<endl;
    						//cout << "rank_coor" << stair[r_s*2] << " "<<stair[r_s*2+1] << " "<<stair[r_e*2] << " "<<stair[r_e*2+1] << endl;
    						rank_inv_tmp[stair[r_s*2]][stair[r_s*2+1]][stair[r_e*2]][stair[r_e*2+1]] +=1;	
    					}
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
    		in_file.close();
    		cout << "time of calculating rank invariance "<<(double)(end_time-start_time)/CLOCKS_PER_SEC <<"s"<<endl;
 		}
    }
}


