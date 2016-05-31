#ifndef __OMP_UTIL_HPP__
#define __OMP_UTIL_HPP__


#define MAT_INC_FACTOR (1.6)


template<typename type>
class matrix{
public:
  matrix(size_t nrow, size_t ncol){
    data = new type[nrow * ncol];
    _ncol = ncol;
    _nrow = nrow;
  }
  ~matrix(){
    delete [] data;
  }
  size_t _nrow;
  size_t _ncol;
  type *data;
  type * operator[](size_t row)
  { return &data[row*_ncol];}
};

template<typename type>
std::vector<type> mat_vec(matrix<type>mat, std::vector<type> g){
	std::vector<type> res(mat._nrow, 0.0);

	for(int i = 0; i < mat._nrow; i++){
		res[i] = 0.0;
		for(int j = 0; j < mat._ncol; j++){
			res[i] += mat[i][j] * g[j];
		}
	}

	return res;
}

template<typename type>
std::vector<type> mat_tran_vec(matrix<type>mat, std::vector<type> g){
	std::vector<type> res(mat._ncol, 0.0);

	for(int i = 0; i < mat._ncol; i++){
		res[i] = 0.0;
		for(int j = 0; j < mat._nrow; j++){
			res[i] += mat[j][i] * g[j];
		}
	}

	return res;
}


//naive mat_mat
template<typename type>
void mat_mat(matrix<type> &mat1, matrix<type> &mat2, matrix<type> &res){
	type sum;
	for(int i = 0; i < mat1._nrow; i++){
		for(int j = 0; j < mat2._ncol; j++){
			sum = 0.0;
			for(int k = 0; k < mat2._nrow; k++){
				sum += mat1[i][k] * mat2[k][j];
			}
			res[i][j] = sum;
		}
	}
}

template<typename type>
void mat_tran_mat(matrix<type> &mat1, matrix<type> &mat2, matrix<type> &res){
	type sum;
	for(int i = 0; i < mat1._ncol; i++){
		for(int j = 0; j < mat2._ncol; j++){
			sum = 0.0;
			for(int k = 0; k < mat2._nrow; k++){
				sum += mat1[k][i] * mat2[k][j];
			}
			res[i][j] = sum;
		}
	}
}

template<typename type>
type l2norm(type *data, size_t len){
	type sum = 0.0;
	for (int i = 0; i < len; ++i)
	{
		sum += data[i] * data[i];
	}
	return sqrt(sum);
}

template<typename type>
void scale(type *data, type scl, size_t len){
	for (int i = 0; i < len; ++i)
	{
		data[i] *= scl;
	}
}

template<typename type>
std::vector<type> vec_plus(
	std::vector<type> a, 
	std::vector<type> b
){
	std::vector<type> res(a.size(), 0.0);

	for (int i = 0; i < a.size(); ++i)
	{
		res[i] = a[i] + b[i];
	}

	return res;
}
#endif