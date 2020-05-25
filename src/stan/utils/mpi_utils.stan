  real[,] swap_index_real(real[,] input) {
    int sizes[2] = dims(input);
    real result[sizes[2],sizes[1]];

    for(i in 1:sizes[1])
      result[:,i] = input[i];

    return result;
  }

  int[,] swap_index_int(int[,] input) {
    int sizes[2] = dims(input);
    int result[sizes[2],sizes[1]];

    for(i in 1:sizes[1])
      result[:,i] = input[i];

    return result;
  }

// converts a ragged matrix into a 2D array with one row per subject
  // with padding
  real[,] ragged2rect_real(int[] M, matrix mdata) {
    int max_row = max(M);
    int num_id = size(M);
    int num_col = cols(mdata);
    real rdata[num_id, max_row * num_col] = rep_array(0.0, num_id, max_row * num_col);
    int ind_start = 1;
    
    for(i in 1:num_id) {
      int ind_end = ind_start + M[i] - 1;
      // column major conversion
      rdata[i, 1:(M[i]*num_col)] = to_array_1d(mdata[ind_start:ind_end]);
      ind_start = ind_end + 1;
    }
    return(rdata);
  }

  int[,] ragged2rect_int(int[] M, int[,] idata) {
    int max_row = max(M);
    int num_id = size(M);
    int num_col = dims(idata)[2];
    int rdata[num_id, max_row * num_col] = rep_array(0, num_id, max_row * num_col);
    int ind_start = 1;
    
    for(i in 1:num_id) {
      int ind_end = ind_start + M[i] - 1;
      // column major conversion
      for(j in 1:num_col) {
        rdata[i, ((j-1)*M[i] + 1):j*M[i]] = idata[ind_start:ind_end,j];
      }
      ind_start = ind_end + 1;
    }
    return(rdata);
  }

// takes as input a id array which must be as long as there are
// columns. Each row is interpreted as a column in the final data.
int[,] long_table2rect_int(int[] id, int[,] row_data) {
  int M[rle_elem_count(id)] = rle_int(id);
  if(size(id) != dims(row_data)[2])
    reject("ID arrays must match in length row_data's 2nd dimension.");
  return ragged2rect_int(M, swap_index_int(row_data));
}
real[,] long_table2rect_real(int[] id, matrix row_data) {
  int M[rle_elem_count(id)] = rle_int(id);
  if(size(id) != cols(row_data))
    reject("ID arrays must match in length row_data's columns.");
  return ragged2rect_real(M, row_data');
}

  real[,] append_real_array_cols(real[,] a, real[,] b) {
    return( to_array_2d(append_col(to_matrix(a), to_matrix(b))) );
  }

real[,] append_real_array_cols_3(real[,] a, real[,] b, real[,] c) {
  return( to_array_2d(append_col(append_col(to_matrix(a), to_matrix(b)), to_matrix(c)) ));
}

real[,] append_real_array_cols_4(real[,] a, real[,] b, real[,] c, real[,] f) {
  return( to_array_2d(append_col(append_col(append_col(to_matrix(a), to_matrix(b)), to_matrix(c)), to_matrix(f) ) ) );
}

  int[,] append_int_array_cols(int[,] a, int[,] b) {

    //return swap_index_int(append_array( swap_index_int(a), swap_index_int(b)));
    int a_row = dims(a)[1];
    int a_col = dims(a)[2];
    int b_row = dims(b)[1];
    int b_col = dims(b)[2];
    int res[a_row,a_col+b_col];

    if(a_row != b_row)
      reject("# of rows do not match. a_row = ", a_row, ", b_row = ", b_row);

    for(i in 1:a_row) {
      res[i, 1:a_col] = a[i];
      res[i, (a_col+1):(a_col+b_col)] = b[i];
    }
        
    return( res );
  }

  int[,] append_int_array_cols_3(int[,] a, int[,] b, int[,] c) {
    return( append_int_array_cols(append_int_array_cols(a, b), c) );
  }
  int[,] append_int_array_cols_4(int[,] a, int[,] b, int[,] c, int[,] d) {
    return( append_int_array_cols(append_int_array_cols(a, b), append_int_array_cols(c, d)) );
  }
  int[,] append_int_array_cols_5(int[,] a, int[,] b, int[,] c, int[,] d, int[,] e) {
    return( append_int_array_cols(append_int_array_cols_4(a, b, c, d), e) );
  }
  int[,] append_int_array_cols_6(int[,] a, int[,] b, int[,] c, int[,] d, int[,] e, int[,] f) {
    return( append_int_array_cols(append_int_array_cols_5(a, b, c, d, e), f) );
  }
  int[,] append_int_array_cols_7(int[,] a, int[,] b, int[,] c, int[,] d, int[,] e, int[,] f, int[,] g) {
    return( append_int_array_cols(append_int_array_cols_6(a, b, c, d, e, f), g) );
  }


int[,] make_2d_int(int[] a) {
    int num_row = size(a);
    int res[num_row,1];
    res[:,1] = a;
    return(res);
  }


matrix append_colwise(vector[] idata);
matrix append_colwise(vector[] idata) {
  int num_cols = size(idata);
  if (num_cols == 1) {
    return to_matrix(idata[1], num_elements(idata[1]), 1, 1);
  } else if (num_cols == 2) {
    return append_col(idata[1], idata[2]);
  }
  return append_col(idata[1], append_colwise(idata[2:num_cols]));
}

  // needed to avoid over-propagation of data to vars due to
  // declaration, fixed in the future
  real[] extract_1d_real(real[] rdata, int offset, int num_row, int c) {
    return(rdata[ (offset + (c-1)*num_row + 1) : (offset + num_row*c)]);
  }

  vector extract_vector(real[] rdata, int offset, int num_row, int c) {
    return(to_vector(extract_1d_real(rdata, offset, num_row, c)));
  }

  int[] extract_1d_int(int[] idata, int offset, int num_row, int c) {
    return(idata[ (offset + (c-1)*num_row + 1) : (offset + num_row*c)]);
  }

// takes as input an array of the form {a,b,c,d} where a,b,c,d are
// real 1d arrays. These are then returned as concatenated 1d array
real[] concatenate_real_array(real[,] elems);
real[] concatenate_real_array(real[,] elems) {
  int num_elems = size(elems);

  if (num_elems == 1)
    return(elems[1]);

  if (num_elems == 2)
    return(to_array_1d(append_row(to_vector(elems[1]), to_vector(elems[2]))));
  return(to_array_1d(append_row(to_vector(elems[1]), to_vector( concatenate_real_array(elems[2:num_elems]) ))));
}


int[] combine_int_array(int[] a, int[] b) {
  int size_a = size(a);
  int size_b = size(b);
  int result[size_a + size_b];
  result[1:size_a] = a;
  result[(size_a+1):(size_a + size_b)] = b;
  return(result);
}
int[] concatenate_int_array(int[,] elems);
int[] concatenate_int_array(int[,] elems) {
  int num_elems = size(elems);
  int num_elem_first = size(elems[1]);

  if (num_elems == 1)
    return(elems[1]);

  if (num_elems == 2) {
    return(combine_int_array(elems[1], elems[2]));
  }

  return(combine_int_array(elems[1], concatenate_int_array(elems[2:num_elems])));
}

// turn the input data columwise into a 2D array, sizes must line up
int[,] to_array_2d_colwise_int(int[] idata, int num_rows, int num_cols) {
  int adata[num_rows,num_cols];
  if(size(idata) != num_rows * num_cols)
    reject("Input data size ", size(idata), " not compatible with ", num_rows, "x", num_cols, " output");
  for(i in 1:num_cols) 
    adata[:,i] = idata[(i-1) * num_rows + 1 : i * num_rows];
  return adata;
}
real[,] to_array_2d_colwise_real(real[] idata, int num_rows, int num_cols) {
  real adata[num_rows,num_cols];
  if(size(idata) != num_rows * num_cols)
    reject("Input data size ", size(idata), " not compatible with ", num_rows, "x", num_cols, " output");
  for(i in 1:num_cols) 
    adata[:,i] = idata[(i-1) * num_rows + 1 : i * num_rows];
  return adata;
}

//  int[] append_array_3_int(int[] a, int[] b, int[] c) {
//    return(append_array(a, append_array(b, c)));
// }

/*
  // forward declare mpi_function
  vector mpi_function(vector mu, vector eta, real[] x_r, int[] x_i);

  vector map_rect_stan(vector theta, vector[] Eta, int[] M, real[,] x_r, int[,] x_i) {
    int num_out = sum(M);
    vector[num_out] out;
    int num_jobs = size(Eta);
    int ind_start = 1;

    for(i in 1:num_jobs) {
      int ind_end = ind_start + M[i] - 1;
      out[ind_start:ind_end] = mpi_function(theta, Eta[i], x_r[i], x_i[i]);
      ind_start = ind_end + 1;
    }

    return(out);
  }

  vector maprect_stan(vector theta, vector[] Eta, int[] M, real[,] x_r, int[,] x_i) {
    int num_out = sum(M);
    vector[num_out] out;
    int num_jobs = size(Eta);
    int ind_start = 1;

    for(i in 1:num_jobs) {
      int ind_end = ind_start + M[i] - 1;
      out[ind_start:ind_end] = mpi_function(theta, Eta[i], x_r[i], x_i[i]);
      ind_start = ind_end + 1;
    }

    return(out);
  }

*/
