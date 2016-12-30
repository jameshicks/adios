#ifndef LINALG_HPP
#define LINALG_HPP

#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <stdexcept>
#include <initializer_list>
#include <cstring> // for std::memcpy

namespace Linalg {

struct Vectorlike {
  size_t size;
  size_t stride;
  double* data;
  explicit Vectorlike(size_t sz, size_t strd) : size(sz), stride(strd)  {}

  inline double get(size_t idx) const { return data[idx*stride]; }
  inline void set(size_t idx, double v)  { data[idx*stride] = v; }
  inline double operator[](size_t idx)  { return data[idx*stride]; }

  void set_all(double v);
  
  double sum(void);
  void apply_inplace(double (*f)(double));
  
  inline bool operator==(const Vectorlike& other) const {
    if (size != other.size) return false;
    for (size_t i=0; i< size; ++i) { if (get(i) != other.get(i)) return false; }
    return true;
  }
  inline bool operator!=(const Vectorlike& other) const { return !(*this == other); }
  inline void inplace_add(const double v) {
     size_t n = size;
     for (size_t i=0; i < n; ++i) set(i, get(i) + v); 
   }
  inline void inplace_sub(const double v) {
     size_t n = size;
     for (size_t i=0; i < n; ++i) set(i, get(i) - v); 
   }
  inline void inplace_mul(const double v) {
     size_t n = size;
     for (size_t i=0; i < n; ++i) set(i, get(i) * v); 
   }
  inline void inplace_div(const double v) {
     size_t n = size;
     const double recip = 1.0 / v;
     for (size_t i=0; i < n; ++i) set(i, get(i) * recip); 
   }
  // inline void boundscheck(size_t i) const { if (i >= size) throw std::invalid_argument("index oor"); }
  inline bool boundscheck(size_t i) const { return i < size; }
  size_t argmax(void) const;
  size_t argmin(void) const;
  double max(void) const;
  double min(void) const;
  void print(void) const;
  
};

struct VectorView : public Vectorlike {
  VectorView(double* d, size_t sz, size_t strd);

  inline double* begin(void) { return data; }
  inline double* end(void) { return data+size; }

};

struct Vector : public Vectorlike {
  void create_data_array(size_t sz);

  Vector apply(double (*f)(double)) {
    Vector outp(*this);
    outp.apply_inplace(f);
    return outp;
  }
  
  void swap(Vector& other);
  Vector(size_t sz);
  Vector(const std::initializer_list<double> v);
  Vector(const std::vector<double> v);
  Vector(size_t sz, double val); 
  Vector(const VectorView& vv); // move constructor
  Vector(Vector& v); // copy constructor
  Vector(Vector&& v);
  ~Vector(void);
  inline double* begin(void) { return data; }
  inline double* end(void) { return data+size; }

  
  Vector& operator=(const double d);
  Vector& operator=(const Vector v); // copy assignment
  // Vector& operator=(Vector&& v); // move assignment

  Vector operator+(const double d);
  Vector operator-(const double d);
  Vector operator*(const double d);
  Vector operator/(const double d);
  Vector operator+(const Vectorlike& v);
  Vector operator-(const Vectorlike& v);
  Vector& operator/=(const double d);
};



class Matrix {
  public:
  double* data;
  void create_data_array(size_t size);
  inline size_t translate_index(size_t i, size_t j) const { return i * ncol + j; }
  
  size_t nrow;
  size_t ncol;
  
  inline size_t size(void) const { return ncol*nrow; }
  
  // Our usual constructor
  Matrix(void);
  Matrix(size_t nr, size_t nc);
  Matrix(size_t nr, size_t  nc, double v);
  Matrix(std::initializer_list<std::initializer_list<double>> m);
  
  // Copy Constructor
  Matrix(const Matrix& other);
  
  // Destructor
  ~Matrix(void) noexcept;

  friend void swap(Matrix& a, Matrix& b) {
    using std::swap;
    swap(a.data, b.data);
    swap(a.ncol, b.ncol);
    swap(a.nrow, b.nrow);
}

  // getters and setters 
  inline double get (size_t i, size_t j) const { return data[i * ncol + j]; }
  inline void set (size_t i, size_t j, double v) { data[i * ncol + j] = v; }
  double sum(void) const;

  // operators
  Matrix& operator=(const double v);
  Matrix& operator=(const Matrix& v);
  Matrix& operator=(Matrix&& v);
  
  bool operator==(const Matrix& other);
  bool operator!=(const Matrix& other);
  
  Matrix operator+(const double d) const;
  Matrix operator-(const double d) const;
  Matrix operator*(const double d) const;
  Matrix operator/(const double d) const;

  Matrix operator+(const Matrix& b) const;  
  Matrix operator-(const Matrix& b) const;

  Matrix apply(double (*f)(double)) const;

  // Misc
  inline bool is_square(void) const { return nrow == ncol; }
  static inline bool conformable (const Matrix& a, const Matrix& b) { return a.ncol == b.ncol && a.nrow == b.nrow; }
  
  // Column getters and setters
  VectorView row_view(size_t rowidx) const;
  VectorView col_view(size_t colidx) const;

  Vector get_row(size_t rowidx) const;
  Vector get_column(size_t colidx) const;
  
  void set_row(const size_t rowidx, const VectorView& vec);
  void set_row(const size_t rowidx, const Vectorlike& vec);
  void set_column(const size_t colidx, const Vectorlike& vec);
  void set_submatrix(const size_t i, const size_t j, const Matrix& m);

  Vector argmax_col(void) const;
  Vector argmin_col(void) const;
  Vector max_col(void) const;
  Vector min_col(void) const;


  Vector argmax_row(void) const;
  Vector argmin_row(void) const;
  Vector max_row(void) const;
  Vector min_row(void) const;


  void print(void) const;
  friend Matrix direct_product(const Matrix& a, const Matrix& b);

};



Matrix diag(const Vectorlike& v);


  inline double dot_product(const Vectorlike & u, const Vectorlike & v) {
    if (u.size != v.size) { throw std::invalid_argument("non-conformable vectors"); }
    double outp = 0.0;
    double val = 0.0;
    size_t n = u.size;
    for (size_t i = 0; i < n; ++i) {
      val = v.get(i) * u.get(i);
      outp += val;
    }
    return outp;
  }


Matrix matrix_product(const Matrix& a, const Matrix& b);
Matrix kronecker_product(const Matrix& a, const Matrix& b);
Matrix direct_product(const Matrix& a, const Matrix& b);


// For preallocated output vectors
void vector_matrix_product(const Vectorlike& v, const Matrix& b, Vectorlike* into);
void matrix_vector_product(const Matrix& a, const Vectorlike& v, Vectorlike* into);


// Allocates a new output vector
Vector vector_matrix_product(const Vectorlike& v, const Matrix& b);
Vector matrix_vector_product(const Matrix& a, const Vectorlike& v);

// A matrix post-multiplied by a diagonal matrix
Matrix dmatrix_matrix_product(const Vectorlike& v, const Matrix& a);

// A matrix post-multiplied by a diagonal matrix
void matrix_dmatrix_product(const Matrix& a, const Vectorlike& v, Matrix* into); // For preallocated outputs
Matrix matrix_dmatrix_product(const Matrix& a, const Vectorlike& v);

// Diagonal matrix pre or post multiplied by a vector
void dmatrix_vector_product(const Vectorlike& d, const Vectorlike& v, Vectorlike* into);
Vector dmatrix_vector_product(const Vectorlike& d, const Vectorlike& v);

}

#endif
