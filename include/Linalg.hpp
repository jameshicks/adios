#ifndef LINALG_HPP
#define LINALG_HPP

#include <iostream>
#include <algorithm>
#include <utility>
#include <vector>
#include <stdexcept>
#include <initializer_list>

namespace Linalg {

struct Vectorlike {
  size_t size;
  explicit Vectorlike(size_t sz) : size(sz) {}
  virtual double get(size_t idx) const =0;
  virtual void set(size_t idx, double v) =0;

  void set_all(double v);
  
  double sum(void);
  void apply_inplace(double (*f)(double));
  
  inline bool operator==(const Vectorlike& other) const {
    if (size != other.size) return false;
    for (size_t i=0; i< size; ++i) { if (get(i) != other.get(i)) return false; }
    return true;
  }
  inline bool operator!=(const Vectorlike& other) const { return !(*this == other); }
  inline void inplace_add(const double v) { for (size_t i=0; i < size; ++i) set(i, get(i) + v); }
  inline void inplace_sub(const double v) { for (size_t i=0; i < size; ++i) set(i, get(i) - v); }
  inline void inplace_mul(const double v) { for (size_t i=0; i < size; ++i) set(i, get(i) * v); }
  inline void inplace_div(const double v) { for (size_t i=0; i < size; ++i) set(i, get(i) / v); }
  inline void boundscheck(size_t i) const { if (i >= size) throw std::invalid_argument("index oor"); }
  // inline void boundscheck(size_t i) const {return;}
  size_t argmax(void) const;
  size_t argmin(void) const;
  double max(void) const;
  double min(void) const;
  void print(void) const;
  // friend Vector column_vec_matrix_prod(const Vectorlike& v, const Matrix& b);
  
};

struct VectorView : public Vectorlike {
  double* const data;
  size_t stride;
  VectorView(double* const d, size_t sz, size_t strd);



  inline double get(size_t idx) const { return data[idx*stride]; }
  inline void set(size_t idx, double v) { data[idx*stride] = v; }
  inline double operator[](size_t idx) { return data[idx*stride]; }
  // inline const double operator[](size_t idx) const { return data[idx*stride]; }
  
  inline double* begin(void) { return data; }
  inline double* end(void) { return data+size; }

};

struct Vector : public Vectorlike {
  double* data;
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
  Vector(const VectorView& vv); // move constructor
  Vector(Vector& v); // copy constructor
  Vector(Vector&& v);
  ~Vector(void);
  inline double* begin(void) { return data; }
  inline double* end(void) { return data+size; }

  
  Vector& operator=(const double d);
  Vector& operator=(const Vector v); // copy assignment
  // Vector& operator=(Vector&& v); // move assignment
  
  inline double get(size_t idx) const { 
    boundscheck(idx);
    return data[idx]; 
  }
  inline void set(size_t idx, const double v) {
  boundscheck(idx);
   data[idx] = v; 
 }
  Vector operator+(const double d);
  Vector operator-(const double d);
  Vector operator*(const double d);
  Vector operator/(const double d);
  Vector operator+(const Vectorlike& v);
  Vector operator-(const Vectorlike& v);
};


struct coord { size_t i; size_t j; };

double dot_product(const Vectorlike& u, const Vectorlike& v);

class Matrix {
  private:
  double* data;
  void create_data_array(size_t size);
  inline size_t translate_index(size_t i, size_t j) const { return i * ncol + j; }
  
  public:
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
  inline double get (size_t i, size_t j) const { return data[translate_index(i,j)]; }
  inline void set (size_t i, size_t j, double v) { data[translate_index(i,j)] = v; }
  double sum(void) const;

  // operators
  Matrix& operator=(const double v);
  Matrix& operator=(const Matrix& v);
  Matrix& operator=(Matrix&& v);
  
  inline double operator[](const coord& c) { return data[translate_index(c.i, c.j)]; }
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

Matrix matrix_product(const Matrix& a, const Matrix& b);
Matrix kronecker_product(const Matrix& a, const Matrix& b);
Matrix direct_product(const Matrix& a, const Matrix& b);
Vector vector_matrix_product(const Vectorlike& v, const Matrix& b);
Vector matrix_vector_product(const Matrix& a, const Vectorlike& v);

// A matrix post-multiplied by a diagonal matrix
Matrix dmatrix_matrix_product(const Vectorlike& v, const Matrix& a);

// A matrix post-multiplied by a diagonal matrix
Matrix matrix_dmatrix_product(const Matrix& a, const Vectorlike& v);

}

#endif