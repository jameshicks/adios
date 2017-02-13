#include "Linalg.hpp"

namespace Linalg
{

double Vectorlike::sum(void)
{
    double outp = 0.0;
    for (size_t i = 0; i < size; ++i) {
        outp += get(i);

    }
    return outp;
}

void Vectorlike::set_all(double v) {
    for (size_t i=0; i<size; ++i) set(i,v);
}

void Vectorlike::apply_inplace(double (*f)(double))
{
    for (size_t i = 0; i < size; ++i) { set(i, f(get(i))); }
}


size_t Vectorlike::argmax(void) const
{
    if (size < 2) { return 0; }
    size_t arg = 0;
    size_t n = size;
    for (size_t i = 1; i < n; ++i) { if (get(i) > get(arg)) { arg = i; } }
    return arg;
}

size_t Vectorlike::argmin(void) const
{
    if (size < 2) { return 0; }
    size_t arg = 0;
    for (size_t i = 1; i < size; ++i) { if (get(i) < get(arg)) { arg = i; } }
    return arg;

}

double Vectorlike::max(void) const { return get(argmax()); }
double Vectorlike::min(void) const { return get(argmin()); }


void Vectorlike::print(void) const
{
    using std::cout; using std::endl;
    cout << "{ ";
    for (size_t i = 0; i < size; ++i) {
        cout << get(i) << ' ';
    }
    cout << "}" << endl;
}


// Vectorview

VectorView::VectorView(double* d, size_t sz, size_t strd) : Vectorlike(sz, strd)

{
    data=d;
}

// Vector
void Vector::create_data_array(size_t sz)
{
    data = (sz == 0) ? nullptr : new double[sz];
}

void Vector::swap(Vector& v)
{
    std::swap(data, v.data);
    std::swap(size, v.size);
}

Vector::Vector(size_t sz) : Vectorlike(sz,1)
{
    create_data_array(sz);
}

Vector::Vector(size_t sz, double val) : Vectorlike(sz, 1) {
    create_data_array(sz);
    set_all(val);
}

Vector::Vector(const std::initializer_list<double> v) : Vectorlike(v.size(), 1)
{
    create_data_array(v.size());
    size_t idx = 0;
    for (auto it = v.begin(); it != v.end(); ++it) {
        data[idx] = *it;
        idx++;
    }
}

Vector::Vector(const std::vector<double> v) : Vectorlike(v.size(), 1)
{
    create_data_array(v.size());
    size_t idx = 0;
    for (auto it : v) {
        data[idx] = it;
        idx++;
    }
}

Vector::Vector(const VectorView& vv) : Vectorlike(vv.size,1)
{
    create_data_array(vv.size);
    for (size_t i = 0; i < size; ++i) { set(i, vv.get(i)); }
}

Vector::Vector(Vector& v) : Vectorlike(v.size, 1)
{
    create_data_array(v.size);
    std::copy(v.data, v.data + v.size, data);
}

Vector::Vector(Vector&& v) : Vectorlike(v.size, 1)
{
    // using std::swap;
    std::swap(data, v.data);
    size = v.size;
}

Vector::~Vector(void) { delete[] data; }

Vector& Vector::operator=(const double d)
{
    for (size_t i = 0; i < size; ++i) {
        data[i] = d;
    }
    return *this;
}

Vector& Vector::operator=(Vector v)
{
    if (this != &v) { swap(v); }
    return *this;
}

// Vector& Vector::operator=(Vector&& v) {
//     delete[] data;
//     data = v.data;
// }


Vector Vector::operator+(const double d)
{
    Vector outp(*this);
    for (size_t i = 0; i < size; ++i) { outp.set(i, get(i) + d); }
    return outp;
}

Vector Vector::operator-(const double d)
{
    Vector outp(*this);
    for (size_t i = 0; i < size; ++i) { outp.set(i, get(i) - d); }
    return outp;
}

Vector Vector::operator*(const double d)
{
    Vector outp(*this);
    for (size_t i = 0; i < size; ++i) { outp.set(i, get(i) * d); }
    return outp;
}

Vector Vector::operator/(const double d)
{
    Vector outp(*this);
    for (size_t i = 0; i < size; ++i) { outp.set(i, get(i) / d); }
    return outp;
}

Vector& Vector::operator/=(const double d) {
    size_t n = size;
    for (size_t i = 0; i < n; ++i) {
        data[i] /= d;
    }
    return *this;
}

Vector Vector::operator+(const Vectorlike& v)
{
    if (size != v.size) {throw std::invalid_argument("Nonconformable vectors"); }
    Vector outp(*this);
    for (size_t i = 0; i < size; ++i) { outp.set(i, get(i) + v.get(i)); }
    return outp;
}

Vector Vector::operator-(const Vectorlike& v)
{
    if (size != v.size) {throw std::invalid_argument("Nonconformable vectors"); }
    Vector outp(*this);
    for (size_t i = 0; i < size; ++i) { outp.set(i, get(i) - v.get(i)); }
    return outp;
}

// Matrix
void Matrix::create_data_array(size_t size)
{
    data = (size == 0) ? nullptr : new double[size];
}
Matrix::Matrix(void) { nrow = 0; ncol = 0; data = nullptr;}
Matrix::Matrix(size_t nr, size_t nc) : nrow(nr), ncol(nc)
{
    create_data_array(nr * nc);
}

Matrix::Matrix(size_t nr, size_t  nc, double v) : Matrix(nr, nc) {
    for (size_t i=0; i<size(); ++i) { data[i] = v; }
}

Matrix::Matrix(std::initializer_list<std::initializer_list<double>> m)
{
    nrow = m.size();
    ncol = m.begin()->size();
    create_data_array(nrow * ncol);

    size_t rowidx = 0;
    for (auto row : m) {
        if (row.size() != ncol) { throw std::invalid_argument("bad matrix literal"); }

        size_t colidx = 0;
        for (auto cell : row) {
            data[translate_index(rowidx, colidx)] = cell;
            colidx++;
        }

        rowidx++;
    }
}
Matrix::Matrix(const Matrix& other)
{
    ncol = other.ncol;
    nrow = other.nrow;

    create_data_array(other.size());
    std::copy(other.data, other.data + other.size(), data);
}

Matrix::~Matrix(void) noexcept { if (data) { delete[] data; } }


Matrix& Matrix::operator=(const double v)
{
    for (size_t i = 0; i < size(); ++i) {
        data[i] = v;
    }
    return *this;
}

Matrix& Matrix::operator=(const Matrix& other)
{
    if (this != &other) {
        double* newdata = new double[other.size()];
        std::copy(other.data, other.data + other.size(), newdata);
        delete[] data;
        data = newdata;
        ncol = other.ncol;
        nrow = other.nrow;
    }
    return *this;
}

Matrix& Matrix::operator=(Matrix&& o) {
    delete[] data;
    data = o.data;
    o.data = nullptr;
    ncol = o.ncol;
    nrow = o.nrow;
    return *this;
}

bool Matrix::operator==(const Matrix& other)
{
    if (!conformable(*this, other)) { return false; }

    for (size_t rawidx = 0; rawidx < size(); ++rawidx) {
        if (data[rawidx] != other.data[rawidx]) { return false; }
    }
    return true;
}


bool Matrix::operator!=(const Matrix& other) { return !(*this == other); }

Matrix Matrix::operator+(const double d) const
{
    Matrix outp(*this);
    for (size_t i = 0; i < size(); ++i) { outp.data[i] += d; }
    return outp;
}

Matrix Matrix::operator*(const double d) const
{
    Matrix outp(*this);
    for (size_t i = 0; i < size(); ++i) { outp.data[i] *= d; }
    return outp;
}

Matrix Matrix::operator-(const double d) const
{
    Matrix outp(*this);
    for (size_t i = 0; i < size(); ++i) { outp.data[i] -= d; }
    return outp;
}

Matrix Matrix::operator/(const double d) const
{
    Matrix outp(*this);
    for (size_t i = 0; i < size(); ++i) { outp.data[i] /= d; }
    return outp;
}

Matrix Matrix::operator+(const Matrix& b) const
{
    Matrix outp(*this);
    for (size_t i = 0; i < size(); ++i) { outp.data[i] += b.data[i]; }
    return outp;
}

Matrix Matrix::operator-(const Matrix& b) const
{
    Matrix outp(*this);
    for (size_t i = 0; i < size(); ++i) { outp.data[i] -= b.data[i]; }
    return outp;
}

Matrix Matrix::apply(double (*f)(double)) const {
    Matrix omat = Matrix(*this);
    for (size_t i=0; i<omat.size(); ++i) { omat.data[i] = f(data[i]); }
    return omat; 
}
VectorView Matrix::row_view(size_t rowidx) const
{
    if (rowidx >= nrow) { throw std::out_of_range("bad index"); }
    return VectorView(data + (rowidx * ncol), ncol, 1);
}
VectorView Matrix::col_view(size_t colidx) const
{
    if (colidx >= ncol) { throw std::out_of_range("bad index"); }
    return VectorView(data + colidx, nrow, ncol);
}



Vector Matrix::get_row(size_t rowidx) const { 
    if (rowidx >= nrow) { throw std::out_of_range("bad index"); }
    return Vector(row_view(rowidx)); 
}
Vector Matrix::get_column(size_t colidx) const { 
if (colidx >= ncol) { throw std::out_of_range("bad index"); }
    return Vector(col_view(colidx)); 
}

void Matrix::set_row(const size_t rowidx, const VectorView& vec)
{
    if (rowidx >= nrow) { throw std::out_of_range("bad index"); }
    if (vec.size != ncol) { throw std::invalid_argument("vector not conformable to matrix"); }
    if (vec.data == data + (ncol * rowidx)) { 
        // Are we setting something to itself? We can just no-op that.
        return; 
    }

    if (vec.stride == 1) {
        // use the stl copy if we can, since its more optimized than the 
        // the for loop we'd have to use otherwise.
        std::copy(vec.data, vec.data + vec.size, data + (ncol * rowidx));
    } else {
        for (size_t idx = 0; idx < ncol; ++idx) {
            set(rowidx, idx, vec.get(idx));

        }
    }
}

void Matrix::set_row(const size_t rowidx, const Vectorlike& vec)
{
    if (rowidx >= nrow) { throw std::out_of_range("bad index"); }
    if (vec.size != ncol) { throw std::invalid_argument("vector not conformable to matrix"); }
    for (size_t idx = 0; idx < ncol; ++idx) {
        set(rowidx, idx, vec.get(idx));

    }
}

void Matrix::set_column(const size_t colidx, const Vectorlike& vec)
{
    if (vec.size != nrow) { throw std::invalid_argument("vector not conformable to matrix"); }
    if (colidx >= ncol) { throw std::out_of_range("bad index"); }
    for (size_t rowidx = 0; rowidx < nrow; ++rowidx) {
        set(rowidx, colidx, vec.get(rowidx));
    }
}

void Matrix::set_submatrix(const size_t i, const size_t j, const Matrix & m)
{
    for (size_t ridx = 0; ridx < m.nrow; ++ridx) {
        auto rv = m.row_view(ridx);
        std::copy(rv.data, rv.data + rv.size, data + translate_index(i + ridx, j));
    }
}

void Matrix::print(void) const
{
    std::cout << '{' << '\n';
    for (size_t rowidx = 0; rowidx < nrow; ++rowidx) {
        std::cout << "\t{ ";
        Vector row = row_view(rowidx);
        for (size_t cellidx = 0; cellidx < row.size; ++cellidx) {
            std::cout << row.data[cellidx] << '\t';

        }
        std::cout << " }\n";

    }
    std::cout << "}" << std::endl;
}

double Matrix::sum(void) const
{
    double o = 0.0;
    for (size_t i = 0; i < size(); ++i) {
        o += data[i];
    }
    return o;
}

Vector Matrix::argmax_col(void) const
{
    Vector outp(ncol);
    for (size_t i = 0; i < ncol; ++i) {
        outp.set(i, col_view(i).argmax());
    }
    return outp;
}
Vector Matrix::argmin_col(void) const
{
    Vector outp(ncol);
    for (size_t i = 0; i < ncol; ++i) {
        outp.set(i, col_view(i).argmin());
    }
    return outp;
}

Vector Matrix::max_col(void) const
{
    Vector outp(ncol);
    for (size_t i = 0; i < ncol; ++i) {
        outp.set(i, col_view(i).max());
    }
    return outp;
}

Vector Matrix::min_col(void) const
{
    Vector outp(ncol);
    for (size_t i = 0; i < ncol; ++i) {
        outp.set(i, col_view(i).min());
    }
    return outp;
}


Vector Matrix::argmax_row(void) const
{
    Vector outp(nrow);
    for (size_t i = 0; i < nrow; ++i) {
        outp.set(i, row_view(i).argmax());
    }
    return outp;
}
Vector Matrix::argmin_row(void) const
{
    Vector outp(nrow);
    for (size_t i = 0; i < nrow; ++i) {
        outp.set(i, row_view(i).argmin());
    }
    return outp;
}

Vector Matrix::max_row(void) const
{
    Vector outp(nrow);
    for (size_t i = 0; i < nrow; ++i) {
        outp.set(i, row_view(i).max());
    }
    return outp;
}

Vector Matrix::min_row(void) const
{
    Vector outp(nrow);
    for (size_t i = 0; i < nrow; ++i) {
        outp.set(i, row_view(i).min());
    }
    return outp;
}


Matrix diag(const Vectorlike& v) {
    size_t n = v.size;
    Matrix d(n, n, 0.0);
    for (size_t i = 0; i < n; ++i) { d.set( i, i, v.get(i)); }
    return d;
}

// Products
Matrix matrix_product(const Matrix& a, const Matrix& b)
{
    //Naive implementation of matrix multiplication

    if (a.ncol != b.nrow) {
        throw std::invalid_argument("Nonconformable arrays");
    }

    Matrix outp = Matrix(a.nrow, b.ncol);
    for (size_t out_i = 0; out_i < a.nrow; ++out_i) {
        for (size_t out_j = 0; out_j < b.ncol; ++out_j) {
            outp.set(out_i, out_j, dot_product(a.row_view(out_i), b.col_view(out_j)));
        }
    }
    return outp;
}

Matrix kronecker_product(const Matrix& a, const Matrix& b)
{
    Matrix prod = Matrix(a.nrow * b.nrow, a.ncol * b.ncol);
    for (size_t i = 0; i < a.nrow; ++i) {
        for (size_t j = 0; j < a.ncol; ++j) {
            Matrix block = b * a.get(i, j);
            prod.set_submatrix(i * a.nrow, j * a.ncol, block);
        }
    }
    return prod;
}

Matrix direct_product(const Matrix& a, const Matrix& b)
{
    if (!Matrix::conformable(a, b)) { throw std::invalid_argument("Nonconformable matrices"); }
    Matrix outp(a);
    for (size_t i = 0; i < a.size(); ++i) { outp.data[i] *= b.data[i]; }
    return outp;
}



void vector_matrix_product(const Vectorlike& v, const Matrix& b, Vectorlike* into) {
    for (size_t idx = 0; idx < v.size; ++idx) {
        into->set(idx, dot_product(v, b.col_view(idx)));
    }
}

void matrix_vector_product(const Matrix& a, const Vectorlike& v, Vectorlike* into) {
    // if (!(v.size == a.nrow == into->size)) { throw std::invalid_argument("Vector and Matrix do not conform"); }

    for (size_t idx = 0; idx < v.size; ++idx) {
        into->set(idx, dot_product(a.row_view(idx), v));
    }

}

Vector vector_matrix_product(const Vectorlike& v, const Matrix& b) {
    if (v.size != b.ncol) { throw std::invalid_argument("Vector and Matrix do not conform"); }
    Vector outp(v.size);

    vector_matrix_product(v, b, &outp);
    return outp;
}

Vector matrix_vector_product(const Matrix& a, const Vectorlike& v) {
    if (v.size != a.nrow) { throw std::invalid_argument("Vector and Matrix do not conform"); }
    Vector outp(v.size);

    matrix_vector_product(a, v, &outp);

    return outp;
}



Matrix dmatrix_matrix_product(const Vectorlike& v, const Matrix& a) {
    // Post-multiplying by a diagonal matrix multiplies row i by
    // the ith diagonal element
    if (!a.is_square() && a.nrow == v.size) { throw std::invalid_argument("Nonconformable operands"); }
    Matrix z(a);
    for (size_t i=0; i < v.size; ++i) {
        auto row = z.row_view(i);
        row.inplace_mul(v.get(i));
    }
    return z;
}


void matrix_dmatrix_product(const Matrix& a, const Vectorlike& v, Matrix* into) {
    for (size_t col = 0; col < a.ncol; col++) {
        double val = v.get(col);
        for (size_t i = 0; i < a.nrow; i++) {
            (into->data)[(i*a.ncol) + col] = a.data[(i*a.ncol)+col] * val; 
        }
    }    
   
}

Matrix matrix_dmatrix_product(const Matrix& a, const Vectorlike& v) {
    if (!a.is_square() && a.nrow == v.size) { throw std::invalid_argument("Nonconformable operands"); }
    Matrix z(a);
    matrix_dmatrix_product(a, v, &z);
    return z; 
}

void dmatrix_vector_product(const Vectorlike& d, const Vectorlike& v, Vectorlike* into) {
    double* dest = into->data;
    for (size_t i = 0; i < d.size; i++) {
        dest[i*into->stride] = d.data[i*d.stride] * v.data[i*v.stride];
    }
}

Vector dmatrix_vector_product(const Vectorlike& d, const Vectorlike& v) {
    if (d.size != v.size) { throw std::invalid_argument("Nonconformable operands"); } 
    Vector outp(d.size);

    dmatrix_vector_product(d,v, &outp);

    return outp;
}

}

