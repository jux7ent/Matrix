#include <iostream>
#include <cstring>
#include <stdexcept>
#include <cmath>

using std::cout;
using std::cin;
using std::endl;

class MatrixWrongSizeError : std::logic_error {
public:
    MatrixWrongSizeError();
};

MatrixWrongSizeError::MatrixWrongSizeError() : std::logic_error("Error") {}

class MatrixIndexError : std::logic_error {
public:
    MatrixIndexError();
};

MatrixIndexError::MatrixIndexError() : std::logic_error("Error") {}

class MatrixIsDegenerateError : std::logic_error {
public:
    MatrixIsDegenerateError();
};

MatrixIsDegenerateError::MatrixIsDegenerateError() : logic_error("Error") {}

template <typename T>
class Matrix {

    template <typename C>
    friend std::istream &operator >> (std::istream &, Matrix<C> &);

    template <typename C>
    friend std::ostream &operator << (std::ostream &, const Matrix<C> &);

    template <typename C>
    friend Matrix<C> operator*(const C, const Matrix<C> &);

protected:

    int row_;

    int column_;

    T **value_;

    void set_size(const int rows, const int columns);

public:
    Matrix();

    Matrix(const int rows, const int columns);

    Matrix(const Matrix<T> &);

    ~Matrix();

    int getRowsNumber() const { return row_; }

    int getColumnsNumber() const { return column_; }

    T operator()(const int row, const int column) const;

    virtual Matrix<T> &operator=(const Matrix<T> &);

    Matrix<T> operator+(const Matrix<T> &) const;

    Matrix<T> operator-(const Matrix<T> &) const;

    virtual Matrix<T> operator*(const Matrix<T> &) const;

    Matrix<T> operator*(const T &) const;

    Matrix<T> operator/(const T &) const;

    Matrix<T> &operator+=(const Matrix<T> &);

    Matrix<T> &operator-=(const Matrix<T> &);

    virtual Matrix<T> &operator*=(const Matrix<T> &);

    virtual Matrix<T> &operator*=(const T);

    virtual Matrix<T> &operator/=(const T);

    Matrix<T> getTransposed() const;

    virtual Matrix<T> &transpose();
};

template <typename T>
void Matrix<T>::set_size(const int rows, const int columns) {
    this->row_ = rows;
    this->column_ = columns;
    this->value_ = new T*[rows];
    for (int i = 0; i < rows; ++i) {
        this->value_[i] = new T[columns];
        for (int j = 0; j < columns; ++j)
            this->value_[i][j] = 0;
    }
}

template <typename T>
Matrix<T>::Matrix() {
    set_size(1, 1);
}

template <typename T>
Matrix<T>::Matrix(const int rows, const int columns) {
    set_size(rows, columns);
}

template <typename T>
Matrix<T>::Matrix(const Matrix<T> &that) {
    this->row_ = that.row_;
    this->column_ = that.column_;
    this->value_ = new T*[this->row_];
    for (int i = 0; i < this->row_; ++i) {
        this->value_[i] = new T[this->column_];
        for (int j = 0; j < this->column_; ++j)
            this->value_[i][j] = that.value_[i][j];
    }
}

template <typename T>
Matrix<T>::~Matrix() {
    for (int i = 0; i < row_; ++i)
        delete[] value_[i];
    delete[] value_;
}

template <typename T>
T Matrix<T>::operator()(const int row, const int column) const{
    if (row >= row_ || column >= column_)
        throw MatrixIndexError();
    return value_[row][column];
}

template <typename C>
std::istream &operator >> (std::istream &in, Matrix<C> &m) {
    for (int i = 0; i < m.row_; ++i)
        for (int j = 0; j < m.column_; ++j)
            in >> m.value_[i][j];
    return in;
}

template <typename C>
std::ostream &operator << (std::ostream &out, const Matrix<C> &m) {
    for (int i = 0; i < m.row_; ++i) {
        for (int j = 0; j < m.column_; ++j)
            out << m.value_[i][j] << ' ';
        out << '\n';
    }
    return out;
}

template <typename T>
Matrix<T> &Matrix<T>::operator=(const Matrix<T> &that) {
    if (this != &that) {
        this->row_ = that.row_;
        this->column_ = that.column_;
        this->value_ = new T*[this->row_];
        for (int i = 0; i < this->row_; ++i) {
            this->value_[i] = new T[this->column_];
            for (int j = 0; j < this->column_; ++j)
                this->value_[i][j] = that.value_[i][j];
        }
    }
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::operator+(const Matrix<T> &that) const {
    if (this->row_ != that.row_ || this->column_ != that.column_)
        throw MatrixWrongSizeError();
    Matrix<T> result(this->row_, this->column_);
    for (int i = 0; i < this->row_; ++i)
        for (int j = 0; j < this->column_; ++j)
            result.value_[i][j] = this->value_[i][j] + that.value_[i][j];
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator-(const Matrix<T> &that) const {
    if (this->row_ != that.row_ || this->column_ != that.column_)
        throw MatrixWrongSizeError();
    Matrix<T> result(this->row_, this->column_);
    for (int i = 0; i < this->row_; ++i)
        for (int j = 0; j < this->column_; ++j)
            result.value_[i][j] = this->value_[i][j] - that.value_[i][j];
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const Matrix<T> &that) const {
    if (this->column_ != that.row_)
        throw MatrixWrongSizeError();
    Matrix<T> result(this->row_, that.column_);
    for (int i = 0; i < this->row_; ++i)
        for (int j = 0; j < that.column_; ++j)
            for (int k = 0; k < this->column_; ++k)
                result.value_[i][j] += this->value_[i][k] * that.value_[k][j];
    return result;
}

template <typename T>
Matrix<T> Matrix<T>::operator*(const T &t) const {
    Matrix<T> result = *this;
    for (int i = 0; i < row_; ++i)
        for (int j = 0; j < column_; ++j)
            result.value_[i][j] *= t;
    return result;
}

template <typename C>
Matrix<C> operator*(const C t, const Matrix<C> &m) {
    return m * t;
}

template <typename T>
Matrix<T> Matrix<T>::operator/(const T &t) const {
    Matrix<T> result = *this;
    for (int i = 0; i < row_; ++i)
        for (int j = 0; j < column_; ++j)
            result.value_[i][j] /= t;
    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::operator+=(const Matrix<T> &that) {
    *this = *this + that;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator-=(const Matrix<T> &that) {
    *this = *this - that;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator*=(const Matrix<T> &that) {
    *this = *this * that;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator*=(const T t) {
    *this = *this * t;
    return *this;
}

template <typename T>
Matrix<T> &Matrix<T>::operator/=(const T t) {
    *this = *this / t;
    return *this;
}

template <typename T>
Matrix<T> Matrix<T>::getTransposed() const{
    Matrix<T> result(this->column_, this->row_);
    for (int i = 0; i < this->column_; ++i)
        for (int j = 0; j < this->row_; ++j)
            result.value_[i][j] = this->value_[j][i];
    return result;
}

template <typename T>
Matrix<T> &Matrix<T>::transpose() {
    *this = this->getTransposed();
    return *this;
}

template <typename T>
class SquareMatrix : public Matrix<T> {

    template<typename C>
    friend SquareMatrix<C> operator*(const C, const SquareMatrix<C> &);

private:

    SquareMatrix<T> my_minor(const int without_row, const int without_col) const;

public:

    SquareMatrix(const int);

    SquareMatrix(const Matrix<T> &);

    SquareMatrix(const SquareMatrix<T> &);

    int getSize() const;

    SquareMatrix<T> operator+(const SquareMatrix<T> &) const;

    SquareMatrix<T> operator-(const SquareMatrix<T> &) const;

    SquareMatrix<T> operator*(const SquareMatrix<T> &) const;

    Matrix<T> operator*(const Matrix<T> &) const;

    SquareMatrix<T> &operator*=(const T);

    SquareMatrix<T> &operator/=(const T);

    SquareMatrix<T> operator*(const T) const;

    SquareMatrix<T> operator/(const T) const;

    SquareMatrix<T> operator=(const SquareMatrix<T> &that);

    SquareMatrix<T> &operator+=(const SquareMatrix<T> &);

    SquareMatrix<T> &operator-=(const SquareMatrix<T> &);

    SquareMatrix<T> &operator*=(const SquareMatrix<T> &);

    Matrix<T> &operator*=(const Matrix<T> &);

    SquareMatrix<T> getTransposed() const;

    SquareMatrix<T> &transpose();

    T getDeterminant() const;

    SquareMatrix<T> getInverse() const;

    SquareMatrix<T> &invert();

    T getTrace() const;
};

template <typename T>
SquareMatrix<T> SquareMatrix<T>::my_minor(const int without_row, const int without_column) const {
    SquareMatrix<T> result(getSize() - 1);
    for (int row = 0; row < getSize(); ++row)
        for (int column = 0; column < getSize(); ++column)
            if (row != without_row && column != without_column) {
                int i = row;
                int j = column;
                if (row > without_row)
                    i = row - 1;
                if (column > without_column)
                    j = column - 1;
                result.value_[i][j] = Matrix<T>::value_[row][column];
            }
    return result;
}

template <typename T>
SquareMatrix<T>::SquareMatrix(const int size) {
    Matrix<T>::set_size(size, size);
}

template <typename T>
SquareMatrix<T>::SquareMatrix(const SquareMatrix<T> &that) {
    this->row_ = that.row_;
    this->column_ = that.column_;
    this->value_ = new T*[this->row_];
    for (int i = 0; i < this->row_; ++i) {
        this->value_[i] = new T[this->column_];
        for (int j = 0; j < this->column_; ++j)
            this->value_[i][j] = that.value_[i][j];
    }
}

template <typename T>
SquareMatrix<T>::SquareMatrix(const Matrix<T> &that) {
    if (that.getRowsNumber() != that.getColumnsNumber())
        throw MatrixWrongSizeError();
    this->row_ = that.getRowsNumber();
    this->column_ = that.getColumnsNumber();
    this->value_ = new T*[this->row_];
    for (int i = 0; i < this->row_; ++i) {
        this->value_[i] = new T[this->column_];
        for (int j = 0; j < this->column_; ++j)
            this->value_[i][j] = that(i, j);
    }
}

template <typename T>
int SquareMatrix<T>::getSize() const{
    return Matrix<T>::getRowsNumber();
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator+(const SquareMatrix<T> &that) const {
    return Matrix<T>::operator+(that);
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator-(const SquareMatrix<T> &that) const {
    return Matrix<T>::operator-(that);
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const SquareMatrix<T> &that) const {
    return Matrix<T>::operator*(that);
}

template <typename T>
Matrix<T> SquareMatrix<T>::operator*(const Matrix<T> &that) const {
    return Matrix<T>::operator*(that);
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator*(const T t) const {
    return Matrix<T>::operator*(t);
}

template<typename C>
SquareMatrix<C> operator*(const C t, const SquareMatrix<C> &m) {
    return m * t;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator/(const T t) const {
    return Matrix<T>::operator/(t);
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::operator=(const SquareMatrix<T> &that) {
    if (this != &that) {
        this->row_ = that.row_;
        this->column_ = that.column_;
        this->value_ = new T*[this->row_];
        for (int i = 0; i < this->row_; ++i) {
            this->value_[i] = new T[this->column_];
            for (int j = 0; j < this->column_; ++j)
                this->value_[i][j] = that.value_[i][j];
        }
    }
    return *this;
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::operator+=(const SquareMatrix<T> &that) {
    *this = *this + that;
    return *this;
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::operator-=(const SquareMatrix<T> &that) {
    *this = *this - that;
    return *this;
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::operator*=(const SquareMatrix<T> &that) {
    *this = *this * that;
    return *this;
}

template <typename T>
Matrix<T> &SquareMatrix<T>::operator*=(const Matrix<T> &that) {
    *this = *this * that;
    return *this;
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::operator*=(const T t) {
    *this = *this * t;
    return *this;
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::operator/=(const T t) {
    *this = *this / t;
    return *this;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::getTransposed() const {
    return Matrix<T>::getTransposed();
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::transpose() {
    *this = getTransposed();
    return *this;
}

template <typename T>
T SquareMatrix<T>::getDeterminant() const {
    SquareMatrix<T> copy = *this;
    T det = T(1);
    bool flag = true;
    for (int column = 0; flag && column < getSize(); ++column) {
        int early_row = column;
        int row = column;
        while (row < copy.getSize() && copy.value_[row][column] == T(0)) {
            ++row;
        }
        if (row == copy.getSize()) {
            flag = false;
        }
        else {
            if (row != early_row) {
                for (int j = 0; j < getSize(); ++j)
                    std::swap(copy.value_[row][j], copy.value_[early_row][j]);
                det *= -1;
            }
            for (int i = early_row + 1; i < copy.getSize(); ++i) {
                T t = copy.value_[i][column] / copy.value_[early_row][column];
                for (int j = column; j < getSize(); ++j)
                    copy.value_[i][j] -= copy.value_[early_row][j] * t;
            }
        }
        det *= copy.value_[early_row][column];
    }
    if (!flag) {
        det = T(0);
    }
    return det;
}

template <typename T>
SquareMatrix<T> SquareMatrix<T>::getInverse() const {
    if (getDeterminant() == T(0))
        throw MatrixIsDegenerateError();
    SquareMatrix<T> result(getSize());
    T det = getDeterminant();
    for (int i = 0; i < getSize(); ++i)
        for (int j = 0; j < getSize(); ++j) {
            result.value_[i][j] = my_minor(i, j).getDeterminant() / det;
            if ((i + j) % 2)
                result.value_[i][j] *= -1;
        }
    result.transpose();
    return result;
}

template <typename T>
SquareMatrix<T> &SquareMatrix<T>::invert() {
    *this = getInverse();
    return *this;
}

template <typename T>
T SquareMatrix<T>::getTrace() const {
    T sum = T(0);
    for (int i = 0; i < getSize(); ++i)
        sum += Matrix<T>::value_[i][i];
    return sum;
}

class RationalDivisionByZero : std::logic_error {
public:
    RationalDivisionByZero();
};

RationalDivisionByZero::RationalDivisionByZero(): std::logic_error("error") {}


class Rational {

    friend std::ostream &operator<<(std::ostream &os, const Rational &a);

    friend std::istream &operator>>(std::istream &is, Rational &a);

private:

    long long numerator_, denominator_;

    long long gcd(long long a, long long b) const;

    void normalize();

    bool isNegative() const;

public:

    Rational();

    Rational(const long long &numerator, const long long &denominator = 1);

    Rational(const Rational &copy);

    const Rational operator-() const;

    const Rational operator+() const;

    Rational &operator+=(const Rational &a);

    Rational &operator-=(const Rational &a);

    Rational &operator*=(const Rational &a);

    Rational &operator/=(const Rational &a);

    Rational &operator=(const Rational &a);

    Rational &operator++();

    Rational &operator--();

    const Rational operator++(int);

    const Rational operator--(int);

    int compare(const Rational &a) const;

    long long getNumerator() const;

    long long getDenominator() const;

};

const Rational operator+(const Rational &a, const Rational &b);

const Rational operator-(const Rational &a, const Rational &b);

const Rational operator*(const Rational &a, const Rational &b);

const Rational operator/(const Rational &a, const Rational &b);

bool operator>(const Rational &a, const Rational &b);

bool operator<(const Rational &a, const Rational &b);

bool operator>=(const Rational &a, const Rational &b);

bool operator<=(const Rational &a, const Rational &b);

bool operator==(const Rational &a, const Rational &b);

bool operator!=(const Rational &a, const Rational &b);


Rational::Rational() : numerator_(0), denominator_(1) {}

Rational::Rational(const long long &numerator, const long long &denominator) : numerator_(numerator),
                                                                               denominator_(denominator) {
    if (denominator == 0) {
        throw RationalDivisionByZero();
    }

    normalize();
}

Rational::Rational(const Rational &copy) : numerator_(copy.numerator_),
                                           denominator_(copy.denominator_) {}

long long Rational::gcd(long long a, long long b) const {
    if (a < 0)
        a = -a;
    if (b < 0)
        b = -b;

    long long temp;

    if (a < b) {
        temp = a;
        a = b;
        b = temp;
    }

    while (b) {
        a %= b;
        temp = a;
        a = b;
        b = temp;
    }

    return a;
}

void Rational::normalize() {
    long long gcd = Rational::gcd(numerator_, denominator_);

    numerator_ /= gcd;
    denominator_ /= gcd;

    if (denominator_ < 0) {
        denominator_ = -denominator_;
        numerator_ = -numerator_;
    }
}

bool Rational::isNegative() const {
    return numerator_ < 0;
}

const Rational Rational::operator-() const {
    Rational result = *this;
    result.numerator_ = -result.numerator_;
    return result;
}

const Rational Rational::operator+() const {
    return *this;
}

Rational &Rational::operator+=(const Rational &a) {
    numerator_ *= a.denominator_;
    numerator_ += a.numerator_ * denominator_;
    denominator_ *= a.denominator_;

    normalize();

    return *this;
}

Rational &Rational::operator-=(const Rational &a) {
    *this += -a;
    return *this;
}

Rational &Rational::operator*=(const Rational &a) {
    numerator_ *= a.numerator_;
    denominator_ *= a.denominator_;

    normalize();

    return *this;
}

Rational &Rational::operator/=(const Rational &a) {
    if (a.numerator_ == 0) {
        throw RationalDivisionByZero();
    }

    numerator_ *= a.denominator_;
    denominator_ *= a.numerator_;

    normalize();

    return *this;
}

Rational &Rational::operator=(const Rational &a) {
    numerator_ = a.numerator_;
    denominator_ = a.denominator_;

    return *this;
}

long long Rational::getNumerator() const {
    return numerator_;
}

long long Rational::getDenominator() const {
    return denominator_;
}

int Rational::compare(const Rational &a) const {
    if (isNegative() != a.isNegative())
        return a.isNegative() ? 1 : -1;

    long long numerator1 = numerator_ * a.denominator_;
    long long numerator2 = a.numerator_ * denominator_;

    if (numerator1 == numerator2)
        return 0;
    else if (numerator1 < numerator2)
        return -1;
    else
        return 1;
}

std::ostream &operator<<(std::ostream &os, const Rational &a) {
    if (a.getNumerator() != 0) {
        if (a.getDenominator() != 1) {
            os << a.getNumerator() << '/' << a.getDenominator();
        }
        else {
            os << a.getNumerator();
        }
    }
    else {
        os << 0;
    }
    return os;
}

std::istream &operator>>(std::istream &is, Rational &a) {
    char string[30];
    is >> string;
    int numerator = 0;
    int denomenator = 0;
    int sign = 1;
    while (string[denomenator] != '/' && string[denomenator] != 0) {
        if (string[denomenator] == '-') {
            sign = -1;
        }
        else {
            numerator = numerator * 10 + string[denomenator] - '0';
        }
        denomenator++;
    }
    numerator *= sign;
    int q = 0;
    if (string[denomenator] == 0) {
        q = 1;
    }
    else {
        denomenator++;
        while (string[denomenator] != 0) {
            q = q * 10 + string[denomenator] - '0';
            denomenator++;
        }
    }
    a.numerator_ = numerator;
    a.denominator_ = q;
    return is;
}

Rational &Rational::operator++() {
    *this += 1;
    return *this;
}

Rational &Rational::operator--() {
    *this -= 1;
    return *this;
}

const Rational Rational::operator++(int) {
    Rational a = *this;
    *this += 1;
    return a;
}

const Rational Rational::operator--(int) {
    Rational a = *this;
    *this -= 1;
    return a;
}

const Rational operator+(const Rational &a, const Rational &b) {
    Rational result = a;
    result += b;

    return result;
}

const Rational operator-(const Rational &a, const Rational &b) {
    Rational result = a;
    result -= b;

    return result;
}

const Rational operator*(const Rational &a, const Rational &b) {
    Rational result = a;
    result *= b;

    return result;
}

const Rational operator/(const Rational &a, const Rational &b) {
    Rational result = a;
    result /= b;

    return result;
}

bool operator>(const Rational &a, const Rational &b) {
    return a.compare(b) > 0;
}

bool operator<(const Rational &a, const Rational &b) {
    return a.compare(b) < 0;
}

bool operator>=(const Rational &a, const Rational &b) {
    return a.compare(b) >= 0;
}

bool operator<=(const Rational &a, const Rational &b) {
    return a.compare(b) <= 0;
}

bool operator==(const Rational &a, const Rational &b) {
    return a.compare(b) == 0;
}

bool operator!=(const Rational &a, const Rational &b) {
    return a.compare(b) != 0;
}

int main() {
    int m, n, p, q;
    cin >> m >> n >> p >> q;

    Matrix<int> A(m, n), B(p, q);
    cin >> A >> B;

    A = A;
    try {
        cout << A + B * 2 - m * A << endl;
        cout << (A -= B += A *= 2) << endl;
        cout << (((A -= B) += A) *= 2) << endl;
    }
    catch (const MatrixWrongSizeError&) {
        cout << "A and B are of different size." << endl;
    }
    B = A;
    cout << B << endl;

    Rational r;
    cin >> r;
    Matrix<Rational> C(m, n), D(p, q);
    cin >> C >> D;
    try {
        cout << C * D << endl;
        cout << (C *= D) << endl;
        cout << C << endl;
    }
    catch (const MatrixWrongSizeError&) {
        cout << "C and D have not appropriate sizes for multiplication." << endl;
    }
    cout << C.getTransposed() * (r * C) << endl;
    cout << C.transpose() << endl;
    cout << C << endl;

    SquareMatrix<Rational> S(m);
    cin >> S;
    SquareMatrix<Rational> P(S);
    const SquareMatrix<Rational>& rS = S;
    cout << rS.getSize() << ' ' << rS.getDeterminant() << ' ' << rS.getTrace() << endl;
    cout << (S = S) * (S + rS) << endl;
    cout << (S *= S) << endl;
    C.transpose();
    cout << rS * C << endl;
    cout << S << endl;
    S = P;
    cout << (Rational(1, 2) * S).getDeterminant() << endl;
    try {
        cout << rS(0, 0) << endl;
        (S(0, 0) *= 2) /= 2;
        cout << rS(0, m) << endl;
    }
    catch (const MatrixIndexError&) {
        cout << "Index out of range." << endl;
    }
    cout << rS.getTransposed() << endl;
    try {
        cout << rS.getInverse() << endl;
        cout << S.invert().getTransposed().getDeterminant() << endl;
        cout << S << endl;
    }
    catch (const MatrixIsDegenerateError&) {
        cout << "Cannot inverse S." << endl;
    }

    return 0;
}
