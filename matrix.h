//
// Created by Алексей on 06.06.2024.
//
// #define MATRIX_SQUARE_MATRIX_IMPLEMENTED

#ifndef MATRIX_H
#define MATRIX_H
#include <iostream>

class MatrixOutOfRange : public std::exception {
  const char *what() const noexcept override {
    return "MatrixOutOfRange";
  }
};

class DivisionByZero : public std::exception {
  const char *what() const noexcept override {
    return "DivisionByZero";
  }
};

template <class T, size_t Rows, size_t Columns>
class Matrix {
 public:
  T matrix[Rows][Columns];

  size_t RowsNumber() const {
    return Rows;
  }

  size_t ColumnsNumber() const {
    return Columns;
  }

  T &operator()(const size_t &row, const size_t &column) {
    return matrix[row][column];
  }

  const T &operator()(const size_t &row, const size_t &column) const {
    return matrix[row][column];
  }

  const T &At(size_t row, size_t column) const;

  T &At(size_t row, size_t column);

  Matrix operator*(T number) const;

  Matrix &operator*=(T number);

  // friend Matrix operator*(T number, const Matrix &matrix);

  Matrix operator/(T number) const;

  Matrix &operator/=(T number);

  Matrix operator+(const Matrix &matr) const;
  Matrix &operator+=(const Matrix &matr);

  Matrix operator-(const Matrix &matr) const;
  Matrix &operator-=(const Matrix &matr);

  template <size_t Columns2>
  Matrix<T, Rows, Columns2> operator*(const Matrix<T, Columns, Columns2> &matrix) const;

  Matrix &operator*=(const Matrix<T, Columns, Columns> &matrix);

  bool operator==(const Matrix &) const;
  bool operator!=(const Matrix &) const;
};

template <class T, size_t Rows, size_t Columns>
Matrix<T, Columns, Rows> GetTransposed(const Matrix<T, Rows, Columns> &matr) {
  Matrix<T, Columns, Rows> matrix = {};
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      matrix(j, i) = matr(i, j);
    }
  }
  return matrix;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> Matrix<T, Rows, Columns>::operator+(const Matrix<T, Rows, Columns> &matr) const {
  Matrix<T, Rows, Columns> result = {};
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      result(i, j) = this->operator()(i, j) + matr(i, j);
    }
  }
  return result;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> &Matrix<T, Rows, Columns>::operator+=(const Matrix<T, Rows, Columns> &matr) {
  Matrix<T, Rows, Columns> result = {};
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      this->operator()(i, j) += matr(i, j);
    }
  }
  return *this;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> &Matrix<T, Rows, Columns>::operator-=(const Matrix<T, Rows, Columns> &matr) {
  Matrix<T, Rows, Columns> result = {};
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      this->operator()(i, j) -= matr(i, j);
    }
  }
  return *this;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> Matrix<T, Rows, Columns>::operator-(const Matrix<T, Rows, Columns> &matr) const {
  return *this + (matr * (-1));
}

template <class T, size_t Rows, size_t Columns>
const T &Matrix<T, Rows, Columns>::At(size_t row, size_t column) const {
  if (row > Rows || column > Columns) {
    throw MatrixOutOfRange{};
  }
  return this->operator()(row, column);
}

template <class T, size_t Rows, size_t Columns>
T &Matrix<T, Rows, Columns>::At(size_t row, size_t column) {
  if (row > Rows || column > Columns) {
    throw MatrixOutOfRange{};
  }
  return this->operator()(row, column);
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> Matrix<T, Rows, Columns>::operator*(T number) const {
  Matrix matr = {};
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      matr(i, j) = this->operator()(i, j) * number;
    }
  }
  return matr;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> &Matrix<T, Rows, Columns>::operator*=(T number) {
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      this->operator()(i, j) *= number;
    }
  }
  return *this;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> Matrix<T, Rows, Columns>::operator/(T number) const {
  if (number == 0) {
    throw DivisionByZero{};
  }
  Matrix matr;
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      matr(i, j) = this->operator()(i, j) / number;
    }
  }
  return matr;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> &Matrix<T, Rows, Columns>::operator/=(T number) {
  if (number == 0) {
    throw DivisionByZero{};
  }
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      this->operator()(i, j) /= number;
    }
  }
  return *this;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> operator*(T number, const Matrix<T, Rows, Columns> &matrix) {
  return matrix * number;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> operator/(T number, Matrix<T, Rows, Columns> matrix) {
  return matrix / number;
}

template <class T, size_t Rows, size_t Columns>
template <size_t Columns2>
Matrix<T, Rows, Columns2> Matrix<T, Rows, Columns>::operator*(const Matrix<T, Columns, Columns2> &matrix1) const {
  Matrix<T, Rows, Columns2> result_matrix = {};
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns2; ++j) {
      for (size_t k = 0; k < Columns; ++k) {
        result_matrix(i, j) += this->operator()(i, k) * matrix1(k, j);
      }
    }
  }
  return result_matrix;
}

template <class T, size_t Rows, size_t Columns>
Matrix<T, Rows, Columns> &Matrix<T, Rows, Columns>::operator*=(const Matrix<T, Columns, Columns> &matrix1) {
  return *this = *this * matrix1;
}

template <class T, size_t Rows, size_t Columns>
bool Matrix<T, Rows, Columns>::operator==(const Matrix<T, Rows, Columns> &matrix) const {
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      if (this->operator()(i, j) != matrix(i, j)) {
        return false;
      }
    }
  }
  return true;
}

template <class T, size_t Rows, size_t Columns>
bool Matrix<T, Rows, Columns>::operator!=(const Matrix<T, Rows, Columns> &matrix) const {
  return !(*this == matrix);
}

template <class T, size_t Rows, size_t Columns>
std::ostream &operator<<(std ::ostream &stream, const Matrix<T, Rows, Columns> &matrix) {
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      stream << matrix(i, j);
      if (j != Columns - 1) {
        stream << " ";
      }
    }
    stream << "\n";
  }
  return stream;
}

template <class T, size_t Rows, size_t Columns>
std::istream &operator>>(std ::istream &stream, Matrix<T, Rows, Columns> &matrix) {
  for (size_t i = 0; i < Rows; ++i) {
    for (size_t j = 0; j < Columns; ++j) {
      stream >> matrix(i, j);
    }
  }
  return stream;
}
#endif  // MATRIX_H
