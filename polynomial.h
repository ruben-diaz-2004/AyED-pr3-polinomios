// AUTOR: Rubén Díaz Marrero
// FECHA: 14/03/2023
// EMAIL: alu0101552613@ull.edu.es
// VERSION: 1.0
// ASIGNATURA: Algoritmos y Estructuras de Datos
// PRÁCTICA Nº: 3
// ESTILO: Google C++ Style Guide
// COMENTARIOS:
// 

#ifndef POLYNOMIAL_H_
#define POLYNOMIAL_H_

#include <iostream>
#include <math.h>  // fabs, pow

#include "vector_t.h"
#include "sparse_vector_t.h"

// Clase para polinomios basados en vectores densos de doubles
class Polynomial : public vector_t<double> {
 public:
  // constructores
  Polynomial(const int n = 0) : vector_t<double>(n) {};
  Polynomial(const Polynomial& pol)
      : vector_t<double>(pol) {}; // constructor de copia

  // destructor
  ~Polynomial() {};

  // E/S
  void Write(std::ostream& = std::cout, const double eps = EPS) const;
  
  // operaciones
  double Eval(const double) const;
  bool IsEqual(const Polynomial&, const double = EPS) const;
 };


// Clase para polinomios basados en vectores dispersos
class SparsePolynomial : public sparse_vector_t {
 public:
  // constructores
  SparsePolynomial(const int n = 0) : sparse_vector_t(n) {};
  SparsePolynomial(const Polynomial& pol) : sparse_vector_t(pol) {};
  SparsePolynomial(const SparsePolynomial&);  // constructor de copia

  // destructor
  ~SparsePolynomial() {};

  // E/S
  void Write(std::ostream& = std::cout) const;
  
  // operaciones
  double Eval(const double) const;
  bool IsEqual(const SparsePolynomial&, const double = EPS) const;
  bool IsEqual(const Polynomial&, const double = EPS) const;
};

// E/S
void Polynomial::Write(std::ostream& os, const double eps) const {
  os << get_size() << ": [ ";
  bool first{true};
  for (int i{0}; i < get_size(); i++)
    if (IsNotZero(at(i), eps)) {
      os << (!first ? " + " : "") << at(i)
	 << (i > 1 ? " x^" : (i == 1) ? " x" : "");
      if (i > 1)
	os << i;
      first = false;
    }
  os << " ]" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const Polynomial& p) {
  p.Write(os);
  return os;
}

// Operaciones con polinomios

// Evaluación de un polinomio representado por vector denso
double Polynomial::Eval(const double x) const {
  double result{0.0};
  for (int i{0}; i < get_size(); ++i) {
  result += get_val(i) * pow(x, i);
  }
  return result;
}

// Comparación si son iguales dos polinomios representados por vectores densos
bool Polynomial::IsEqual(const Polynomial& pol, const double eps) const {
  bool differents = false;
  if (get_size() < pol.get_size()) for (int i{get_size()}; i < pol.get_size(); ++i) if (pol.get_val(i) != 0) return false;
  else if (pol.get_size() < get_size()) for (int i{pol.get_size()}; i < get_size(); ++i) if (get_val(i) != 0) return false;
  for (int i{0}; i < get_size() && i < pol.get_size(); ++i) if (!(fabs(Eval(i) - pol.Eval(i)) < eps)) differents = true;
  return !differents;
}

// constructor de copia
SparsePolynomial::SparsePolynomial(const SparsePolynomial& spol) {
  *this = spol;   // se invoca directamente al operator=
}

// E/S
void SparsePolynomial::Write(std::ostream& os) const {
  os << get_n() << "(" << get_nz() << "): [ ";
  bool first{true};
  for (int i{0}; i < get_nz(); i++) {
    int inx{at(i).get_inx()};
    os << (!first ? " + " : "") << at(i).get_val()
       << (inx > 1 ? " x^" : (inx == 1) ? " x" : "");
    if (inx > 1)
      os << inx;
    first = false;
  }
  os << " ]" << std::endl;
}

std::ostream& operator<<(std::ostream& os, const SparsePolynomial& p) {
  p.Write(os);
  return os;
}

// Operaciones con polinomios

// Evaluación de un polinomio representado por vector disperso
double SparsePolynomial::Eval(const double x) const {
  double result{0.0};
  for (int i{0}; i < get_nz(); ++i) {
   result += at(i).get_val() * pow(x, at(i).get_inx());
  }
  // poner el código aquí
  return result;
}
  
// Comparación si son iguales dos polinomios representados por vectores dispersos
bool SparsePolynomial::IsEqual(const SparsePolynomial& spol
			       , const double eps) const {
  bool differents = false;
  if (get_nz() < spol.get_nz()) for (int i{get_nz()}; i < spol.get_nz(); ++i) if (spol.at(i).get_val() != 0) return false;
  else if (spol.get_nz() < get_nz()) for (int i{spol.get_nz()}; i < get_nz(); ++i) if (at(i).get_val() != 0) return false;
  for (int i{0}; i < get_nz() && spol.get_nz(); ++i) if (!(fabs(Eval(i) - spol.Eval(i)) < eps)) differents = true;
  return !differents;
}

// Comparación si son iguales dos polinomios representados por
// vector disperso y vector denso
bool SparsePolynomial::IsEqual(const Polynomial& pol, const double eps) const {
  bool differents = false;
  if (get_nz() > pol.get_size()) return false;
  if (get_n() < pol.get_size()) for (int i{get_n()}; i < pol.get_size(); ++i) if (pol.get_val(i) != 0) return false;
  int sparce_position{0};
  int pol_position{0};
  while (sparce_position < get_nz() && pol_position < pol.get_size()) {
    if (at(sparce_position).get_inx() == pol_position) {
      if (at(sparce_position).get_val() != pol.get_val(pol_position)) return false;
      ++sparce_position;
    } else if (pol_position < at(sparce_position).get_inx()) if (pol.get_val(pol_position) != 0) return false;
    ++pol_position;
  }
  return !differents;
}


#endif  // POLYNOMIAL_H_
