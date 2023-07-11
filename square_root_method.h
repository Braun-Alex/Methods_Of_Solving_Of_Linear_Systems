#pragma once
#include "filling.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
double sgn(double value) noexcept;
void fill_S(const Eigen::MatrixXd &A, Eigen::MatrixXd &S, Eigen::MatrixXd &D, int i, int j) noexcept;
void fill_D(const Eigen::MatrixXd &A, Eigen::MatrixXd &S, Eigen::MatrixXd &D, int i) noexcept;
void Square_Root_Method(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, int count_of_symbols) noexcept;