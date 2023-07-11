#pragma once
#include "filling.h"
#include <iostream>
#include <iomanip>
#include <Eigen/Dense>
void Jacobi_Method(Eigen::MatrixXd &A, Eigen::VectorXd &b, int count_of_symbols, double accuracy) noexcept;