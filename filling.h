#pragma once
#include <Eigen/Dense>
double fill_element_of_matrix_A(int i, int j) noexcept;
void fill_matrix_A(Eigen::MatrixXd &A) noexcept;
void fill_vector_b(Eigen::VectorXd &b) noexcept;
double compute_norm(const Eigen::MatrixXd &matrix) noexcept;