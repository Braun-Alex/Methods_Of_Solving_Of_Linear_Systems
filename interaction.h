#pragma once
#include <Eigen/Dense>
class Interaction
{
public:
Interaction(int entered_count_of_symbols, double entered_accuracy) noexcept;
void Set_N() noexcept;
void Initialize_And_Print_Matrix_A(const std::function<void(Eigen::MatrixXd&)> &fill_matrix) noexcept;
void Initialize_And_Print_Vector_b(const std::function<void(Eigen::VectorXd&)> &fill_vector) noexcept;
void Enter_Numeric_Method() noexcept;
void Execute_Numeric_Method() noexcept;
void Ask_For_Continue() noexcept;
std::string Get_Answer() const noexcept;
private:
int N, count_of_symbols;
double accuracy;
std::string symbol;
Eigen::MatrixXd A;
Eigen::VectorXd b;
};