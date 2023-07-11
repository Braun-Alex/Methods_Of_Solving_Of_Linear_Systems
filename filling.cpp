#include "filling.h"
double fill_element_of_matrix_A(int i, int j) noexcept
{
    return 140.0/std::pow((3.0+0.5*static_cast<double>(i+1)*static_cast<double>(j+1)), 2) ;
}
void fill_matrix_A(Eigen::MatrixXd &A) noexcept
{
    auto n=A.rows();
    for (int i=0; i<n; i++)
    {
        for (int j=0; j<n; j++)
            A(i, j)=fill_element_of_matrix_A(i, j);
    }
}
void fill_vector_b(Eigen::VectorXd &b) noexcept
{
    auto n=b.rows();
    for (int i=0; i<n; i++) b(i)=i+1;
}
double compute_norm(const Eigen::MatrixXd &matrix) noexcept
{
    auto count_of_rows=matrix.rows(), count_of_columns=matrix.cols();
    Eigen::VectorXd container=Eigen::VectorXd::Zero(count_of_rows, 1);
    for (int i=0; i<count_of_rows; i++)
        for (int j=0; j<count_of_columns; j++)
            container(i)+=std::abs(matrix(i, j));
    return container.maxCoeff();
}