#include "square_root_method.h"
double sgn(double value) noexcept
{
    if (value<0.0) return -1.0;
    if (value==0.0) return 0.0;
    return 1.0;
}
void fill_S(const Eigen::MatrixXd &A, Eigen::MatrixXd &S, Eigen::MatrixXd &D, int i, int j) noexcept
{
    if (i>j) return;
    if (i==0&j==0) S(0, 0)=std::sqrt(std::abs(A(0, 0)));
    else if (i==j)
    {
        double sum=0.0;
        for (int p=0; p<=i-1; p++) sum+=std::pow(S(p, i), 2)*D(p, p);
        S(i, i)=std::sqrt(std::abs(A(i, i)-sum));
    }
    else
    {
        double sum=0.0;
        for (int p=0; p<=i-1; p++) sum+=S(p, i)*D(p, p)*S(p, j);
        S(i, j)=(A(i, j)-sum)/(D(i, i)*S(i, i));
    }
}
void fill_D(const Eigen::MatrixXd &A, Eigen::MatrixXd &S, Eigen::MatrixXd &D, int i) noexcept
{
    if (i==0) D(i, i)=sgn(A(0, 0));
    else
    {
        double sum=0.0;
        for (int p=0; p<=i-1; p++) sum+=std::pow(S(p, i), 2)*D(p, p);
        D(i, i)=sgn(A(i, i)-sum);
    }
}
void Square_Root_Method(const Eigen::MatrixXd &A, const Eigen::VectorXd &b, int count_of_symbols)
noexcept
{
    std::cout<<std::setw(count_of_symbols)<<std::setfill('*')<<'\n'
    <<"You have chose square root method. So, transpose matrix A to check whether we "
    <<"may use square root method at all"<<'\n'
    <<"Matrix A is represented as:"<<'\n'<<A<<'\n';
    auto transposed_A=A.transpose();
    std::cout<<"Transposed matrix A is represented as:"<<'\n'<<transposed_A<<'\n';
    if (A!=transposed_A)
    {
        std::cout<<"Matrix A and transposed matrix A are not equal. Terminate"<<'\n';
        return;
    }
    std::cout<<"Matrices are equal. So, now we need to build matrices S and D. Do it right now"<<'\n';
    auto n=A.rows();
    Eigen::MatrixXd S=Eigen::MatrixXd::Zero(n, n);
    Eigen::MatrixXd D=Eigen::MatrixXd::Zero(n, n);
    for (int i=0; i<n; i++)
    {
        fill_D(A, S, D, i);
        for (int j=0; j<n; j++) fill_S(A, S, D, i, j);
    }
    std::cout<<"Matrix S is represented as:"<<'\n'<<S<<'\n'
    <<"Matrix D is represented as:"<<'\n'<<D<<'\n'
    <<"So, now we need to transpose matrix S. Transposed matrix S is represented as:"<<'\n';
    Eigen::MatrixXd transposed_S=S.transpose();
    std::cout<<transposed_S<<'\n'
    <<"And now multiply transposed matrix S and matrix D. The given matrix is represented as:"<<'\n';
    Eigen::MatrixXd M=transposed_S*D;
    std::cout<<M<<'\n'
    <<"Let M be multiplied matrix. Find vector y from matrix equation M * y = b. So, vector y is "
    <<"represented as:"<<'\n';
    Eigen::VectorXd y=M.colPivHouseholderQr().solve(b);
    std::cout<<y<<'\n'
    <<"So, final stage. From matrix equation equation S * x = y find necessary vector x. The root x is "
    <<"represented as:"<<'\n';
    Eigen::VectorXd x=S.colPivHouseholderQr().solve(y);
    std::cout<<x<<'\n'
    <<std::setw(count_of_symbols)<<std::setfill('-')<<'\n'
    <<"Let r be (A * x) - b. Compute r. So, r is represented as:"<<'\n';
    Eigen::VectorXd r=A*x-b;
    std::cout<<r<<'\n'
    <<std::setw(count_of_symbols)<<std::setfill('-')<<'\n';
    double K=compute_norm(A), L=compute_norm(A.inverse()), U=K*L;
    std::cout<<"Find condition number of matrix A using the norm. So, "<<'\n'
    <<"the norm of matrix A: "<<K<<". Let it be K"<<'\n'
    <<"the norm of inverse matrix A: "<<L<<". Let it be L"<<'\n'
    <<"Let U be condition number. Then U = K * L. Finally, U equals: "<<U<<'\n'
    <<std::setw(count_of_symbols)<<std::setfill('-')<<'\n'
    <<"Find determinant of matrix A using matrices S and D. So, the determinant equals: "<<'\n';
    double determinant=1.0;
    for (int i=0; i<n; i++) determinant*=D(i, i)*std::pow(S(i, i), 2);
    std::cout<<determinant<<'\n'
    <<std::setw(count_of_symbols)<<std::setfill('-')<<'\n'
    <<"Find inverse matrix A. Do it right now. So, "
    <<"inverse matrix A is represented as:"<<'\n';
    Eigen::MatrixXd inverse_A=Eigen::MatrixXd::Zero(n, n);
    for (int k=0; k<n; k++)
    {
        Eigen::VectorXd e=Eigen::VectorXd::Zero(n, 1);
        e(k)=1;
        y=M.colPivHouseholderQr().solve(e);
        Eigen::VectorXd v=S.colPivHouseholderQr().solve(y);
        for (int l=0; l<n; l++) inverse_A(k, l)=v(l);
    }
    std::cout<<inverse_A<<'\n'
    <<"Multiply inverse matrix A and matrix A. The given matrix is represented as:"<<'\n'
    <<inverse_A*A<<'\n'<<std::setw(count_of_symbols)<<std::setfill('*')<<'\n';
}