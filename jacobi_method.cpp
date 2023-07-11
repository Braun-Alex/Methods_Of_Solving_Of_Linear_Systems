#include "jacobi_method.h"
void Jacobi_Method(Eigen::MatrixXd &A, Eigen::VectorXd &b, int  count_of_symbols, double accuracy)
noexcept
{
    std::cout<<std::setw(count_of_symbols)<<std::setfill('*')<<'\n'
    <<"You have chose Jacobi method. So, take a look at matrix A to check whether we "
    <<"may use square Jacobi method at all"<<'\n'
    <<"Matrix A is represented as:"<<'\n'<<A<<'\n';
    auto n=A.rows();
    double sum;
    bool converging=true;
    for (int i=0; i<n; i++)
    {
        sum=0.0;
        for (int j=0; j<n; j++) sum+=std::abs(A(i, j));
        sum-=std::abs(A(i, i));
        if (std::abs(A(i, i))<sum)
        {
            std::cout<<"Matrix A is not convergent for Jacobi Method. So, change it"<<'\n';
            converging=false;
            break;
        }
    }
    if (!converging)
    {
        for (int i=0; i<n; i++) A(i, i)+=50;
        std::cout<<"Changed matrix A is already convergent for Jacobi Method and is represented as:"
        <<'\n'<<A<<'\n';
    }
    else std::cout<<"Matrix A is convergent for Jacobi Method"<<'\n';
    std::cout<<"Now compute matrices D, L and U. Afterwards, compute matrix "
    <<"T = -(inverse_D) * (L + U) and vector C = inverse_D * b"<<'\n';
    Eigen::MatrixXd D=Eigen::MatrixXd::Zero(n, n);
    for (int i=0; i<n; i++) D(i, i)=A(i, i);
    std::cout<<"Matrix D is represented as:"<<'\n'<<D<<'\n';
    Eigen::MatrixXd L=Eigen::MatrixXd::Zero(n, n);
    for (int j=0; j<n; j++)
        for (int i=j+1; i<n; ++i) L(i, j)=A(i, j);
    std::cout<<"Matrix L is represented as:"<<'\n'<<L<<'\n';
    Eigen::MatrixXd U=Eigen::MatrixXd::Zero(n, n);
    for (int j=n-1; j>0; j--)
        for (int i=j-1; i>=0; i--) U(i, j)=A(i, j);
    std::cout<<"Matrix U is represented as:"<<'\n'<<U<<'\n';
    Eigen::MatrixXd T=Eigen::MatrixXd::Zero(n, n), inverse_D=D.inverse();
    T=(-1)*inverse_D*(L+U);
    std::cout<<"Matrix T is represented as:"<<'\n'<<T<<'\n';
    Eigen::VectorXd C=Eigen::MatrixXd::Zero(n, 1);
    C=inverse_D*b;
    std::cout<<"Vector C is represented as:"<<'\n'<<C<<'\n';
    Eigen::VectorXd previous_iteration_x=Eigen::VectorXd::Zero(n, 1);
    Eigen::VectorXd iteration_x=Eigen::VectorXd::Zero(n, 1);
    std::cout<<"Let iteration root be zero-vector, so it is represented at zero iteration as:"
    <<'\n'<<iteration_x<<'\n'
    <<"Start iteration process. It is described as x(k+1) = T * x(k) + C"<<'\n';
    size_t iteration_count=0;
    do
    {
        ++iteration_count;
        previous_iteration_x=std::move(iteration_x);
        iteration_x=T*previous_iteration_x+C;
        std::cout<<"Iteration root at "<<iteration_count<<" iteration:"<<'\n'<<iteration_x<<'\n';
    }
    while (compute_norm(iteration_x-previous_iteration_x)>=accuracy||
    compute_norm(A*iteration_x-b)>=accuracy);
    std::cout<<"The accuracy has been achieved. Iteration process has been stopped. So, the "
    <<"root is represented as:"<<'\n'<<iteration_x<<'\n'
    <<std::setw(count_of_symbols)<<std::setfill('-')<<'\n'
    <<"Let r be (A * x) - b. Compute r. So, r is represented as:"<<'\n';
    Eigen::VectorXd r=A*iteration_x-b;
    std::cout<<r<<'\n'
    <<std::setw(count_of_symbols)<<std::setfill('-')<<'\n'
    <<"Totally needed "<<iteration_count<<" iterations"<<'\n'
    <<std::setw(count_of_symbols)<<std::setfill('*')<<'\n';
}