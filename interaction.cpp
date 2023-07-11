#include "interaction.h"
#include "square_root_method.h"
#include "jacobi_method.h"
Interaction::Interaction(int entered_count_of_symbols,
double entered_accuracy) noexcept: count_of_symbols(entered_count_of_symbols),
accuracy(entered_accuracy) {}
void Interaction::Set_N() noexcept
{
    std::string count;
    bool check=false;
    std::cout<<"Enter the size N of matrix NxN, where natural N is between 7 and 15:"<<'\n';
    do
    {
        std::getline(std::cin, count);
        if (!count.empty()&&count.size()<=2&&((count.size()==1&&count[0]>=55&&count[0]<=57)
        ||(count.size()==2&&count[0]==49&&count[1]>=48&&count[1]<=53))) check=true;
        else std::cout<<"Natural N must be between 7 and 15! "
        <<"Please, enter it again:"<<'\n';
    }
    while (!check);
    if (count.size()==1) N=count[0]-48;
    else N=count[1]-38;
}
void Interaction::Initialize_And_Print_Matrix_A(const std::function<void(Eigen::MatrixXd&)> &fill_matrix)
noexcept
{
A=Eigen::MatrixXd::Zero(N, N);
fill_matrix(A);
std::cout<<"Matrix A is represented as:"<<'\n'<<A<<'\n';
}
void Interaction::Initialize_And_Print_Vector_b(const std::function<void(Eigen::VectorXd&)> &fill_vector)
noexcept
{
b=Eigen::VectorXd::Zero(N, 1);
fill_vector(b);
std::cout<<"Vector b is represented as:"<<'\n'<<b<<'\n';
}
void Interaction::Enter_Numeric_Method() noexcept
{
std::cout<<"Enter the numeric method: square root method (S) or Jacobi method (J):"<<'\n';
bool check=false;
do
{
    std::getline(std::cin, symbol);
    if (symbol=="S"||symbol=="J") check=true;
    else std::cout<<"Entered incorrect data! Please, enter S or J:"<<'\n';
}
while (!check);
}
void Interaction::Execute_Numeric_Method() noexcept
{
if (symbol=="S") Square_Root_Method(A, b, count_of_symbols);
else Jacobi_Method(A, b, count_of_symbols, accuracy);
}
void Interaction::Ask_For_Continue() noexcept
{
std::cout<<R"(Do you want to continue computing? 'Y' means "yes", another means "no":)"<<'\n';
std::getline(std::cin, symbol);
}
std::string Interaction::Get_Answer() const noexcept
{
return symbol;
}