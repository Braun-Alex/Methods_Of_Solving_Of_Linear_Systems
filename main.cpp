#include "interaction.h"
#include "filling.h"
void Compute(Interaction &interaction)
{
interaction.Set_N();
interaction.Initialize_And_Print_Matrix_A(fill_matrix_A);
interaction.Initialize_And_Print_Vector_b(fill_vector_b);
interaction.Enter_Numeric_Method();
interaction.Execute_Numeric_Method();
interaction.Ask_For_Continue();
if (interaction.Get_Answer()=="Y") Compute(interaction);
}
int main()
{
Interaction interaction(120, 0.0001);
Compute(interaction);
return 0;
}