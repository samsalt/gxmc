#include "gxcInitiate.h"
#include "gxcProblem.h"
int main()
{

  LogInitiate();
  gxcProblem Problem;
  Problem.TimeIntegration();

  return 0;
}