#include "functions.hpp"

int main(int argc, char* argv[])
{
  if(argc != 3) {std::cerr << "Enter (1) lattice size, (2) Gamma.\n"; std::exit(1);}
  size = atoi(argv[1]);
  G = atof(argv[2]);

  MatrixXd H0 = construct_h0_2d();
  Eigen::VectorXd eivals_2d = Eigenvalues(H0);

  std::ofstream outfile_dos("2d_dos.dat");

  for(double omega = -5*t; omega < 5*t; omega += 0.01)
  {
    outfile_dos << omega << " " << dos(omega, eivals_2d) << endl;
  }
  
  return 0;
}