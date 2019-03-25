#include "functions.hpp"

int main(int argc, char* argv[])
{
  if(argc != 3) {std::cerr << "Enter (1) lattice size, (2) Gamma.\n"; std::exit(1);}
  size = atoi(argv[1]);
  G = atof(argv[2]);

  MatrixXd H0 = construct_h0();
  Eigen::VectorXd eivals_1d = Eigenvalues(H0);

  std::ofstream outfile_dos("1d_dos.dat");

  for(double omega = -3*t; omega < 3*t; omega += 0.01)
  {
    outfile_dos << omega << " " << dos(omega, eivals_1d) << endl;
  }
  
  return 0;
}