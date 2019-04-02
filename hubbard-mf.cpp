#include "functions.hpp"

int main(int argc, char* argv[])
{
  if(argc!=5){std::cerr << "Enter (1) Lattice size, (2) U, (3) temperature, (4) max steps.\n"; exit(1);}
  size = atoi(argv[1]);
  U = atof(argv[2]);
  double temperature = atof(argv[3]);
  int MAX_STEPS = atoi(argv[4]);
  G = 8.0/(size*size);

  VectorXd last_eivals = VectorXd::Zero(2*size*size);
  VectorXd n_up_avg = VectorXd::Constant(size*size, 0.49);
  VectorXd n_down_avg = VectorXd::Constant(size*size,0.51);

  MatrixXd H0 = construct_h0_2d();
  std::vector < std::pair<double,double> > DoS;

  for(double omega = -4*t-U; omega< 4*t+U; omega += 0.1)
  {
    DoS.push_back(std::make_pair(omega,0));
  }

  int step = 0;
  for(step = 0; step < MAX_STEPS; step++)
  {
    cout << n_up_avg.unaryExpr(&filter_d).transpose() << endl << n_down_avg.unaryExpr(&filter_d).transpose() << endl;

    MatrixXd Hmf = H0 + H_int(n_down_avg, n_up_avg);
    std::pair<MatrixXcd, VectorXd> mf_spectrum = Eigenspectrum(Hmf);
    double mu = get_mu(temperature, mf_spectrum.second);
    cout << "mu = " << mu << endl;
    VectorXd eivals_minus_mu = mf_spectrum.second - VectorXd::Constant(mf_spectrum.second.size(),mu);

    for(int i=0; i<n_up_avg.size(); i++) n_up_avg(i) = get_n_sigma(eivals_minus_mu, mf_spectrum.first, temperature, i);
    for(int i=0; i<n_up_avg.size(); i++) n_down_avg(i) = get_n_sigma(eivals_minus_mu, mf_spectrum.first, temperature, i+size*size);
    
    if((mf_spectrum.second-last_eivals).norm() < 1e-10)
    {
      cout << "converged after " << step << " steps." << endl;
      break;
    } 
    
    for(auto i:DoS) 
    {
      i.second += dos(i.first, eivals_minus_mu);
    }

    last_eivals = mf_spectrum.second;
    cout << endl;
  }
  
  // for(auto i:DoS) {i.second = i.second/double(step);}

  std::ofstream  fout("mf_dos.dat");
  for(auto i:DoS)
  {
    fout << i.first << " " << i.second << endl;
  }

  // cout << endl << n_up_avg.unaryExpr(&filter_d).transpose() << endl << n_down_avg.unaryExpr(&filter_d).transpose() << endl;


}