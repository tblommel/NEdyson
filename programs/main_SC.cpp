#include "params.h"
#include "molNEsim.h"
#include "gfmol/repn.h"
#include <iostream>
#include "gfmol/utils.h"

int main(const int argc, char *const *const argv)
{
  using namespace NEdyson;

  int nt = 128;
  int ntau = 256;
  int nao = 2;
  int k = 5;
  double beta = 5;

  gfmol::Mode mode = gfmol::Mode::GF2;

  dyson dyson_sol(nt, ntau, nao, k, mode);

  GREEN G(nt, ntau, nao, -1);
  GREEN Sigma(nt, ntau, nao, -1);
  ZTensor<3> hmf;
  
  h5e::File repr_f(argv[1]);
  gfmol::ChebyshevRepr frepr(repr_f, "/fermi", gfmol::Stats::Fermi, beta);
  gfmol::ChebyshevRepr brepr(repr_f, "/bose", gfmol::Stats::Bose, beta);

  h5e::File input_f(argv[2]);
  gfmol::HartreeFock hf(input_f);
  gfmol::Simulation<gfmol::ChebyshevRepr> MatSim(hf, frepr, brepr, mode, 0.);

  MatSim.run(2000, 1e-11, nullptr);

  int nL = MatSim.frepr().nl();
  DMatrix Trans(ntau+1, nL);
  for(int t=0; t<=ntau; t++){
    double x = dyson_sol.Convolution().collocation().x_i()(t);
    for(int l=0; l<nL; l++){
      Trans(t,l) = boost::math::chebyshev_t(l,x);
    }
  }
  int nao2 = nao*nao;
  for(int i=0; i<nao2; i++){
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(G.matptr(0)+i, ntau+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(MatSim.gl().data()+i, nL, Eigen::InnerStride<>(nao2));
    Eigen::Map<ZColVector, 0, Eigen::InnerStride<> >(Sigma.matptr(0)+i, ntau+1, Eigen::InnerStride<>(nao2)) = Trans *
      Eigen::Map<const DColVector, 0, Eigen::InnerStride<> >(MatSim.sigmal().data()+i, nL, Eigen::InnerStride<>(nao2));
  }

  h5e::File output_f(argv[3], h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
  h5e::dump(output_f, "G/mat", ZColVectorMap(G.matptr(0), (ntau+1)*nao*nao));
  h5e::dump(output_f, "Sigma/mat", ZColVectorMap(Sigma.matptr(0), (ntau+1)*nao*nao));
  h5e::dump(output_f, "it0B", DColVectorConstMap(dyson_sol.Convolution().collocation().x_i().data(), ntau+1));





/*
  gfmol::HartreeFock hf(input_f);

  if (!hf.transformed())
    throw std::invalid_argument("Input Hartree-Fock should be in orthogonal basis!");

  if (p.repr == "cheb") {
    // load repr
    std::cout << "Loading Chebyshev repr at " << p.repr_file << std::endl;
    h5e::File repr_f(p.repr_file);
    gfmol::ChebyshevRepr frepr(repr_f, "/fermi", gfmol::Stats::Fermi, p.beta);
    gfmol::ChebyshevRepr brepr(repr_f, "/bose", gfmol::Stats::Bose, p.beta);

    // Make the sim
    std::unique_ptr<SimulationBase> p_sim =
        make_simulation(hf, frepr, brepr, p);

    // Run the sim
    p_sim->run();

    // Save the simulation
    if(p.boolOutput) {
      h5e::File fout(p.output, h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
      p_sim->save(fout, "");
    }
    if(p.boolPumpProbe && p.boolOutputPP) {
      h5e::File fout(p.outputPP, h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);
      p_sim->save_PP(fout, "");
    }
  } 
  else if(p.repr == "ir") {
    throw std::runtime_error("ir basis for NEdyson not yet implemented");
  }
*/
  return 0;
}
