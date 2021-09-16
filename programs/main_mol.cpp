#include "params.h"
#include "molNEsim.h"
#include "gfmol/repn.h"
#include <iostream>
#include "gfmol/utils.h"

int main(const int argc, char *const *const argv)
{
  using namespace NEdyson;
  NEdyson::Params p = NEdyson::parse_args(argc, argv);
  std::cout<<"hello"<<std::endl;
  std::cout<<p<<std::endl;

  h5e::File input_f(p.hf_input);
  std::cout<<"hello1"<<std::endl;
  gfmol::HartreeFock hf(input_f);
  std::cout<<"hello2"<<std::endl;

  if (!hf.transformed())
    throw std::invalid_argument("Input Hartree-Fock should be in orthogonal basis!");

  std::cout<<"hello3"<<std::endl;
  if (p.repr == "cheb") {
    // load repr
    std::cout << "Loading Chebyshev repr at " << p.repr_file << std::endl;
    h5e::File repr_f(p.repr_file);
    std::cout<<"hello4"<<std::endl;
    gfmol::ChebyshevRepr frepr(repr_f, "/fermi", gfmol::Stats::Fermi, p.beta);
    std::cout<<"hello5"<<std::endl;
    gfmol::ChebyshevRepr brepr(repr_f, "/bose", gfmol::Stats::Bose, p.beta);
    std::cout<<"hello6"<<std::endl;

    // Make the sim
    std::unique_ptr<SimulationBase> p_sim =
        make_simulation(hf, frepr, brepr, p);
    std::cout<<"hello7"<<std::endl;

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

  return 0;
}
