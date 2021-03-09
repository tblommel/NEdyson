#include "params.h"
#include "molNEsim.h"
#include "gfmol/repn.h"
#include <iostream>
#include "gfmol/utils.h"

int main(const int argc, char *const *const argv)
{
  using namespace NEdyson;
  NEdyson::Params p = NEdyson::parse_args(argc, argv);

  std::cout<<p<<std::endl;

  h5e::File input_f(p.hf_input);
  gfmol::HartreeFock hf(input_f);

  if (!hf.transformed())
    throw std::invalid_argument("Input Hartree-Fock should be in orthogonal basis!");

  gfmol::Mode mode;
  if (p.mode == "GF2")
    mode = gfmol::Mode::GF2;
  else if (p.mode == "GW")
    throw std::runtime_error("GW solver not implemented");

  h5e::File fout(p.output, h5e::File::Overwrite | h5e::File::ReadWrite | h5e::File::Create);

  if (p.repr == "cheb") {
    // load repr
    std::cout << "Loading Chebyshev repr at " << p.repr_file << std::endl;
    h5e::File repr_f(p.repr_file);
    gfmol::ChebyshevRepr frepr(repr_f, "/fermi", gfmol::Stats::Fermi, p.beta);
    gfmol::ChebyshevRepr brepr(repr_f, "/bose", gfmol::Stats::Bose, p.beta);

    // Make the sim
    std::unique_ptr<SimulationBase> p_sim =
        make_simulation(p.tti, p.unrestricted, p.decomposed, hf, frepr, brepr, mode, p.nt, p.ntau, p.k, p.dt, p.maxiter, p.etol, p.BootMaxIter, p.BootMaxErr, p.CorrSteps, p.damping, p.decomp_prec, p.hf, p.boolPumpProbe, p.PumpProbe_file, p.hf_input, p.lPumpProbe, p.nPumpProbe);

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
