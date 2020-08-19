//
// Created by katherlee on 2020-04-20.
//

#ifndef _PARAMS_H_
#define _PARAMS_H_

#include <args.hxx>
#include <iniparser.hpp>
#include <iostream>

namespace NEdyson {

struct Params {
  std::string mode;
  bool tti;
  bool unrestricted;
  bool decomposed;
  std::string repr;
  std::string repr_file;
  double decomp_prec;
  double etol;
  int maxiter;
  double beta;
  double damping;
  std::string hf_input;
  std::string output;
  bool checkpoint;
  std::string checkpoint_file;
  bool boolOutput;

  int nt;
  int ntau;
  int k;
  double dt;
  int BootMaxIter;
  double BootMaxErr;
  int CorrSteps;

  int nw;
  double wmax;

  void validate() const
  {
    if (mode != "GW" && mode != "GF2")
      throw std::runtime_error("parameter mode: invalid mode " + mode);

    if (repr != "cheb" && repr != "ir")
      throw std::runtime_error("parameter repr: invalid repr " + repr);

    if (maxiter < 0)
      throw std::runtime_error("parameter maxiter: invalid value " + std::to_string(maxiter));

    if (beta <= 0)
      throw std::runtime_error("parameter beta: invalid value " + std::to_string(beta));

    if (damping < 0 || damping > 1)
      throw std::runtime_error("parameter damping: invalid value " + std::to_string(damping));
    
    if (nt < 0)
      throw std::runtime_error("parameter nt: invalid value " + std::to_string(nt));
      
    if (ntau <= 0)
      throw std::runtime_error("parameter ntau: invalid value " + std::to_string(ntau));

    if (k <= 0)
      throw std::runtime_error("parameter k: invalid value " + std::to_string(k));

    if(nt < k || ntau < k)
      throw std::runtime_error("timesteps less than integrator accuracy: k " + std::to_string(k) + ", nt " + std::to_string(nt) + ", ntau " + std::to_string(ntau));

    if (dt <= 0)
      throw std::runtime_error("parameter dt: invalid value " + std::to_string(dt));
    
    if (BootMaxIter <= 0)
      throw std::runtime_error("parameter BootMaxIter: invalid value " + std::to_string(BootMaxIter));

    if (BootMaxErr <= 0)
      throw std::runtime_error("parameter BootMaxErr: invalid value " + std::to_string(BootMaxErr));

    if (CorrSteps <= 0)
      throw std::runtime_error("parameter CorrSteps: invalid value " + std::to_string(CorrSteps));

    if (nw <= 0 || nw%2 ==0)
      throw std::runtime_error("parameter nw: invalid value " + std::to_string(nw));

    if (wmax <= 0) 
      throw std::runtime_error("parameter wmax: invalid value " + std::to_string(wmax));
  }
};

inline std::ostream &operator<<(std::ostream &os, const Params &p)
{
  os << std::boolalpha;
  os << "Parameters:" << std::endl;
  os << "    tti:             " << p.tti <<std::endl;
  os << "    mode:            " << p.mode << std::endl;
  os << "    unrestricted:    " << p.unrestricted << std::endl;
  os << "    decomposed:      " << p.decomposed << std::endl;
  os << "    repr:            " << p.repr << std::endl;
  os << "    repr_file:       " << p.repr_file << std::endl;
  os << "    decomp_prec:     " << p.decomp_prec << std::endl;
  os << "    etol:            " << p.etol << std::endl;
  os << "    maxiter:         " << p.maxiter << std::endl;
  os << "    beta:            " << p.beta << std::endl;
  os << "    damping:         " << p.damping << std::endl;
  os << "    hf_input:        " << p.hf_input << std::endl;
  os << "    boolOutput:      " << p.boolOutput << std::endl;
  os << "    output:          " << p.output << std::endl;
  os << "    checkpoint:      " << p.checkpoint << std::endl;
  os << "    checkpoint_file: " << p.checkpoint_file << std::endl;
  os << "    nt:              " << p.nt << std::endl;
  os << "    ntau:            " << p.ntau << std::endl;
  os << "    k:               " << p.k << std::endl;
  os << "    dt:              " << p.dt << std::endl;
  os << "    BootMaxIter:     " << p.BootMaxIter << std::endl;
  os << "    BootMaxErr:      " << p.BootMaxErr << std::endl;
  os << "    CorrSteps:       " << p.CorrSteps << std::endl;
  os << "    nw:              " << p.nw << std::endl;
  os << "    wmax:            " << p.wmax << std::endl;
  os << std::noboolalpha;
  return os;
}

// From Sergei
template <typename T>
T extract_value(const INI::File &f, args::ValueFlag<T> &parameter)
{
  return parameter ? parameter.Get() : f.GetValue(parameter.Name(), parameter.GetDefault()).template Get<T>();
}

inline bool extract_value(const INI::File &f, args::Flag &parameter)
{
  return parameter ? parameter.Get() : f.GetValue(parameter.Name(), false).Get<bool>();
}

inline Params parse_args(const int argc, char *const *const argv)
{
  args::ArgumentParser parser("NEdyson: Non-Equilibrium Green's function solver for molecules");
  args::Positional<std::string> inifile(parser, "inifile", "The parameter file");

  // Argument list
  args::ValueFlag<std::string> mode(parser, "mode", "GF2 or GW", {"mode"}, "GF2");
  args::Flag tti(parser, "tti", "use time translationally invariant functions if set", {"tti"});
  args::Flag unrestricted(parser, "unrestricted", "use unrestricted spin if set", {"unrestricted"});
  args::Flag decomposed(parser, "decomposed", "use unrestricted spin if set", {"decomposed"});
  args::ValueFlag<std::string> repr(parser, "repr", "type of representation: cheb or ir", {"repr"}, "cheb");
  args::ValueFlag<std::string> repr_file(parser, "repr-file", "path to representation file", {"repr-file"});
  args::ValueFlag<double> decomp_prec(parser, "decomp-prec", "decomposition precision", {"decomp-prec"}, 1e-12);
  args::ValueFlag<double> etol(parser, "etol", "energy convergence tolarence", {"etol"}, 1e-8);
  args::ValueFlag<int> maxiter(parser, "maxiter", "maximum number of iterations", {"maxiter"}, 50);
  args::ValueFlag<double> beta(parser, "beta", "inverse temperature", {"beta"}, 100);
  args::ValueFlag<double> damping(parser, "damping", "damping factor", {"damping"});
  args::ValueFlag<std::string> hf_input(parser, "hf-input", "input hdf5 file", {"hf-input"});
  args::ValueFlag<std::string> output(parser, "output", "path to output HDF5 file", {"output"}, "gfmol-output.h5");
  args::Flag checkpoint(parser, "checkpoint", "use checkpoint file", {"checkpoint"});
  args::ValueFlag<std::string> checkpoint_file(parser,
                                               "checkpoint-file",
                                               "path to checkpoint HDF5 file",
                                               {"checkpoint-file"},
                                               "gfmol-chk.h5");
  args::Flag boolOutput(parser, "boolOutput", "output results", {"boolOutput"});
  args::ValueFlag<int> nt(parser, "nt", "number of real time points", {"nt"});
  args::ValueFlag<int> ntau(parser, "ntau", "number of imaginary time points", {"ntau"});
  args::ValueFlag<int> k(parser, "k", "calculus accuracy order", {"k"});
  args::ValueFlag<double> dt(parser, "dt", "real timestep length", {"dt"});
  args::ValueFlag<int> BootMaxIter(parser, "BootMaxIter", "Maximum number of bootstrap iterations", {"BootMaxIter"});
  args::ValueFlag<double> BootMaxErr(parser, "BootMaxErr", "Maximum bootstrap error", {"BootMaxErr"});
  args::ValueFlag<int> CorrSteps(parser, "CorrSteps", "Number of predictor-corrector iterations in timestepping", {"CorrSteps"});
  args::ValueFlag<int> nw(parser, "nw", "Number of frequencies sampled in spectral function", {"nw"});
  args::ValueFlag<double> wmax(parser, "wmax", "absolute value of the bounds of the sampled frequency space", {"wmax"});

  args::HelpFlag help(parser, "help", "Display this help menu", {'h', "help"});
  try {
    parser.ParseCLI(argc, argv);
  } catch (const args::Help &) {
    std::cout << parser;
    exit(0);
  }

  INI::File param_file;
  param_file.Load(inifile.Get(), true);

  Params p{extract_value(param_file, mode),
           extract_value(param_file, tti),
           extract_value(param_file, unrestricted),
           extract_value(param_file, decomposed),
           extract_value(param_file, repr),
           extract_value(param_file, repr_file),
           extract_value(param_file, decomp_prec),
           extract_value(param_file, etol),
           extract_value(param_file, maxiter),
           extract_value(param_file, beta),
           extract_value(param_file, damping),
           extract_value(param_file, hf_input),
           extract_value(param_file, output),
           extract_value(param_file, checkpoint),
           extract_value(param_file, checkpoint_file),
           extract_value(param_file, boolOutput),
           extract_value(param_file, nt),
           extract_value(param_file, ntau),
           extract_value(param_file, k),
           extract_value(param_file, dt),
           extract_value(param_file, BootMaxIter),
           extract_value(param_file, BootMaxErr),
           extract_value(param_file, CorrSteps),
           extract_value(param_file, nw),
           extract_value(param_file, wmax) };
  try {
    p.validate();
  } catch (const std::runtime_error &e) {
    std::cerr << "Error encountered! " << e.what() << std::endl;
    std::cout << parser;
    throw e;
  }
  return p;
}

} // namespace NEdyson

#endif //_PARAMS_H_
