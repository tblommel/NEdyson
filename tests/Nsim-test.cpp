#include "params.h"
#include "molNEsim.h"
#include "gfmol/repn.h"
#include "gfmol/utils.h"
#include "tests-def.h"

using namespace NEdyson;

TEST_CASE("Nsim Tests") {
  int nt = 100, ntau = 100, k = 5, nw = 401, BootMaxIter = 70, MatMaxIter = 10, CorrSteps = 5;
  gfmol::Mode mode = gfmol::Mode::GF2;
  double dt = 0.01, wmax = 2, MatMaxErr = 1e-8, BootMaxErr = 1e-10, beta = 50;

  std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
  std::ofstream   stream("/dev/null");
  std::cout.rdbuf(stream.rdbuf()); // redirect 'cout' to a 'fout'

  // Read in actual values
  h5e::File check_file(std::string(TEST_DATA_DIR) + "/Nsim.h5");
  GREEN Gcheck;
  SPECT Acheck;
  Gcheck.read_from_file(check_file, "/G");
  Acheck.read_from_file(check_file, "/A");

  // Read in hf
  h5e::File input_f(std::string(TEST_DATA_DIR) + "/hf.hdf5");
  gfmol::HartreeFock hf(input_f);

  // Read in repn
  h5e::File repr_f(std::string(TEST_DATA_DIR) + "/repn.hdf5");
  gfmol::ChebyshevRepr frepr(repr_f, "/fermi", gfmol::Stats::Fermi, beta);
  gfmol::ChebyshevRepr brepr(repr_f, "/bose", gfmol::Stats::Bose, beta);
  
  // Make the sim
  std::unique_ptr<Simulation<gfmol::ChebyshevRepr>> p_sim = std::unique_ptr<Simulation<gfmol::ChebyshevRepr>>(new Simulation<gfmol::ChebyshevRepr>(hf, frepr, brepr, nt, ntau, k, dt, nw, wmax, MatMaxIter, MatMaxErr, BootMaxIter, BootMaxErr, CorrSteps, mode, 0));

  // Run the sim
  p_sim->run();


  // Check the results
  ZColVectorMap GMVec = ZColVectorMap(p_sim->G.matptr(0), (ntau+1)*Gcheck.element_size());
  ZColVectorMap GMcVec = ZColVectorMap(Gcheck.matptr(0), (ntau+1)*Gcheck.element_size());
  REQUIRE((GMVec-GMcVec).norm() < 1e-12);

  ZColVectorMap GRVec = ZColVectorMap(p_sim->G.retptr(0,0), (nt+1)*(nt+2)/2*Gcheck.element_size());
  ZColVectorMap GRcVec = ZColVectorMap(Gcheck.retptr(0,0), (nt+1)*(nt+2)/2*Gcheck.element_size());
  REQUIRE((GRVec-GRcVec).norm() < 1e-12);

  ZColVectorMap GLVec = ZColVectorMap(p_sim->G.lesptr(0,0), (nt+1)*(nt+2)/2*Gcheck.element_size());
  ZColVectorMap GLcVec = ZColVectorMap(Gcheck.lesptr(0,0), (nt+1)*(nt+2)/2*Gcheck.element_size());
  REQUIRE((GLVec-GLcVec).norm() < 1e-12);

  ZColVectorMap GTVVec = ZColVectorMap(p_sim->G.tvptr(0,0), (nt+1)*(ntau+1)*Gcheck.element_size());
  ZColVectorMap GTVcVec = ZColVectorMap(Gcheck.tvptr(0,0), (nt+1)*(ntau+1)*Gcheck.element_size());
  REQUIRE((GTVVec-GTVcVec).norm() < 1e-12);
  
  DColVectorMap AVec = DColVectorMap(p_sim->A.ptr(0,0), (Acheck.nt()+1)*Acheck.nw()*Acheck.es());
  DColVectorMap AcVec = DColVectorMap(Acheck.ptr(0,0), (Acheck.nt()+1)*Acheck.nw()*Acheck.es());
  REQUIRE((AVec-AcVec).norm() < 1e-12);
 
  std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
}
