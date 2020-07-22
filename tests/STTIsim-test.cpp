#include "params.h"
#include "molNEsim.h"
#include "gfmol/repn.h"
#include "gfmol/utils.h"
#include "tests-def.h"

using namespace NEdyson;

TEST_CASE("STTIsim Tests") {
  int nt = 100, ntau = 100, k = 5, nw = 401, BootMaxIter = 70, MatMaxIter = 10, CorrSteps = 5;
  gfmol::Mode mode = gfmol::Mode::GF2;
  double dt = 0.01, wmax = 2, MatMaxErr = 1e-8, BootMaxErr = 1e-10, beta = 50;
  int es = 4;

  std::streambuf* cout_sbuf = std::cout.rdbuf(); // save original sbuf
  std::ofstream   stream("/dev/null");
  std::cout.rdbuf(stream.rdbuf()); // redirect 'cout' to a 'fout'

  // Read in actual values
  h5e::File check_file(std::string(TEST_DATA_DIR) + "/STTIsim.h5");
  TTI_GREEN Gdowncheck;
  TTI_GREEN Gupcheck;
  TTI_SPECT Acheck;
  Gupcheck.read_from_file(check_file, "/Gup");
  Acheck.read_from_file(check_file, "/Aup");
  Gdowncheck.read_from_file(check_file, "/Gdown");

  // Read in hf
  h5e::File input_f(std::string(TEST_DATA_DIR) + "/hf-spin.hdf5");
  gfmol::HartreeFock hf(input_f);

  // Read in repn
  h5e::File repr_f(std::string(TEST_DATA_DIR) + "/repn.hdf5");
  gfmol::ChebyshevRepr frepr(repr_f, "/fermi", gfmol::Stats::Fermi, beta);
  gfmol::ChebyshevRepr brepr(repr_f, "/bose", gfmol::Stats::Bose, beta);
  
  // Make the sim
  std::unique_ptr<tti_SpinSimulation<gfmol::ChebyshevRepr>> p_sim = std::unique_ptr<tti_SpinSimulation<gfmol::ChebyshevRepr>>(new tti_SpinSimulation<gfmol::ChebyshevRepr>(hf, frepr, brepr, nt, ntau, k, dt, nw, wmax, MatMaxIter, MatMaxErr, BootMaxIter, BootMaxErr, CorrSteps, mode, 0));

  // Run the sim
  p_sim->run();


  // Check the results
  ZColVectorMap GuMVec = ZColVectorMap(p_sim->G[0].get().matptr(0), (ntau+1)*es);
  ZColVectorMap GuMcVec = ZColVectorMap(Gupcheck.matptr(0), (ntau+1)*es);
  REQUIRE((GuMVec-GuMcVec).norm() < 1e-12);

  ZColVectorMap GuRVec = ZColVectorMap(p_sim->G[0].get().retptr(0), (nt+1)*es);
  ZColVectorMap GuRcVec = ZColVectorMap(Gupcheck.retptr(0), (nt+1)*es);
  REQUIRE((GuRVec-GuRcVec).norm() < 1e-12);

  ZColVectorMap GuLVec = ZColVectorMap(p_sim->G[0].get().lesptr(0), (nt+1)*es);
  ZColVectorMap GuLcVec = ZColVectorMap(Gupcheck.lesptr(0), (nt+1)*es);
  REQUIRE((GuLVec-GuLcVec).norm() < 1e-12);

  ZColVectorMap GuTVVec = ZColVectorMap(p_sim->G[0].get().tvptr(0,0), (nt+1)*(ntau+1)*es);
  ZColVectorMap GuTVcVec = ZColVectorMap(Gupcheck.tvptr(0,0), (nt+1)*(ntau+1)*es);
  REQUIRE((GuTVVec-GuTVcVec).norm() < 1e-12);
  
  DColVectorMap AuVec = DColVectorMap(p_sim->A.ptr(0), Acheck.nw()*es);
  DColVectorMap AucVec = DColVectorMap(Acheck.ptr(0), Acheck.nw()*es);
  REQUIRE((AuVec-AucVec).norm() < 1e-12);
 
  // Check the results
  ZColVectorMap GdMVec = ZColVectorMap(p_sim->G[1].get().matptr(0), (ntau+1)*es);
  ZColVectorMap GdMcVec = ZColVectorMap(Gdowncheck.matptr(0), (ntau+1)*es);
  REQUIRE((GdMVec-GdMcVec).norm() < 1e-12);

  ZColVectorMap GdRVec = ZColVectorMap(p_sim->G[1].get().retptr(0), (nt+1)*es);
  ZColVectorMap GdRcVec = ZColVectorMap(Gdowncheck.retptr(0), (nt+1)*es);
  REQUIRE((GdRVec-GdRcVec).norm() < 1e-12);

  ZColVectorMap GdLVec = ZColVectorMap(p_sim->G[1].get().lesptr(0), (nt+1)*es);
  ZColVectorMap GdLcVec = ZColVectorMap(Gdowncheck.lesptr(0), (nt+1)*es);
  REQUIRE((GdLVec-GdLcVec).norm() < 1e-12);

  ZColVectorMap GdTVVec = ZColVectorMap(p_sim->G[1].get().tvptr(0,0), (nt+1)*(ntau+1)*es);
  ZColVectorMap GdTVcVec = ZColVectorMap(Gdowncheck.tvptr(0,0), (nt+1)*(ntau+1)*es);
  REQUIRE((GdTVVec-GdTVcVec).norm() < 1e-12);
  
  std::cout.rdbuf(cout_sbuf); // restore the original stream buffer
}
