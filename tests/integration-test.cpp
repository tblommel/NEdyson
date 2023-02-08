#include "tests-def.h"
#include "integration.h"
#include "utils.h"

TEST_CASE("Polynomial Calculus Coefficients")
{
  INTEG I(5);

  h5e::File file(std::string(TEST_DATA_DIR)+"/calccoeff.h5");

  Eigen::VectorXd PolyInterp = h5e::load<Eigen::VectorXd>(file,"/PInterp");
  Eigen::VectorXd PolyDiff = h5e::load<Eigen::VectorXd>(file,"/PDiff");
  Eigen::VectorXd PolyInteg = h5e::load<Eigen::VectorXd>(file,"/PInteg");
  Eigen::VectorXd BackDiff = h5e::load<Eigen::VectorXd>(file,"/BD");
  Eigen::VectorXd Start = h5e::load<Eigen::VectorXd>(file,"/Start");
  Eigen::VectorXd Sigma = h5e::load<Eigen::VectorXd>(file,"/Sigma");
  Eigen::VectorXd Rcorr = h5e::load<Eigen::VectorXd>(file,"/R");
  Eigen::VectorXd omega = h5e::load<Eigen::VectorXd>(file,"/w");

  for(int i=0;i<=5;i++){
    for(int j=0;j<=5;j++){
      PolyInterp(i*6+j)-=I.poly_interp(i,j);
    }
  }
  REQUIRE(PolyInterp.norm()<1e-13);
  
  for(int i=0;i<=5;i++){
    for(int j=0;j<=5;j++){
      PolyDiff(i*6+j)-=I.poly_diff(i,j);
    }
  }
  REQUIRE(PolyDiff.norm()<1e-12);
  
  for(int i=0;i<=5;i++){
    for(int j=0;j<=5;j++){
      for(int k=0;k<=5;k++){
        PolyInteg(i*36+j*6+k)-=I.poly_integ(i,j,k);
      }
    }
  }
  REQUIRE(PolyInteg.norm()<1e-12);

  for(int i=0;i<=6;i++){
    BackDiff(i)-=I.bd_weights(i);
  }
  REQUIRE(BackDiff.norm()<1e-12);

  for(int i=0;i<=5;i++){
    for(int j=0;j<=5;j++){
      Start(i*6+j)-=I.gregory_weights(i,j);
    }
  }
  REQUIRE(Start.norm()<1e-12);

  int count=0;
  for(int i=6;i<=11;i++){
    for(int j=0;j<=i;j++){
      Sigma(count)-=I.gregory_weights(i,j);
      count++;
    }
  }
  REQUIRE(Sigma.norm()<1e-12);

  for(int i=0;i<=5;i++){
    omega(i)-=I.omega(i);
  }
  REQUIRE(omega.norm()<1e-12);

}
