#ifndef UTILS_DECL
#define UTILS_DECL

#define PI 3.1415926535897932384626433832795028841971693

#include <Eigen/Eigen>

#define INTEG Integration::Integrator
#define CFUNC NEdyson::function
#define GREEN NEdyson::green_func
#define GREEN_TSTP NEdyson::green_func_tstp
#define SPECT NEdyson::spectral

namespace NEdyson {

using cplx = std::complex<double>;

using dvector = Eigen::VectorXd;
using cdmatrix = Eigen::MatrixXcd;
using cdmmap = Eigen::Map<Eigen::Matrix<cplx,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >;

using cdmmap2 = Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >;
using cdmmap3 = Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >;
using cdmmap4 = Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >;

using ccdmmap = Eigen::Map<const Eigen::Matrix<cplx,Eigen::Dynamic,Eigen::Dynamic,Eigen::RowMajor> >;

using ccdmmap2 = Eigen::Map<const Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >;
using ccdmmap3 = Eigen::Map<const Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >;
using ccdmmap4 = Eigen::Map<const Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >;
} // namespace NEdyson

namespace Hubb{
 
  using cplx = std::complex<double>;

  using cdmatrix = Eigen::MatrixXcd; 

} // namespace Hubb

namespace Integration{
 
  using cplx = std::complex<double>;

  using cdmatrix = Eigen::MatrixXcd; 

} // namespace Integration

#endif // UTILS_DECL
