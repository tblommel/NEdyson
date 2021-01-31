#ifndef UTILS_DECL
#define UTILS_DECL

#define PI 3.1415926535897932384626433832795028841971693

#include <Eigen/Eigen>
#include <alps/numeric/tensors.hpp>
#include <highfive/H5Easy.hpp>

#define INTEG Integration::Integrator

#define GREEN NEdyson::green_func
#define TTI_GREEN NEdyson::tti_green_func

#define DYSON NEdyson::dyson
#define TTI_DYSON NEdyson::tti_dyson


namespace h5 = HighFive;
namespace h5e= H5Easy;

namespace NEdyson {

using cplx = std::complex<double>;

// Tensor Aliases
template <typename T, size_t N, class C>
using TensorBase = alps::numerics::detail::tensor_base<T, N, C>;

template <typename T, size_t N>
using Tensor = alps::numerics::tensor<T, N>;

template <typename T, size_t N>
using TensorView = alps::numerics::tensor_view<T, N>;

template <size_t N>
using DTensor = alps::numerics::tensor<double, N>;

template <size_t N>
using ZTensor = alps::numerics::tensor<std::complex<double>, N>;

template <size_t N>
using DTensorView = alps::numerics::tensor_view<double, N>;

template <size_t N>
using ZTensorView = alps::numerics::tensor_view<std::complex<double>, N>;

// Matrix Aliases
template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
template <typename T>
using MatrixMap = Eigen::Map<Matrix<T>>;
template <typename T>
using MatrixConstMap = Eigen::Map<const Matrix<T>>;

template <typename T>
using RowVector = Eigen::Matrix<T, 1, Eigen::Dynamic, Eigen::RowMajor>;
template <typename T>
using RowVectorMap = Eigen::Map<RowVector<T>>;
template <typename T>
using RowVectorConstMap = Eigen::Map<const RowVector<T>>;

template <typename T>
using ColVector = Eigen::Matrix<T, Eigen::Dynamic, 1, Eigen::ColMajor>;
template <typename T>
using ColVectorMap = Eigen::Map<ColVector<T>>;
template <typename T>
using ColVectorConstMap = Eigen::Map<const ColVector<T>>;


using DMatrix = Matrix<double>;
using ZMatrix = Matrix<cplx>;
using DColVector = ColVector<double>;
using ZColVector = ColVector<std::complex<double>>;
using DRowVector = RowVector<double>;
using ZRowVector = RowVector<std::complex<double>>;

using DMatrixMap = MatrixMap<double>;
using ZMatrixMap = MatrixMap<std::complex<double>>;
using ZFMatrixMap = MatrixMap<std::complex<float>>;
using ZLMatrixMap = MatrixMap<std::complex<long double>>;
using DColVectorMap = ColVectorMap<double>;
using ZColVectorMap = ColVectorMap<std::complex<double>>;
using ZFColVectorMap = ColVectorMap<std::complex<float>>;
using ZLColVectorMap = ColVectorMap<std::complex<long double>>;
using DRowVectorMap = RowVectorMap<double>;
using ZRowVectorMap = RowVectorMap<std::complex<double>>;
using ZFRowVectorMap = RowVectorMap<std::complex<float>>;
using ZLRowVectorMap = RowVectorMap<std::complex<long double>>;

using DMatrixConstMap = MatrixConstMap<double>;
using FMatrixConstMap = MatrixConstMap<float>;
using LMatrixConstMap = MatrixConstMap<long double>;
using ZMatrixConstMap = MatrixConstMap<std::complex<double>>;
using DColVectorConstMap = ColVectorConstMap<double>;
using ZColVectorConstMap = ColVectorConstMap<std::complex<double>>;
using DRowVectorConstMap = RowVectorConstMap<double>;
using ZRowVectorConstMap = RowVectorConstMap<std::complex<double>>;

using ZMatrixMap2 = Eigen::Map<Eigen::Matrix<cplx,2,2,Eigen::RowMajor> >;
using ZMatrixMap3 = Eigen::Map<Eigen::Matrix<cplx,3,3,Eigen::RowMajor> >;
using ZMatrixMap4 = Eigen::Map<Eigen::Matrix<cplx,4,4,Eigen::RowMajor> >;



} // namespace NEdyson

namespace Integration{
 
using cplx = std::complex<double>;

template <typename T>
using Matrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
  using ZMatrix = Matrix<cplx>;

} // namespace Integration

#endif // UTILS_DECL
