#ifndef OMP_INTEGRALS_DECL
#define OMP_INTEGRALS_DECL

#include <vector>
#include "greens.h"
#include "integration.h"
#include "utils.h"

namespace NEdyson{

void incr_convolution_ret(int tstp, const std::vector<bool> &mask_ret, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I , double dt);

void incr_convolution_tv(int tstp, const std::vector<bool> &mask_tv, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I, double beta, double dt);

void incr_convolution_les(int tstp, const std::vector<bool> &mask_les, GREEN &C, const GREEN &A, const GREEN &Acc, const GREEN &B, const GREEN &Bcc, const INTEG &I, double beta, double dt);

}//namespace
#endif
