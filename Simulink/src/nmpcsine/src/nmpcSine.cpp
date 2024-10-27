//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
//
// File: nmpcSine.cpp
//
// Code generated for Simulink model 'nmpcSine'.
//
// Model version                  : 1.102
// Simulink Coder version         : 24.2 (R2024b) 21-Jun-2024
// C/C++ source code generated on : Sun Oct 27 09:42:42 2024
//
// Target selection: ert.tlc
// Embedded hardware selection: Intel->x86-64 (Linux 64)
// Code generation objectives: Unspecified
// Validation result: Not run
//
#include "nmpcSine.h"
#include "nmpcSine_types.h"
#include "rtwtypes.h"

extern "C"
{

#include "rt_nonfinite.h"

}

#include <string.h>
#include <math.h>
#include <emmintrin.h>
#include "rmw/qos_profiles.h"
#include <stddef.h>
#include "nmpcSine_private.h"
#include <float.h>

// Named constants for Chart: '<Root>/Control Logic'
const uint8_T nmpcSine_IN_Control = 1U;
const uint8_T nmpcSine_IN_Stop = 2U;
const uint8_T nmpcSine_IN_Wait = 3U;
static void rate_scheduler(RT_MODEL_nmpcSine_T *const nmpcSine_M);
int32_T div_nde_s32_floor(int32_T numerator, int32_T denominator)
{
  return (((numerator < 0) != (denominator < 0)) && (numerator % denominator !=
           0) ? -1 : 0) + numerator / denominator;
}

//
//         This function updates active task flag for each subrate.
//         The function is called at model base rate, hence the
//         generated code self-manages all its subrates.
//
static void rate_scheduler(RT_MODEL_nmpcSine_T *const nmpcSine_M)
{
  // Compute which subrates run during the next base time step.  Subrates
  //  are an integer multiple of the base rate counter.  Therefore, the subtask
  //  counter is reset when it reaches its limit (zero means run).

  (nmpcSine_M->Timing.TaskCounters.TID[1])++;
  if ((nmpcSine_M->Timing.TaskCounters.TID[1]) > 3) {// Sample time: [0.4s, 0.0s] 
    nmpcSine_M->Timing.TaskCounters.TID[1] = 0;
  }

  (nmpcSine_M->Timing.TaskCounters.TID[2])++;
  if ((nmpcSine_M->Timing.TaskCounters.TID[2]) > 9) {// Sample time: [1.0s, 0.0s] 
    nmpcSine_M->Timing.TaskCounters.TID[2] = 0;
  }
}

real_T rt_remd_snf(real_T u0, real_T u1)
{
  real_T y;
  if (rtIsNaN(u0) || rtIsNaN(u1) || rtIsInf(u0)) {
    y = (rtNaN);
  } else if (rtIsInf(u1)) {
    y = u0;
  } else {
    real_T q;
    if (u1 < 0.0) {
      q = ceil(u1);
    } else {
      q = floor(u1);
    }

    if ((u1 != 0.0) && (u1 != q)) {
      q = fabs(u0 / u1);
      if (!(fabs(q - floor(q + 0.5)) > DBL_EPSILON * q)) {
        y = 0.0 * u0;
      } else {
        y = fmod(u0, u1);
      }
    } else {
      y = fmod(u0, u1);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S3>/LatLonToXY'
void nmpcSine::nmpcSine_sind(real_T *x)
{
  real_T absx;
  real_T b_x;
  int8_T n;
  if (rtIsInf(*x) || rtIsNaN(*x)) {
    *x = (rtNaN);
  } else {
    b_x = rt_remd_snf(*x, 360.0);
    absx = fabs(b_x);
    if (absx > 180.0) {
      if (b_x > 0.0) {
        b_x -= 360.0;
      } else {
        b_x += 360.0;
      }

      absx = fabs(b_x);
    }

    if (absx <= 45.0) {
      b_x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (b_x > 0.0) {
        b_x = (b_x - 90.0) * 0.017453292519943295;
        n = 1;
      } else {
        b_x = (b_x + 90.0) * 0.017453292519943295;
        n = -1;
      }
    } else if (b_x > 0.0) {
      b_x = (b_x - 180.0) * 0.017453292519943295;
      n = 2;
    } else {
      b_x = (b_x + 180.0) * 0.017453292519943295;
      n = -2;
    }

    switch (n) {
     case 0:
      *x = sin(b_x);
      break;

     case 1:
      *x = cos(b_x);
      break;

     case -1:
      *x = -cos(b_x);
      break;

     default:
      *x = -sin(b_x);
      break;
    }
  }
}

// Function for MATLAB Function: '<S3>/LatLonToXY'
void nmpcSine::nmpcSine_cosd(real_T *x)
{
  real_T absx;
  real_T b_x;
  int8_T n;
  if (rtIsInf(*x) || rtIsNaN(*x)) {
    *x = (rtNaN);
  } else {
    b_x = rt_remd_snf(*x, 360.0);
    absx = fabs(b_x);
    if (absx > 180.0) {
      if (b_x > 0.0) {
        b_x -= 360.0;
      } else {
        b_x += 360.0;
      }

      absx = fabs(b_x);
    }

    if (absx <= 45.0) {
      b_x *= 0.017453292519943295;
      n = 0;
    } else if (absx <= 135.0) {
      if (b_x > 0.0) {
        b_x = (b_x - 90.0) * 0.017453292519943295;
        n = 1;
      } else {
        b_x = (b_x + 90.0) * 0.017453292519943295;
        n = -1;
      }
    } else if (b_x > 0.0) {
      b_x = (b_x - 180.0) * 0.017453292519943295;
      n = 2;
    } else {
      b_x = (b_x + 180.0) * 0.017453292519943295;
      n = -2;
    }

    switch (n) {
     case 0:
      *x = cos(b_x);
      break;

     case 1:
      *x = -sin(b_x);
      break;

     case -1:
      *x = sin(b_x);
      break;

     default:
      *x = -cos(b_x);
      break;
    }
  }
}

real_T nmpcSine::nmpcSine_xnrm2(int32_T n, const real_T x[16], int32_T ix0)
{
  real_T y;

  // Start for MATLABSystem: '<S18>/MATLAB System'
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (int32_T k = ix0; k < kend; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  // End of Start for MATLABSystem: '<S18>/MATLAB System'
  return y;
}

real_T nmpcSine::nmpcSine_rt_hypotd_snf(real_T u0, real_T u1)
{
  real_T y;
  nmpcSine_B.a_d = fabs(u0);
  nmpcSine_B.b_g4 = fabs(u1);
  if (nmpcSine_B.a_d < nmpcSine_B.b_g4) {
    nmpcSine_B.a_d /= nmpcSine_B.b_g4;
    y = sqrt(nmpcSine_B.a_d * nmpcSine_B.a_d + 1.0) * nmpcSine_B.b_g4;
  } else if (nmpcSine_B.a_d > nmpcSine_B.b_g4) {
    nmpcSine_B.b_g4 /= nmpcSine_B.a_d;
    y = sqrt(nmpcSine_B.b_g4 * nmpcSine_B.b_g4 + 1.0) * nmpcSine_B.a_d;
  } else if (rtIsNaN(nmpcSine_B.b_g4)) {
    y = (rtNaN);
  } else {
    y = nmpcSine_B.a_d * 1.4142135623730951;
  }

  return y;
}

void nmpcSine::nmpcSine_qr(const real_T A[16], real_T Q[16], real_T R[4])
{
  __m128d tmp;
  real_T b_tau[2];
  real_T work[2];
  real_T atmp;
  real_T beta1;
  real_T c;
  int32_T d;
  int32_T exitg1;
  int32_T i;
  int32_T iac;
  int32_T iaii;
  int32_T ii;
  int32_T jA;
  int32_T knt;
  int32_T lastc;
  int32_T lastv;
  int32_T scalarLB;
  boolean_T exitg2;

  // Start for MATLABSystem: '<S18>/MATLAB System'
  b_tau[0] = 0.0;
  b_tau[1] = 0.0;
  memcpy(&Q[0], &A[0], sizeof(real_T) << 4U);
  work[0] = 0.0;
  work[1] = 0.0;
  for (iaii = 0; iaii < 2; iaii++) {
    ii = (iaii << 3) + iaii;
    i = ii + 2;
    atmp = Q[ii];
    b_tau[iaii] = 0.0;
    beta1 = nmpcSine_xnrm2(7 - iaii, Q, ii + 2);
    if (beta1 != 0.0) {
      c = Q[ii];
      beta1 = nmpcSine_rt_hypotd_snf(c, beta1);
      if (c >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          scalarLB = ii - iaii;
          lastc = (((((scalarLB - ii) + 7) / 2) << 1) + ii) + 2;
          jA = lastc - 2;
          for (lastv = i; lastv <= jA; lastv += 2) {
            tmp = _mm_loadu_pd(&Q[lastv - 1]);
            _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (lastv = lastc; lastv <= scalarLB + 8; lastv++) {
            Q[lastv - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          atmp *= 9.9792015476736E+291;
        } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt + 1 < 20));

        beta1 = nmpcSine_rt_hypotd_snf(atmp, nmpcSine_xnrm2(7 - iaii, Q, ii + 2));
        if (atmp >= 0.0) {
          beta1 = -beta1;
        }

        b_tau[iaii] = (beta1 - atmp) / beta1;
        atmp = 1.0 / (atmp - beta1);
        for (lastv = i; lastv <= jA; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = lastc; lastv <= scalarLB + 8; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        for (lastv = 0; lastv <= knt; lastv++) {
          beta1 *= 1.0020841800044864E-292;
        }

        atmp = beta1;
      } else {
        b_tau[iaii] = (beta1 - c) / beta1;
        atmp = 1.0 / (c - beta1);
        knt = ii - iaii;
        scalarLB = (((((knt - ii) + 7) / 2) << 1) + ii) + 2;
        lastc = scalarLB - 2;
        for (lastv = i; lastv <= lastc; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = scalarLB; lastv <= knt + 8; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        atmp = beta1;
      }
    }

    Q[ii] = atmp;
    if (iaii + 1 < 2) {
      Q[ii] = 1.0;
      scalarLB = ii + 9;
      if (b_tau[iaii] != 0.0) {
        lastv = 8 - iaii;
        i = ii - iaii;
        while ((lastv > 0) && (Q[i + 7] == 0.0)) {
          lastv--;
          i--;
        }

        knt = 1 - iaii;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          jA = ii + 9;
          do {
            exitg1 = 0;
            if (jA <= (ii + lastv) + 8) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt = 0;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            work[0] = 0.0;
          }

          knt = ((lastc << 3) + ii) + 9;
          for (iac = scalarLB; iac <= knt; iac += 8) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[(ii + jA) - iac] * Q[jA - 1];
            }

            jA = ((iac - ii) - 9) >> 3;
            work[jA] += c;
          }
        }

        if (!(-b_tau[iaii] == 0.0)) {
          jA = ii;
          for (iac = 0; iac <= lastc; iac++) {
            if (work[0] != 0.0) {
              c = work[0] * -b_tau[iaii];
              knt = jA + 9;
              scalarLB = (lastv + jA) + 8;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((ii + d) - jA) - 9] * c;
              }
            }

            jA += 8;
          }
        }
      }

      Q[ii] = atmp;
    }
  }

  for (ii = 0; ii < 2; ii++) {
    for (iaii = 0; iaii <= ii; iaii++) {
      // Start for MATLABSystem: '<S18>/MATLAB System'
      R[iaii + (ii << 1)] = Q[(ii << 3) + iaii];
    }

    if (ii + 2 <= 2) {
      R[(ii << 1) + 1] = 0.0;
    }

    // Start for MATLABSystem: '<S18>/MATLAB System'
    work[ii] = 0.0;
  }

  // Start for MATLABSystem: '<S18>/MATLAB System'
  for (i = 1; i >= 0; i--) {
    iaii = ((i << 3) + i) + 8;
    if (i + 1 < 2) {
      Q[iaii - 8] = 1.0;
      scalarLB = iaii + 1;
      if (b_tau[i] != 0.0) {
        lastv = 8 - i;
        ii = (iaii - i) - 1;
        while ((lastv > 0) && (Q[ii] == 0.0)) {
          lastv--;
          ii--;
        }

        knt = 1 - i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          jA = iaii + 1;
          do {
            exitg1 = 0;
            if (jA <= iaii + lastv) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt = 0;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            work[0] = 0.0;
          }

          knt = ((lastc << 3) + iaii) + 1;
          for (iac = scalarLB; iac <= knt; iac += 8) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[((iaii + jA) - iac) - 8] * Q[jA - 1];
            }

            jA = ((iac - iaii) - 1) >> 3;
            work[jA] += c;
          }
        }

        if (!(-b_tau[i] == 0.0)) {
          jA = iaii;
          for (iac = 0; iac <= lastc; iac++) {
            if (work[0] != 0.0) {
              c = work[0] * -b_tau[i];
              knt = jA + 1;
              scalarLB = lastv + jA;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((iaii + d) - jA) - 9] * c;
              }
            }

            jA += 8;
          }
        }
      }
    }

    knt = iaii - i;
    scalarLB = (((((knt - iaii) + 7) / 2) << 1) + iaii) - 6;
    lastc = scalarLB - 2;
    for (lastv = iaii - 6; lastv <= lastc; lastv += 2) {
      tmp = _mm_loadu_pd(&Q[lastv - 1]);
      _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(-b_tau[i])));
    }

    for (lastv = scalarLB; lastv <= knt; lastv++) {
      Q[lastv - 1] *= -b_tau[i];
    }

    Q[iaii - 8] = 1.0 - b_tau[i];
    if (i - 1 >= 0) {
      Q[iaii - 9] = 0.0;
    }
  }
}

void nmpcSine::nmpcSine_trisolve(const real_T A[4], real_T B[12])
{
  // Start for MATLABSystem: '<S18>/MATLAB System'
  for (int32_T b_j = 0; b_j < 6; b_j++) {
    int32_T jBcol;
    jBcol = (b_j << 1) - 1;
    for (int32_T b_k = 0; b_k < 2; b_k++) {
      real_T B_0;
      int32_T B_tmp;
      int32_T k;
      int32_T kAcol;
      k = b_k + 1;
      kAcol = (b_k << 1) - 1;
      B_tmp = (b_k + jBcol) + 1;
      B_0 = B[B_tmp];
      if (B_0 != 0.0) {
        B[B_tmp] = B_0 / A[(b_k + kAcol) + 1];
        for (int32_T i = k + 1; i < 3; i++) {
          B[jBcol + 2] -= A[kAcol + 2] * B[B_tmp];
        }
      }
    }
  }

  // End of Start for MATLABSystem: '<S18>/MATLAB System'
}

void nmpcSine::nmpcSine_trisolve_m(const real_T A[4], real_T B[12])
{
  // Start for MATLABSystem: '<S18>/MATLAB System'
  for (int32_T b_j = 0; b_j < 6; b_j++) {
    int32_T jBcol;
    jBcol = b_j << 1;
    for (int32_T k = 1; k >= 0; k--) {
      real_T tmp;
      int32_T kAcol;
      int32_T tmp_0;
      kAcol = k << 1;
      tmp_0 = k + jBcol;
      tmp = B[tmp_0];
      if (tmp != 0.0) {
        B[tmp_0] = tmp / A[k + kAcol];
        for (int32_T b_i = 0; b_i < k; b_i++) {
          B[jBcol] -= B[tmp_0] * A[kAcol];
        }
      }
    }
  }

  // End of Start for MATLABSystem: '<S18>/MATLAB System'
}

real_T nmpcSine::nmpcSine_xnrm2_m(int32_T n, const real_T x[48], int32_T ix0)
{
  real_T y;

  // Start for MATLABSystem: '<S18>/MATLAB System'
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (int32_T k = ix0; k < kend; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  // End of Start for MATLABSystem: '<S18>/MATLAB System'
  return y;
}

void nmpcSine::nmpcSine_qr_m(const real_T A[48], real_T Q[48], real_T R[36])
{
  __m128d tmp;
  real_T b_tau[6];
  real_T work[6];
  real_T atmp;
  real_T beta1;
  real_T c;
  int32_T d;
  int32_T exitg1;
  int32_T i;
  int32_T iac;
  int32_T iaii;
  int32_T ii;
  int32_T jA;
  int32_T knt;
  int32_T lastc;
  int32_T lastv;
  int32_T scalarLB;
  boolean_T exitg2;
  for (i = 0; i < 6; i++) {
    // Start for MATLABSystem: '<S18>/MATLAB System'
    b_tau[i] = 0.0;
  }

  // Start for MATLABSystem: '<S18>/MATLAB System'
  memcpy(&Q[0], &A[0], 48U * sizeof(real_T));
  for (i = 0; i < 6; i++) {
    // Start for MATLABSystem: '<S18>/MATLAB System'
    work[i] = 0.0;
  }

  // Start for MATLABSystem: '<S18>/MATLAB System'
  for (iaii = 0; iaii < 6; iaii++) {
    ii = (iaii << 3) + iaii;
    i = ii + 2;
    atmp = Q[ii];
    b_tau[iaii] = 0.0;
    beta1 = nmpcSine_xnrm2_m(7 - iaii, Q, ii + 2);
    if (beta1 != 0.0) {
      c = Q[ii];
      beta1 = nmpcSine_rt_hypotd_snf(c, beta1);
      if (c >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          scalarLB = ii - iaii;
          lastc = (((((scalarLB - ii) + 7) / 2) << 1) + ii) + 2;
          jA = lastc - 2;
          for (lastv = i; lastv <= jA; lastv += 2) {
            tmp = _mm_loadu_pd(&Q[lastv - 1]);
            _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (lastv = lastc; lastv <= scalarLB + 8; lastv++) {
            Q[lastv - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          atmp *= 9.9792015476736E+291;
        } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt + 1 < 20));

        beta1 = nmpcSine_rt_hypotd_snf(atmp, nmpcSine_xnrm2_m(7 - iaii, Q, ii +
          2));
        if (atmp >= 0.0) {
          beta1 = -beta1;
        }

        b_tau[iaii] = (beta1 - atmp) / beta1;
        atmp = 1.0 / (atmp - beta1);
        for (lastv = i; lastv <= jA; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = lastc; lastv <= scalarLB + 8; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        for (lastv = 0; lastv <= knt; lastv++) {
          beta1 *= 1.0020841800044864E-292;
        }

        atmp = beta1;
      } else {
        b_tau[iaii] = (beta1 - c) / beta1;
        atmp = 1.0 / (c - beta1);
        knt = ii - iaii;
        scalarLB = (((((knt - ii) + 7) / 2) << 1) + ii) + 2;
        lastc = scalarLB - 2;
        for (lastv = i; lastv <= lastc; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = scalarLB; lastv <= knt + 8; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        atmp = beta1;
      }
    }

    Q[ii] = atmp;
    if (iaii + 1 < 6) {
      Q[ii] = 1.0;
      scalarLB = ii + 9;
      if (b_tau[iaii] != 0.0) {
        lastv = 8 - iaii;
        i = ii - iaii;
        while ((lastv > 0) && (Q[i + 7] == 0.0)) {
          lastv--;
          i--;
        }

        knt = 5 - iaii;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          lastc = ((knt - 1) << 3) + ii;
          jA = lastc + 9;
          do {
            exitg1 = 0;
            if (jA <= (lastc + lastv) + 8) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            memset(&work[0], 0, static_cast<uint32_T>(lastc + 1) * sizeof(real_T));
          }

          knt = ((lastc << 3) + ii) + 9;
          for (iac = scalarLB; iac <= knt; iac += 8) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[(ii + jA) - iac] * Q[jA - 1];
            }

            jA = ((iac - ii) - 9) >> 3;
            work[jA] += c;
          }
        }

        if (!(-b_tau[iaii] == 0.0)) {
          jA = ii;
          for (iac = 0; iac <= lastc; iac++) {
            c = work[iac];
            if (c != 0.0) {
              c *= -b_tau[iaii];
              knt = jA + 9;
              scalarLB = (lastv + jA) + 8;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((ii + d) - jA) - 9] * c;
              }
            }

            jA += 8;
          }
        }
      }

      Q[ii] = atmp;
    }
  }

  for (ii = 0; ii < 6; ii++) {
    for (iaii = 0; iaii <= ii; iaii++) {
      // Start for MATLABSystem: '<S18>/MATLAB System'
      R[iaii + 6 * ii] = Q[(ii << 3) + iaii];
    }

    for (iaii = ii + 2; iaii < 7; iaii++) {
      R[(iaii + 6 * ii) - 1] = 0.0;
    }

    // Start for MATLABSystem: '<S18>/MATLAB System'
    work[ii] = 0.0;
  }

  // Start for MATLABSystem: '<S18>/MATLAB System'
  for (i = 5; i >= 0; i--) {
    iaii = ((i << 3) + i) + 8;
    if (i + 1 < 6) {
      Q[iaii - 8] = 1.0;
      scalarLB = iaii + 1;
      if (b_tau[i] != 0.0) {
        lastv = 8 - i;
        ii = (iaii - i) - 1;
        while ((lastv > 0) && (Q[ii] == 0.0)) {
          lastv--;
          ii--;
        }

        knt = 5 - i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          lastc = ((knt - 1) << 3) + iaii;
          jA = lastc + 1;
          do {
            exitg1 = 0;
            if (jA <= lastc + lastv) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            memset(&work[0], 0, static_cast<uint32_T>(lastc + 1) * sizeof(real_T));
          }

          knt = ((lastc << 3) + iaii) + 1;
          for (iac = scalarLB; iac <= knt; iac += 8) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[((iaii + jA) - iac) - 8] * Q[jA - 1];
            }

            jA = ((iac - iaii) - 1) >> 3;
            work[jA] += c;
          }
        }

        if (!(-b_tau[i] == 0.0)) {
          jA = iaii;
          for (iac = 0; iac <= lastc; iac++) {
            c = work[iac];
            if (c != 0.0) {
              c *= -b_tau[i];
              knt = jA + 1;
              scalarLB = lastv + jA;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((iaii + d) - jA) - 9] * c;
              }
            }

            jA += 8;
          }
        }
      }
    }

    knt = iaii - i;
    scalarLB = (((((knt - iaii) + 7) / 2) << 1) + iaii) - 6;
    lastc = scalarLB - 2;
    for (lastv = iaii - 6; lastv <= lastc; lastv += 2) {
      tmp = _mm_loadu_pd(&Q[lastv - 1]);
      _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(-b_tau[i])));
    }

    for (lastv = scalarLB; lastv <= knt; lastv++) {
      Q[lastv - 1] *= -b_tau[i];
    }

    Q[iaii - 8] = 1.0 - b_tau[i];
    for (ii = 0; ii < i; ii++) {
      Q[(iaii - ii) - 9] = 0.0;
    }
  }
}

real_T nmpcSine::nmpcSine_mod(real_T x, real_T y)
{
  real_T r;

  // Start for MATLABSystem: '<S19>/MATLAB System'
  r = x;
  if (y == 0.0) {
    if (x == 0.0) {
      r = y;
    }
  } else if (rtIsNaN(x) || rtIsNaN(y) || rtIsInf(x)) {
    r = (rtNaN);
  } else if (x == 0.0) {
    r = 0.0 / y;
  } else {
    boolean_T rEQ0;
    r = fmod(x, y);
    rEQ0 = (r == 0.0);
    if ((!rEQ0) && (y > floor(y))) {
      real_T q;
      q = fabs(x / y);
      rEQ0 = !(fabs(q - floor(q + 0.5)) > 2.2204460492503131E-16 * q);
    }

    if (rEQ0) {
      r = y * 0.0;
    } else if (((r < 0.0) && (!(y < 0.0))) || ((!(r < 0.0)) && (y < 0.0))) {
      r += y;
    }
  }

  // End of Start for MATLABSystem: '<S19>/MATLAB System'
  return r;
}

void nmpcSine::nmpcSine_trisolve_m3(real_T A, real_T B[6])
{
  // Start for MATLABSystem: '<S19>/MATLAB System' incorporates:
  //   MATLABSystem: '<S20>/MATLAB System'

  for (int32_T b_j = 0; b_j < 6; b_j++) {
    real_T B_0;
    B_0 = B[b_j];
    if (B_0 != 0.0) {
      B[b_j] = B_0 / A;
    }
  }

  // End of Start for MATLABSystem: '<S19>/MATLAB System'
}

real_T nmpcSine::nmpcSine_xnrm2_m3(int32_T n, const real_T x[42], int32_T ix0)
{
  real_T y;
  y = 0.0;

  // Start for MATLABSystem: '<S19>/MATLAB System'
  if (n == 1) {
    y = fabs(x[ix0 - 1]);
  } else {
    int32_T kend;
    nmpcSine_B.scale_f = 3.3121686421112381E-170;
    kend = ix0 + n;
    for (int32_T k = ix0; k < kend; k++) {
      nmpcSine_B.absxk_h = fabs(x[k - 1]);
      if (nmpcSine_B.absxk_h > nmpcSine_B.scale_f) {
        nmpcSine_B.t_m = nmpcSine_B.scale_f / nmpcSine_B.absxk_h;
        y = y * nmpcSine_B.t_m * nmpcSine_B.t_m + 1.0;
        nmpcSine_B.scale_f = nmpcSine_B.absxk_h;
      } else {
        nmpcSine_B.t_m = nmpcSine_B.absxk_h / nmpcSine_B.scale_f;
        y += nmpcSine_B.t_m * nmpcSine_B.t_m;
      }
    }

    y = nmpcSine_B.scale_f * sqrt(y);
  }

  // End of Start for MATLABSystem: '<S19>/MATLAB System'
  return y;
}

void nmpcSine::EKFCorrector_correctStateAndSqr(const real_T x[6], const real_T
  S[36], real_T residue, const real_T Pxy[6], real_T Sy, const real_T H[6],
  real_T Rsqrt, real_T b_x[6], real_T b_S[36])
{
  __m128d tmp;
  int32_T exitg1;
  int32_T i;
  int32_T vectorUB_tmp;
  boolean_T exitg2;
  for (nmpcSine_B.ii = 0; nmpcSine_B.ii < 6; nmpcSine_B.ii++) {
    // Start for MATLABSystem: '<S19>/MATLAB System'
    nmpcSine_B.work[nmpcSine_B.ii] = Pxy[nmpcSine_B.ii];
  }

  // Start for MATLABSystem: '<S19>/MATLAB System'
  nmpcSine_trisolve_m3(Sy, nmpcSine_B.work);
  for (nmpcSine_B.ii = 0; nmpcSine_B.ii < 6; nmpcSine_B.ii++) {
    // Start for MATLABSystem: '<S19>/MATLAB System'
    nmpcSine_B.K_a[nmpcSine_B.ii] = nmpcSine_B.work[nmpcSine_B.ii];
  }

  // Start for MATLABSystem: '<S19>/MATLAB System'
  nmpcSine_trisolve_m3(Sy, nmpcSine_B.K_a);
  for (i = 0; i <= 4; i += 2) {
    tmp = _mm_loadu_pd(&nmpcSine_B.K_a[i]);
    _mm_storeu_pd(&b_x[i], _mm_add_pd(_mm_mul_pd(tmp, _mm_set1_pd(residue)),
      _mm_loadu_pd(&x[i])));
    _mm_storeu_pd(&nmpcSine_B.work[i], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
  }

  for (nmpcSine_B.b_i_e = 0; nmpcSine_B.b_i_e < 6; nmpcSine_B.b_i_e++) {
    for (nmpcSine_B.ii = 0; nmpcSine_B.ii <= 4; nmpcSine_B.ii += 2) {
      tmp = _mm_loadu_pd(&nmpcSine_B.work[nmpcSine_B.ii]);

      // Start for MATLABSystem: '<S19>/MATLAB System'
      _mm_storeu_pd(&nmpcSine_B.A_j[nmpcSine_B.ii + 6 * nmpcSine_B.b_i_e],
                    _mm_mul_pd(tmp, _mm_set1_pd(H[nmpcSine_B.b_i_e])));
    }
  }

  for (nmpcSine_B.b_i_e = 0; nmpcSine_B.b_i_e < 6; nmpcSine_B.b_i_e++) {
    // Start for MATLABSystem: '<S19>/MATLAB System'
    nmpcSine_B.ii = 6 * nmpcSine_B.b_i_e + nmpcSine_B.b_i_e;
    nmpcSine_B.A_j[nmpcSine_B.ii]++;
  }

  for (nmpcSine_B.ii = 0; nmpcSine_B.ii < 6; nmpcSine_B.ii++) {
    // Start for MATLABSystem: '<S19>/MATLAB System'
    nmpcSine_B.b_tau[nmpcSine_B.ii] = 0.0;
    for (nmpcSine_B.b_i_e = 0; nmpcSine_B.b_i_e < 6; nmpcSine_B.b_i_e++) {
      // Start for MATLABSystem: '<S19>/MATLAB System'
      nmpcSine_B.aoffset_b = nmpcSine_B.b_i_e * 6 - 1;
      nmpcSine_B.s = 0.0;
      for (i = 0; i < 6; i++) {
        // Start for MATLABSystem: '<S19>/MATLAB System'
        nmpcSine_B.s += S[(nmpcSine_B.aoffset_b + i) + 1] * nmpcSine_B.A_j[i * 6
          + nmpcSine_B.ii];
      }

      // Start for MATLABSystem: '<S19>/MATLAB System'
      nmpcSine_B.c_A[nmpcSine_B.b_i_e + 7 * nmpcSine_B.ii] = nmpcSine_B.s;
    }

    // Start for MATLABSystem: '<S19>/MATLAB System'
    nmpcSine_B.c_A[7 * nmpcSine_B.ii + 6] = nmpcSine_B.K_a[nmpcSine_B.ii] *
      Rsqrt;
    nmpcSine_B.work[nmpcSine_B.ii] = 0.0;
  }

  for (nmpcSine_B.b_i_e = 0; nmpcSine_B.b_i_e < 6; nmpcSine_B.b_i_e++) {
    // Start for MATLABSystem: '<S19>/MATLAB System'
    nmpcSine_B.ii = nmpcSine_B.b_i_e * 7 + nmpcSine_B.b_i_e;
    nmpcSine_B.aoffset_b = nmpcSine_B.ii + 2;
    nmpcSine_B.s = nmpcSine_B.c_A[nmpcSine_B.ii];
    nmpcSine_B.b_tau[nmpcSine_B.b_i_e] = 0.0;

    // Start for MATLABSystem: '<S19>/MATLAB System'
    nmpcSine_B.beta1 = nmpcSine_xnrm2_m3(6 - nmpcSine_B.b_i_e, nmpcSine_B.c_A,
      nmpcSine_B.ii + 2);
    if (nmpcSine_B.beta1 != 0.0) {
      // Start for MATLABSystem: '<S19>/MATLAB System'
      nmpcSine_B.c_cu = nmpcSine_B.c_A[nmpcSine_B.ii];
      nmpcSine_B.beta1 = nmpcSine_rt_hypotd_snf(nmpcSine_B.c_cu,
        nmpcSine_B.beta1);
      if (nmpcSine_B.c_cu >= 0.0) {
        nmpcSine_B.beta1 = -nmpcSine_B.beta1;
      }

      if (fabs(nmpcSine_B.beta1) < 1.0020841800044864E-292) {
        nmpcSine_B.knt = -1;
        do {
          nmpcSine_B.knt++;
          nmpcSine_B.lastc = nmpcSine_B.ii - nmpcSine_B.b_i_e;
          nmpcSine_B.jA = (((((nmpcSine_B.lastc - nmpcSine_B.ii) + 6) / 2) << 1)
                           + nmpcSine_B.ii) + 2;
          vectorUB_tmp = nmpcSine_B.jA - 2;
          for (i = nmpcSine_B.aoffset_b; i <= vectorUB_tmp; i += 2) {
            tmp = _mm_loadu_pd(&nmpcSine_B.c_A[i - 1]);
            _mm_storeu_pd(&nmpcSine_B.c_A[i - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (i = nmpcSine_B.jA; i <= nmpcSine_B.lastc + 7; i++) {
            nmpcSine_B.c_A[i - 1] *= 9.9792015476736E+291;
          }

          nmpcSine_B.beta1 *= 9.9792015476736E+291;
          nmpcSine_B.s *= 9.9792015476736E+291;
        } while ((fabs(nmpcSine_B.beta1) < 1.0020841800044864E-292) &&
                 (nmpcSine_B.knt + 1 < 20));

        nmpcSine_B.beta1 = nmpcSine_rt_hypotd_snf(nmpcSine_B.s,
          nmpcSine_xnrm2_m3(6 - nmpcSine_B.b_i_e, nmpcSine_B.c_A, nmpcSine_B.ii
                            + 2));
        if (nmpcSine_B.s >= 0.0) {
          nmpcSine_B.beta1 = -nmpcSine_B.beta1;
        }

        nmpcSine_B.b_tau[nmpcSine_B.b_i_e] = (nmpcSine_B.beta1 - nmpcSine_B.s) /
          nmpcSine_B.beta1;
        nmpcSine_B.s = 1.0 / (nmpcSine_B.s - nmpcSine_B.beta1);
        for (i = nmpcSine_B.aoffset_b; i <= vectorUB_tmp; i += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.c_A[i - 1]);
          _mm_storeu_pd(&nmpcSine_B.c_A[i - 1], _mm_mul_pd(tmp, _mm_set1_pd
            (nmpcSine_B.s)));
        }

        for (i = nmpcSine_B.jA; i <= nmpcSine_B.lastc + 7; i++) {
          nmpcSine_B.c_A[i - 1] *= nmpcSine_B.s;
        }

        for (i = 0; i <= nmpcSine_B.knt; i++) {
          nmpcSine_B.beta1 *= 1.0020841800044864E-292;
        }

        nmpcSine_B.s = nmpcSine_B.beta1;
      } else {
        nmpcSine_B.b_tau[nmpcSine_B.b_i_e] = (nmpcSine_B.beta1 - nmpcSine_B.c_cu)
          / nmpcSine_B.beta1;
        nmpcSine_B.s = 1.0 / (nmpcSine_B.c_cu - nmpcSine_B.beta1);
        nmpcSine_B.knt = nmpcSine_B.ii - nmpcSine_B.b_i_e;
        nmpcSine_B.lastc = (((((nmpcSine_B.knt - nmpcSine_B.ii) + 6) / 2) << 1)
                            + nmpcSine_B.ii) + 2;
        nmpcSine_B.jA = nmpcSine_B.lastc - 2;
        for (i = nmpcSine_B.aoffset_b; i <= nmpcSine_B.jA; i += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.c_A[i - 1]);
          _mm_storeu_pd(&nmpcSine_B.c_A[i - 1], _mm_mul_pd(tmp, _mm_set1_pd
            (nmpcSine_B.s)));
        }

        for (i = nmpcSine_B.lastc; i <= nmpcSine_B.knt + 7; i++) {
          nmpcSine_B.c_A[i - 1] *= nmpcSine_B.s;
        }

        nmpcSine_B.s = nmpcSine_B.beta1;
      }
    }

    nmpcSine_B.c_A[nmpcSine_B.ii] = nmpcSine_B.s;

    // Start for MATLABSystem: '<S19>/MATLAB System'
    if (nmpcSine_B.b_i_e + 1 < 6) {
      nmpcSine_B.c_A[nmpcSine_B.ii] = 1.0;
      nmpcSine_B.knt = nmpcSine_B.ii + 8;
      if (nmpcSine_B.b_tau[nmpcSine_B.b_i_e] != 0.0) {
        nmpcSine_B.aoffset_b = 7 - nmpcSine_B.b_i_e;
        i = nmpcSine_B.ii - nmpcSine_B.b_i_e;
        while ((nmpcSine_B.aoffset_b > 0) && (nmpcSine_B.c_A[i + 6] == 0.0)) {
          nmpcSine_B.aoffset_b--;
          i--;
        }

        nmpcSine_B.lastc = 5 - nmpcSine_B.b_i_e;
        exitg2 = false;
        while ((!exitg2) && (nmpcSine_B.lastc > 0)) {
          i = (nmpcSine_B.lastc - 1) * 7 + nmpcSine_B.ii;
          nmpcSine_B.jA = i + 8;
          do {
            exitg1 = 0;
            if (nmpcSine_B.jA <= (i + nmpcSine_B.aoffset_b) + 7) {
              if (nmpcSine_B.c_A[nmpcSine_B.jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                nmpcSine_B.jA++;
              }
            } else {
              nmpcSine_B.lastc--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }
      } else {
        nmpcSine_B.aoffset_b = 0;
        nmpcSine_B.lastc = 0;
      }

      if (nmpcSine_B.aoffset_b > 0) {
        if (nmpcSine_B.lastc != 0) {
          memset(&nmpcSine_B.work[0], 0, static_cast<uint32_T>(nmpcSine_B.lastc)
                 * sizeof(real_T));
          i = ((nmpcSine_B.lastc - 1) * 7 + nmpcSine_B.ii) + 8;
          for (vectorUB_tmp = nmpcSine_B.knt; vectorUB_tmp <= i; vectorUB_tmp +=
               7) {
            nmpcSine_B.c_cu = 0.0;
            nmpcSine_B.e_a = vectorUB_tmp + nmpcSine_B.aoffset_b;
            for (nmpcSine_B.jA = vectorUB_tmp; nmpcSine_B.jA < nmpcSine_B.e_a;
                 nmpcSine_B.jA++) {
              nmpcSine_B.c_cu += nmpcSine_B.c_A[(nmpcSine_B.ii + nmpcSine_B.jA)
                - vectorUB_tmp] * nmpcSine_B.c_A[nmpcSine_B.jA - 1];
            }

            nmpcSine_B.jA = div_nde_s32_floor((vectorUB_tmp - nmpcSine_B.ii) - 8,
              7);
            nmpcSine_B.work[nmpcSine_B.jA] += nmpcSine_B.c_cu;
          }
        }

        if (!(-nmpcSine_B.b_tau[nmpcSine_B.b_i_e] == 0.0)) {
          nmpcSine_B.jA = nmpcSine_B.ii;
          nmpcSine_B.knt = nmpcSine_B.lastc - 1;
          for (nmpcSine_B.lastc = 0; nmpcSine_B.lastc <= nmpcSine_B.knt;
               nmpcSine_B.lastc++) {
            nmpcSine_B.c_cu = nmpcSine_B.work[nmpcSine_B.lastc];
            if (nmpcSine_B.c_cu != 0.0) {
              nmpcSine_B.c_cu *= -nmpcSine_B.b_tau[nmpcSine_B.b_i_e];
              vectorUB_tmp = nmpcSine_B.jA + 8;
              i = (nmpcSine_B.aoffset_b + nmpcSine_B.jA) + 7;
              for (nmpcSine_B.e_a = vectorUB_tmp; nmpcSine_B.e_a <= i;
                   nmpcSine_B.e_a++) {
                nmpcSine_B.c_A[nmpcSine_B.e_a - 1] += nmpcSine_B.c_A
                  [((nmpcSine_B.ii + nmpcSine_B.e_a) - nmpcSine_B.jA) - 8] *
                  nmpcSine_B.c_cu;
              }
            }

            nmpcSine_B.jA += 7;
          }
        }
      }

      nmpcSine_B.c_A[nmpcSine_B.ii] = nmpcSine_B.s;
    }
  }

  for (nmpcSine_B.ii = 0; nmpcSine_B.ii < 6; nmpcSine_B.ii++) {
    nmpcSine_B.knt = nmpcSine_B.ii + 1;

    // Start for MATLABSystem: '<S19>/MATLAB System'
    memcpy(&nmpcSine_B.A_j[nmpcSine_B.ii * 6], &nmpcSine_B.c_A[nmpcSine_B.ii * 7],
           static_cast<uint32_T>(nmpcSine_B.knt) * sizeof(real_T));
    for (nmpcSine_B.b_i_e = nmpcSine_B.ii + 2; nmpcSine_B.b_i_e < 7;
         nmpcSine_B.b_i_e++) {
      nmpcSine_B.A_j[(nmpcSine_B.b_i_e + 6 * nmpcSine_B.ii) - 1] = 0.0;
    }
  }

  for (nmpcSine_B.b_i_e = 0; nmpcSine_B.b_i_e < 6; nmpcSine_B.b_i_e++) {
    for (nmpcSine_B.ii = 0; nmpcSine_B.ii < 6; nmpcSine_B.ii++) {
      // Start for MATLABSystem: '<S19>/MATLAB System'
      b_S[nmpcSine_B.ii + 6 * nmpcSine_B.b_i_e] = nmpcSine_B.A_j[6 *
        nmpcSine_B.ii + nmpcSine_B.b_i_e];
    }
  }
}

real_T nmpcSine::nmpcSine_xnrm2_m32(int32_T n, const real_T x[7], int32_T ix0)
{
  real_T y;

  // Start for MATLABSystem: '<S20>/MATLAB System'
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      int32_T kend;
      nmpcSine_B.scale_m = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (int32_T k = ix0; k < kend; k++) {
        nmpcSine_B.absxk_m = fabs(x[k - 1]);
        if (nmpcSine_B.absxk_m > nmpcSine_B.scale_m) {
          nmpcSine_B.t_j = nmpcSine_B.scale_m / nmpcSine_B.absxk_m;
          y = y * nmpcSine_B.t_j * nmpcSine_B.t_j + 1.0;
          nmpcSine_B.scale_m = nmpcSine_B.absxk_m;
        } else {
          nmpcSine_B.t_j = nmpcSine_B.absxk_m / nmpcSine_B.scale_m;
          y += nmpcSine_B.t_j * nmpcSine_B.t_j;
        }
      }

      y = nmpcSine_B.scale_m * sqrt(y);
    }
  }

  // End of Start for MATLABSystem: '<S20>/MATLAB System'
  return y;
}

void nmpcSine::EKFCorrectorAdditive_getMeasure(real_T Rs, const real_T x[6],
  const real_T S[36], real_T *zEstimated, real_T Pxy[6], real_T *Sy, real_T
  dHdx[6], real_T *Rsqrt)
{
  __m128d tmp;
  int32_T i;
  int32_T knt;

  // Start for MATLABSystem: '<S20>/MATLAB System'
  //  IMU Measurement Function
  //  Angular velocity is measured
  *Rsqrt = Rs;

  //  IMU Measurement Function
  *zEstimated = x[5];

  //  Angular velocity is measured
  for (nmpcSine_B.b_j_ew = 0; nmpcSine_B.b_j_ew < 6; nmpcSine_B.b_j_ew++) {
    nmpcSine_B.epsilon_j = 1.4901161193847656E-8 * fabs(x[nmpcSine_B.b_j_ew]);
    if ((nmpcSine_B.epsilon_j <= 1.4901161193847656E-8) || rtIsNaN
        (nmpcSine_B.epsilon_j)) {
      nmpcSine_B.epsilon_j = 1.4901161193847656E-8;
    }

    //  IMU Measurement Function
    //  Angular velocity is measured
    for (i = 0; i < 6; i++) {
      nmpcSine_B.imvec_i[i] = x[i];
      nmpcSine_B.beta1_e = 0.0;
      for (knt = 0; knt < 6; knt++) {
        nmpcSine_B.beta1_e += S[6 * knt + nmpcSine_B.b_j_ew] * S[6 * knt + i];
      }

      nmpcSine_B.S[nmpcSine_B.b_j_ew + 6 * i] = nmpcSine_B.beta1_e;
    }

    nmpcSine_B.imvec_i[nmpcSine_B.b_j_ew] = x[nmpcSine_B.b_j_ew] +
      nmpcSine_B.epsilon_j;
    dHdx[nmpcSine_B.b_j_ew] = (nmpcSine_B.imvec_i[5] - x[5]) /
      nmpcSine_B.epsilon_j;
  }

  for (nmpcSine_B.b_j_ew = 0; nmpcSine_B.b_j_ew < 6; nmpcSine_B.b_j_ew++) {
    // Start for MATLABSystem: '<S20>/MATLAB System'
    i = nmpcSine_B.b_j_ew * 6 - 1;
    nmpcSine_B.epsilon_j = 0.0;
    nmpcSine_B.beta1_e = 0.0;

    // Start for MATLABSystem: '<S20>/MATLAB System'
    for (knt = 0; knt < 6; knt++) {
      _mm_storeu_pd(&nmpcSine_B.dv9[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd(S[(i +
        knt) + 1], nmpcSine_B.S[6 * knt + nmpcSine_B.b_j_ew]), _mm_set1_pd
        (dHdx[knt])), _mm_set_pd(nmpcSine_B.epsilon_j, nmpcSine_B.beta1_e)));
      nmpcSine_B.beta1_e = nmpcSine_B.dv9[0];
      nmpcSine_B.epsilon_j = nmpcSine_B.dv9[1];
    }

    Pxy[nmpcSine_B.b_j_ew] = nmpcSine_B.beta1_e;
    nmpcSine_B.M[nmpcSine_B.b_j_ew] = nmpcSine_B.epsilon_j;
  }

  // Start for MATLABSystem: '<S20>/MATLAB System'
  nmpcSine_B.M[6] = Rs;
  for (nmpcSine_B.b_j_ew = 0; nmpcSine_B.b_j_ew < 1; nmpcSine_B.b_j_ew++) {
    nmpcSine_B.epsilon_j = nmpcSine_B.M[0];
    nmpcSine_B.beta1_e = nmpcSine_xnrm2_m32(6, nmpcSine_B.M, 2);
    if (nmpcSine_B.beta1_e != 0.0) {
      nmpcSine_B.beta1_e = nmpcSine_rt_hypotd_snf(nmpcSine_B.M[0],
        nmpcSine_B.beta1_e);
      if (nmpcSine_B.M[0] >= 0.0) {
        nmpcSine_B.beta1_e = -nmpcSine_B.beta1_e;
      }

      if (fabs(nmpcSine_B.beta1_e) < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          for (i = 2; i <= 6; i += 2) {
            tmp = _mm_loadu_pd(&nmpcSine_B.M[i - 1]);
            _mm_storeu_pd(&nmpcSine_B.M[i - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          nmpcSine_B.beta1_e *= 9.9792015476736E+291;
          nmpcSine_B.epsilon_j *= 9.9792015476736E+291;
        } while ((fabs(nmpcSine_B.beta1_e) < 1.0020841800044864E-292) && (knt +
                  1 < 20));

        nmpcSine_B.beta1_e = nmpcSine_rt_hypotd_snf(nmpcSine_B.epsilon_j,
          nmpcSine_xnrm2_m32(6, nmpcSine_B.M, 2));
        if (nmpcSine_B.epsilon_j >= 0.0) {
          nmpcSine_B.beta1_e = -nmpcSine_B.beta1_e;
        }

        nmpcSine_B.epsilon_j = 1.0 / (nmpcSine_B.epsilon_j - nmpcSine_B.beta1_e);
        for (i = 2; i <= 6; i += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.M[i - 1]);
          _mm_storeu_pd(&nmpcSine_B.M[i - 1], _mm_mul_pd(tmp, _mm_set1_pd
            (nmpcSine_B.epsilon_j)));
        }

        for (i = 0; i <= knt; i++) {
          nmpcSine_B.beta1_e *= 1.0020841800044864E-292;
        }

        nmpcSine_B.epsilon_j = nmpcSine_B.beta1_e;
      } else {
        nmpcSine_B.epsilon_j = 1.0 / (nmpcSine_B.M[0] - nmpcSine_B.beta1_e);
        for (i = 2; i <= 6; i += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.M[i - 1]);
          _mm_storeu_pd(&nmpcSine_B.M[i - 1], _mm_mul_pd(tmp, _mm_set1_pd
            (nmpcSine_B.epsilon_j)));
        }

        nmpcSine_B.epsilon_j = nmpcSine_B.beta1_e;
      }
    }

    nmpcSine_B.M[0] = nmpcSine_B.epsilon_j;
  }

  *Sy = nmpcSine_B.M[0];
}

real_T nmpcSine::nmpcSine_xnrm2_m32u(int32_T n, const real_T x[42], int32_T ix0)
{
  real_T y;

  // Start for MATLABSystem: '<S20>/MATLAB System'
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (int32_T k = ix0; k < kend; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  // End of Start for MATLABSystem: '<S20>/MATLAB System'
  return y;
}

void nmpcSine::nmpcSine_qr_m3(const real_T A[42], real_T Q[42], real_T R[36])
{
  __m128d tmp;
  real_T b_tau[6];
  real_T work[6];
  real_T atmp;
  real_T beta1;
  real_T c;
  int32_T d;
  int32_T exitg1;
  int32_T i;
  int32_T iac;
  int32_T iaii;
  int32_T ii;
  int32_T jA;
  int32_T knt;
  int32_T lastc;
  int32_T lastv;
  int32_T scalarLB;
  boolean_T exitg2;
  for (i = 0; i < 6; i++) {
    // Start for MATLABSystem: '<S20>/MATLAB System'
    b_tau[i] = 0.0;
  }

  // Start for MATLABSystem: '<S20>/MATLAB System'
  memcpy(&Q[0], &A[0], 42U * sizeof(real_T));
  for (i = 0; i < 6; i++) {
    // Start for MATLABSystem: '<S20>/MATLAB System'
    work[i] = 0.0;
  }

  // Start for MATLABSystem: '<S20>/MATLAB System'
  for (iaii = 0; iaii < 6; iaii++) {
    ii = iaii * 7 + iaii;
    i = ii + 2;
    atmp = Q[ii];
    b_tau[iaii] = 0.0;
    beta1 = nmpcSine_xnrm2_m32u(6 - iaii, Q, ii + 2);
    if (beta1 != 0.0) {
      c = Q[ii];
      beta1 = nmpcSine_rt_hypotd_snf(c, beta1);
      if (c >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          scalarLB = ii - iaii;
          lastc = (((((scalarLB - ii) + 6) / 2) << 1) + ii) + 2;
          jA = lastc - 2;
          for (lastv = i; lastv <= jA; lastv += 2) {
            tmp = _mm_loadu_pd(&Q[lastv - 1]);
            _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (lastv = lastc; lastv <= scalarLB + 7; lastv++) {
            Q[lastv - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          atmp *= 9.9792015476736E+291;
        } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt + 1 < 20));

        beta1 = nmpcSine_rt_hypotd_snf(atmp, nmpcSine_xnrm2_m32u(6 - iaii, Q, ii
          + 2));
        if (atmp >= 0.0) {
          beta1 = -beta1;
        }

        b_tau[iaii] = (beta1 - atmp) / beta1;
        atmp = 1.0 / (atmp - beta1);
        for (lastv = i; lastv <= jA; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = lastc; lastv <= scalarLB + 7; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        for (lastv = 0; lastv <= knt; lastv++) {
          beta1 *= 1.0020841800044864E-292;
        }

        atmp = beta1;
      } else {
        b_tau[iaii] = (beta1 - c) / beta1;
        atmp = 1.0 / (c - beta1);
        knt = ii - iaii;
        scalarLB = (((((knt - ii) + 6) / 2) << 1) + ii) + 2;
        lastc = scalarLB - 2;
        for (lastv = i; lastv <= lastc; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = scalarLB; lastv <= knt + 7; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        atmp = beta1;
      }
    }

    Q[ii] = atmp;
    if (iaii + 1 < 6) {
      Q[ii] = 1.0;
      scalarLB = ii + 8;
      if (b_tau[iaii] != 0.0) {
        lastv = 7 - iaii;
        i = ii - iaii;
        while ((lastv > 0) && (Q[i + 6] == 0.0)) {
          lastv--;
          i--;
        }

        knt = 5 - iaii;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          lastc = (knt - 1) * 7 + ii;
          jA = lastc + 8;
          do {
            exitg1 = 0;
            if (jA <= (lastc + lastv) + 7) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            memset(&work[0], 0, static_cast<uint32_T>(lastc + 1) * sizeof(real_T));
          }

          knt = (7 * lastc + ii) + 8;
          for (iac = scalarLB; iac <= knt; iac += 7) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[(ii + jA) - iac] * Q[jA - 1];
            }

            jA = div_nde_s32_floor((iac - ii) - 8, 7);
            work[jA] += c;
          }
        }

        if (!(-b_tau[iaii] == 0.0)) {
          jA = ii;
          for (iac = 0; iac <= lastc; iac++) {
            c = work[iac];
            if (c != 0.0) {
              c *= -b_tau[iaii];
              knt = jA + 8;
              scalarLB = (lastv + jA) + 7;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((ii + d) - jA) - 8] * c;
              }
            }

            jA += 7;
          }
        }
      }

      Q[ii] = atmp;
    }
  }

  for (ii = 0; ii < 6; ii++) {
    for (iaii = 0; iaii <= ii; iaii++) {
      // Start for MATLABSystem: '<S20>/MATLAB System'
      R[iaii + 6 * ii] = Q[7 * ii + iaii];
    }

    for (iaii = ii + 2; iaii < 7; iaii++) {
      R[(iaii + 6 * ii) - 1] = 0.0;
    }

    // Start for MATLABSystem: '<S20>/MATLAB System'
    work[ii] = 0.0;
  }

  // Start for MATLABSystem: '<S20>/MATLAB System'
  for (i = 5; i >= 0; i--) {
    iaii = (i * 7 + i) + 7;
    if (i + 1 < 6) {
      Q[iaii - 7] = 1.0;
      scalarLB = iaii + 1;
      if (b_tau[i] != 0.0) {
        lastv = 7 - i;
        ii = (iaii - i) - 1;
        while ((lastv > 0) && (Q[ii] == 0.0)) {
          lastv--;
          ii--;
        }

        knt = 5 - i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          lastc = (knt - 1) * 7 + iaii;
          jA = lastc + 1;
          do {
            exitg1 = 0;
            if (jA <= lastc + lastv) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            memset(&work[0], 0, static_cast<uint32_T>(lastc + 1) * sizeof(real_T));
          }

          knt = (7 * lastc + iaii) + 1;
          for (iac = scalarLB; iac <= knt; iac += 7) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[((iaii + jA) - iac) - 7] * Q[jA - 1];
            }

            jA = div_nde_s32_floor((iac - iaii) - 1, 7);
            work[jA] += c;
          }
        }

        if (!(-b_tau[i] == 0.0)) {
          jA = iaii;
          for (iac = 0; iac <= lastc; iac++) {
            c = work[iac];
            if (c != 0.0) {
              c *= -b_tau[i];
              knt = jA + 1;
              scalarLB = lastv + jA;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((iaii + d) - jA) - 8] * c;
              }
            }

            jA += 7;
          }
        }
      }
    }

    knt = iaii - i;
    scalarLB = (((((knt - iaii) + 6) / 2) << 1) + iaii) - 5;
    lastc = scalarLB - 2;
    for (lastv = iaii - 5; lastv <= lastc; lastv += 2) {
      tmp = _mm_loadu_pd(&Q[lastv - 1]);
      _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(-b_tau[i])));
    }

    for (lastv = scalarLB; lastv <= knt; lastv++) {
      Q[lastv - 1] *= -b_tau[i];
    }

    Q[iaii - 7] = 1.0 - b_tau[i];
    for (ii = 0; ii < i; ii++) {
      Q[(iaii - ii) - 8] = 0.0;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_getXUe(const real_T z[125], const real_T x[6], real_T X
  [126], real_T U[84], real_T *e)
{
  static const int8_T y[160] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
    1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
    0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };

  memset(&X[0], 0, 126U * sizeof(real_T));
  memset(&U[0], 0, 84U * sizeof(real_T));
  memset(&nmpcSine_B.Umv[0], 0, 42U * sizeof(real_T));
  nmpcSine_B.z_ch = z[121];
  nmpcSine_B.z_a = z[120];
  nmpcSine_B.z_da = z[122];
  nmpcSine_B.z_af = z[123];
  for (nmpcSine_B.i3 = 0; nmpcSine_B.i3 <= 38; nmpcSine_B.i3 += 2) {
    _mm_storeu_pd(&nmpcSine_B.y_j[nmpcSine_B.i3], _mm_add_pd(_mm_add_pd
      (_mm_add_pd(_mm_mul_pd(_mm_set_pd(static_cast<real_T>(y[nmpcSine_B.i3 + 41]),
      static_cast<real_T>(y[nmpcSine_B.i3 + 40])), _mm_set1_pd(nmpcSine_B.z_ch)),
                  _mm_mul_pd(_mm_set_pd(static_cast<real_T>(y[nmpcSine_B.i3 + 1]),
      static_cast<real_T>(y[nmpcSine_B.i3])), _mm_set1_pd(nmpcSine_B.z_a))),
       _mm_mul_pd(_mm_set_pd(static_cast<real_T>(y[nmpcSine_B.i3 + 81]),
      static_cast<real_T>(y[nmpcSine_B.i3 + 80])), _mm_set1_pd(nmpcSine_B.z_da))),
      _mm_mul_pd(_mm_set_pd(static_cast<real_T>(y[nmpcSine_B.i3 + 121]),
      static_cast<real_T>(y[nmpcSine_B.i3 + 120])), _mm_set1_pd(nmpcSine_B.z_af))));
  }

  for (nmpcSine_B.i3 = 0; nmpcSine_B.i3 < 20; nmpcSine_B.i3++) {
    nmpcSine_B.Umv_tmp = nmpcSine_B.i3 << 1;
    nmpcSine_B.Umv[nmpcSine_B.i3] = nmpcSine_B.y_j[nmpcSine_B.Umv_tmp];
    nmpcSine_B.Umv[nmpcSine_B.i3 + 21] = nmpcSine_B.y_j[nmpcSine_B.Umv_tmp + 1];
  }

  *e = z[124];
  memcpy(&nmpcSine_B.z_d[0], &z[0], 120U * sizeof(real_T));
  for (nmpcSine_B.i3 = 0; nmpcSine_B.i3 < 6; nmpcSine_B.i3++) {
    for (nmpcSine_B.Umv_tmp = 0; nmpcSine_B.Umv_tmp < 20; nmpcSine_B.Umv_tmp++)
    {
      X[(nmpcSine_B.Umv_tmp + 21 * nmpcSine_B.i3) + 1] = nmpcSine_B.z_d[6 *
        nmpcSine_B.Umv_tmp + nmpcSine_B.i3];
    }

    X[21 * nmpcSine_B.i3] = x[nmpcSine_B.i3];
  }

  nmpcSine_B.Umv[20] = nmpcSine_B.Umv[19];
  nmpcSine_B.Umv[41] = nmpcSine_B.Umv[40];
  for (nmpcSine_B.i3 = 0; nmpcSine_B.i3 < 21; nmpcSine_B.i3++) {
    U[nmpcSine_B.i3] = nmpcSine_B.Umv[nmpcSine_B.i3];
    U[nmpcSine_B.i3 + 21] = nmpcSine_B.Umv[nmpcSine_B.i3 + 21];
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_costFcn(const real_T X[126], const real_T U[84], const
  real_T data_References[80], real_T Q, real_T R, real_T Qt)
{
  __m128d tmp;
  __m128d tmp_0;
  nmpcSine_B.data_References_nb = data_References[0];
  nmpcSine_B.data_References_i = data_References[20];
  for (nmpcSine_B.k_f = 0; nmpcSine_B.k_f <= 18; nmpcSine_B.k_f += 2) {
    _mm_storeu_pd(&nmpcSine_B.err[nmpcSine_B.k_f], _mm_sub_pd(_mm_set1_pd
      (nmpcSine_B.data_References_nb), _mm_loadu_pd(&X[nmpcSine_B.k_f + 1])));
    _mm_storeu_pd(&nmpcSine_B.err[nmpcSine_B.k_f + 20], _mm_sub_pd(_mm_set1_pd
      (nmpcSine_B.data_References_i), _mm_loadu_pd(&X[nmpcSine_B.k_f + 22])));
  }

  for (nmpcSine_B.k_f = 0; nmpcSine_B.k_f <= 16; nmpcSine_B.k_f += 2) {
    tmp = _mm_loadu_pd(&nmpcSine_B.err[nmpcSine_B.k_f]);
    tmp_0 = _mm_loadu_pd(&nmpcSine_B.err[nmpcSine_B.k_f + 20]);
    _mm_storeu_pd(&nmpcSine_B.y_n[nmpcSine_B.k_f], _mm_add_pd(_mm_mul_pd(tmp,
      tmp), _mm_mul_pd(tmp_0, tmp_0)));
  }

  for (nmpcSine_B.k_f = 18; nmpcSine_B.k_f < 19; nmpcSine_B.k_f++) {
    nmpcSine_B.data_References_nb = nmpcSine_B.err[nmpcSine_B.k_f];
    nmpcSine_B.data_References_i = nmpcSine_B.err[nmpcSine_B.k_f + 20];
    nmpcSine_B.y_n[nmpcSine_B.k_f] = nmpcSine_B.data_References_nb *
      nmpcSine_B.data_References_nb + nmpcSine_B.data_References_i *
      nmpcSine_B.data_References_i;
  }

  nmpcSine_B.data_References_nb = nmpcSine_B.y_n[0];
  for (nmpcSine_B.k_f = 0; nmpcSine_B.k_f < 18; nmpcSine_B.k_f++) {
    nmpcSine_B.data_References_nb += nmpcSine_B.y_n[nmpcSine_B.k_f + 1];
  }

  for (nmpcSine_B.k_f = 0; nmpcSine_B.k_f <= 16; nmpcSine_B.k_f += 2) {
    tmp = _mm_loadu_pd(&U[nmpcSine_B.k_f]);
    tmp_0 = _mm_loadu_pd(&U[nmpcSine_B.k_f + 21]);
    _mm_storeu_pd(&nmpcSine_B.y_n[nmpcSine_B.k_f], _mm_add_pd(_mm_mul_pd(tmp,
      tmp), _mm_mul_pd(tmp_0, tmp_0)));
  }

  for (nmpcSine_B.k_f = 18; nmpcSine_B.k_f < 19; nmpcSine_B.k_f++) {
    nmpcSine_B.data_References_i = U[nmpcSine_B.k_f];
    nmpcSine_B.U_m = U[nmpcSine_B.k_f + 21];
    nmpcSine_B.y_n[nmpcSine_B.k_f] = nmpcSine_B.data_References_i *
      nmpcSine_B.data_References_i + nmpcSine_B.U_m * nmpcSine_B.U_m;
  }

  nmpcSine_B.data_References_i = nmpcSine_B.y_n[0];
  for (nmpcSine_B.k_f = 0; nmpcSine_B.k_f < 18; nmpcSine_B.k_f++) {
    nmpcSine_B.data_References_i += nmpcSine_B.y_n[nmpcSine_B.k_f + 1];
  }

  return (nmpcSine_B.err[19] * nmpcSine_B.err[19] + nmpcSine_B.err[39] *
          nmpcSine_B.err[39]) * Qt + (Q * nmpcSine_B.data_References_nb + R *
    nmpcSine_B.data_References_i);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_mtimes(const real_T A_data[], const int32_T A_size[2],
  real_T C_data[], int32_T C_size[2])
{
  static const int8_T b[160] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
    1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
    0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };

  nmpcSine_B.mc = A_size[0] - 1;
  C_size[0] = A_size[0];
  C_size[1] = 4;
  for (nmpcSine_B.j = 0; nmpcSine_B.j < 4; nmpcSine_B.j++) {
    nmpcSine_B.coffset_f = (nmpcSine_B.mc + 1) * nmpcSine_B.j;
    nmpcSine_B.boffset = nmpcSine_B.j * 40;
    if (nmpcSine_B.mc >= 0) {
      memset(&C_data[nmpcSine_B.coffset_f], 0, static_cast<uint32_T>
             (nmpcSine_B.mc + 1) * sizeof(real_T));
    }

    for (nmpcSine_B.i_i = 0; nmpcSine_B.i_i < 40; nmpcSine_B.i_i++) {
      nmpcSine_B.aoffset_n = nmpcSine_B.i_i * A_size[0];
      nmpcSine_B.bkj = b[nmpcSine_B.boffset + nmpcSine_B.i_i];
      nmpcSine_B.scalarLB_i = ((nmpcSine_B.mc + 1) / 2) << 1;
      nmpcSine_B.vectorUB_k = nmpcSine_B.scalarLB_i - 2;
      for (nmpcSine_B.b_i_l = 0; nmpcSine_B.b_i_l <= nmpcSine_B.vectorUB_k;
           nmpcSine_B.b_i_l += 2) {
        __m128d tmp;
        nmpcSine_B.i7 = nmpcSine_B.coffset_f + nmpcSine_B.b_i_l;
        tmp = _mm_loadu_pd(&C_data[nmpcSine_B.i7]);
        _mm_storeu_pd(&C_data[nmpcSine_B.i7], _mm_add_pd(_mm_mul_pd(_mm_loadu_pd
          (&A_data[nmpcSine_B.aoffset_n + nmpcSine_B.b_i_l]), _mm_set1_pd(
          static_cast<real_T>(nmpcSine_B.bkj))), tmp));
      }

      for (nmpcSine_B.b_i_l = nmpcSine_B.scalarLB_i; nmpcSine_B.b_i_l <=
           nmpcSine_B.mc; nmpcSine_B.b_i_l++) {
        nmpcSine_B.i7 = nmpcSine_B.coffset_f + nmpcSine_B.b_i_l;
        C_data[nmpcSine_B.i7] += A_data[nmpcSine_B.aoffset_n + nmpcSine_B.b_i_l]
          * static_cast<real_T>(nmpcSine_B.bkj);
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_getUBounds(const real_T runtimedata_lastMV[2], const
  real_T runtimedata_MVMin[40], const real_T runtimedata_MVMax[40], const real_T
  runtimedata_MVRateMin[40], const real_T runtimedata_MVRateMax[40], real_T
  A_data[], int32_T A_size[2], real_T Bu_data[], int32_T Bu_size[1])
{
  int32_T Au_tmp;
  int32_T Au_tmp_0;
  int32_T Auf_data_tmp;
  int32_T Auf_data_tmp_0;
  boolean_T exitg1;
  memset(&nmpcSine_B.Au[0], 0, 6400U * sizeof(int8_T));
  memset(&nmpcSine_B.b_Bu[0], 0, 160U * sizeof(real_T));
  memset(&nmpcSine_B.x[0], 0, 160U * sizeof(boolean_T));
  nmpcSine_B.ic_idx_0_d = 1;
  nmpcSine_B.ic_idx_1_e = 2;
  for (nmpcSine_B.i_a = 0; nmpcSine_B.i_a < 20; nmpcSine_B.i_a++) {
    nmpcSine_B.runtimedata_MVRateMin = runtimedata_MVRateMin[nmpcSine_B.i_a];
    nmpcSine_B.x[nmpcSine_B.ic_idx_0_d - 1] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVRateMin)) && (!rtIsNaN
      (nmpcSine_B.runtimedata_MVRateMin)));
    nmpcSine_B.runtimedata_MVRateMin_k = runtimedata_MVRateMin[nmpcSine_B.i_a +
      20];
    nmpcSine_B.x[nmpcSine_B.ic_idx_1_e - 1] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVRateMin_k)) && (!rtIsNaN
      (nmpcSine_B.runtimedata_MVRateMin_k)));
    nmpcSine_B.runtimedata_MVRateMax = runtimedata_MVRateMax[nmpcSine_B.i_a];
    nmpcSine_B.x[nmpcSine_B.ic_idx_0_d + 1] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVRateMax)) && (!rtIsNaN
      (nmpcSine_B.runtimedata_MVRateMax)));
    nmpcSine_B.runtimedata_MVRateMax_i = runtimedata_MVRateMax[nmpcSine_B.i_a +
      20];
    nmpcSine_B.x[nmpcSine_B.ic_idx_1_e + 1] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVRateMax_i)) && (!rtIsNaN
      (nmpcSine_B.runtimedata_MVRateMax_i)));
    nmpcSine_B.runtimedata_MVMin = runtimedata_MVMin[nmpcSine_B.i_a];
    nmpcSine_B.x[nmpcSine_B.ic_idx_0_d + 3] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVMin)) && (!rtIsNaN(nmpcSine_B.runtimedata_MVMin)));
    nmpcSine_B.runtimedata_MVMin_o = runtimedata_MVMin[nmpcSine_B.i_a + 20];
    nmpcSine_B.x[nmpcSine_B.ic_idx_1_e + 3] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVMin_o)) && (!rtIsNaN
      (nmpcSine_B.runtimedata_MVMin_o)));
    nmpcSine_B.runtimedata_MVMax = runtimedata_MVMax[nmpcSine_B.i_a];
    nmpcSine_B.x[nmpcSine_B.ic_idx_0_d + 5] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVMax)) && (!rtIsNaN(nmpcSine_B.runtimedata_MVMax)));
    nmpcSine_B.runtimedata_MVMax_m = runtimedata_MVMax[nmpcSine_B.i_a + 20];
    nmpcSine_B.x[nmpcSine_B.ic_idx_1_e + 5] = ((!rtIsInf
      (nmpcSine_B.runtimedata_MVMax_m)) && (!rtIsNaN
      (nmpcSine_B.runtimedata_MVMax_m)));
    Au_tmp = 320 * nmpcSine_B.i_a + nmpcSine_B.ic_idx_0_d;
    nmpcSine_B.Au[Au_tmp - 1] = -1;
    Au_tmp_0 = 320 * nmpcSine_B.i_a + nmpcSine_B.ic_idx_1_e;
    nmpcSine_B.Au[Au_tmp_0 - 1] = 0;
    nmpcSine_B.Au[Au_tmp + 159] = 0;
    nmpcSine_B.Au[Au_tmp_0 + 159] = -1;
    nmpcSine_B.Au[Au_tmp + 1] = 1;
    nmpcSine_B.Au[Au_tmp_0 + 1] = 0;
    nmpcSine_B.Au[Au_tmp + 161] = 0;
    nmpcSine_B.Au[Au_tmp_0 + 161] = 1;
    nmpcSine_B.Au[Au_tmp + 3] = -1;
    nmpcSine_B.Au[Au_tmp_0 + 3] = 0;
    nmpcSine_B.Au[Au_tmp + 163] = 0;
    nmpcSine_B.Au[Au_tmp_0 + 163] = -1;
    nmpcSine_B.Au[Au_tmp + 5] = 1;
    nmpcSine_B.Au[Au_tmp_0 + 5] = 0;
    nmpcSine_B.Au[Au_tmp + 165] = 0;
    nmpcSine_B.Au[Au_tmp_0 + 165] = 1;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_0_d - 1] =
      -nmpcSine_B.runtimedata_MVRateMin;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e - 1] =
      -nmpcSine_B.runtimedata_MVRateMin_k;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_0_d + 1] =
      nmpcSine_B.runtimedata_MVRateMax;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e + 1] =
      nmpcSine_B.runtimedata_MVRateMax_i;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_0_d + 3] = -nmpcSine_B.runtimedata_MVMin;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e + 3] = -nmpcSine_B.runtimedata_MVMin_o;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_0_d + 5] = nmpcSine_B.runtimedata_MVMax;
    nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e + 5] = nmpcSine_B.runtimedata_MVMax_m;
    if (nmpcSine_B.i_a + 1 == 1) {
      nmpcSine_B.runtimedata_MVRateMin = nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e -
        1] - runtimedata_lastMV[1];
      nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_0_d - 1] -= runtimedata_lastMV[0];
      nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e - 1] =
        nmpcSine_B.runtimedata_MVRateMin;
      nmpcSine_B.runtimedata_MVRateMin = nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e +
        1] + runtimedata_lastMV[1];
      nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_0_d + 1] += runtimedata_lastMV[0];
      nmpcSine_B.b_Bu[nmpcSine_B.ic_idx_1_e + 1] =
        nmpcSine_B.runtimedata_MVRateMin;
    } else {
      Au_tmp = (nmpcSine_B.i_a - 1) * 320;
      Au_tmp_0 = nmpcSine_B.ic_idx_0_d + Au_tmp;
      nmpcSine_B.Au[Au_tmp_0 - 1] = 1;
      Au_tmp += nmpcSine_B.ic_idx_1_e;
      nmpcSine_B.Au[Au_tmp - 1] = 0;
      nmpcSine_B.Au[Au_tmp_0 + 159] = 0;
      nmpcSine_B.Au[Au_tmp + 159] = 1;
      nmpcSine_B.Au[Au_tmp_0 + 1] = -1;
      nmpcSine_B.Au[Au_tmp + 1] = 0;
      nmpcSine_B.Au[Au_tmp_0 + 161] = 0;
      nmpcSine_B.Au[Au_tmp + 161] = -1;
    }

    nmpcSine_B.ic_idx_0_d += 8;
    nmpcSine_B.ic_idx_1_e += 8;
  }

  nmpcSine_B.i_a = 0;
  nmpcSine_B.ic_idx_0_d = 0;
  exitg1 = false;
  while ((!exitg1) && (nmpcSine_B.ic_idx_0_d < 160)) {
    if (nmpcSine_B.x[nmpcSine_B.ic_idx_0_d]) {
      nmpcSine_B.i_a++;
      nmpcSine_B.ii_data[nmpcSine_B.i_a - 1] = static_cast<uint8_T>
        (nmpcSine_B.ic_idx_0_d + 1);
      if (nmpcSine_B.i_a >= 160) {
        exitg1 = true;
      } else {
        nmpcSine_B.ic_idx_0_d++;
      }
    } else {
      nmpcSine_B.ic_idx_0_d++;
    }
  }

  if (nmpcSine_B.i_a < 1) {
    nmpcSine_B.i_a = 0;
  }

  if (nmpcSine_B.i_a > 0) {
    Bu_size[0] = nmpcSine_B.i_a;
    nmpcSine_B.ic_idx_0_d = (nmpcSine_B.i_a / 2) << 1;
    nmpcSine_B.ic_idx_1_e = nmpcSine_B.ic_idx_0_d - 2;
    for (Au_tmp = 0; Au_tmp <= nmpcSine_B.ic_idx_1_e; Au_tmp += 2) {
      Bu_data[Au_tmp] = nmpcSine_B.b_Bu[nmpcSine_B.ii_data[Au_tmp] - 1];
      Bu_data[Au_tmp + 1] = nmpcSine_B.b_Bu[nmpcSine_B.ii_data[Au_tmp + 1] - 1];
    }

    for (Au_tmp = nmpcSine_B.ic_idx_0_d; Au_tmp < nmpcSine_B.i_a; Au_tmp++) {
      Bu_data[Au_tmp] = nmpcSine_B.b_Bu[nmpcSine_B.ii_data[Au_tmp] - 1];
    }

    for (Au_tmp = 0; Au_tmp < 2; Au_tmp++) {
      for (Au_tmp_0 = 0; Au_tmp_0 < 20; Au_tmp_0++) {
        for (nmpcSine_B.b_i_f = 0; nmpcSine_B.b_i_f <= nmpcSine_B.ic_idx_1_e;
             nmpcSine_B.b_i_f += 2) {
          Auf_data_tmp = nmpcSine_B.i_a * Au_tmp;
          Auf_data_tmp_0 = (nmpcSine_B.i_a << 1) * Au_tmp_0;
          nmpcSine_B.Auf_data[(nmpcSine_B.b_i_f + Auf_data_tmp) + Auf_data_tmp_0]
            = nmpcSine_B.Au[((160 * Au_tmp + nmpcSine_B.ii_data[nmpcSine_B.b_i_f])
                             + 320 * Au_tmp_0) - 1];
          nmpcSine_B.Auf_data[((nmpcSine_B.b_i_f + Auf_data_tmp) +
                               Auf_data_tmp_0) + 1] = nmpcSine_B.Au[((160 *
            Au_tmp + nmpcSine_B.ii_data[nmpcSine_B.b_i_f + 1]) + 320 * Au_tmp_0)
            - 1];
        }

        for (nmpcSine_B.b_i_f = nmpcSine_B.ic_idx_0_d; nmpcSine_B.b_i_f <
             nmpcSine_B.i_a; nmpcSine_B.b_i_f++) {
          nmpcSine_B.Auf_data[(nmpcSine_B.b_i_f + nmpcSine_B.i_a * Au_tmp) +
            (nmpcSine_B.i_a << 1) * Au_tmp_0] = nmpcSine_B.Au[((160 * Au_tmp +
            nmpcSine_B.ii_data[nmpcSine_B.b_i_f]) + 320 * Au_tmp_0) - 1];
        }
      }
    }

    nmpcSine_B.b_i[0] = nmpcSine_B.i_a;
    nmpcSine_B.b_i[1] = 40;
    nmpcSine_mtimes(nmpcSine_B.Auf_data, nmpcSine_B.b_i, nmpcSine_B.tmp_data_g,
                    nmpcSine_B.tmp_size_f);
    A_size[0] = nmpcSine_B.i_a;
    A_size[1] = 125;
    nmpcSine_B.ic_idx_0_d = nmpcSine_B.i_a * 120;
    memset(&A_data[0], 0, static_cast<uint32_T>(nmpcSine_B.ic_idx_0_d) * sizeof
           (real_T));
    nmpcSine_B.ic_idx_0_d = nmpcSine_B.i_a << 2;
    for (Au_tmp = 0; Au_tmp < nmpcSine_B.ic_idx_0_d; Au_tmp++) {
      A_data[Au_tmp + nmpcSine_B.i_a * 120] = nmpcSine_B.tmp_data_g[Au_tmp];
    }

    memset(&A_data[nmpcSine_B.i_a * 120 + nmpcSine_B.ic_idx_0_d], 0,
           static_cast<uint32_T>((((nmpcSine_B.i_a + nmpcSine_B.i_a * 120) +
              nmpcSine_B.ic_idx_0_d) - nmpcSine_B.i_a * 120) -
            nmpcSine_B.ic_idx_0_d) * sizeof(real_T));
  } else {
    Bu_size[0] = 0;
    A_size[0] = 0;
    A_size[1] = 161;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_eye(real_T b_I[36])
{
  memset(&b_I[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.k_d = 0; nmpcSine_B.k_d < 6; nmpcSine_B.k_d++) {
    b_I[nmpcSine_B.k_d + 6 * nmpcSine_B.k_d] = 1.0;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nm_stateTransitionFcnJacobianAD(const real_T inputVariables[10],
  const real_T extraParams[67], real_T obj[6], real_T grad[60])
{
  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  __m128d tmp_2;
  __m128d tmp_3;
  __m128d tmp_4;
  __m128d tmp_5;
  __m128d tmp_6;
  __m128d tmp_7;
  __m128d tmp_8;
  nmpcSine_B.arg6_tmp = sin(inputVariables[2]);
  nmpcSine_B.arg10_tmp = cos(inputVariables[2]);
  tmp_6 = _mm_set_pd(nmpcSine_B.arg6_tmp, nmpcSine_B.arg10_tmp);
  _mm_storeu_pd(&obj[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd(inputVariables[4],
    -inputVariables[4]), _mm_set_pd(nmpcSine_B.arg10_tmp, nmpcSine_B.arg6_tmp)),
    _mm_mul_pd(tmp_6, _mm_set1_pd(inputVariables[3]))));
  obj[2] = inputVariables[5];
  nmpcSine_B.arg45 = extraParams[3] * inputVariables[4];
  obj[3] = (((((nmpcSine_B.arg10_tmp * inputVariables[8] + nmpcSine_B.arg6_tmp *
                inputVariables[9]) - extraParams[0] * inputVariables[3]) +
              extraParams[1] * inputVariables[6]) + extraParams[2] *
             inputVariables[7]) + nmpcSine_B.arg45 * inputVariables[5] *
            extraParams[4]) + inputVariables[5] * inputVariables[5] *
    extraParams[5] * extraParams[6] * extraParams[7];
  nmpcSine_B.arg96 = inputVariables[5] * extraParams[16] * extraParams[17];
  tmp_8 = _mm_set1_pd(inputVariables[5]);
  _mm_storeu_pd(&nmpcSine_B.dv6[0], _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(tmp_8,
    _mm_set_pd(extraParams[27], extraParams[24])), _mm_set_pd(extraParams[28],
    extraParams[25])), _mm_set_pd(extraParams[29], extraParams[26])));
  nmpcSine_B.arg121 = nmpcSine_B.dv6[0];
  nmpcSine_B.arg132 = nmpcSine_B.dv6[1];
  tmp_7 = _mm_set1_pd(inputVariables[4]);
  _mm_storeu_pd(&nmpcSine_B.dv6[0], _mm_mul_pd(_mm_mul_pd(tmp_7, _mm_set_pd
    (extraParams[33], extraParams[30])), _mm_set_pd(extraParams[34],
    extraParams[31])));
  nmpcSine_B.arg140 = nmpcSine_B.dv6[0];
  nmpcSine_B.arg150 = nmpcSine_B.dv6[1];
  _mm_storeu_pd(&nmpcSine_B.dv6[0], _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(_mm_set_pd
    (extraParams[39], extraParams[36]), tmp_6), _mm_set_pd(extraParams[40],
    extraParams[37])), _mm_set_pd(extraParams[41], extraParams[38])));
  nmpcSine_B.arg163 = nmpcSine_B.dv6[0];
  nmpcSine_B.arg176 = nmpcSine_B.dv6[1];
  obj[4] = ((((((((((((nmpcSine_B.arg10_tmp * inputVariables[9] -
                       nmpcSine_B.arg6_tmp * inputVariables[8]) -
                      inputVariables[4] * extraParams[8] * extraParams[9]) -
                     inputVariables[4] * extraParams[10] * extraParams[11] *
                     extraParams[12]) + inputVariables[5] * extraParams[13] *
                    extraParams[14] * extraParams[15]) - nmpcSine_B.arg96 *
                   inputVariables[3]) - inputVariables[6] * extraParams[18] *
                  extraParams[19] * extraParams[20]) + inputVariables[7] *
                 extraParams[21] * extraParams[22] * extraParams[23]) -
                nmpcSine_B.arg121 * inputVariables[3]) + nmpcSine_B.arg132 *
               inputVariables[3]) - nmpcSine_B.arg140 * inputVariables[3] *
              extraParams[32]) + nmpcSine_B.arg150 * inputVariables[3] *
             extraParams[35]) + nmpcSine_B.dv6[0] * inputVariables[9]) -
    nmpcSine_B.dv6[1] * inputVariables[8];
  _mm_storeu_pd(&nmpcSine_B.dv6[0], _mm_mul_pd(_mm_mul_pd(tmp_7, _mm_set_pd
    (extraParams[53], extraParams[51])), _mm_set_pd(extraParams[54],
    extraParams[52])));
  nmpcSine_B.arg212 = nmpcSine_B.dv6[0];
  nmpcSine_B.arg222 = nmpcSine_B.dv6[1];
  _mm_storeu_pd(&nmpcSine_B.dv6[0], _mm_mul_pd(_mm_mul_pd(tmp_8, _mm_set_pd
    (extraParams[58], extraParams[55])), _mm_set_pd(extraParams[59],
    extraParams[56])));
  nmpcSine_B.arg232 = nmpcSine_B.dv6[0];
  nmpcSine_B.arg244 = nmpcSine_B.dv6[1];
  _mm_storeu_pd(&nmpcSine_B.dv6[0], _mm_mul_pd(_mm_mul_pd(_mm_set_pd
    (extraParams[64], extraParams[61]), tmp_6), _mm_set_pd(extraParams[65],
    extraParams[62])));
  nmpcSine_B.arg258 = nmpcSine_B.dv6[0];
  nmpcSine_B.arg272 = nmpcSine_B.dv6[1];
  obj[5] = ((((((((-(inputVariables[5] * extraParams[42] * extraParams[43]) +
                   inputVariables[6] * extraParams[44] * extraParams[45]) -
                  inputVariables[7] * extraParams[46] * extraParams[47]) +
                 inputVariables[4] * extraParams[48] * extraParams[49] *
                 extraParams[50]) + nmpcSine_B.arg212 * inputVariables[3]) -
               nmpcSine_B.arg222 * inputVariables[3]) + nmpcSine_B.arg232 *
              inputVariables[3] * extraParams[57]) - nmpcSine_B.arg244 *
             inputVariables[3] * extraParams[60]) - nmpcSine_B.dv6[0] *
            inputVariables[9] * extraParams[63]) + nmpcSine_B.dv6[1] *
    inputVariables[8] * extraParams[66];
  memset(&nmpcSine_B.arg478[0], 0, 60U * sizeof(real_T));
  nmpcSine_eye(nmpcSine_B.arg278);
  nmpcSine_B.extraParams = extraParams[66];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg281[nmpcSine_B.i1] = nmpcSine_B.arg278[6 * nmpcSine_B.i1 + 5] *
      nmpcSine_B.extraParams;
  }

  memset(&nmpcSine_B.arg283[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 4; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_mul_pd(_mm_loadu_pd(&nmpcSine_B.arg281[nmpcSine_B.i1]),
                       _mm_set1_pd(nmpcSine_B.arg272));
    _mm_storeu_pd(&nmpcSine_B.dv6[0], tmp_6);
    nmpcSine_B.arg283[(nmpcSine_B.i1 << 2) + 2] = nmpcSine_B.dv6[0];
    nmpcSine_B.arg283[((nmpcSine_B.i1 + 1) << 2) + 2] = nmpcSine_B.dv6[1];
  }

  memset(&nmpcSine_B.arg287[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[8];
  nmpcSine_B.extraParams = extraParams[65];
  nmpcSine_B.extraParams_c = extraParams[64];
  nmpcSine_B.extraParams_k = extraParams[63];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg287[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg272 * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.arg10_tmp;
    nmpcSine_B.arg288_tmp_tmp = -nmpcSine_B.arg278[6 * nmpcSine_B.i1 + 5];
    nmpcSine_B.arg288_tmp[nmpcSine_B.i1] = nmpcSine_B.arg288_tmp_tmp;
    nmpcSine_B.arg281[nmpcSine_B.i1] = nmpcSine_B.arg288_tmp_tmp *
      nmpcSine_B.extraParams_k;
  }

  memset(&nmpcSine_B.arg292[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 4; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_mul_pd(_mm_loadu_pd(&nmpcSine_B.arg281[nmpcSine_B.i1]),
                       _mm_set1_pd(nmpcSine_B.arg258));
    _mm_storeu_pd(&nmpcSine_B.dv6[0], tmp_6);
    nmpcSine_B.arg292[(nmpcSine_B.i1 << 2) + 3] = nmpcSine_B.dv6[0];
    nmpcSine_B.arg292[((nmpcSine_B.i1 + 1) << 2) + 3] = nmpcSine_B.dv6[1];
  }

  nmpcSine_B.arg258 = sin(-inputVariables[2]);
  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[9];
  nmpcSine_B.extraParams = extraParams[62];
  nmpcSine_B.extraParams_c = extraParams[61];
  nmpcSine_B.extraParams_k = extraParams[60];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg272 * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.arg258;
    nmpcSine_B.arg281[nmpcSine_B.i1] = nmpcSine_B.arg288_tmp[nmpcSine_B.i1] *
      nmpcSine_B.extraParams_k;
  }

  memset(&nmpcSine_B.arg301[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg301[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg244;
  }

  memset(&nmpcSine_B.arg305[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[59];
  nmpcSine_B.extraParams_c = extraParams[58];
  nmpcSine_B.extraParams_k = extraParams[57];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg305_tmp = 6 * nmpcSine_B.i1 + 5;
    nmpcSine_B.arg305[nmpcSine_B.arg305_tmp] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg272 * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c;
    nmpcSine_B.arg281[nmpcSine_B.i1] = nmpcSine_B.arg278[nmpcSine_B.arg305_tmp] *
      nmpcSine_B.extraParams_k;
  }

  memset(&nmpcSine_B.arg309[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg309[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg232;
  }

  memset(&nmpcSine_B.arg313[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[56];
  nmpcSine_B.extraParams_c = extraParams[55];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg313[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg272 * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c;
  }

  memset(&nmpcSine_B.arg317[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg317[6 * nmpcSine_B.i1 + 3] =
      nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.arg222;
  }

  memset(&nmpcSine_B.arg321[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[54];
  nmpcSine_B.extraParams_c = extraParams[53];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg321[6 * nmpcSine_B.i1 + 4] =
      nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.arg272 *
      nmpcSine_B.extraParams * nmpcSine_B.extraParams_c;
  }

  memset(&nmpcSine_B.arg324[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg324[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 5] * nmpcSine_B.arg212;
  }

  memset(&nmpcSine_B.arg328[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[52];
  nmpcSine_B.extraParams_c = extraParams[51];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg328[6 * nmpcSine_B.i1 + 4] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 5] * nmpcSine_B.arg272 * nmpcSine_B.extraParams *
      nmpcSine_B.extraParams_c;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    tmp_7 = _mm_loadu_pd(&nmpcSine_B.arg301[nmpcSine_B.i1]);
    tmp = _mm_loadu_pd(&nmpcSine_B.arg305[nmpcSine_B.i1]);
    tmp_0 = _mm_loadu_pd(&nmpcSine_B.arg309[nmpcSine_B.i1]);
    tmp_1 = _mm_loadu_pd(&nmpcSine_B.arg313[nmpcSine_B.i1]);
    tmp_2 = _mm_loadu_pd(&nmpcSine_B.arg317[nmpcSine_B.i1]);
    tmp_3 = _mm_loadu_pd(&nmpcSine_B.arg321[nmpcSine_B.i1]);
    tmp_4 = _mm_loadu_pd(&nmpcSine_B.arg324[nmpcSine_B.i1]);
    tmp_5 = _mm_loadu_pd(&nmpcSine_B.arg328[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(_mm_add_pd
      (_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd
      (_mm_add_pd(tmp_6, tmp_8), tmp_7), tmp), tmp_0), tmp_1), tmp_2), tmp_3),
       tmp_4), tmp_5));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.extraParams = extraParams[50];
  nmpcSine_B.extraParams_c = extraParams[49];
  nmpcSine_B.extraParams_k = extraParams[48];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 5] * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.extraParams_k;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg337[0], 0, 24U * sizeof(real_T));
  nmpcSine_B.extraParams = extraParams[47];
  nmpcSine_B.extraParams_c = extraParams[46];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 4; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_mul_pd(_mm_mul_pd(_mm_loadu_pd
      (&nmpcSine_B.arg288_tmp[nmpcSine_B.i1]), _mm_set1_pd
      (nmpcSine_B.extraParams)), _mm_set1_pd(nmpcSine_B.extraParams_c));
    _mm_storeu_pd(&nmpcSine_B.dv6[0], tmp_6);
    nmpcSine_B.arg337[(nmpcSine_B.i1 << 2) + 1] = nmpcSine_B.dv6[0];
    nmpcSine_B.arg337[((nmpcSine_B.i1 + 1) << 2) + 1] = nmpcSine_B.dv6[1];
  }

  memset(&nmpcSine_B.arg340[0], 0, 24U * sizeof(real_T));
  nmpcSine_B.extraParams = extraParams[45];
  nmpcSine_B.extraParams_c = extraParams[44];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg340[nmpcSine_B.i1 << 2] = nmpcSine_B.arg278[6 * nmpcSine_B.i1
      + 5] * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c;
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.extraParams = extraParams[43];
  nmpcSine_B.extraParams_c = extraParams[42];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg305_tmp = 6 * nmpcSine_B.i1 + 5;
    nmpcSine_B.arg296[nmpcSine_B.arg305_tmp] =
      -nmpcSine_B.arg278[nmpcSine_B.arg305_tmp] * nmpcSine_B.extraParams *
      nmpcSine_B.extraParams_c;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg288_tmp[nmpcSine_B.i1] = -nmpcSine_B.arg278[6 * nmpcSine_B.i1
      + 4];
  }

  memset(&nmpcSine_B.arg349[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 4; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_mul_pd(_mm_loadu_pd(&nmpcSine_B.arg288_tmp[nmpcSine_B.i1]),
                       _mm_set1_pd(nmpcSine_B.arg176));
    _mm_storeu_pd(&nmpcSine_B.dv6[0], tmp_6);
    nmpcSine_B.arg349[(nmpcSine_B.i1 << 2) + 2] = nmpcSine_B.dv6[0];
    nmpcSine_B.arg349[((nmpcSine_B.i1 + 1) << 2) + 2] = nmpcSine_B.dv6[1];
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[8];
  nmpcSine_B.extraParams = extraParams[41];
  nmpcSine_B.extraParams_c = extraParams[40];
  nmpcSine_B.extraParams_k = extraParams[39];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg176 = nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.arg272;
    nmpcSine_B.arg296_tmp[nmpcSine_B.i1] = nmpcSine_B.arg176;
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg176 *
      nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.extraParams_k * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg357[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg357[(nmpcSine_B.i1 << 2) + 3] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 4] * nmpcSine_B.arg163;
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[9];
  nmpcSine_B.extraParams = extraParams[38];
  nmpcSine_B.extraParams_c = extraParams[37];
  nmpcSine_B.extraParams_k = extraParams[36];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg176 = nmpcSine_B.arg278[6 * nmpcSine_B.i1 + 4] *
      nmpcSine_B.arg272;
    nmpcSine_B.arg296_tmp_p[nmpcSine_B.i1] = nmpcSine_B.arg176;
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg176 *
      nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.extraParams_k * nmpcSine_B.arg258;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  nmpcSine_B.extraParams = extraParams[35];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg281[nmpcSine_B.i1] = nmpcSine_B.arg278[6 * nmpcSine_B.i1 + 4] *
      nmpcSine_B.extraParams;
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg150;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[34];
  nmpcSine_B.extraParams_c = extraParams[33];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg272 * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  nmpcSine_B.extraParams = extraParams[32];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 4; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg288_tmp[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg281[nmpcSine_B.i1], _mm_mul_pd(tmp_6,
      _mm_set1_pd(nmpcSine_B.extraParams)));
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg140;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[31];
  nmpcSine_B.extraParams_c = extraParams[30];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg272 * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 4] * nmpcSine_B.arg132;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[29];
  nmpcSine_B.extraParams_c = extraParams[28];
  nmpcSine_B.extraParams_k = extraParams[27];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 4] * nmpcSine_B.arg272 * nmpcSine_B.extraParams *
      nmpcSine_B.extraParams_c * nmpcSine_B.extraParams_k;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 3] =
      nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.arg121;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[3];
  nmpcSine_B.extraParams = extraParams[26];
  nmpcSine_B.extraParams_c = extraParams[25];
  nmpcSine_B.extraParams_k = extraParams[24];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg121 = nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.arg272;
    nmpcSine_B.arg281[nmpcSine_B.i1] = nmpcSine_B.arg121;
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg121 *
      nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.extraParams_k;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg396[0], 0, 24U * sizeof(real_T));
  nmpcSine_B.extraParams = extraParams[23];
  nmpcSine_B.extraParams_c = extraParams[22];
  nmpcSine_B.extraParams_k = extraParams[21];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg396[(nmpcSine_B.i1 << 2) + 1] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 4] * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.extraParams_k;
  }

  memset(&nmpcSine_B.arg401[0], 0, 24U * sizeof(real_T));
  nmpcSine_B.extraParams = extraParams[20];
  nmpcSine_B.extraParams_c = extraParams[19];
  nmpcSine_B.extraParams_k = extraParams[18];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 4; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_mul_pd(_mm_mul_pd(_mm_mul_pd(_mm_loadu_pd
      (&nmpcSine_B.arg288_tmp[nmpcSine_B.i1]), _mm_set1_pd
      (nmpcSine_B.extraParams)), _mm_set1_pd(nmpcSine_B.extraParams_c)),
                       _mm_set1_pd(nmpcSine_B.extraParams_k));
    _mm_storeu_pd(&nmpcSine_B.dv6[0], tmp_6);
    nmpcSine_B.arg401[nmpcSine_B.i1 << 2] = nmpcSine_B.dv6[0];
    nmpcSine_B.arg401[(nmpcSine_B.i1 + 1) << 2] = nmpcSine_B.dv6[1];
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 3] =
      nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.arg96;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.extraParams = extraParams[17];
  nmpcSine_B.extraParams_c = extraParams[16];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.extraParams * nmpcSine_B.extraParams_c;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.extraParams = extraParams[15];
  nmpcSine_B.extraParams_c = extraParams[14];
  nmpcSine_B.extraParams_k = extraParams[13];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 4] * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.extraParams_k;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.extraParams = extraParams[12];
  nmpcSine_B.extraParams_c = extraParams[11];
  nmpcSine_B.extraParams_k = extraParams[10];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] =
      nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.extraParams *
      nmpcSine_B.extraParams_c * nmpcSine_B.extraParams_k;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.extraParams = extraParams[9];
  nmpcSine_B.extraParams_c = extraParams[8];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] =
      nmpcSine_B.arg288_tmp[nmpcSine_B.i1] * nmpcSine_B.extraParams *
      nmpcSine_B.extraParams_c;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg425[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 4; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_mul_pd(_mm_loadu_pd(&nmpcSine_B.arg288_tmp[nmpcSine_B.i1]),
                       _mm_set1_pd(nmpcSine_B.arg6_tmp));
    _mm_storeu_pd(&nmpcSine_B.dv6[0], tmp_6);
    nmpcSine_B.arg425[(nmpcSine_B.i1 << 2) + 2] = nmpcSine_B.dv6[0];
    nmpcSine_B.arg425[((nmpcSine_B.i1 + 1) << 2) + 2] = nmpcSine_B.dv6[1];
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] =
      nmpcSine_B.arg296_tmp[nmpcSine_B.i1] * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg430[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg430[(nmpcSine_B.i1 << 2) + 3] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 4] * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 22; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg292[nmpcSine_B.i1]);
    tmp_7 = _mm_loadu_pd(&nmpcSine_B.arg337[nmpcSine_B.i1]);
    tmp = _mm_loadu_pd(&nmpcSine_B.arg340[nmpcSine_B.i1]);
    tmp_0 = _mm_loadu_pd(&nmpcSine_B.arg349[nmpcSine_B.i1]);
    tmp_1 = _mm_loadu_pd(&nmpcSine_B.arg357[nmpcSine_B.i1]);
    tmp_2 = _mm_loadu_pd(&nmpcSine_B.arg396[nmpcSine_B.i1]);
    tmp_3 = _mm_loadu_pd(&nmpcSine_B.arg401[nmpcSine_B.i1]);
    tmp_4 = _mm_loadu_pd(&nmpcSine_B.arg425[nmpcSine_B.i1]);
    tmp_5 = _mm_loadu_pd(&nmpcSine_B.arg430[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1], _mm_add_pd(_mm_add_pd
      (_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd
      (_mm_add_pd(tmp_6, tmp_8), tmp_7), tmp), tmp_0), tmp_1), tmp_2), tmp_3),
       tmp_4), tmp_5));
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] =
      nmpcSine_B.arg296_tmp_p[nmpcSine_B.i1] * nmpcSine_B.arg258;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.extraParams = extraParams[7];
  nmpcSine_B.extraParams_c = extraParams[6];
  nmpcSine_B.extraParams_k = extraParams[5];
  nmpcSine_B.arg272 = inputVariables[5];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 3] * nmpcSine_B.extraParams * nmpcSine_B.extraParams_c *
      nmpcSine_B.extraParams_k * 2.0 * nmpcSine_B.arg272;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  nmpcSine_B.extraParams = extraParams[4];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg281[nmpcSine_B.i1] = nmpcSine_B.arg278[6 * nmpcSine_B.i1 + 3] *
      nmpcSine_B.extraParams;
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg45;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[5];
  nmpcSine_B.extraParams = extraParams[3];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] = nmpcSine_B.arg281[nmpcSine_B.i1] *
      nmpcSine_B.arg272 * nmpcSine_B.extraParams;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg292[0], 0, 24U * sizeof(real_T));
  nmpcSine_B.extraParams = extraParams[2];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg292[(nmpcSine_B.i1 << 2) + 1] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 3] * nmpcSine_B.extraParams;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 22; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg292[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg292[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.extraParams = extraParams[1];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg292[nmpcSine_B.i1 << 2] = nmpcSine_B.arg278[6 * nmpcSine_B.i1
      + 3] * nmpcSine_B.extraParams;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 22; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg292[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.extraParams = extraParams[0];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg305_tmp = 6 * nmpcSine_B.i1 + 3;
    nmpcSine_B.arg296[nmpcSine_B.arg305_tmp] =
      -nmpcSine_B.arg278[nmpcSine_B.arg305_tmp] * nmpcSine_B.extraParams;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg292[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg292[(nmpcSine_B.i1 << 2) + 3] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 3] * nmpcSine_B.arg6_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 22; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg292[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg283[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[9];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 3] * nmpcSine_B.arg272 * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
  }

  memset(&nmpcSine_B.arg292[0], 0, 24U * sizeof(real_T));
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg305_tmp = nmpcSine_B.i1 << 2;
    nmpcSine_B.arg45 = nmpcSine_B.arg278[6 * nmpcSine_B.i1 + 3] *
      nmpcSine_B.arg10_tmp;
    nmpcSine_B.arg292[nmpcSine_B.arg305_tmp + 2] = nmpcSine_B.arg45;
    tmp_6 = _mm_add_pd(_mm_loadu_pd(&nmpcSine_B.arg283[nmpcSine_B.arg305_tmp]),
                       _mm_loadu_pd(&nmpcSine_B.arg292[nmpcSine_B.arg305_tmp]));
    _mm_storeu_pd(&nmpcSine_B.arg478[10 * nmpcSine_B.i1 + 6], tmp_6);
    nmpcSine_B.arg478[10 * nmpcSine_B.i1 + 8] =
      nmpcSine_B.arg283[nmpcSine_B.arg305_tmp + 2] + nmpcSine_B.arg45;
    nmpcSine_B.arg478[10 * nmpcSine_B.i1 + 9] =
      nmpcSine_B.arg283[nmpcSine_B.arg305_tmp + 3] +
      nmpcSine_B.arg292[nmpcSine_B.arg305_tmp + 3];
  }

  memset(&nmpcSine_B.arg296[0], 0, 36U * sizeof(real_T));
  nmpcSine_B.arg272 = inputVariables[8];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 3] * nmpcSine_B.arg272 * nmpcSine_B.arg258;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 5] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 2];
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 1] * nmpcSine_B.arg6_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[3];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 1] * nmpcSine_B.arg272 * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 1] * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[4];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1 + 1] * nmpcSine_B.arg272 * nmpcSine_B.arg258;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 3] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1] * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = inputVariables[3];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1] * nmpcSine_B.arg272 * nmpcSine_B.arg258;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  nmpcSine_B.arg272 = -inputVariables[4];
  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 2] = nmpcSine_B.arg278[6 *
      nmpcSine_B.i1] * nmpcSine_B.arg272 * nmpcSine_B.arg10_tmp;
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 <= 34; nmpcSine_B.i1 += 2) {
    tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1]);
    tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1]);
    _mm_storeu_pd(&nmpcSine_B.arg287[nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    _mm_storeu_pd(&nmpcSine_B.arg296[nmpcSine_B.i1], _mm_set1_pd(0.0));
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    nmpcSine_B.arg296[6 * nmpcSine_B.i1 + 4] = -(nmpcSine_B.arg278[6 *
      nmpcSine_B.i1] * nmpcSine_B.arg6_tmp);
    for (nmpcSine_B.arg305_tmp = 0; nmpcSine_B.arg305_tmp <= 4;
         nmpcSine_B.arg305_tmp += 2) {
      nmpcSine_B.i2 = 6 * nmpcSine_B.i1 + nmpcSine_B.arg305_tmp;
      tmp_6 = _mm_loadu_pd(&nmpcSine_B.arg287[nmpcSine_B.i2]);
      tmp_8 = _mm_loadu_pd(&nmpcSine_B.arg296[nmpcSine_B.i2]);
      _mm_storeu_pd(&nmpcSine_B.arg478[nmpcSine_B.arg305_tmp + 10 *
                    nmpcSine_B.i1], _mm_add_pd(tmp_6, tmp_8));
    }
  }

  for (nmpcSine_B.i1 = 0; nmpcSine_B.i1 < 6; nmpcSine_B.i1++) {
    for (nmpcSine_B.arg305_tmp = 0; nmpcSine_B.arg305_tmp < 10;
         nmpcSine_B.arg305_tmp++) {
      grad[nmpcSine_B.i1 + 6 * nmpcSine_B.arg305_tmp] = nmpcSine_B.arg478[10 *
        nmpcSine_B.i1 + nmpcSine_B.arg305_tmp];
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_stateEvolution(const real_T X[126], const real_T U[84],
  real_T c[120], real_T J[15000])
{
  __m128d tmp;
  static const real_T d[67] = { 0.80982825268897241, 0.079415032681564857,
    0.079415032681564857, 63.0144902620752, 0.022897834317090127,
    63.0144902620752, 0.022897834317090127, 0.699999982041001, 74.8832173268989,
    0.015869365852854363, 74.8832173268989, 0.48999997485740171,
    0.03333333342744111, 125.558050313511, 0.03333333342744111,
    0.699999982041001, 43.6722524126937, 0.015869365852854363,
    0.11560777874744173, 0.699999982041001, 0.415, 0.11560777874744173,
    0.699999982041001, 0.415, 43.6722524126937, 0.48999997485740171,
    0.03333333342744111, 63.0144902620752, 0.48999997485740171,
    0.03333333342744111, 43.6722524126937, 0.03333333342744111,
    0.699999982041001, 63.0144902620752, 0.03333333342744111, 0.699999982041001,
    63.0144902620752, 0.48999997485740171, 0.03333333342744111, 63.0144902620752,
    0.48999997485740171, 0.03333333342744111, 125.558050313511,
    0.03333333342744111, 0.11560777874744173, 0.415, 0.11560777874744173, 0.415,
    74.8832173268989, 0.03333333342744111, 0.699999982041001, 43.6722524126937,
    0.03333333342744111, 63.0144902620752, 0.03333333342744111, 43.6722524126937,
    0.03333333342744111, 0.699999982041001, 63.0144902620752,
    0.03333333342744111, 0.699999982041001, 63.0144902620752,
    0.03333333342744111, 0.699999982041001, 63.0144902620752,
    0.03333333342744111, 0.699999982041001 };

  static const int8_T b[160] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
    1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0,
    0, 0, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1,
    0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1, 0, 1 };

  memset(&nmpcSine_B.Jx_c[0], 0, 14400U * sizeof(real_T));
  memset(&nmpcSine_B.Jmv[0], 0, 4800U * sizeof(real_T));
  memset(&c[0], 0, 120U * sizeof(real_T));
  for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 6; nmpcSine_B.Jx_tmp++) {
    nmpcSine_B.ic[nmpcSine_B.Jx_tmp] = static_cast<real_T>(nmpcSine_B.Jx_tmp) +
      1.0;
  }

  for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 21; nmpcSine_B.Jx_tmp++) {
    nmpcSine_B.b_U_tmp = nmpcSine_B.Jx_tmp << 2;
    nmpcSine_B.b_U[nmpcSine_B.b_U_tmp] = U[nmpcSine_B.Jx_tmp];
    nmpcSine_B.b_U[nmpcSine_B.b_U_tmp + 1] = U[nmpcSine_B.Jx_tmp + 21];
    nmpcSine_B.b_U[nmpcSine_B.b_U_tmp + 2] = U[nmpcSine_B.Jx_tmp + 42];
    nmpcSine_B.b_U[nmpcSine_B.b_U_tmp + 3] = U[nmpcSine_B.Jx_tmp + 63];
  }

  for (nmpcSine_B.i_c = 0; nmpcSine_B.i_c < 6; nmpcSine_B.i_c++) {
    for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 21; nmpcSine_B.Jx_tmp++) {
      nmpcSine_B.b_X[nmpcSine_B.i_c + 6 * nmpcSine_B.Jx_tmp] = X[21 *
        nmpcSine_B.i_c + nmpcSine_B.Jx_tmp];
    }
  }

  for (nmpcSine_B.i_c = 0; nmpcSine_B.i_c < 20; nmpcSine_B.i_c++) {
    nmpcSine_B.b_X_mv = nmpcSine_B.b_X[6 * nmpcSine_B.i_c + 2];
    nmpcSine_B.t2 = cos(nmpcSine_B.b_X_mv);
    nmpcSine_B.t3 = sin(nmpcSine_B.b_X_mv);
    if (!nmpcSine_DW.ADdata_not_empty) {
      memcpy(&nmpcSine_DW.ADdata.constants[0], &d[0], 67U * sizeof(real_T));
      nmpcSine_DW.ADdata_not_empty = true;
    }

    for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 6; nmpcSine_B.Jx_tmp++) {
      nmpcSine_B.b_X_m[nmpcSine_B.Jx_tmp] = nmpcSine_B.b_X[6 * nmpcSine_B.i_c +
        nmpcSine_B.Jx_tmp];
    }

    nmpcSine_B.b_U_tmp = nmpcSine_B.i_c << 2;
    nmpcSine_B.b_U_f = nmpcSine_B.b_U[nmpcSine_B.b_U_tmp];
    nmpcSine_B.b_X_m[6] = nmpcSine_B.b_U_f;
    nmpcSine_B.b_U_p = nmpcSine_B.b_U[nmpcSine_B.b_U_tmp + 1];
    nmpcSine_B.b_X_m[7] = nmpcSine_B.b_U_p;
    nmpcSine_B.b_U_e = nmpcSine_B.b_U[nmpcSine_B.b_U_tmp + 2];
    nmpcSine_B.b_X_m[8] = nmpcSine_B.b_U_e;
    nmpcSine_B.b_U_o = nmpcSine_B.b_U[nmpcSine_B.b_U_tmp + 3];
    nmpcSine_B.b_X_m[9] = nmpcSine_B.b_U_o;
    nm_stateTransitionFcnJacobianAD(nmpcSine_B.b_X_m,
      nmpcSine_DW.ADdata.constants, nmpcSine_B.a__4, nmpcSine_B.b_J);
    nmpcSine_B.b_U_tmp = (nmpcSine_B.i_c + 1) * 6;
    nmpcSine_B.b_X_mv = nmpcSine_B.b_X[nmpcSine_B.b_U_tmp + 2];
    nmpcSine_B.b_t2 = cos(nmpcSine_B.b_X_mv);
    nmpcSine_B.b_t3 = sin(nmpcSine_B.b_X_mv);
    nmpcSine_B.b_X_m[6] = nmpcSine_B.b_U_f;
    nmpcSine_B.b_X_m[7] = nmpcSine_B.b_U_p;
    nmpcSine_B.b_X_m[8] = nmpcSine_B.b_U_e;
    nmpcSine_B.b_X_m[9] = nmpcSine_B.b_U_o;
    nmpcSine_B.b_X_mv = nmpcSine_B.b_X[6 * nmpcSine_B.i_c + 4];
    nmpcSine_B.b_X_mj = nmpcSine_B.b_X[6 * nmpcSine_B.i_c + 3];
    nmpcSine_B.a__4[0] = -nmpcSine_B.b_X_mv * nmpcSine_B.t3 + nmpcSine_B.t2 *
      nmpcSine_B.b_X_mj;
    nmpcSine_B.a__4[1] = nmpcSine_B.b_X_mv * nmpcSine_B.t2 + nmpcSine_B.t3 *
      nmpcSine_B.b_X_mj;
    nmpcSine_B.b_X_c = nmpcSine_B.b_X[6 * nmpcSine_B.i_c + 5];
    nmpcSine_B.a__4[2] = nmpcSine_B.b_X_c;
    nmpcSine_B.a__4[3] = (((((nmpcSine_B.t2 * nmpcSine_B.b_U_e + nmpcSine_B.t3 *
      nmpcSine_B.b_U_o) - 0.80982825268897241 * nmpcSine_B.b_X_mj) +
      0.079415032681564857 * nmpcSine_B.b_U_f) + 0.079415032681564857 *
      nmpcSine_B.b_U_p) + 63.0144902620752 * nmpcSine_B.b_X_mv *
                          nmpcSine_B.b_X_c * 0.022897834317090127) +
      nmpcSine_B.b_X_c * nmpcSine_B.b_X_c * 63.0144902620752 *
      0.022897834317090127 * 0.699999982041001;
    nmpcSine_B.b_X_tmp = 125.558050313511 * nmpcSine_B.b_X_c *
      0.03333333342744111;
    nmpcSine_B.b_X_tmp_h = 43.6722524126937 * nmpcSine_B.b_X_mv *
      0.03333333342744111 * nmpcSine_B.b_X_mj;
    nmpcSine_B.b_X_tmp_l = 63.0144902620752 * nmpcSine_B.b_X_mv *
      0.03333333342744111 * nmpcSine_B.b_X_mj;
    nmpcSine_B.b_X_tmp_h2 = 0.11560777874744173 * nmpcSine_B.b_U_f *
      0.699999982041001 * 0.415;
    nmpcSine_B.b_X_tmp_m = 0.11560777874744173 * nmpcSine_B.b_U_p *
      0.699999982041001 * 0.415;
    nmpcSine_B.a__4[4] = ((((((((((((nmpcSine_B.t2 * nmpcSine_B.b_U_o -
      nmpcSine_B.t3 * nmpcSine_B.b_U_e) - 74.8832173268989 * nmpcSine_B.b_X_mv *
      0.015869365852854363) - 74.8832173268989 * nmpcSine_B.b_X_mv *
      0.48999997485740171 * 0.03333333342744111) + nmpcSine_B.b_X_tmp *
      0.699999982041001) - 43.6722524126937 * nmpcSine_B.b_X_c *
      0.015869365852854363 * nmpcSine_B.b_X_mj) - nmpcSine_B.b_X_tmp_h2) +
      nmpcSine_B.b_X_tmp_m) - 43.6722524126937 * nmpcSine_B.b_X_c *
      0.48999997485740171 * 0.03333333342744111 * nmpcSine_B.b_X_mj) +
      63.0144902620752 * nmpcSine_B.b_X_c * 0.48999997485740171 *
      0.03333333342744111 * nmpcSine_B.b_X_mj) - nmpcSine_B.b_X_tmp_h *
      0.699999982041001) + nmpcSine_B.b_X_tmp_l * 0.699999982041001) +
                          63.0144902620752 * nmpcSine_B.t2 * 0.48999997485740171
                          * 0.03333333342744111 * nmpcSine_B.b_U_o) -
      63.0144902620752 * nmpcSine_B.t3 * 0.48999997485740171 *
      0.03333333342744111 * nmpcSine_B.b_U_e;
    nmpcSine_B.b_X_tmp_mc = 0.11560777874744173 * nmpcSine_B.b_U_f * 0.415;
    nmpcSine_B.b_X_tmp_h3 = 0.11560777874744173 * nmpcSine_B.b_U_p * 0.415;
    nmpcSine_B.a__4[5] = ((((((74.8832173268989 * nmpcSine_B.b_X_mv *
      0.03333333342744111 * 0.699999982041001 + ((-nmpcSine_B.b_X_tmp +
      nmpcSine_B.b_X_tmp_mc) - nmpcSine_B.b_X_tmp_h3)) + nmpcSine_B.b_X_tmp_h) -
      nmpcSine_B.b_X_tmp_l) + 43.6722524126937 * nmpcSine_B.b_X_c *
      0.03333333342744111 * nmpcSine_B.b_X_mj * 0.699999982041001) -
      63.0144902620752 * nmpcSine_B.b_X_c * 0.03333333342744111 *
      nmpcSine_B.b_X_mj * 0.699999982041001) - 63.0144902620752 * nmpcSine_B.t2 *
                          0.03333333342744111 * nmpcSine_B.b_U_o *
                          0.699999982041001) + 63.0144902620752 * nmpcSine_B.t3 *
      0.03333333342744111 * nmpcSine_B.b_U_e * 0.699999982041001;
    nmpcSine_B.b_X_mv = nmpcSine_B.b_X[nmpcSine_B.b_U_tmp + 4];
    nmpcSine_B.b_X_mj = nmpcSine_B.b_X[nmpcSine_B.b_U_tmp + 3];
    nmpcSine_B.b_X_p[0] = -nmpcSine_B.b_X_mv * nmpcSine_B.b_t3 + nmpcSine_B.b_t2
      * nmpcSine_B.b_X_mj;
    nmpcSine_B.b_X_p[1] = nmpcSine_B.b_X_mv * nmpcSine_B.b_t2 + nmpcSine_B.b_t3 *
      nmpcSine_B.b_X_mj;
    nmpcSine_B.b_X_c = nmpcSine_B.b_X[nmpcSine_B.b_U_tmp + 5];
    nmpcSine_B.b_X_p[2] = nmpcSine_B.b_X_c;
    nmpcSine_B.b_X_p[3] = (((((nmpcSine_B.b_t2 * nmpcSine_B.b_U_e +
      nmpcSine_B.b_t3 * nmpcSine_B.b_U_o) - 0.80982825268897241 *
      nmpcSine_B.b_X_mj) + 0.079415032681564857 * nmpcSine_B.b_U_f) +
      0.079415032681564857 * nmpcSine_B.b_U_p) + 63.0144902620752 *
      nmpcSine_B.b_X_mv * nmpcSine_B.b_X_c * 0.022897834317090127) +
      nmpcSine_B.b_X_c * nmpcSine_B.b_X_c * 63.0144902620752 *
      0.022897834317090127 * 0.699999982041001;
    nmpcSine_B.b_X_tmp = 125.558050313511 * nmpcSine_B.b_X_c *
      0.03333333342744111;
    nmpcSine_B.b_X_tmp_h = 43.6722524126937 * nmpcSine_B.b_X_mv *
      0.03333333342744111 * nmpcSine_B.b_X_mj;
    nmpcSine_B.b_X_tmp_l = 63.0144902620752 * nmpcSine_B.b_X_mv *
      0.03333333342744111 * nmpcSine_B.b_X_mj;
    nmpcSine_B.b_X_p[4] = ((((((((((((nmpcSine_B.b_t2 * nmpcSine_B.b_U_o -
      nmpcSine_B.b_t3 * nmpcSine_B.b_U_e) - 74.8832173268989 * nmpcSine_B.b_X_mv
      * 0.015869365852854363) - 74.8832173268989 * nmpcSine_B.b_X_mv *
      0.48999997485740171 * 0.03333333342744111) + nmpcSine_B.b_X_tmp *
      0.699999982041001) - 43.6722524126937 * nmpcSine_B.b_X_c *
      0.015869365852854363 * nmpcSine_B.b_X_mj) - nmpcSine_B.b_X_tmp_h2) +
      nmpcSine_B.b_X_tmp_m) - 43.6722524126937 * nmpcSine_B.b_X_c *
      0.48999997485740171 * 0.03333333342744111 * nmpcSine_B.b_X_mj) +
      63.0144902620752 * nmpcSine_B.b_X_c * 0.48999997485740171 *
      0.03333333342744111 * nmpcSine_B.b_X_mj) - nmpcSine_B.b_X_tmp_h *
      0.699999982041001) + nmpcSine_B.b_X_tmp_l * 0.699999982041001) +
      63.0144902620752 * nmpcSine_B.b_t2 * 0.48999997485740171 *
      0.03333333342744111 * nmpcSine_B.b_U_o) - 63.0144902620752 *
      nmpcSine_B.b_t3 * 0.48999997485740171 * 0.03333333342744111 *
      nmpcSine_B.b_U_e;
    nmpcSine_B.b_X_p[5] = ((((((74.8832173268989 * nmpcSine_B.b_X_mv *
      0.03333333342744111 * 0.699999982041001 + ((-nmpcSine_B.b_X_tmp +
      nmpcSine_B.b_X_tmp_mc) - nmpcSine_B.b_X_tmp_h3)) + nmpcSine_B.b_X_tmp_h) -
      nmpcSine_B.b_X_tmp_l) + 43.6722524126937 * nmpcSine_B.b_X_c *
      0.03333333342744111 * nmpcSine_B.b_X_mj * 0.699999982041001) -
      63.0144902620752 * nmpcSine_B.b_X_c * 0.03333333342744111 *
      nmpcSine_B.b_X_mj * 0.699999982041001) - 63.0144902620752 *
      nmpcSine_B.b_t2 * 0.03333333342744111 * nmpcSine_B.b_U_o *
      0.699999982041001) + 63.0144902620752 * nmpcSine_B.b_t3 *
      0.03333333342744111 * nmpcSine_B.b_U_e * 0.699999982041001;
    for (nmpcSine_B.k_h = 0; nmpcSine_B.k_h < 6; nmpcSine_B.k_h++) {
      nmpcSine_B.b_X_mv = nmpcSine_B.b_X[nmpcSine_B.b_U_tmp + nmpcSine_B.k_h];
      nmpcSine_B.b_X_m[nmpcSine_B.k_h] = nmpcSine_B.b_X_mv;
      nmpcSine_B.c_tmp = static_cast<int32_T>(nmpcSine_B.ic[nmpcSine_B.k_h]);
      c[nmpcSine_B.c_tmp - 1] = (nmpcSine_B.b_X[6 * nmpcSine_B.i_c +
        nmpcSine_B.k_h] + (nmpcSine_B.a__4[nmpcSine_B.k_h] +
                           nmpcSine_B.b_X_p[nmpcSine_B.k_h]) * 0.05) -
        nmpcSine_B.b_X_mv;
      if (nmpcSine_B.i_c + 1 > 1) {
        for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 6; nmpcSine_B.Jx_tmp++)
        {
          nmpcSine_B.Jx_c[((static_cast<int32_T>(nmpcSine_B.ic[nmpcSine_B.Jx_tmp])
                            + 120 * nmpcSine_B.k_h) + 720 * (nmpcSine_B.i_c - 1))
            - 1] = nmpcSine_B.b_J[6 * nmpcSine_B.k_h + nmpcSine_B.Jx_tmp] * 0.05;
        }

        nmpcSine_B.Jx_tmp = ((120 * nmpcSine_B.k_h + nmpcSine_B.c_tmp) +
                             (nmpcSine_B.i_c - 1) * 720) - 1;
        nmpcSine_B.Jx_c[nmpcSine_B.Jx_tmp]++;
      }
    }

    nm_stateTransitionFcnJacobianAD(nmpcSine_B.b_X_m,
      nmpcSine_DW.ADdata.constants, nmpcSine_B.a__4, nmpcSine_B.c_J);
    for (nmpcSine_B.b_U_tmp = 0; nmpcSine_B.b_U_tmp < 6; nmpcSine_B.b_U_tmp++) {
      for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 6; nmpcSine_B.Jx_tmp++) {
        nmpcSine_B.Jx_c[((static_cast<int32_T>(nmpcSine_B.ic[nmpcSine_B.Jx_tmp])
                          + 120 * nmpcSine_B.b_U_tmp) + 720 * nmpcSine_B.i_c) -
          1] = nmpcSine_B.c_J[6 * nmpcSine_B.b_U_tmp + nmpcSine_B.Jx_tmp] * 0.05;
      }

      nmpcSine_B.Jx_tmp = ((120 * nmpcSine_B.b_U_tmp + static_cast<int32_T>
                            (nmpcSine_B.ic[nmpcSine_B.b_U_tmp])) + 720 *
                           nmpcSine_B.i_c) - 1;
      nmpcSine_B.Jx_c[nmpcSine_B.Jx_tmp]--;
    }

    for (nmpcSine_B.b_U_tmp = 0; nmpcSine_B.b_U_tmp < 2; nmpcSine_B.b_U_tmp++) {
      for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 6; nmpcSine_B.Jx_tmp++) {
        nmpcSine_B.k_h = (nmpcSine_B.b_U_tmp + 6) * 6 + nmpcSine_B.Jx_tmp;
        nmpcSine_B.Jmv[((static_cast<int32_T>(nmpcSine_B.ic[nmpcSine_B.Jx_tmp])
                         + 120 * nmpcSine_B.b_U_tmp) + 240 * nmpcSine_B.i_c) - 1]
          = (nmpcSine_B.b_J[nmpcSine_B.k_h] + nmpcSine_B.c_J[nmpcSine_B.k_h]) *
          0.05;
      }
    }

    for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp <= 4; nmpcSine_B.Jx_tmp += 2)
    {
      tmp = _mm_loadu_pd(&nmpcSine_B.ic[nmpcSine_B.Jx_tmp]);
      _mm_storeu_pd(&nmpcSine_B.ic[nmpcSine_B.Jx_tmp], _mm_add_pd(tmp,
        _mm_set1_pd(6.0)));
    }
  }

  for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 4; nmpcSine_B.Jx_tmp++) {
    for (nmpcSine_B.i_c = 0; nmpcSine_B.i_c < 120; nmpcSine_B.i_c++) {
      nmpcSine_B.b_X_mv = 0.0;
      for (nmpcSine_B.b_U_tmp = 0; nmpcSine_B.b_U_tmp < 40; nmpcSine_B.b_U_tmp++)
      {
        nmpcSine_B.b_X_mv += nmpcSine_B.Jmv[120 * nmpcSine_B.b_U_tmp +
          nmpcSine_B.i_c] * static_cast<real_T>(b[40 * nmpcSine_B.Jx_tmp +
          nmpcSine_B.b_U_tmp]);
      }

      nmpcSine_B.Jmv_n[nmpcSine_B.Jx_tmp + (nmpcSine_B.i_c << 2)] =
        nmpcSine_B.b_X_mv;
    }
  }

  for (nmpcSine_B.Jx_tmp = 0; nmpcSine_B.Jx_tmp < 120; nmpcSine_B.Jx_tmp++) {
    for (nmpcSine_B.i_c = 0; nmpcSine_B.i_c < 120; nmpcSine_B.i_c++) {
      J[nmpcSine_B.i_c + 125 * nmpcSine_B.Jx_tmp] = nmpcSine_B.Jx_c[120 *
        nmpcSine_B.i_c + nmpcSine_B.Jx_tmp];
    }

    nmpcSine_B.i_c = nmpcSine_B.Jx_tmp << 2;
    J[125 * nmpcSine_B.Jx_tmp + 120] = nmpcSine_B.Jmv_n[nmpcSine_B.i_c];
    J[125 * nmpcSine_B.Jx_tmp + 121] = nmpcSine_B.Jmv_n[nmpcSine_B.i_c + 1];
    J[125 * nmpcSine_B.Jx_tmp + 122] = nmpcSine_B.Jmv_n[nmpcSine_B.i_c + 2];
    J[125 * nmpcSine_B.Jx_tmp + 123] = nmpcSine_B.Jmv_n[nmpcSine_B.i_c + 3];
    J[125 * nmpcSine_B.Jx_tmp + 124] = 0.0;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_all(const boolean_T x[80], boolean_T y[4])
{
  int32_T i2;
  y[0] = true;
  y[1] = true;
  y[2] = true;
  y[3] = true;
  i2 = 1;
  for (int32_T i = 0; i < 4; i++) {
    int32_T a;
    int32_T ix;
    boolean_T exitg1;
    a = i2 + 19;
    ix = i2;
    i2 += 20;
    exitg1 = false;
    while ((!exitg1) && (ix <= a)) {
      if (!x[ix - 1]) {
        y[i] = false;
        exitg1 = true;
      } else {
        ix++;
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
boolean_T nmpcSine::nmpcSine_any(const boolean_T x[8])
{
  int32_T k;
  boolean_T exitg1;
  boolean_T y;
  y = false;
  k = 0;
  exitg1 = false;
  while ((!exitg1) && (k <= 7)) {
    if (x[k]) {
      y = true;
      exitg1 = true;
    } else {
      k++;
    }
  }

  return y;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_reformJacobian(const real_T Jx_data[], const int32_T
  Jx_size[3], const real_T Jmv_data[], const real_T Je_data[], const int32_T
  Je_size[1], real_T Jc_data[], int32_T Jc_size[2])
{
  if (Jx_size[0] == 0) {
    Jc_size[0] = 0;
    Jc_size[1] = 0;
  } else {
    nmpcSine_B.Jx_e[0] = static_cast<uint8_T>(Jx_size[0]);
    nmpcSine_B.loop_ub_f = Jx_size[0];
    for (nmpcSine_B.i6 = 0; nmpcSine_B.i6 < nmpcSine_B.loop_ub_f; nmpcSine_B.i6
         ++) {
      for (nmpcSine_B.i5 = 0; nmpcSine_B.i5 < 120; nmpcSine_B.i5++) {
        nmpcSine_B.varargin_1_data_m[nmpcSine_B.i5 + 120 * nmpcSine_B.i6] =
          Jx_data[nmpcSine_B.Jx_e[0] * nmpcSine_B.i5 + nmpcSine_B.i6];
      }
    }

    nmpcSine_B.Jx_f[0] = Jx_size[0];
    nmpcSine_B.Jx_f[1] = 40;
    nmpcSine_mtimes(Jmv_data, nmpcSine_B.Jx_f, nmpcSine_B.tmp_data_m,
                    nmpcSine_B.tmp_size);
    nmpcSine_B.loop_ub_cm = nmpcSine_B.tmp_size[0];
    for (nmpcSine_B.i6 = 0; nmpcSine_B.i6 < nmpcSine_B.loop_ub_cm; nmpcSine_B.i6
         ++) {
      nmpcSine_B.varargin_2_data[nmpcSine_B.i6 << 2] =
        nmpcSine_B.tmp_data_m[nmpcSine_B.i6];
      nmpcSine_B.varargin_2_data[1 + (nmpcSine_B.i6 << 2)] =
        nmpcSine_B.tmp_data_m[nmpcSine_B.i6 + nmpcSine_B.tmp_size[0]];
      nmpcSine_B.varargin_2_data[2 + (nmpcSine_B.i6 << 2)] =
        nmpcSine_B.tmp_data_m[(nmpcSine_B.tmp_size[0] << 1) + nmpcSine_B.i6];
      nmpcSine_B.varargin_2_data[3 + (nmpcSine_B.i6 << 2)] =
        nmpcSine_B.tmp_data_m[nmpcSine_B.tmp_size[0] * 3 + nmpcSine_B.i6];
    }

    nmpcSine_B.Jx_e[0] = 120U;
    if (nmpcSine_B.tmp_size[0] != 0) {
      nmpcSine_B.varargin_2[0] = 4U;
    } else {
      nmpcSine_B.varargin_2[0] = 0U;
    }

    nmpcSine_B.loop_ub_cm = Je_size[0];
    if (nmpcSine_B.loop_ub_cm - 1 >= 0) {
      memcpy(&nmpcSine_B.Je_data[0], &Je_data[0], static_cast<uint32_T>
             (nmpcSine_B.loop_ub_cm) * sizeof(real_T));
    }

    nmpcSine_B.Je_n[0] = (Je_size[0] != 0);
    Jc_size[0] = (nmpcSine_B.varargin_2[0] + nmpcSine_B.Je_n[0]) + 120;
    Jc_size[1] = Jx_size[0];
    nmpcSine_B.loop_ub_cm = nmpcSine_B.varargin_2[0];
    nmpcSine_B.loop_ub_jk = nmpcSine_B.Je_n[0];
    for (nmpcSine_B.i6 = 0; nmpcSine_B.i6 < nmpcSine_B.loop_ub_f; nmpcSine_B.i6
         ++) {
      for (nmpcSine_B.i5 = 0; nmpcSine_B.i5 < 120; nmpcSine_B.i5++) {
        Jc_data[nmpcSine_B.i5 + Jc_size[0] * nmpcSine_B.i6] =
          nmpcSine_B.varargin_1_data_m[nmpcSine_B.Jx_e[0] * nmpcSine_B.i6 +
          nmpcSine_B.i5];
      }

      for (nmpcSine_B.i5 = 0; nmpcSine_B.i5 < nmpcSine_B.loop_ub_cm;
           nmpcSine_B.i5++) {
        Jc_data[(nmpcSine_B.i5 + Jc_size[0] * nmpcSine_B.i6) + 120] =
          nmpcSine_B.varargin_2_data[nmpcSine_B.varargin_2[0] * nmpcSine_B.i6 +
          nmpcSine_B.i5];
      }

      if (nmpcSine_B.loop_ub_jk - 1 >= 0) {
        Jc_data[(nmpcSine_B.varargin_2[0] + Jc_size[0] * nmpcSine_B.i6) + 120] =
          nmpcSine_B.Je_data[nmpcSine_B.Je_n[0] * nmpcSine_B.i6];
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_outputBounds(const real_T runtimedata_OutputMin[80],
  const real_T runtimedata_OutputMax[80], const real_T X[126], real_T e, real_T
  c_data[], int32_T c_size[2], real_T Jc_data[], int32_T Jc_size[2])
{
  static const int8_T val[24] = { -1, 0, 0, 0, 0, -1, 0, 0, 0, 0, -1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, -1 };

  static const int8_T d[24] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1 };

  boolean_T exitg1;
  for (nmpcSine_B.c_k = 0; nmpcSine_B.c_k < 80; nmpcSine_B.c_k++) {
    nmpcSine_B.bv1[nmpcSine_B.c_k] = rtIsInf
      (runtimedata_OutputMin[nmpcSine_B.c_k]);
  }

  nmpcSine_all(nmpcSine_B.bv1, nmpcSine_B.x_e);
  nmpcSine_B.y_e = true;
  nmpcSine_B.c_k = 0;
  exitg1 = false;
  while ((!exitg1) && (nmpcSine_B.c_k < 4)) {
    if (!nmpcSine_B.x_e[nmpcSine_B.c_k]) {
      nmpcSine_B.y_e = false;
      exitg1 = true;
    } else {
      nmpcSine_B.c_k++;
    }
  }

  if (nmpcSine_B.y_e) {
    for (nmpcSine_B.c_k = 0; nmpcSine_B.c_k < 80; nmpcSine_B.c_k++) {
      nmpcSine_B.bv1[nmpcSine_B.c_k] = rtIsInf
        (runtimedata_OutputMax[nmpcSine_B.c_k]);
    }

    nmpcSine_all(nmpcSine_B.bv1, nmpcSine_B.x_e);
    nmpcSine_B.y_e = true;
    nmpcSine_B.c_k = 0;
    exitg1 = false;
    while ((!exitg1) && (nmpcSine_B.c_k < 4)) {
      if (!nmpcSine_B.x_e[nmpcSine_B.c_k]) {
        nmpcSine_B.y_e = false;
        exitg1 = true;
      } else {
        nmpcSine_B.c_k++;
      }
    }
  } else {
    nmpcSine_B.y_e = false;
  }

  if (nmpcSine_B.y_e) {
    c_size[0] = 0;
    c_size[1] = 0;
    Jc_size[0] = 0;
    Jc_size[1] = 0;
  } else {
    for (nmpcSine_B.i_o = 0; nmpcSine_B.i_o < 160; nmpcSine_B.i_o++) {
      nmpcSine_B.b_c[nmpcSine_B.i_o] = 0.0;
      nmpcSine_B.icf_i[nmpcSine_B.i_o] = true;
    }

    memset(&nmpcSine_B.Jx[0], 0, 19200U * sizeof(real_T));
    memset(&nmpcSine_B.Je[0], 0, 160U * sizeof(int8_T));
    nmpcSine_B.ic_idx_0 = 1.0;
    nmpcSine_B.ic_idx_1 = 2.0;
    nmpcSine_B.ic_idx_2 = 3.0;
    nmpcSine_B.ic_idx_3 = 4.0;
    for (nmpcSine_B.i_o = 0; nmpcSine_B.i_o < 20; nmpcSine_B.i_o++) {
      nmpcSine_B.runtimedata_OutputMin = runtimedata_OutputMin[nmpcSine_B.i_o];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_0) - 1] =
        ((!rtIsInf(nmpcSine_B.runtimedata_OutputMin)) && (!rtIsNaN
          (nmpcSine_B.runtimedata_OutputMin)));
      nmpcSine_B.runtimedata_OutputMin_l = runtimedata_OutputMin[nmpcSine_B.i_o
        + 20];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_1) - 1] =
        ((!rtIsInf(nmpcSine_B.runtimedata_OutputMin_l)) && (!rtIsNaN
          (nmpcSine_B.runtimedata_OutputMin_l)));
      nmpcSine_B.runtimedata_OutputMin_p = runtimedata_OutputMin[nmpcSine_B.i_o
        + 40];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_2) - 1] =
        ((!rtIsInf(nmpcSine_B.runtimedata_OutputMin_p)) && (!rtIsNaN
          (nmpcSine_B.runtimedata_OutputMin_p)));
      nmpcSine_B.runtimedata_OutputMin_pt = runtimedata_OutputMin[nmpcSine_B.i_o
        + 60];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_3) - 1] =
        ((!rtIsInf(nmpcSine_B.runtimedata_OutputMin_pt)) && (!rtIsNaN
          (nmpcSine_B.runtimedata_OutputMin_pt)));
      nmpcSine_B.runtimedata_OutputMax = runtimedata_OutputMax[nmpcSine_B.i_o];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_0 + 4.0) - 1] = ((
        !rtIsInf(nmpcSine_B.runtimedata_OutputMax)) && (!rtIsNaN
        (nmpcSine_B.runtimedata_OutputMax)));
      nmpcSine_B.icf_tmp_o[0] = static_cast<int32_T>(nmpcSine_B.ic_idx_0) - 1;
      nmpcSine_B.icf_tmp_o[4] = static_cast<int32_T>(nmpcSine_B.ic_idx_0 + 4.0)
        - 1;
      nmpcSine_B.runtimedata_OutputMax_f = runtimedata_OutputMax[nmpcSine_B.i_o
        + 20];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_1 + 4.0) - 1] = ((
        !rtIsInf(nmpcSine_B.runtimedata_OutputMax_f)) && (!rtIsNaN
        (nmpcSine_B.runtimedata_OutputMax_f)));
      nmpcSine_B.icf_tmp_o[1] = static_cast<int32_T>(nmpcSine_B.ic_idx_1) - 1;
      nmpcSine_B.icf_tmp_o[5] = static_cast<int32_T>(nmpcSine_B.ic_idx_1 + 4.0)
        - 1;
      nmpcSine_B.runtimedata_OutputMax_i = runtimedata_OutputMax[nmpcSine_B.i_o
        + 40];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_2 + 4.0) - 1] = ((
        !rtIsInf(nmpcSine_B.runtimedata_OutputMax_i)) && (!rtIsNaN
        (nmpcSine_B.runtimedata_OutputMax_i)));
      nmpcSine_B.icf_tmp_o[2] = static_cast<int32_T>(nmpcSine_B.ic_idx_2) - 1;
      nmpcSine_B.icf_tmp_o[6] = static_cast<int32_T>(nmpcSine_B.ic_idx_2 + 4.0)
        - 1;
      nmpcSine_B.runtimedata_OutputMax_o = runtimedata_OutputMax[nmpcSine_B.i_o
        + 60];
      nmpcSine_B.icf_i[static_cast<int32_T>(nmpcSine_B.ic_idx_3 + 4.0) - 1] = ((
        !rtIsInf(nmpcSine_B.runtimedata_OutputMax_o)) && (!rtIsNaN
        (nmpcSine_B.runtimedata_OutputMax_o)));
      nmpcSine_B.icf_tmp_o[3] = static_cast<int32_T>(nmpcSine_B.ic_idx_3) - 1;
      nmpcSine_B.icf_tmp_o[7] = static_cast<int32_T>(nmpcSine_B.ic_idx_3 + 4.0)
        - 1;
      for (nmpcSine_B.c_k = 0; nmpcSine_B.c_k < 8; nmpcSine_B.c_k++) {
        nmpcSine_B.icf_ip[nmpcSine_B.c_k] =
          nmpcSine_B.icf_i[nmpcSine_B.icf_tmp_o[nmpcSine_B.c_k]];
      }

      if (nmpcSine_any(nmpcSine_B.icf_ip)) {
        nmpcSine_B.yk_idx_0 = X[nmpcSine_B.i_o + 1];
        nmpcSine_B.yk_idx_1 = X[nmpcSine_B.i_o + 22];
        nmpcSine_B.yk_idx_2 = X[nmpcSine_B.i_o + 43];
        nmpcSine_B.yk_idx_3 = X[nmpcSine_B.i_o + 106];
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_0) - 1] =
          (nmpcSine_B.runtimedata_OutputMin - e) - nmpcSine_B.yk_idx_0;
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_1) - 1] =
          (nmpcSine_B.runtimedata_OutputMin_l - e) - nmpcSine_B.yk_idx_1;
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_2) - 1] =
          (nmpcSine_B.runtimedata_OutputMin_p - e) - nmpcSine_B.yk_idx_2;
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_3) - 1] =
          (nmpcSine_B.runtimedata_OutputMin_pt - e) - nmpcSine_B.yk_idx_3;
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_0 + 4.0) - 1] =
          (nmpcSine_B.yk_idx_0 - nmpcSine_B.runtimedata_OutputMax) - e;
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_1 + 4.0) - 1] =
          (nmpcSine_B.yk_idx_1 - nmpcSine_B.runtimedata_OutputMax_f) - e;
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_2 + 4.0) - 1] =
          (nmpcSine_B.yk_idx_2 - nmpcSine_B.runtimedata_OutputMax_i) - e;
        nmpcSine_B.b_c[static_cast<int32_T>(nmpcSine_B.ic_idx_3 + 4.0) - 1] =
          (nmpcSine_B.yk_idx_3 - nmpcSine_B.runtimedata_OutputMax_o) - e;
        for (nmpcSine_B.k_pi = 0; nmpcSine_B.k_pi < 6; nmpcSine_B.k_pi++) {
          nmpcSine_B.c_k = nmpcSine_B.k_pi << 2;
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_0) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            val[nmpcSine_B.c_k];
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_1) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            val[nmpcSine_B.c_k + 1];
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_2) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            val[nmpcSine_B.c_k + 2];
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_3) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            val[nmpcSine_B.c_k + 3];
        }

        for (nmpcSine_B.k_pi = 0; nmpcSine_B.k_pi < 6; nmpcSine_B.k_pi++) {
          nmpcSine_B.c_k = nmpcSine_B.k_pi << 2;
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_0 + 4.0) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            d[nmpcSine_B.c_k];
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_1 + 4.0) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            d[nmpcSine_B.c_k + 1];
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_2 + 4.0) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            d[nmpcSine_B.c_k + 2];
          nmpcSine_B.Jx[((static_cast<int32_T>(nmpcSine_B.ic_idx_3 + 4.0) + 160 *
                          nmpcSine_B.k_pi) + 960 * nmpcSine_B.i_o) - 1] =
            d[nmpcSine_B.c_k + 3];
        }

        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_0) - 1] = -1;
        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_1) - 1] = -1;
        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_2) - 1] = -1;
        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_3) - 1] = -1;
        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_0 + 4.0) - 1] = -1;
        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_1 + 4.0) - 1] = -1;
        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_2 + 4.0) - 1] = -1;
        nmpcSine_B.Je[static_cast<int32_T>(nmpcSine_B.ic_idx_3 + 4.0) - 1] = -1;
      }

      nmpcSine_B.ic_idx_0 += 8.0;
      nmpcSine_B.ic_idx_1 += 8.0;
      nmpcSine_B.ic_idx_2 += 8.0;
      nmpcSine_B.ic_idx_3 += 8.0;
    }

    nmpcSine_B.c_k = 0;
    for (nmpcSine_B.i_o = 0; nmpcSine_B.i_o < 160; nmpcSine_B.i_o++) {
      if (nmpcSine_B.icf_i[nmpcSine_B.i_o]) {
        nmpcSine_B.c_k++;
      }
    }

    nmpcSine_B.tmp_size_idx_0 = nmpcSine_B.c_k;
    nmpcSine_B.c_k = 0;
    for (nmpcSine_B.i_o = 0; nmpcSine_B.i_o < 160; nmpcSine_B.i_o++) {
      if (nmpcSine_B.icf_i[nmpcSine_B.i_o]) {
        nmpcSine_B.tmp_data_o[nmpcSine_B.c_k] = static_cast<uint8_T>
          (nmpcSine_B.i_o);
        nmpcSine_B.c_k++;
      }
    }

    c_size[0] = nmpcSine_B.tmp_size_idx_0;
    c_size[1] = 1;
    for (nmpcSine_B.c_k = 0; nmpcSine_B.c_k < nmpcSine_B.tmp_size_idx_0;
         nmpcSine_B.c_k++) {
      c_data[nmpcSine_B.c_k] =
        nmpcSine_B.b_c[nmpcSine_B.tmp_data_o[nmpcSine_B.c_k]];
    }

    nmpcSine_B.Jx_size[0] = nmpcSine_B.tmp_size_idx_0;
    nmpcSine_B.Jx_size[1] = 6;
    nmpcSine_B.Jx_size[2] = 20;
    for (nmpcSine_B.c_k = 0; nmpcSine_B.c_k < 20; nmpcSine_B.c_k++) {
      for (nmpcSine_B.k_pi = 0; nmpcSine_B.k_pi < 6; nmpcSine_B.k_pi++) {
        nmpcSine_B.i_o = (nmpcSine_B.tmp_size_idx_0 / 2) << 1;
        nmpcSine_B.vectorUB_l = nmpcSine_B.i_o - 2;
        for (nmpcSine_B.i4 = 0; nmpcSine_B.i4 <= nmpcSine_B.vectorUB_l;
             nmpcSine_B.i4 += 2) {
          nmpcSine_B.Jx_data_tmp = nmpcSine_B.tmp_size_idx_0 * nmpcSine_B.k_pi;
          nmpcSine_B.Jx_data_tmp_k = nmpcSine_B.tmp_size_idx_0 * 6 *
            nmpcSine_B.c_k;
          nmpcSine_B.Jx_data[(nmpcSine_B.i4 + nmpcSine_B.Jx_data_tmp) +
            nmpcSine_B.Jx_data_tmp_k] = nmpcSine_B.Jx[(160 * nmpcSine_B.k_pi +
            nmpcSine_B.tmp_data_o[nmpcSine_B.i4]) + 960 * nmpcSine_B.c_k];
          nmpcSine_B.Jx_data[((nmpcSine_B.i4 + nmpcSine_B.Jx_data_tmp) +
                              nmpcSine_B.Jx_data_tmp_k) + 1] = nmpcSine_B.Jx
            [(160 * nmpcSine_B.k_pi + nmpcSine_B.tmp_data_o[nmpcSine_B.i4 + 1])
            + 960 * nmpcSine_B.c_k];
        }

        for (nmpcSine_B.i4 = nmpcSine_B.i_o; nmpcSine_B.i4 <
             nmpcSine_B.tmp_size_idx_0; nmpcSine_B.i4++) {
          nmpcSine_B.Jx_data[(nmpcSine_B.i4 + nmpcSine_B.tmp_size_idx_0 *
                              nmpcSine_B.k_pi) + nmpcSine_B.tmp_size_idx_0 * 6 *
            nmpcSine_B.c_k] = nmpcSine_B.Jx[(160 * nmpcSine_B.k_pi +
            nmpcSine_B.tmp_data_o[nmpcSine_B.i4]) + 960 * nmpcSine_B.c_k];
        }
      }
    }

    nmpcSine_B.k_pi = (nmpcSine_B.tmp_size_idx_0 << 1) * 20;
    if (nmpcSine_B.k_pi - 1 >= 0) {
      memset(&nmpcSine_B.tmp_data[0], 0, static_cast<uint32_T>(nmpcSine_B.k_pi) *
             sizeof(real_T));
    }

    nmpcSine_B.Je_size[0] = nmpcSine_B.tmp_size_idx_0;
    for (nmpcSine_B.c_k = 0; nmpcSine_B.c_k < nmpcSine_B.tmp_size_idx_0;
         nmpcSine_B.c_k++) {
      nmpcSine_B.b_c[nmpcSine_B.c_k] =
        nmpcSine_B.Je[nmpcSine_B.tmp_data_o[nmpcSine_B.c_k]];
    }

    nmpcSine_reformJacobian(nmpcSine_B.Jx_data, nmpcSine_B.Jx_size,
      nmpcSine_B.tmp_data, nmpcSine_B.b_c, nmpcSine_B.Je_size, Jc_data, Jc_size);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_c4_mpclib_anonFcn2(const real_T runtimedata_x[6], const
  real_T runtimedata_OutputMin[80], const real_T runtimedata_OutputMax[80],
  const real_T z[125], real_T varargout_1_data[], int32_T varargout_1_size[2],
  real_T varargout_2[120], real_T varargout_3_data[], int32_T varargout_3_size[2],
  real_T varargout_4[15000])
{
  nmpcSine_getXUe(z, runtimedata_x, nmpcSine_B.X, nmpcSine_B.U, &nmpcSine_B.e_l);
  nmpcSine_stateEvolution(nmpcSine_B.X, nmpcSine_B.U, varargout_2, varargout_4);
  nmpcSine_outputBounds(runtimedata_OutputMin, runtimedata_OutputMax,
                        nmpcSine_B.X, nmpcSine_B.e_l,
                        nmpcSine_B.varargin_1_data_l, nmpcSine_B.varargin_1_size,
                        nmpcSine_B.b_varargin_1_data,
                        nmpcSine_B.b_varargin_1_size);
  nmpcSine_B.sizes_idx_1_tmp = ((nmpcSine_B.varargin_1_size[0] != 0) &&
    (nmpcSine_B.varargin_1_size[1] != 0));
  if (!nmpcSine_B.sizes_idx_1_tmp) {
    nmpcSine_B.sizes[0] = static_cast<uint8_T>(nmpcSine_B.varargin_1_size[0]);
  } else if (nmpcSine_B.sizes_idx_1_tmp) {
    nmpcSine_B.sizes[0] = static_cast<uint8_T>(nmpcSine_B.varargin_1_size[0]);
  } else {
    nmpcSine_B.sizes[0] = 0U;
  }

  varargout_1_size[0] = nmpcSine_B.sizes[0];
  varargout_1_size[1] = nmpcSine_B.sizes_idx_1_tmp;
  nmpcSine_B.loop_ub_o = nmpcSine_B.sizes_idx_1_tmp;
  for (nmpcSine_B.i_d = 0; nmpcSine_B.i_d < nmpcSine_B.loop_ub_o; nmpcSine_B.i_d
       ++) {
    nmpcSine_B.loop_ub_j = nmpcSine_B.sizes[0];
    if (nmpcSine_B.loop_ub_j - 1 >= 0) {
      memcpy(&varargout_1_data[0], &nmpcSine_B.varargin_1_data_l[0],
             static_cast<uint32_T>(nmpcSine_B.loop_ub_j) * sizeof(real_T));
    }
  }

  nmpcSine_B.sizes_idx_1_tmp = ((nmpcSine_B.b_varargin_1_size[0] != 0) &&
    (nmpcSine_B.b_varargin_1_size[1] != 0));
  if (nmpcSine_B.sizes_idx_1_tmp) {
    nmpcSine_B.sizes_idx_0 = static_cast<int8_T>(nmpcSine_B.b_varargin_1_size[0]);
  } else {
    nmpcSine_B.sizes_idx_0 = 0;
  }

  varargout_3_size[0] = nmpcSine_B.sizes_idx_0;
  if (nmpcSine_B.sizes_idx_0 == 0) {
    varargout_3_size[1] = nmpcSine_B.b_varargin_1_size[1];
    nmpcSine_B.u = static_cast<uint8_T>(nmpcSine_B.b_varargin_1_size[1]);
  } else if (nmpcSine_B.sizes_idx_1_tmp) {
    varargout_3_size[1] = nmpcSine_B.b_varargin_1_size[1];
    nmpcSine_B.u = static_cast<uint8_T>(nmpcSine_B.b_varargin_1_size[1]);
  } else {
    varargout_3_size[1] = 0;
    nmpcSine_B.u = 0U;
  }

  nmpcSine_B.loop_ub_o = nmpcSine_B.sizes_idx_0 * nmpcSine_B.u;
  if (nmpcSine_B.loop_ub_o - 1 >= 0) {
    memcpy(&varargout_3_data[0], &nmpcSine_B.b_varargin_1_data[0],
           static_cast<uint32_T>(nmpcSine_B.loop_ub_o) * sizeof(real_T));
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factoryConstruct(int32_T nVarMax, int32_T mConstrMax,
  int32_T mIneq, int32_T mNonlinIneq, s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *obj)
{
  obj->nVarMax = nVarMax;
  obj->mNonlinIneq = mNonlinIneq;
  obj->mNonlinEq = 120;
  obj->mIneq = mIneq;
  obj->mEq = 120;
  obj->iNonIneq0 = (mIneq - mNonlinIneq) + 1;
  obj->iNonEq0 = 1;
  obj->sqpFval = 0.0;
  obj->sqpFval_old = 0.0;
  obj->cIneq.size[0] = mIneq;
  obj->cIneq_old.size[0] = mIneq;
  obj->grad.size[0] = nVarMax;
  obj->grad_old.size[0] = nVarMax;
  obj->FunctionEvaluations = 0;
  obj->sqpIterations = 0;
  obj->sqpExitFlag = 0;
  obj->lambdasqp.size[0] = mConstrMax;
  if (mConstrMax - 1 >= 0) {
    memset(&obj->lambdasqp.data[0], 0, static_cast<uint32_T>(mConstrMax) *
           sizeof(real_T));
  }

  obj->lambdaStopTest.size[0] = mConstrMax;
  obj->lambdaStopTestPrev.size[0] = mConstrMax;
  obj->steplength = 1.0;
  obj->delta_x.size[0] = nVarMax;
  if (nVarMax - 1 >= 0) {
    memset(&obj->delta_x.data[0], 0, static_cast<uint32_T>(nVarMax) * sizeof
           (real_T));
  }

  obj->socDirection.size[0] = nVarMax;
  obj->workingset_old.size[0] = mConstrMax;
  if (mNonlinIneq > 0) {
    obj->JacCineqTrans_old.size[0] = nVarMax;
    obj->JacCineqTrans_old.size[1] = mNonlinIneq;
  } else {
    obj->JacCineqTrans_old.size[0] = 0;
    obj->JacCineqTrans_old.size[1] = 0;
  }

  obj->JacCeqTrans_old.size[0] = nVarMax;
  obj->JacCeqTrans_old.size[1] = 120;
  obj->gradLag.size[0] = nVarMax;
  obj->delta_gradLag.size[0] = nVarMax;
  obj->xstar.size[0] = nVarMax;
  obj->fstar = 0.0;
  obj->firstorderopt = 0.0;
  obj->lambda.size[0] = mConstrMax;
  if (mConstrMax - 1 >= 0) {
    memset(&obj->lambda.data[0], 0, static_cast<uint32_T>(mConstrMax) * sizeof
           (real_T));
  }

  obj->state = 0;
  obj->maxConstr = 0.0;
  obj->iterations = 0;
  obj->searchDir.size[0] = nVarMax;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factoryConstruct_h53b(int32_T MaxVars, int32_T
  obj_grad_size[1], int32_T obj_Hx_size[1], boolean_T *obj_hasLinear, int32_T
  *obj_nvar, int32_T *obj_maxVar, real_T *obj_beta, real_T *obj_rho, int32_T
  *obj_objtype, int32_T *obj_prev_objtype, int32_T *obj_prev_nvar, boolean_T
  *obj_prev_hasLinear, real_T *obj_gammaScalar)
{
  obj_grad_size[0] = MaxVars;
  obj_Hx_size[0] = MaxVars - 1;
  *obj_hasLinear = false;
  *obj_nvar = 0;
  *obj_maxVar = MaxVars;
  *obj_beta = 0.0;
  *obj_rho = 0.0;
  *obj_objtype = 3;
  *obj_prev_objtype = 3;
  *obj_prev_nvar = 0;
  *obj_prev_hasLinear = false;
  *obj_gammaScalar = 0.0;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factoryConstruct_h53bm(int32_T mIneqMax, int32_T nVarMax,
  int32_T mConstrMax, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj)
{
  obj->mConstr = 0;
  obj->mConstrOrig = 0;
  obj->mConstrMax = mConstrMax;
  obj->nVar = 125;
  obj->nVarOrig = 125;
  obj->nVarMax = nVarMax;
  obj->ldA = nVarMax;
  obj->Aineq.size[0] = mIneqMax * nVarMax;
  obj->bineq.size[0] = mIneqMax;
  obj->Aeq.size[0] = 120 * nVarMax;
  obj->lb.size[0] = nVarMax;
  obj->ub.size[0] = nVarMax;
  obj->indexLB.size[0] = nVarMax;
  obj->indexUB.size[0] = nVarMax;
  obj->indexFixed.size[0] = nVarMax;
  obj->mEqRemoved = 0;
  obj->ATwset.size[0] = nVarMax * mConstrMax;
  obj->bwset.size[0] = mConstrMax;
  obj->nActiveConstr = 0;
  obj->maxConstrWorkspace.size[0] = mConstrMax;
  for (int32_T i = 0; i < 5; i++) {
    obj->sizes[i] = 0;
    obj->sizesNormal[i] = 0;
    obj->sizesPhaseOne[i] = 0;
    obj->sizesRegularized[i] = 0;
    obj->sizesRegPhaseOne[i] = 0;
  }

  for (int32_T i = 0; i < 6; i++) {
    obj->isActiveIdx[i] = 0;
    obj->isActiveIdxNormal[i] = 0;
    obj->isActiveIdxPhaseOne[i] = 0;
    obj->isActiveIdxRegularized[i] = 0;
    obj->isActiveIdxRegPhaseOne[i] = 0;
  }

  obj->isActiveConstr.size[0] = mConstrMax;
  obj->Wid.size[0] = mConstrMax;
  obj->Wlocalidx.size[0] = mConstrMax;
  for (int32_T i = 0; i < 5; i++) {
    obj->nWConstr[i] = 0;
  }

  obj->probType = 3;
  obj->SLACK0 = 1.0E-5;
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_costFcn_d(const real_T X[126], const real_T U[84],
  const real_T data_References[80], real_T Q, real_T R, real_T Qt)
{
  __m128d tmp;
  __m128d tmp_0;
  nmpcSine_B.data_References = data_References[0];
  nmpcSine_B.data_References_n = data_References[20];
  for (nmpcSine_B.k_j = 0; nmpcSine_B.k_j <= 18; nmpcSine_B.k_j += 2) {
    _mm_storeu_pd(&nmpcSine_B.err_data[nmpcSine_B.k_j], _mm_sub_pd(_mm_set1_pd
      (nmpcSine_B.data_References), _mm_loadu_pd(&X[nmpcSine_B.k_j + 1])));
    _mm_storeu_pd(&nmpcSine_B.err_data[nmpcSine_B.k_j + 20], _mm_sub_pd
                  (_mm_set1_pd(nmpcSine_B.data_References_n), _mm_loadu_pd
                   (&X[nmpcSine_B.k_j + 22])));
  }

  for (nmpcSine_B.k_j = 0; nmpcSine_B.k_j <= 16; nmpcSine_B.k_j += 2) {
    tmp = _mm_loadu_pd(&nmpcSine_B.err_data[nmpcSine_B.k_j]);
    tmp_0 = _mm_loadu_pd(&nmpcSine_B.err_data[nmpcSine_B.k_j + 20]);
    _mm_storeu_pd(&nmpcSine_B.x_data[nmpcSine_B.k_j], _mm_add_pd(_mm_mul_pd(tmp,
      tmp), _mm_mul_pd(tmp_0, tmp_0)));
  }

  for (nmpcSine_B.k_j = 18; nmpcSine_B.k_j < 19; nmpcSine_B.k_j++) {
    nmpcSine_B.data_References_n = nmpcSine_B.err_data[nmpcSine_B.k_j];
    nmpcSine_B.varargin_1 = nmpcSine_B.err_data[nmpcSine_B.k_j + 20];
    nmpcSine_B.x_data[nmpcSine_B.k_j] = nmpcSine_B.data_References_n *
      nmpcSine_B.data_References_n + nmpcSine_B.varargin_1 *
      nmpcSine_B.varargin_1;
  }

  nmpcSine_B.data_References = nmpcSine_B.x_data[0];
  for (nmpcSine_B.k_j = 0; nmpcSine_B.k_j < 18; nmpcSine_B.k_j++) {
    nmpcSine_B.data_References += nmpcSine_B.x_data[nmpcSine_B.k_j + 1];
  }

  for (nmpcSine_B.k_j = 0; nmpcSine_B.k_j <= 16; nmpcSine_B.k_j += 2) {
    tmp = _mm_loadu_pd(&U[nmpcSine_B.k_j]);
    tmp_0 = _mm_loadu_pd(&U[nmpcSine_B.k_j + 21]);
    _mm_storeu_pd(&nmpcSine_B.x_data[nmpcSine_B.k_j], _mm_add_pd(_mm_mul_pd(tmp,
      tmp), _mm_mul_pd(tmp_0, tmp_0)));
  }

  for (nmpcSine_B.k_j = 18; nmpcSine_B.k_j < 19; nmpcSine_B.k_j++) {
    nmpcSine_B.data_References_n = U[nmpcSine_B.k_j];
    nmpcSine_B.varargin_1 = U[nmpcSine_B.k_j + 21];
    nmpcSine_B.x_data[nmpcSine_B.k_j] = nmpcSine_B.data_References_n *
      nmpcSine_B.data_References_n + nmpcSine_B.varargin_1 *
      nmpcSine_B.varargin_1;
  }

  nmpcSine_B.data_References_n = nmpcSine_B.x_data[0];
  for (nmpcSine_B.k_j = 0; nmpcSine_B.k_j < 18; nmpcSine_B.k_j++) {
    nmpcSine_B.data_References_n += nmpcSine_B.x_data[nmpcSine_B.k_j + 1];
  }

  nmpcSine_B.varargin_1 = nmpcSine_B.err_data[19];
  nmpcSine_B.J_tmp = nmpcSine_B.err_data[39];
  return (nmpcSine_B.varargin_1 * nmpcSine_B.varargin_1 + nmpcSine_B.J_tmp *
          nmpcSine_B.J_tmp) * Qt + (Q * nmpcSine_B.data_References + R *
    nmpcSine_B.data_References_n);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::computeObjectiveAndUserGradient(const
  s_I4XPpWw7d7shktLagLlNtD_nmpc_T *obj_next_next_next_next_next_ne, const real_T
  x[125], real_T grad_workspace_data[], real_T *fval, int32_T *status)
{
  __m128d tmp;
  __m128d tmp_0;
  __m128d tmp_1;
  static const int8_T b[160] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0,
    0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0,
    0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0,
    1, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 1 };

  nmpcSine_getXUe(x, obj_next_next_next_next_next_ne->runtimedata.x,
                  nmpcSine_B.X_d, nmpcSine_B.U_l, &nmpcSine_B.e_pb);
  *fval = nmpcSine_costFcn_d(nmpcSine_B.X_d, nmpcSine_B.U_l,
    obj_next_next_next_next_next_ne->userdata.References,
    obj_next_next_next_next_next_ne->runtimedata.Parameters[0],
    obj_next_next_next_next_next_ne->runtimedata.Parameters[1],
    obj_next_next_next_next_next_ne->runtimedata.Parameters[2]);
  nmpcSine_B.e_pb = obj_next_next_next_next_next_ne->runtimedata.Parameters[0];
  nmpcSine_B.R_m = obj_next_next_next_next_next_ne->runtimedata.Parameters[1];
  nmpcSine_B.obj_next_next_next_next_next_ne =
    obj_next_next_next_next_next_ne->userdata.References[0];
  nmpcSine_B.obj_next_next_next_next_next__o =
    obj_next_next_next_next_next_ne->userdata.References[20];
  for (nmpcSine_B.idx_current = 0; nmpcSine_B.idx_current <= 18;
       nmpcSine_B.idx_current += 2) {
    tmp_0 = _mm_loadu_pd(&nmpcSine_B.X_d[nmpcSine_B.idx_current + 1]);
    _mm_storeu_pd(&nmpcSine_B.err_x_data[nmpcSine_B.idx_current], _mm_sub_pd
                  (_mm_set1_pd(nmpcSine_B.obj_next_next_next_next_next_ne),
                   tmp_0));
    tmp_0 = _mm_loadu_pd(&nmpcSine_B.X_d[nmpcSine_B.idx_current + 22]);
    _mm_storeu_pd(&nmpcSine_B.err_y_data[nmpcSine_B.idx_current], _mm_sub_pd
                  (_mm_set1_pd(nmpcSine_B.obj_next_next_next_next_next__o),
                   tmp_0));
  }

  memset(&nmpcSine_B.Gfxu_data_l[0], 0, 120U * sizeof(real_T));
  memset(&nmpcSine_B.Gfuu_data[0], 0, 40U * sizeof(real_T));
  tmp_0 = _mm_set1_pd(-2.0);
  _mm_storeu_pd(&nmpcSine_B.dv7[0], _mm_mul_pd(_mm_mul_pd(tmp_0, _mm_set1_pd
    (obj_next_next_next_next_next_ne->runtimedata.Parameters[2])), _mm_set_pd
    (nmpcSine_B.err_y_data[19], nmpcSine_B.err_x_data[19])));
  nmpcSine_B.Gfxu_data_l[19] = nmpcSine_B.dv7[0];
  nmpcSine_B.Gfxu_data_l[39] = nmpcSine_B.dv7[1];
  for (nmpcSine_B.idx_current = 0; nmpcSine_B.idx_current <= 16;
       nmpcSine_B.idx_current += 2) {
    tmp = _mm_loadu_pd(&nmpcSine_B.err_x_data[nmpcSine_B.idx_current]);
    tmp_1 = _mm_set1_pd(-2.0 * nmpcSine_B.e_pb);
    _mm_storeu_pd(&nmpcSine_B.Gfxu_data_l[nmpcSine_B.idx_current], _mm_mul_pd
                  (tmp_1, tmp));
    tmp = _mm_loadu_pd(&nmpcSine_B.err_y_data[nmpcSine_B.idx_current]);
    nmpcSine_B.Gfxu_data_tmp = nmpcSine_B.idx_current + 20;
    _mm_storeu_pd(&nmpcSine_B.Gfxu_data_l[nmpcSine_B.Gfxu_data_tmp], _mm_mul_pd
                  (tmp_1, tmp));
    tmp = _mm_loadu_pd(&nmpcSine_B.U_l[nmpcSine_B.idx_current]);
    tmp_1 = _mm_set1_pd(2.0 * nmpcSine_B.R_m);
    _mm_storeu_pd(&nmpcSine_B.Gfuu_data[nmpcSine_B.idx_current], _mm_mul_pd
                  (tmp_1, tmp));
    tmp = _mm_loadu_pd(&nmpcSine_B.U_l[nmpcSine_B.idx_current + 21]);
    _mm_storeu_pd(&nmpcSine_B.Gfuu_data[nmpcSine_B.Gfxu_data_tmp], _mm_mul_pd
                  (tmp_1, tmp));
  }

  for (nmpcSine_B.idx_current = 18; nmpcSine_B.idx_current < 19;
       nmpcSine_B.idx_current++) {
    _mm_storeu_pd(&nmpcSine_B.dv7[0], _mm_mul_pd(_mm_mul_pd(tmp_0, _mm_set1_pd
      (nmpcSine_B.e_pb)), _mm_set_pd
      (nmpcSine_B.err_y_data[nmpcSine_B.idx_current],
       nmpcSine_B.err_x_data[nmpcSine_B.idx_current])));
    nmpcSine_B.Gfxu_data_l[nmpcSine_B.idx_current] = nmpcSine_B.dv7[0];
    nmpcSine_B.Gfxu_data_tmp = nmpcSine_B.idx_current + 20;
    nmpcSine_B.Gfxu_data_l[nmpcSine_B.Gfxu_data_tmp] = nmpcSine_B.dv7[1];
    _mm_storeu_pd(&nmpcSine_B.dv7[0], _mm_mul_pd(_mm_mul_pd(_mm_set1_pd(2.0),
      _mm_set1_pd(nmpcSine_B.R_m)), _mm_set_pd
      (nmpcSine_B.U_l[nmpcSine_B.idx_current + 21],
       nmpcSine_B.U_l[nmpcSine_B.idx_current])));
    nmpcSine_B.Gfuu_data[nmpcSine_B.idx_current] = nmpcSine_B.dv7[0];
    nmpcSine_B.Gfuu_data[nmpcSine_B.Gfxu_data_tmp] = nmpcSine_B.dv7[1];
  }

  for (nmpcSine_B.idx_current = 0; nmpcSine_B.idx_current < 20;
       nmpcSine_B.idx_current++) {
    for (nmpcSine_B.Gfxu_data_tmp = 0; nmpcSine_B.Gfxu_data_tmp < 6;
         nmpcSine_B.Gfxu_data_tmp++) {
      nmpcSine_B.Gfxu_data_o[nmpcSine_B.Gfxu_data_tmp + 6 *
        nmpcSine_B.idx_current] = nmpcSine_B.Gfxu_data_l[20 *
        nmpcSine_B.Gfxu_data_tmp + nmpcSine_B.idx_current];
    }

    nmpcSine_B.Gfuu_data_f[nmpcSine_B.idx_current << 1] =
      nmpcSine_B.Gfuu_data[nmpcSine_B.idx_current];
    nmpcSine_B.Gfuu_data_f[1 + (nmpcSine_B.idx_current << 1)] =
      nmpcSine_B.Gfuu_data[nmpcSine_B.idx_current + 20];
  }

  for (nmpcSine_B.idx_current = 0; nmpcSine_B.idx_current < 4;
       nmpcSine_B.idx_current++) {
    nmpcSine_B.e_pb = 0.0;
    for (nmpcSine_B.Gfxu_data_tmp = 0; nmpcSine_B.Gfxu_data_tmp < 40;
         nmpcSine_B.Gfxu_data_tmp++) {
      nmpcSine_B.e_pb += static_cast<real_T>(b[(nmpcSine_B.Gfxu_data_tmp << 2) +
        nmpcSine_B.idx_current]) *
        nmpcSine_B.Gfuu_data_f[nmpcSine_B.Gfxu_data_tmp];
    }

    nmpcSine_B.b[nmpcSine_B.idx_current] = nmpcSine_B.e_pb;
  }

  memcpy(&nmpcSine_B.Gfxu_data[0], &nmpcSine_B.Gfxu_data_o[0], 120U * sizeof
         (real_T));
  nmpcSine_B.Gfxu_data[120] = nmpcSine_B.b[0];
  nmpcSine_B.Gfxu_data[121] = nmpcSine_B.b[1];
  nmpcSine_B.Gfxu_data[122] = nmpcSine_B.b[2];
  nmpcSine_B.Gfxu_data[123] = nmpcSine_B.b[3];
  nmpcSine_B.Gfxu_data[124] = 0.0;
  memcpy(&grad_workspace_data[0], &nmpcSine_B.Gfxu_data[0], 125U * sizeof(real_T));
  *status = 1;
  nmpcSine_B.allFinite_n = rtIsNaN(*fval);
  if (rtIsInf(*fval) || nmpcSine_B.allFinite_n) {
    if (nmpcSine_B.allFinite_n) {
      *status = -3;
    } else if (*fval < 0.0) {
      *status = -1;
    } else {
      *status = -2;
    }
  } else {
    nmpcSine_B.allFinite_n = true;
    nmpcSine_B.idx_current = -1;
    while (nmpcSine_B.allFinite_n && (nmpcSine_B.idx_current + 2 <= 125)) {
      nmpcSine_B.e_pb = grad_workspace_data[nmpcSine_B.idx_current + 1];
      nmpcSine_B.allFinite_n = ((!rtIsInf(nmpcSine_B.e_pb)) && (!rtIsNaN
        (nmpcSine_B.e_pb)));
      nmpcSine_B.idx_current++;
    }

    if (!nmpcSine_B.allFinite_n) {
      if (rtIsNaN(grad_workspace_data[nmpcSine_B.idx_current])) {
        *status = -3;
      } else if (grad_workspace_data[nmpcSine_B.idx_current] < 0.0) {
        *status = -1;
      } else {
        *status = -2;
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
int32_T nmpcSine::nmpcSine_checkVectorNonFinite(int32_T N, const real_T
  vec_data[], int32_T iv0)
{
  int32_T idx_current;
  int32_T idx_end;
  int32_T status;
  boolean_T allFinite;
  status = 1;
  allFinite = true;
  idx_current = iv0 - 2;
  idx_end = (iv0 + N) - 1;
  while (allFinite && (idx_current + 2 <= idx_end)) {
    real_T allFinite_tmp;
    allFinite_tmp = vec_data[idx_current + 1];
    allFinite = ((!rtIsInf(allFinite_tmp)) && (!rtIsNaN(allFinite_tmp)));
    idx_current++;
  }

  if (!allFinite) {
    if (rtIsNaN(vec_data[idx_current])) {
      status = -3;
    } else if (vec_data[idx_current] < 0.0) {
      status = -1;
    } else {
      status = -2;
    }
  }

  return status;
}

// Function for MATLAB Function: '<S24>/NLMPC'
int32_T nmpcSine::nmpcSine_checkVectorNonFinite_n(const real_T vec[120])
{
  int32_T idx_current;
  int32_T status;
  boolean_T allFinite;
  status = 1;
  allFinite = true;
  idx_current = -1;
  while (allFinite && (idx_current + 2 <= 120)) {
    real_T allFinite_tmp;
    allFinite_tmp = vec[idx_current + 1];
    allFinite = ((!rtIsInf(allFinite_tmp)) && (!rtIsNaN(allFinite_tmp)));
    idx_current++;
  }

  if (!allFinite) {
    if (rtIsNaN(vec[idx_current])) {
      status = -3;
    } else if (vec[idx_current] < 0.0) {
      status = -1;
    } else {
      status = -2;
    }
  }

  return status;
}

// Function for MATLAB Function: '<S24>/NLMPC'
int32_T nmpcSine::computeConstraintsAndUserJacobi(int32_T
  obj_next_next_next_next_next_b_, const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
  *obj_next_next_next_next_next_ne, const real_T x[125], real_T
  Cineq_workspace_data[], int32_T ineq0, real_T Ceq_workspace[120], real_T
  JacIneqTrans_workspace_data[], int32_T iJI_col, int32_T ldJI, real_T
  JacEqTrans_workspace_data[], int32_T ldJE)
{
  int32_T status;
  if (obj_next_next_next_next_next_b_ > 0) {
    nmpcSine_c4_mpclib_anonFcn2(obj_next_next_next_next_next_ne->x,
      obj_next_next_next_next_next_ne->OutputMin,
      obj_next_next_next_next_next_ne->OutputMax, x, nmpcSine_B.a__3_data,
      nmpcSine_B.a__3_size, nmpcSine_B.b_x, nmpcSine_B.a__4_data,
      nmpcSine_B.a__4_size, nmpcSine_B.JacEqTrans_tmp);
    nmpcSine_B.col = static_cast<uint8_T>(obj_next_next_next_next_next_b_);
    for (nmpcSine_B.row = 0; nmpcSine_B.row < nmpcSine_B.col; nmpcSine_B.row++)
    {
      Cineq_workspace_data[(ineq0 + nmpcSine_B.row) - 1] =
        nmpcSine_B.a__3_data[nmpcSine_B.row];
    }

    memcpy(&Ceq_workspace[0], &nmpcSine_B.b_x[0], 120U * sizeof(real_T));
    nmpcSine_B.col_end = nmpcSine_B.a__4_size[0];
    for (nmpcSine_B.row = 0; nmpcSine_B.row < nmpcSine_B.col_end; nmpcSine_B.row
         ++) {
      nmpcSine_B.idx_mat = nmpcSine_B.a__4_size[1];
      for (nmpcSine_B.col = 0; nmpcSine_B.col < nmpcSine_B.idx_mat;
           nmpcSine_B.col++) {
        JacIneqTrans_workspace_data[nmpcSine_B.row + ldJI * ((iJI_col +
          nmpcSine_B.col) - 1)] = nmpcSine_B.a__4_data[nmpcSine_B.a__4_size[0] *
          nmpcSine_B.col + nmpcSine_B.row];
      }
    }

    for (nmpcSine_B.row = 0; nmpcSine_B.row < 125; nmpcSine_B.row++) {
      for (nmpcSine_B.col = 0; nmpcSine_B.col < 120; nmpcSine_B.col++) {
        JacEqTrans_workspace_data[nmpcSine_B.row + ldJE * nmpcSine_B.col] =
          nmpcSine_B.JacEqTrans_tmp[125 * nmpcSine_B.col + nmpcSine_B.row];
      }
    }
  } else {
    nmpcSine_c4_mpclib_anonFcn2(obj_next_next_next_next_next_ne->x,
      obj_next_next_next_next_next_ne->OutputMin,
      obj_next_next_next_next_next_ne->OutputMax, x, nmpcSine_B.a__3_data,
      nmpcSine_B.a__3_size, nmpcSine_B.b_x, nmpcSine_B.a__4_data,
      nmpcSine_B.a__4_size, nmpcSine_B.JacEqTrans_tmp);
    memcpy(&Ceq_workspace[0], &nmpcSine_B.b_x[0], 120U * sizeof(real_T));
    for (nmpcSine_B.row = 0; nmpcSine_B.row < 125; nmpcSine_B.row++) {
      for (nmpcSine_B.col = 0; nmpcSine_B.col < 120; nmpcSine_B.col++) {
        JacEqTrans_workspace_data[nmpcSine_B.row + ldJE * nmpcSine_B.col] =
          nmpcSine_B.JacEqTrans_tmp[125 * nmpcSine_B.col + nmpcSine_B.row];
      }
    }
  }

  status = nmpcSine_checkVectorNonFinite(obj_next_next_next_next_next_b_,
    Cineq_workspace_data, ineq0);
  if (status == 1) {
    status = nmpcSine_checkVectorNonFinite_n(Ceq_workspace);
    if (status == 1) {
      nmpcSine_B.allFinite = true;
      nmpcSine_B.row = -1;
      nmpcSine_B.col = iJI_col;
      nmpcSine_B.col_end = (iJI_col + obj_next_next_next_next_next_b_) - 1;
      while (nmpcSine_B.allFinite && (nmpcSine_B.col <= nmpcSine_B.col_end)) {
        nmpcSine_B.row = -1;
        while (nmpcSine_B.allFinite && (nmpcSine_B.row + 2 <= 125)) {
          nmpcSine_B.idx_mat = ((nmpcSine_B.col - 1) * ldJI + nmpcSine_B.row) +
            1;
          nmpcSine_B.allFinite = ((!rtIsInf
            (JacIneqTrans_workspace_data[nmpcSine_B.idx_mat])) && (!rtIsNaN
            (JacIneqTrans_workspace_data[nmpcSine_B.idx_mat])));
          nmpcSine_B.row++;
        }

        nmpcSine_B.col++;
      }

      if (!nmpcSine_B.allFinite) {
        nmpcSine_B.idx_mat = (nmpcSine_B.col - 2) * ldJI + nmpcSine_B.row;
        if (rtIsNaN(JacIneqTrans_workspace_data[nmpcSine_B.idx_mat])) {
          status = -3;
        } else if (JacIneqTrans_workspace_data[nmpcSine_B.idx_mat] < 0.0) {
          status = -1;
        } else {
          status = -2;
        }
      } else {
        nmpcSine_B.allFinite = true;
        nmpcSine_B.row = -1;
        nmpcSine_B.col = -1;
        while (nmpcSine_B.allFinite && (nmpcSine_B.col + 2 <= 120)) {
          nmpcSine_B.row = -1;
          while (nmpcSine_B.allFinite && (nmpcSine_B.row + 2 <= 125)) {
            nmpcSine_B.col_end = ((nmpcSine_B.col + 1) * ldJE + nmpcSine_B.row)
              + 1;
            nmpcSine_B.allFinite = ((!rtIsInf
              (JacEqTrans_workspace_data[nmpcSine_B.col_end])) && (!rtIsNaN
              (JacEqTrans_workspace_data[nmpcSine_B.col_end])));
            nmpcSine_B.row++;
          }

          nmpcSine_B.col++;
        }

        if (!nmpcSine_B.allFinite) {
          nmpcSine_B.col_end = ldJE * nmpcSine_B.col + nmpcSine_B.row;
          if (rtIsNaN(JacEqTrans_workspace_data[nmpcSine_B.col_end])) {
            status = -3;
          } else if (JacEqTrans_workspace_data[nmpcSine_B.col_end] < 0.0) {
            status = -1;
          } else {
            status = -2;
          }
        }
      }
    }
  }

  return status;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::evalObjAndConstrAndDerivatives(int32_T
  obj_next_next_next_next_next_b_, const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
  *obj_next_next_next_next_next_ne, const s_I4XPpWw7d7shktLagLlNtD_nmpc_T
  *obj_next_next_next_next_next__0, const real_T x[125], real_T
  grad_workspace_data[], real_T Cineq_workspace_data[], int32_T ineq0, real_T
  Ceq_workspace[120], real_T JacIneqTrans_workspace_data[], int32_T iJI_col,
  int32_T ldJI, real_T JacEqTrans_workspace_data[], int32_T ldJE, real_T *fval,
  int32_T *status)
{
  computeObjectiveAndUserGradient(obj_next_next_next_next_next__0, x,
    grad_workspace_data, fval, status);
  if (*status == 1) {
    *status = computeConstraintsAndUserJacobi(obj_next_next_next_next_next_b_,
      obj_next_next_next_next_next_ne, x, Cineq_workspace_data, ineq0,
      Ceq_workspace, JacIneqTrans_workspace_data, iJI_col, ldJI,
      JacEqTrans_workspace_data, ldJE);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSin_modifyOverheadPhaseOne_(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *obj)
{
  int32_T d;
  int32_T idxEq;
  idxEq = static_cast<uint16_T>(obj->sizes[0]);
  for (int32_T idx = 0; idx < idxEq; idx++) {
    obj->ATwset.data[(obj->nVar + obj->ldA * idx) - 1] = 0.0;
  }

  for (int32_T idx = 0; idx < 120; idx++) {
    idxEq = (obj->ldA * idx + obj->nVar) - 1;
    obj->Aeq.data[idxEq] = 0.0;
    obj->ATwset.data[idxEq + obj->ldA * (obj->isActiveIdx[1] - 1)] = 0.0;
  }

  idxEq = static_cast<uint16_T>(obj->sizes[2]);
  for (int32_T idx = 0; idx < idxEq; idx++) {
    obj->Aineq.data[(obj->nVar + obj->ldA * idx) - 1] = -1.0;
  }

  obj->indexLB.data[obj->sizes[3] - 1] = obj->nVar;
  obj->lb.data[obj->nVar - 1] = 1.0E-5;
  idxEq = obj->isActiveIdx[2];
  d = obj->nActiveConstr;
  for (int32_T idx = idxEq; idx <= d; idx++) {
    obj->ATwset.data[(obj->nVar + obj->ldA * (idx - 1)) - 1] = -1.0;
  }

  idxEq = obj->isActiveIdx[4] - 1;
  if (obj->nWConstr[4] > 0) {
    d = obj->sizesNormal[4];
    for (int32_T idx = d; idx >= 1; idx--) {
      int32_T tmp;
      tmp = idxEq + idx;
      obj->isActiveConstr.data[tmp] = obj->isActiveConstr.data[tmp - 1];
    }
  } else {
    obj->isActiveConstr.data[(obj->isActiveIdx[4] + obj->sizesNormal[4]) - 1] =
      false;
  }

  obj->isActiveConstr.data[obj->isActiveIdx[4] - 1] = false;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_setProblemType(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
  int32_T PROBLEM_TYPE)
{
  int32_T c;
  int32_T colOffsetATw;
  int32_T colOffsetAineq;
  int32_T d;
  int32_T g;
  int32_T idxUpperExisting;
  int32_T offsetEq1;
  int32_T offsetEq2;
  switch (PROBLEM_TYPE) {
   case 3:
    obj->nVar = 125;
    obj->mConstr = obj->mConstrOrig;
    if (obj->nWConstr[4] > 0) {
      idxUpperExisting = obj->isActiveIdx[4] - 1;
      offsetEq1 = static_cast<uint16_T>(obj->sizesNormal[4]);
      for (colOffsetATw = 0; colOffsetATw < offsetEq1; colOffsetATw++) {
        c = idxUpperExisting + colOffsetATw;
        obj->isActiveConstr.data[(obj->isActiveIdxNormal[4] + colOffsetATw) - 1]
          = obj->isActiveConstr.data[c];
        obj->isActiveConstr.data[c] = false;
      }
    }

    for (c = 0; c < 5; c++) {
      obj->sizes[c] = obj->sizesNormal[c];
    }

    for (c = 0; c < 6; c++) {
      obj->isActiveIdx[c] = obj->isActiveIdxNormal[c];
    }
    break;

   case 1:
    obj->nVar = 126;
    obj->mConstr = obj->mConstrOrig + 1;
    for (c = 0; c < 5; c++) {
      obj->sizes[c] = obj->sizesPhaseOne[c];
    }

    nmpcSin_modifyOverheadPhaseOne_(obj);
    for (c = 0; c < 6; c++) {
      obj->isActiveIdx[c] = obj->isActiveIdxPhaseOne[c];
    }
    break;

   case 2:
    obj->nVar = obj->nVarMax - 1;
    obj->mConstr = obj->mConstrMax - 1;
    for (c = 0; c < 5; c++) {
      obj->sizes[c] = obj->sizesRegularized[c];
    }

    if (obj->probType != 4) {
      offsetEq2 = obj->sizes[2] + 245;
      offsetEq1 = obj->sizes[2] + 125;
      c = static_cast<uint16_T>(obj->sizes[0]);
      for (idxUpperExisting = 0; idxUpperExisting < c; idxUpperExisting++) {
        memset(&obj->ATwset.data[obj->ldA * idxUpperExisting + 125], 0,
               static_cast<uint32_T>(obj->nVar - 125) * sizeof(real_T));
      }

      idxUpperExisting = static_cast<uint16_T>(obj->sizes[2]);
      for (colOffsetATw = 0; colOffsetATw < idxUpperExisting; colOffsetATw++) {
        colOffsetAineq = obj->ldA * colOffsetATw - 1;
        for (c = 126; c <= colOffsetATw + 125; c++) {
          obj->Aineq.data[c + colOffsetAineq] = 0.0;
        }

        obj->Aineq.data[(colOffsetATw + colOffsetAineq) + 126] = -1.0;
        d = obj->nVar;
        for (c = colOffsetATw + 127; c <= d; c++) {
          obj->Aineq.data[c + colOffsetAineq] = 0.0;
        }
      }

      for (idxUpperExisting = 0; idxUpperExisting < 120; idxUpperExisting++) {
        colOffsetAineq = obj->ldA * idxUpperExisting - 1;
        colOffsetATw = (obj->isActiveIdx[1] - 1) * obj->ldA + colOffsetAineq;
        if (offsetEq1 >= 126) {
          memset(&obj->Aeq.data[colOffsetAineq + 126], 0, static_cast<uint32_T>
                 (offsetEq1 - 125) * sizeof(real_T));
        }

        if (offsetEq1 >= 126) {
          memset(&obj->ATwset.data[colOffsetATw + 126], 0, static_cast<uint32_T>
                 (offsetEq1 - 125) * sizeof(real_T));
        }

        d = offsetEq2 + idxUpperExisting;
        if (offsetEq2 - 119 <= d - 120) {
          memset(&obj->Aeq.data[(offsetEq2 + colOffsetAineq) + -119], 0,
                 static_cast<uint32_T>(((((d - 120) + colOffsetAineq) -
                    offsetEq2) - colOffsetAineq) + 120) * sizeof(real_T));
        }

        if (offsetEq2 - 119 <= d - 120) {
          memset(&obj->ATwset.data[(offsetEq2 + colOffsetATw) + -119], 0,
                 static_cast<uint32_T>(((((d - 120) + colOffsetATw) - offsetEq2)
                   - colOffsetATw) + 120) * sizeof(real_T));
        }

        c = d + colOffsetAineq;
        obj->Aeq.data[c - 119] = -1.0;
        g = d + colOffsetATw;
        obj->ATwset.data[g - 119] = -1.0;
        if (d - 118 <= offsetEq2) {
          memset(&obj->Aeq.data[(d - 118) + colOffsetAineq], 0,
                 static_cast<uint32_T>((((offsetEq2 + colOffsetAineq) - (d - 118))
                   - colOffsetAineq) + 1) * sizeof(real_T));
        }

        if (d - 118 <= offsetEq2) {
          memset(&obj->ATwset.data[(d - 118) + colOffsetATw], 0,
                 static_cast<uint32_T>((((offsetEq2 + colOffsetATw) - (d - 118))
                   - colOffsetATw) + 1) * sizeof(real_T));
        }

        if (offsetEq2 + 1 <= d) {
          memset(&obj->Aeq.data[(offsetEq2 + colOffsetAineq) + 1], 0,
                 static_cast<uint32_T>(((d + colOffsetAineq) - offsetEq2) -
                  colOffsetAineq) * sizeof(real_T));
        }

        if (offsetEq2 + 1 <= d) {
          memset(&obj->ATwset.data[(offsetEq2 + colOffsetATw) + 1], 0,
                 static_cast<uint32_T>(((d + colOffsetATw) - offsetEq2) -
                  colOffsetATw) * sizeof(real_T));
        }

        obj->Aeq.data[c + 1] = 1.0;
        obj->ATwset.data[g + 1] = 1.0;
        d += 2;
        if (d <= obj->nVar) {
          memset(&obj->Aeq.data[d + colOffsetAineq], 0, static_cast<uint32_T>
                 ((((obj->nVar + colOffsetAineq) - d) - colOffsetAineq) + 1) *
                 sizeof(real_T));
        }

        if (d <= obj->nVar) {
          memset(&obj->ATwset.data[d + colOffsetATw], 0, static_cast<uint32_T>
                 ((((obj->nVar + colOffsetATw) - d) - colOffsetATw) + 1) *
                 sizeof(real_T));
        }
      }

      idxUpperExisting = 125;
      offsetEq1 = obj->sizesNormal[3] + 1;
      offsetEq2 = obj->sizesRegularized[3];
      for (colOffsetATw = offsetEq1; colOffsetATw <= offsetEq2; colOffsetATw++)
      {
        idxUpperExisting++;
        obj->indexLB.data[colOffsetATw - 1] = idxUpperExisting;
      }

      if (obj->nWConstr[4] > 0) {
        idxUpperExisting = static_cast<uint16_T>(obj->sizesRegularized[4]);
        for (colOffsetATw = 0; colOffsetATw < idxUpperExisting; colOffsetATw++)
        {
          obj->isActiveConstr.data[obj->isActiveIdxRegularized[4] + colOffsetATw]
            = obj->isActiveConstr.data[(obj->isActiveIdx[4] + colOffsetATw) - 1];
        }
      }

      idxUpperExisting = obj->isActiveIdx[4];
      offsetEq1 = obj->isActiveIdxRegularized[4];
      if (idxUpperExisting <= offsetEq1 - 1) {
        memset(&obj->isActiveConstr.data[idxUpperExisting + -1], 0, static_cast<
               uint32_T>(offsetEq1 - idxUpperExisting) * sizeof(boolean_T));
      }

      idxUpperExisting = obj->sizes[2] + 365;
      memset(&obj->lb.data[125], 0, static_cast<uint32_T>(idxUpperExisting - 125)
             * sizeof(real_T));
      offsetEq1 = obj->isActiveIdx[2];
      offsetEq2 = obj->nActiveConstr;
      for (idxUpperExisting = offsetEq1; idxUpperExisting <= offsetEq2;
           idxUpperExisting++) {
        colOffsetATw = (idxUpperExisting - 1) * obj->ldA - 1;
        if (obj->Wid.data[idxUpperExisting - 1] == 3) {
          c = obj->Wlocalidx.data[idxUpperExisting - 1];
          if (c + 124 >= 126) {
            memset(&obj->ATwset.data[colOffsetATw + 126], 0,
                   static_cast<uint32_T>((c + 124) - 125) * sizeof(real_T));
          }

          obj->ATwset.data[(c + colOffsetATw) + 125] = -1.0;
          if (c + 126 <= obj->nVar) {
            memset(&obj->ATwset.data[(c + 126) + colOffsetATw], 0,
                   static_cast<uint32_T>((((obj->nVar + colOffsetATw) - (c + 126))
                     - colOffsetATw) + 1) * sizeof(real_T));
          }
        } else {
          memset(&obj->ATwset.data[colOffsetATw + 126], 0, static_cast<uint32_T>
                 (obj->nVar - 125) * sizeof(real_T));
        }
      }
    }

    for (c = 0; c < 6; c++) {
      obj->isActiveIdx[c] = obj->isActiveIdxRegularized[c];
    }
    break;

   default:
    obj->nVar = obj->nVarMax;
    obj->mConstr = obj->mConstrMax;
    for (c = 0; c < 5; c++) {
      obj->sizes[c] = obj->sizesRegPhaseOne[c];
    }

    nmpcSin_modifyOverheadPhaseOne_(obj);
    for (c = 0; c < 6; c++) {
      obj->isActiveIdx[c] = obj->isActiveIdxRegPhaseOne[c];
    }
    break;
  }

  obj->probType = PROBLEM_TYPE;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_initActiveSet(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj)
{
  int32_T colOffsetATw;
  int32_T f;
  int32_T iATw0;
  int32_T iAeq0;
  int32_T idx;
  int32_T idxFillStart;
  nmpcSine_setProblemType(obj, 3);
  idxFillStart = obj->isActiveIdx[2];
  if (idxFillStart <= obj->mConstrMax) {
    memset(&obj->isActiveConstr.data[idxFillStart + -1], 0, static_cast<uint32_T>
           ((obj->mConstrMax - idxFillStart) + 1) * sizeof(boolean_T));
  }

  obj->nWConstr[0] = obj->sizes[0];
  obj->nWConstr[1] = 120;
  obj->nWConstr[2] = 0;
  obj->nWConstr[3] = 0;
  obj->nWConstr[4] = 0;
  obj->nActiveConstr = obj->nWConstr[0] + 120;
  idxFillStart = static_cast<uint16_T>(obj->sizes[0]);
  for (idx = 0; idx < idxFillStart; idx++) {
    obj->Wid.data[idx] = 1;
    obj->Wlocalidx.data[idx] = idx + 1;
    obj->isActiveConstr.data[idx] = true;
    colOffsetATw = obj->ldA * idx;
    iATw0 = static_cast<uint8_T>(obj->indexFixed.data[idx] - 1);
    if (iATw0 - 1 >= 0) {
      memset(&obj->ATwset.data[colOffsetATw], 0, static_cast<uint32_T>(iATw0) *
             sizeof(real_T));
    }

    obj->ATwset.data[(obj->indexFixed.data[idx] + colOffsetATw) - 1] = 1.0;
    iATw0 = obj->indexFixed.data[idx] + 1;
    if (iATw0 <= obj->nVar) {
      memset(&obj->ATwset.data[(iATw0 + colOffsetATw) + -1], 0, static_cast<
             uint32_T>((((obj->nVar + colOffsetATw) - iATw0) - colOffsetATw) + 1)
             * sizeof(real_T));
    }

    obj->bwset.data[idx] = obj->ub.data[obj->indexFixed.data[idx] - 1];
  }

  for (idx = 0; idx < 120; idx++) {
    colOffsetATw = obj->sizes[0] + idx;
    obj->Wid.data[colOffsetATw] = 2;
    obj->Wlocalidx.data[colOffsetATw] = idx + 1;
    obj->isActiveConstr.data[colOffsetATw] = true;
    iAeq0 = obj->ldA * idx;
    iATw0 = obj->ldA * colOffsetATw;
    f = obj->nVar;
    for (idxFillStart = 0; idxFillStart < f; idxFillStart++) {
      obj->ATwset.data[iATw0 + idxFillStart] = obj->Aeq.data[iAeq0 +
        idxFillStart];
    }

    obj->bwset.data[colOffsetATw] = obj->beq[idx];
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factoryConstruct_h5(int32_T maxRows, int32_T maxCols,
  int32_T *obj_ldq, int32_T obj_QR_size[2], real_T obj_Q_data[], int32_T
  obj_Q_size[2], int32_T obj_jpvt_data[], int32_T obj_jpvt_size[1], int32_T
  *obj_mrows, int32_T *obj_ncols, int32_T obj_tau_size[1], int32_T
  *obj_minRowCol, boolean_T *obj_usedPivoting)
{
  int32_T loop_ub_tmp;
  *obj_ldq = maxRows;
  obj_QR_size[0] = maxRows;
  obj_QR_size[1] = maxCols;
  obj_Q_size[0] = maxRows;
  obj_Q_size[1] = maxRows;
  loop_ub_tmp = maxRows * maxRows;
  if (loop_ub_tmp - 1 >= 0) {
    memset(&obj_Q_data[0], 0, static_cast<uint32_T>(loop_ub_tmp) * sizeof(real_T));
  }

  obj_jpvt_size[0] = maxCols;
  if (maxCols - 1 >= 0) {
    memset(&obj_jpvt_data[0], 0, static_cast<uint32_T>(maxCols) * sizeof(int32_T));
  }

  *obj_mrows = 0;
  *obj_ncols = 0;
  if (maxRows <= maxCols) {
    obj_tau_size[0] = maxRows;
  } else {
    obj_tau_size[0] = maxCols;
  }

  *obj_minRowCol = 0;
  *obj_usedPivoting = false;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factoryConstruct_h53(int32_T MaxDims, int32_T
  obj_FMat_size[2], int32_T *obj_ldm, int32_T *obj_ndims, int32_T *obj_info,
  real_T *obj_scaleFactor, boolean_T *obj_ConvexCheck, real_T *obj_regTol_,
  real_T *obj_workspace_, real_T *obj_workspace2_)
{
  obj_FMat_size[0] = MaxDims;
  obj_FMat_size[1] = MaxDims;
  *obj_ldm = MaxDims;
  *obj_ndims = 0;
  *obj_info = 0;
  *obj_scaleFactor = 0.0;
  *obj_ConvexCheck = true;
  *obj_regTol_ = (rtInf);
  *obj_workspace_ = (rtInf);
  *obj_workspace2_ = (rtInf);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_computeGradLag(real_T workspace_data[], int32_T ldA,
  int32_T nVar, const real_T grad_data[], int32_T mIneq, const real_T
  AineqTrans_data[], const real_T AeqTrans_data[], const int32_T
  finiteFixed_data[], int32_T mFixed, const int32_T finiteLB_data[], int32_T mLB,
  const int32_T finiteUB_data[], int32_T mUB, const real_T lambda_data[])
{
  int32_T b;
  int32_T f;
  int32_T finiteFixed;
  int32_T g;
  int32_T iL0;
  int32_T ix;
  memcpy(&workspace_data[0], &grad_data[0], static_cast<uint16_T>(nVar) * sizeof
         (real_T));
  b = static_cast<uint16_T>(mFixed);
  for (iL0 = 0; iL0 < b; iL0++) {
    finiteFixed = finiteFixed_data[iL0];
    workspace_data[finiteFixed - 1] += lambda_data[iL0];
  }

  ix = mFixed;
  f = ldA * 119 + 1;
  for (b = 1; ldA < 0 ? b >= f : b <= f; b += ldA) {
    g = (b + nVar) - 1;
    for (finiteFixed = b; finiteFixed <= g; finiteFixed++) {
      iL0 = finiteFixed - b;
      workspace_data[iL0] += AeqTrans_data[finiteFixed - 1] * lambda_data[ix];
    }

    ix++;
  }

  if (mIneq != 0) {
    ix = mFixed + 120;
    f = (mIneq - 1) * ldA + 1;
    for (b = 1; ldA < 0 ? b >= f : b <= f; b += ldA) {
      g = (b + nVar) - 1;
      for (finiteFixed = b; finiteFixed <= g; finiteFixed++) {
        iL0 = finiteFixed - b;
        workspace_data[iL0] += AineqTrans_data[finiteFixed - 1] * lambda_data[ix];
      }

      ix++;
    }
  }

  iL0 = (mFixed + mIneq) + 120;
  finiteFixed = static_cast<uint16_T>(mLB) - 1;
  for (b = 0; b <= finiteFixed; b++) {
    ix = finiteLB_data[b];
    workspace_data[ix - 1] -= lambda_data[iL0 + b];
  }

  if (static_cast<uint16_T>(mLB) - 1 >= 0) {
    iL0 += static_cast<uint16_T>(mLB);
  }

  finiteFixed = static_cast<uint16_T>(mUB) - 1;
  for (b = 0; b <= finiteFixed; b++) {
    ix = finiteUB_data[b];
    workspace_data[ix - 1] += lambda_data[iL0 + b];
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_computePrimalFeasError(const real_T x[125], int32_T
  mLinIneq, int32_T mNonlinIneq, const real_T cIneq_data[], const real_T cEq[120],
  const int32_T finiteLB_data[], int32_T mLB, const real_T lb[125], const
  int32_T finiteUB_data[], int32_T mUB, const real_T ub[125])
{
  real_T feasError;
  real_T u1;
  int32_T finiteLB;
  int32_T mIneq;
  feasError = 0.0;
  mIneq = mNonlinIneq + mLinIneq;
  for (int32_T idx = 0; idx < 120; idx++) {
    u1 = fabs(cEq[idx]);
    if ((!(feasError >= u1)) && (!rtIsNaN(u1))) {
      feasError = u1;
    }
  }

  for (int32_T idx = 0; idx < mIneq; idx++) {
    u1 = cIneq_data[idx];
    if ((!(feasError >= u1)) && (!rtIsNaN(u1))) {
      feasError = u1;
    }
  }

  mIneq = static_cast<uint16_T>(mLB);
  for (int32_T idx = 0; idx < mIneq; idx++) {
    finiteLB = finiteLB_data[idx];
    u1 = lb[finiteLB - 1] - x[finiteLB - 1];
    if ((!(feasError >= u1)) && (!rtIsNaN(u1))) {
      feasError = u1;
    }
  }

  mIneq = static_cast<uint16_T>(mUB);
  for (int32_T idx = 0; idx < mIneq; idx++) {
    finiteLB = finiteUB_data[idx];
    u1 = x[finiteLB - 1] - ub[finiteLB - 1];
    if ((!(feasError >= u1)) && (!rtIsNaN(u1))) {
      feasError = u1;
    }
  }

  return feasError;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_computeDualFeasError(int32_T nVar, const real_T
  gradLag_data[], boolean_T *gradOK, real_T *val)
{
  int32_T idx;
  boolean_T exitg1;
  *gradOK = true;
  *val = 0.0;
  idx = 0;
  exitg1 = false;
  while ((!exitg1) && (idx <= static_cast<uint16_T>(nVar) - 1)) {
    *gradOK = ((!rtIsInf(gradLag_data[idx])) && (!rtIsNaN(gradLag_data[idx])));
    if (!*gradOK) {
      exitg1 = true;
    } else {
      real_T u1;
      u1 = fabs(gradLag_data[idx]);
      if ((!(*val >= u1)) && (!rtIsNaN(u1))) {
        *val = u1;
      }

      idx++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_test_exit(sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction,
  const s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState, const real_T lb[125], const
  real_T ub[125], boolean_T *Flags_gradOK, boolean_T *Flags_fevalOK, boolean_T
  *Flags_done, boolean_T *Flags_stepAccepted, boolean_T *Flags_failedLineSearch,
  int32_T *Flags_stepType)
{
  real_T s;
  real_T smax;
  int32_T idx_max;
  int32_T k;
  int32_T mLambda;
  int32_T nVar;
  boolean_T isFeasible;
  *Flags_fevalOK = true;
  *Flags_stepAccepted = false;
  *Flags_failedLineSearch = false;
  *Flags_stepType = 1;
  nVar = WorkingSet->nVar;
  mLambda = (((WorkingSet->sizes[0] + WorkingSet->sizes[2]) + WorkingSet->sizes
              [3]) + WorkingSet->sizes[4]) + 119;
  if (mLambda >= 0) {
    memcpy(&TrialState->lambdaStopTest.data[0], &TrialState->lambdasqp.data[0],
           static_cast<uint32_T>(mLambda + 1) * sizeof(real_T));
  }

  nmpcSine_computeGradLag(TrialState->gradLag.data, WorkingSet->ldA,
    WorkingSet->nVar, TrialState->grad.data, WorkingSet->sizes[2],
    WorkingSet->Aineq.data, WorkingSet->Aeq.data, WorkingSet->indexFixed.data,
    WorkingSet->sizes[0], WorkingSet->indexLB.data, WorkingSet->sizes[3],
    WorkingSet->indexUB.data, WorkingSet->sizes[4],
    TrialState->lambdaStopTest.data);
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = fabs(TrialState->grad.data[0]);
      for (k = 2; k <= nVar; k++) {
        s = fabs(TrialState->grad.data[k - 1]);
        if (s > smax) {
          idx_max = k;
          smax = s;
        }
      }
    }
  }

  smax = fabs(TrialState->grad.data[idx_max - 1]);
  if ((smax <= 1.0) || rtIsNaN(smax)) {
    smax = 1.0;
  }

  if (rtIsInf(smax)) {
    smax = 1.0;
  }

  MeritFunction->nlpPrimalFeasError = nmpcSine_computePrimalFeasError
    (TrialState->xstarsqp, WorkingSet->sizes[2] - TrialState->mNonlinIneq,
     TrialState->mNonlinIneq, TrialState->cIneq.data, TrialState->cEq,
     WorkingSet->indexLB.data, WorkingSet->sizes[3], lb,
     WorkingSet->indexUB.data, WorkingSet->sizes[4], ub);
  if ((MeritFunction->nlpPrimalFeasError <= 1.0) || rtIsNaN
      (MeritFunction->nlpPrimalFeasError)) {
    MeritFunction->feasRelativeFactor = 1.0;
  } else {
    MeritFunction->feasRelativeFactor = MeritFunction->nlpPrimalFeasError;
  }

  isFeasible = (MeritFunction->nlpPrimalFeasError <= 0.001 *
                MeritFunction->feasRelativeFactor);
  nmpcSine_computeDualFeasError(WorkingSet->nVar, TrialState->gradLag.data,
    Flags_gradOK, &MeritFunction->nlpDualFeasError);
  if (!*Flags_gradOK) {
    *Flags_done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    MeritFunction->nlpComplError = 0.0;
    if (MeritFunction->nlpDualFeasError >= 0.0) {
      MeritFunction->firstOrderOpt = MeritFunction->nlpDualFeasError;
    } else {
      MeritFunction->firstOrderOpt = 0.0;
    }

    if (mLambda >= 0) {
      memcpy(&TrialState->lambdaStopTestPrev.data[0],
             &TrialState->lambdaStopTest.data[0], static_cast<uint32_T>(mLambda
              + 1) * sizeof(real_T));
    }

    if (isFeasible && (MeritFunction->nlpDualFeasError <= 0.001 * smax)) {
      *Flags_done = true;
      TrialState->sqpExitFlag = 1;
    } else {
      *Flags_done = false;
      if (isFeasible && (TrialState->sqpFval < -1.0E+20)) {
        *Flags_done = true;
        TrialState->sqpExitFlag = -3;
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_saveJacobian(s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *obj,
  int32_T nVar, int32_T mIneq, const real_T JacCineqTrans_data[], int32_T
  ineqCol0, const real_T JacCeqTrans_data[], int32_T ldJ)
{
  int32_T b;
  int32_T iCol;
  int32_T iCol_old;
  int32_T loop_ub_tmp;
  iCol = (ineqCol0 - 1) * ldJ;
  iCol_old = 0;
  b = mIneq - ineqCol0;
  for (int32_T idx_col = 0; idx_col <= b; idx_col++) {
    int32_T c;
    loop_ub_tmp = obj->JacCineqTrans_old.size[0] * obj->JacCineqTrans_old.size[1];
    if (loop_ub_tmp - 1 >= 0) {
      memcpy(&nmpcSine_B.y_data[0], &obj->JacCineqTrans_old.data[0],
             static_cast<uint32_T>(loop_ub_tmp) * sizeof(real_T));
    }

    c = static_cast<uint16_T>(nVar);
    for (int32_T k = 0; k < c; k++) {
      nmpcSine_B.y_data[iCol_old + k] = JacCineqTrans_data[iCol + k];
    }

    if (loop_ub_tmp - 1 >= 0) {
      memcpy(&obj->JacCineqTrans_old.data[0], &nmpcSine_B.y_data[0],
             static_cast<uint32_T>(loop_ub_tmp) * sizeof(real_T));
    }

    iCol += ldJ;
    iCol_old += ldJ;
  }

  iCol = 0;
  iCol_old = 0;
  loop_ub_tmp = static_cast<uint16_T>(nVar);
  for (int32_T idx_col = 0; idx_col < 120; idx_col++) {
    for (b = 0; b < loop_ub_tmp; b++) {
      obj->JacCeqTrans_old.data[iCol_old + b] = JacCeqTrans_data[iCol + b];
    }

    iCol += ldJ;
    iCol_old = iCol;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_computeComplError(const int32_T
  fscales_lineq_constraint_size[1], const int32_T fscales_cineq_constraint_size
  [1], const real_T xCurrent[125], int32_T mIneq, const real_T cIneq_data[],
  const int32_T finiteLB_data[], int32_T mLB, const real_T lb[125], const
  int32_T finiteUB_data[], int32_T mUB, const real_T ub[125], const real_T
  lambda_data[], int32_T iL0)
{
  real_T nlpComplError;
  int32_T ubOffset;
  nlpComplError = 0.0;
  ubOffset = fscales_lineq_constraint_size[0];
  if ((mIneq + mLB) + mUB > 0) {
    real_T lbDelta;
    real_T lbLambda;
    real_T u0;
    int32_T b_idx;
    int32_T c;
    int32_T iLineq0;
    for (iLineq0 = 0; iLineq0 < ubOffset; iLineq0++) {
      lbDelta = cIneq_data[iLineq0];
      lbLambda = lambda_data[(iL0 + iLineq0) - 1];
      u0 = fabs(lbDelta);
      if ((!(u0 <= lbLambda)) && (!rtIsNaN(lbLambda))) {
        u0 = lbLambda;
      }

      lbDelta = fabs(lbDelta * lbLambda);
      if ((lbDelta <= u0) || rtIsNaN(u0)) {
        u0 = lbDelta;
      }

      if ((!(nlpComplError >= u0)) && (!rtIsNaN(u0))) {
        nlpComplError = u0;
      }
    }

    iLineq0 = (iL0 + fscales_lineq_constraint_size[0]) - 1;
    c = fscales_cineq_constraint_size[0];
    for (b_idx = 0; b_idx < c; b_idx++) {
      lbDelta = cIneq_data[ubOffset + b_idx];
      lbLambda = lambda_data[iLineq0 + b_idx];
      u0 = fabs(lbDelta);
      if ((!(u0 <= lbLambda)) && (!rtIsNaN(lbLambda))) {
        u0 = lbLambda;
      }

      lbDelta = fabs(lbDelta * lbLambda);
      if ((lbDelta <= u0) || rtIsNaN(u0)) {
        u0 = lbDelta;
      }

      if ((!(nlpComplError >= u0)) && (!rtIsNaN(u0))) {
        nlpComplError = u0;
      }
    }

    iLineq0 = (iL0 + mIneq) - 1;
    ubOffset = iLineq0 + mLB;
    c = static_cast<uint16_T>(mLB);
    for (b_idx = 0; b_idx < c; b_idx++) {
      int32_T finiteLB;
      finiteLB = finiteLB_data[b_idx];
      lbDelta = xCurrent[finiteLB - 1] - lb[finiteLB - 1];
      lbLambda = lambda_data[iLineq0 + b_idx];
      u0 = fabs(lbDelta);
      if ((!(u0 <= lbLambda)) && (!rtIsNaN(lbLambda))) {
        u0 = lbLambda;
      }

      lbDelta = fabs(lbDelta * lbLambda);
      if ((lbDelta <= u0) || rtIsNaN(u0)) {
        u0 = lbDelta;
      }

      if ((!(nlpComplError >= u0)) && (!rtIsNaN(u0))) {
        nlpComplError = u0;
      }
    }

    iLineq0 = static_cast<uint16_T>(mUB);
    for (c = 0; c < iLineq0; c++) {
      b_idx = finiteUB_data[c];
      lbDelta = ub[b_idx - 1] - xCurrent[b_idx - 1];
      lbLambda = lambda_data[ubOffset + c];
      u0 = fabs(lbDelta);
      if ((!(u0 <= lbLambda)) && (!rtIsNaN(lbLambda))) {
        u0 = lbLambda;
      }

      lbDelta = fabs(lbDelta * lbLambda);
      if ((lbDelta <= u0) || rtIsNaN(u0)) {
        u0 = lbDelta;
      }

      if ((!(nlpComplError >= u0)) && (!rtIsNaN(u0))) {
        nlpComplError = u0;
      }
    }
  }

  return nlpComplError;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_computeGradLag_e(real_T workspace_data[], int32_T ldA,
  int32_T nVar, const real_T grad_data[], int32_T mIneq, const real_T
  AineqTrans_data[], const real_T AeqTrans_data[], const int32_T
  finiteFixed_data[], int32_T mFixed, const int32_T finiteLB_data[], int32_T mLB,
  const int32_T finiteUB_data[], int32_T mUB, const real_T lambda_data[])
{
  int32_T c;
  int32_T f;
  int32_T finiteFixed;
  int32_T g;
  int32_T iL0;
  int32_T ix;
  memcpy(&workspace_data[0], &grad_data[0], static_cast<uint16_T>(nVar) * sizeof
         (real_T));
  c = static_cast<uint16_T>(mFixed);
  for (iL0 = 0; iL0 < c; iL0++) {
    finiteFixed = finiteFixed_data[iL0];
    workspace_data[finiteFixed - 1] += lambda_data[iL0];
  }

  ix = mFixed;
  f = ldA * 119 + 1;
  for (c = 1; ldA < 0 ? c >= f : c <= f; c += ldA) {
    g = (c + nVar) - 1;
    for (finiteFixed = c; finiteFixed <= g; finiteFixed++) {
      iL0 = finiteFixed - c;
      workspace_data[iL0] += AeqTrans_data[finiteFixed - 1] * lambda_data[ix];
    }

    ix++;
  }

  if (mIneq != 0) {
    ix = mFixed + 120;
    f = (mIneq - 1) * ldA + 1;
    for (c = 1; ldA < 0 ? c >= f : c <= f; c += ldA) {
      g = (c + nVar) - 1;
      for (finiteFixed = c; finiteFixed <= g; finiteFixed++) {
        iL0 = finiteFixed - c;
        workspace_data[iL0] += AineqTrans_data[finiteFixed - 1] * lambda_data[ix];
      }

      ix++;
    }
  }

  iL0 = (mFixed + mIneq) + 120;
  finiteFixed = static_cast<uint16_T>(mLB) - 1;
  for (c = 0; c <= finiteFixed; c++) {
    ix = finiteLB_data[c];
    workspace_data[ix - 1] -= lambda_data[iL0 + c];
  }

  if (static_cast<uint16_T>(mLB) - 1 >= 0) {
    iL0 += static_cast<uint16_T>(mLB);
  }

  finiteFixed = static_cast<uint16_T>(mUB) - 1;
  for (c = 0; c <= finiteFixed; c++) {
    ix = finiteUB_data[c];
    workspace_data[ix - 1] += lambda_data[iL0 + c];
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_computeDualFeasError_k(int32_T nVar, const real_T
  gradLag_data[], boolean_T *gradOK, real_T *val)
{
  int32_T idx;
  boolean_T exitg1;
  *gradOK = true;
  *val = 0.0;
  idx = 0;
  exitg1 = false;
  while ((!exitg1) && (idx <= static_cast<uint16_T>(nVar) - 1)) {
    *gradOK = ((!rtIsInf(gradLag_data[idx])) && (!rtIsNaN(gradLag_data[idx])));
    if (!*gradOK) {
      exitg1 = true;
    } else {
      real_T u1;
      u1 = fabs(gradLag_data[idx]);
      if ((!(*val >= u1)) && (!rtIsNaN(u1))) {
        *val = u1;
      }

      idx++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSi_updateWorkingSetForNewQP(const real_T xk[125],
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet, int32_T mIneq, int32_T
  mNonlinIneq, const real_T cIneq_data[], const real_T cEq[120], int32_T mLB,
  const real_T lb[125], int32_T mUB, const real_T ub[125], int32_T mFixed)
{
  real_T tmp[2];
  int32_T i;
  int32_T iEq0;
  int32_T iw0;
  int32_T nVar;
  nVar = WorkingSet->nVar;
  iw0 = WorkingSet->ldA * mFixed;
  iEq0 = 0;
  for (int32_T idx = 0; idx < 120; idx++) {
    real_T WorkingSet_beq;
    WorkingSet_beq = -cEq[idx];
    WorkingSet->beq[idx] = WorkingSet_beq;
    WorkingSet->bwset.data[mFixed + idx] = WorkingSet_beq;
    for (i = 0; i < nVar; i++) {
      WorkingSet->ATwset.data[iw0 + i] = WorkingSet->Aeq.data[iEq0 + i];
    }

    iw0 += WorkingSet->ldA;
    iEq0 += WorkingSet->ldA;
  }

  i = static_cast<uint16_T>(mIneq);
  iw0 = (static_cast<uint16_T>(mIneq) / 2) << 1;
  iEq0 = iw0 - 2;
  for (int32_T idx = 0; idx <= iEq0; idx += 2) {
    _mm_storeu_pd(&WorkingSet->bineq.data[idx], _mm_mul_pd(_mm_loadu_pd
      (&cIneq_data[idx]), _mm_set1_pd(-1.0)));
  }

  for (int32_T idx = iw0; idx < i; idx++) {
    WorkingSet->bineq.data[idx] = -cIneq_data[idx];
  }

  i = static_cast<uint16_T>(mLB);
  for (int32_T idx = 0; idx < i; idx++) {
    WorkingSet->lb.data[WorkingSet->indexLB.data[idx] - 1] = -lb
      [WorkingSet->indexLB.data[idx] - 1] + xk[WorkingSet->indexLB.data[idx] - 1];
  }

  i = static_cast<uint16_T>(mUB);
  for (int32_T idx = 0; idx < i; idx++) {
    WorkingSet->ub.data[WorkingSet->indexUB.data[idx] - 1] = ub
      [WorkingSet->indexUB.data[idx] - 1] - xk[WorkingSet->indexUB.data[idx] - 1];
  }

  i = static_cast<uint16_T>(mFixed);
  for (int32_T idx = 0; idx < i; idx++) {
    _mm_storeu_pd(&tmp[0], _mm_sub_pd(_mm_set1_pd(ub[WorkingSet->
      indexFixed.data[idx] - 1]), _mm_set1_pd(xk[WorkingSet->indexFixed.data[idx]
      - 1])));
    WorkingSet->ub.data[WorkingSet->indexFixed.data[idx] - 1] = tmp[0];
    WorkingSet->bwset.data[idx] = tmp[1];
  }

  if (WorkingSet->nActiveConstr > mFixed + 120) {
    iw0 = WorkingSet->nActiveConstr;
    for (int32_T idx = mFixed + 121; idx <= iw0; idx++) {
      switch (WorkingSet->Wid.data[idx - 1]) {
       case 4:
        WorkingSet->bwset.data[idx - 1] = WorkingSet->lb.data
          [WorkingSet->indexLB.data[WorkingSet->Wlocalidx.data[idx - 1] - 1] - 1];
        break;

       case 5:
        WorkingSet->bwset.data[idx - 1] = WorkingSet->ub.data
          [WorkingSet->indexUB.data[WorkingSet->Wlocalidx.data[idx - 1] - 1] - 1];
        break;

       default:
        {
          i = WorkingSet->Wlocalidx.data[idx - 1];
          WorkingSet->bwset.data[idx - 1] = WorkingSet->bineq.data[i - 1];
          if ((mNonlinIneq > 0) && (i > mIneq - mNonlinIneq)) {
            int32_T g;
            int32_T ix0;
            iEq0 = (idx - 1) * WorkingSet->ldA;
            ix0 = (i - 1) * WorkingSet->ldA;
            g = static_cast<uint16_T>(nVar);
            for (i = 0; i < g; i++) {
              WorkingSet->ATwset.data[iEq0 + i] = WorkingSet->Aineq.data[ix0 + i];
            }
          }
        }
        break;
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xswap(int32_T n, real_T x_data[], int32_T ix0, int32_T
  iy0)
{
  for (int32_T k = 0; k < n; k++) {
    real_T temp;
    int32_T temp_tmp;
    int32_T tmp;
    temp_tmp = (ix0 + k) - 1;
    temp = x_data[temp_tmp];
    tmp = (iy0 + k) - 1;
    x_data[temp_tmp] = x_data[tmp];
    x_data[tmp] = temp;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_xnrm2_n(int32_T n, const real_T x_data[], int32_T ix0)
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x_data[ix0 - 1]);
    } else {
      nmpcSine_B.scale_p = 3.3121686421112381E-170;
      nmpcSine_B.kend = ix0 + n;
      for (nmpcSine_B.k_pl = ix0; nmpcSine_B.k_pl < nmpcSine_B.kend;
           nmpcSine_B.k_pl++) {
        nmpcSine_B.absxk_b = fabs(x_data[nmpcSine_B.k_pl - 1]);
        if (nmpcSine_B.absxk_b > nmpcSine_B.scale_p) {
          nmpcSine_B.t_c = nmpcSine_B.scale_p / nmpcSine_B.absxk_b;
          y = y * nmpcSine_B.t_c * nmpcSine_B.t_c + 1.0;
          nmpcSine_B.scale_p = nmpcSine_B.absxk_b;
        } else {
          nmpcSine_B.t_c = nmpcSine_B.absxk_b / nmpcSine_B.scale_p;
          y += nmpcSine_B.t_c * nmpcSine_B.t_c;
        }
      }

      y = nmpcSine_B.scale_p * sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_xzlarfg(int32_T n, real_T *alpha1, real_T x_data[],
  int32_T ix0)
{
  __m128d tmp;
  real_T tau;
  tau = 0.0;
  if (n > 0) {
    nmpcSine_B.xnorm = nmpcSine_xnrm2_n(n - 1, x_data, ix0);
    if (nmpcSine_B.xnorm != 0.0) {
      nmpcSine_B.xnorm = nmpcSine_rt_hypotd_snf(*alpha1, nmpcSine_B.xnorm);
      if (*alpha1 >= 0.0) {
        nmpcSine_B.xnorm = -nmpcSine_B.xnorm;
      }

      if (fabs(nmpcSine_B.xnorm) < 1.0020841800044864E-292) {
        nmpcSine_B.knt_o = 0;
        nmpcSine_B.scalarLB_c = (ix0 + n) - 2;
        do {
          nmpcSine_B.knt_o++;
          nmpcSine_B.vectorUB_h = ((((nmpcSine_B.scalarLB_c - ix0) + 1) / 2) <<
            1) + ix0;
          nmpcSine_B.vectorUB_tmp_i = nmpcSine_B.vectorUB_h - 2;
          for (nmpcSine_B.d = ix0; nmpcSine_B.d <= nmpcSine_B.vectorUB_tmp_i;
               nmpcSine_B.d += 2) {
            tmp = _mm_loadu_pd(&x_data[nmpcSine_B.d - 1]);
            _mm_storeu_pd(&x_data[nmpcSine_B.d - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (nmpcSine_B.d = nmpcSine_B.vectorUB_h; nmpcSine_B.d <=
               nmpcSine_B.scalarLB_c; nmpcSine_B.d++) {
            x_data[nmpcSine_B.d - 1] *= 9.9792015476736E+291;
          }

          nmpcSine_B.xnorm *= 9.9792015476736E+291;
          *alpha1 *= 9.9792015476736E+291;
        } while ((fabs(nmpcSine_B.xnorm) < 1.0020841800044864E-292) &&
                 (nmpcSine_B.knt_o < 20));

        nmpcSine_B.xnorm = nmpcSine_rt_hypotd_snf(*alpha1, nmpcSine_xnrm2_n(n -
          1, x_data, ix0));
        if (*alpha1 >= 0.0) {
          nmpcSine_B.xnorm = -nmpcSine_B.xnorm;
        }

        tau = (nmpcSine_B.xnorm - *alpha1) / nmpcSine_B.xnorm;
        nmpcSine_B.a = 1.0 / (*alpha1 - nmpcSine_B.xnorm);
        for (nmpcSine_B.d = ix0; nmpcSine_B.d <= nmpcSine_B.vectorUB_tmp_i;
             nmpcSine_B.d += 2) {
          tmp = _mm_loadu_pd(&x_data[nmpcSine_B.d - 1]);
          _mm_storeu_pd(&x_data[nmpcSine_B.d - 1], _mm_mul_pd(tmp, _mm_set1_pd
            (nmpcSine_B.a)));
        }

        for (nmpcSine_B.d = nmpcSine_B.vectorUB_h; nmpcSine_B.d <=
             nmpcSine_B.scalarLB_c; nmpcSine_B.d++) {
          x_data[nmpcSine_B.d - 1] *= nmpcSine_B.a;
        }

        for (nmpcSine_B.d = 0; nmpcSine_B.d < nmpcSine_B.knt_o; nmpcSine_B.d++)
        {
          nmpcSine_B.xnorm *= 1.0020841800044864E-292;
        }

        *alpha1 = nmpcSine_B.xnorm;
      } else {
        tau = (nmpcSine_B.xnorm - *alpha1) / nmpcSine_B.xnorm;
        nmpcSine_B.a = 1.0 / (*alpha1 - nmpcSine_B.xnorm);
        nmpcSine_B.d = (ix0 + n) - 2;
        nmpcSine_B.scalarLB_c = ((((nmpcSine_B.d - ix0) + 1) / 2) << 1) + ix0;
        nmpcSine_B.vectorUB_h = nmpcSine_B.scalarLB_c - 2;
        for (nmpcSine_B.knt_o = ix0; nmpcSine_B.knt_o <= nmpcSine_B.vectorUB_h;
             nmpcSine_B.knt_o += 2) {
          tmp = _mm_loadu_pd(&x_data[nmpcSine_B.knt_o - 1]);
          _mm_storeu_pd(&x_data[nmpcSine_B.knt_o - 1], _mm_mul_pd(tmp,
            _mm_set1_pd(nmpcSine_B.a)));
        }

        for (nmpcSine_B.knt_o = nmpcSine_B.scalarLB_c; nmpcSine_B.knt_o <=
             nmpcSine_B.d; nmpcSine_B.knt_o++) {
          x_data[nmpcSine_B.knt_o - 1] *= nmpcSine_B.a;
        }

        *alpha1 = nmpcSine_B.xnorm;
      }
    }
  }

  return tau;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemv(int32_T m, int32_T n, const real_T A_data[],
  int32_T ia0, int32_T lda, const real_T x_data[], int32_T ix0, real_T y_data[])
{
  if (n != 0) {
    if (n - 1 >= 0) {
      memset(&y_data[0], 0, static_cast<uint32_T>(n) * sizeof(real_T));
    }

    nmpcSine_B.iy_e = 0;
    nmpcSine_B.b_jb = (n - 1) * lda + ia0;
    for (nmpcSine_B.b_iy_j = ia0; lda < 0 ? nmpcSine_B.b_iy_j >= nmpcSine_B.b_jb
         : nmpcSine_B.b_iy_j <= nmpcSine_B.b_jb; nmpcSine_B.b_iy_j += lda) {
      nmpcSine_B.c_n3 = 0.0;
      nmpcSine_B.d_g = (nmpcSine_B.b_iy_j + m) - 1;
      for (nmpcSine_B.ia_o = nmpcSine_B.b_iy_j; nmpcSine_B.ia_o <=
           nmpcSine_B.d_g; nmpcSine_B.ia_o++) {
        nmpcSine_B.c_n3 += x_data[((ix0 + nmpcSine_B.ia_o) - nmpcSine_B.b_iy_j)
          - 1] * A_data[nmpcSine_B.ia_o - 1];
      }

      y_data[nmpcSine_B.iy_e] += nmpcSine_B.c_n3;
      nmpcSine_B.iy_e++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgerc(int32_T m, int32_T n, real_T alpha1, int32_T ix0,
  const real_T y_data[], real_T A_data[], int32_T ia0, int32_T lda)
{
  if (!(alpha1 == 0.0)) {
    nmpcSine_B.jA_k = ia0;
    for (nmpcSine_B.j_p = 0; nmpcSine_B.j_p < n; nmpcSine_B.j_p++) {
      nmpcSine_B.temp_n = y_data[nmpcSine_B.j_p];
      if (nmpcSine_B.temp_n != 0.0) {
        nmpcSine_B.temp_n *= alpha1;
        nmpcSine_B.b_p5 = (m + nmpcSine_B.jA_k) - 1;
        for (nmpcSine_B.ijA_m = nmpcSine_B.jA_k; nmpcSine_B.ijA_m <=
             nmpcSine_B.b_p5; nmpcSine_B.ijA_m++) {
          A_data[nmpcSine_B.ijA_m - 1] += A_data[((ix0 + nmpcSine_B.ijA_m) -
            nmpcSine_B.jA_k) - 1] * nmpcSine_B.temp_n;
        }
      }

      nmpcSine_B.jA_k += lda;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xzlarf(int32_T m, int32_T n, int32_T iv0, real_T tau,
  real_T C_data[], int32_T ic0, int32_T ldc, real_T work_data[])
{
  int32_T exitg1;
  boolean_T exitg2;
  if (tau != 0.0) {
    nmpcSine_B.lastv = m;
    nmpcSine_B.lastc_g = iv0 + m;
    while ((nmpcSine_B.lastv > 0) && (C_data[nmpcSine_B.lastc_g - 2] == 0.0)) {
      nmpcSine_B.lastv--;
      nmpcSine_B.lastc_g--;
    }

    nmpcSine_B.lastc_g = n;
    exitg2 = false;
    while ((!exitg2) && (nmpcSine_B.lastc_g > 0)) {
      nmpcSine_B.coltop = (nmpcSine_B.lastc_g - 1) * ldc + ic0;
      nmpcSine_B.ia_j = nmpcSine_B.coltop;
      do {
        exitg1 = 0;
        if (nmpcSine_B.ia_j <= (nmpcSine_B.coltop + nmpcSine_B.lastv) - 1) {
          if (C_data[nmpcSine_B.ia_j - 1] != 0.0) {
            exitg1 = 1;
          } else {
            nmpcSine_B.ia_j++;
          }
        } else {
          nmpcSine_B.lastc_g--;
          exitg1 = 2;
        }
      } while (exitg1 == 0);

      if (exitg1 == 1) {
        exitg2 = true;
      }
    }
  } else {
    nmpcSine_B.lastv = 0;
    nmpcSine_B.lastc_g = 0;
  }

  if (nmpcSine_B.lastv > 0) {
    nmpcSine_xgemv(nmpcSine_B.lastv, nmpcSine_B.lastc_g, C_data, ic0, ldc,
                   C_data, iv0, work_data);
    nmpcSine_xgerc(nmpcSine_B.lastv, nmpcSine_B.lastc_g, -tau, iv0, work_data,
                   C_data, ic0, ldc);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_qrf(real_T A_data[], const int32_T A_size[2], int32_T m,
  int32_T n, int32_T nfxd, real_T tau_data[])
{
  nmpcSine_B.lda = A_size[0];
  nmpcSine_B.loop_ub_j0 = A_size[1];
  if (nmpcSine_B.loop_ub_j0 - 1 >= 0) {
    memset(&nmpcSine_B.work_data[0], 0, static_cast<uint32_T>
           (nmpcSine_B.loop_ub_j0) * sizeof(real_T));
  }

  nmpcSine_B.loop_ub_j0 = static_cast<uint16_T>(nfxd);
  for (nmpcSine_B.i_e0 = 0; nmpcSine_B.i_e0 < nmpcSine_B.loop_ub_j0;
       nmpcSine_B.i_e0++) {
    nmpcSine_B.ii_o = nmpcSine_B.i_e0 * nmpcSine_B.lda + nmpcSine_B.i_e0;
    nmpcSine_B.mmi_i = m - nmpcSine_B.i_e0;
    if (nmpcSine_B.i_e0 + 1 < m) {
      nmpcSine_B.b_atmp = A_data[nmpcSine_B.ii_o];
      nmpcSine_B.tau = nmpcSine_xzlarfg(nmpcSine_B.mmi_i, &nmpcSine_B.b_atmp,
        A_data, nmpcSine_B.ii_o + 2);
      tau_data[nmpcSine_B.i_e0] = nmpcSine_B.tau;
      A_data[nmpcSine_B.ii_o] = nmpcSine_B.b_atmp;
    } else {
      nmpcSine_B.tau = 0.0;
      tau_data[nmpcSine_B.i_e0] = 0.0;
    }

    if (nmpcSine_B.i_e0 + 1 < n) {
      nmpcSine_B.b_atmp = A_data[nmpcSine_B.ii_o];
      A_data[nmpcSine_B.ii_o] = 1.0;
      nmpcSine_xzlarf(nmpcSine_B.mmi_i, (n - nmpcSine_B.i_e0) - 1,
                      nmpcSine_B.ii_o + 1, nmpcSine_B.tau, A_data,
                      (nmpcSine_B.ii_o + nmpcSine_B.lda) + 1, nmpcSine_B.lda,
                      nmpcSine_B.work_data);
      A_data[nmpcSine_B.ii_o] = nmpcSine_B.b_atmp;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_qrpf(real_T A_data[], const int32_T A_size[2], int32_T m,
  int32_T n, int32_T nfxd, real_T tau_data[], int32_T jpvt_data[])
{
  nmpcSine_B.ma = A_size[0];
  if (m <= n) {
    nmpcSine_B.minmn_h = m;
  } else {
    nmpcSine_B.minmn_h = n;
  }

  nmpcSine_B.mmi = A_size[1];
  if (nmpcSine_B.mmi - 1 >= 0) {
    memset(&nmpcSine_B.work_data_b[0], 0, static_cast<uint32_T>(nmpcSine_B.mmi) *
           sizeof(real_T));
  }

  if (nmpcSine_B.mmi - 1 >= 0) {
    memset(&nmpcSine_B.vn1_data[0], 0, static_cast<uint32_T>(nmpcSine_B.mmi) *
           sizeof(real_T));
  }

  if (nmpcSine_B.mmi - 1 >= 0) {
    memset(&nmpcSine_B.vn2_data[0], 0, static_cast<uint32_T>(nmpcSine_B.mmi) *
           sizeof(real_T));
  }

  for (nmpcSine_B.j_j = nfxd + 1; nmpcSine_B.j_j <= n; nmpcSine_B.j_j++) {
    nmpcSine_B.smax = nmpcSine_xnrm2_n(m - nfxd, A_data, ((nmpcSine_B.j_j - 1) *
      nmpcSine_B.ma + nfxd) + 1);
    nmpcSine_B.vn1_data[nmpcSine_B.j_j - 1] = nmpcSine_B.smax;
    nmpcSine_B.vn2_data[nmpcSine_B.j_j - 1] = nmpcSine_B.smax;
  }

  for (nmpcSine_B.j_j = nfxd + 1; nmpcSine_B.j_j <= nmpcSine_B.minmn_h;
       nmpcSine_B.j_j++) {
    nmpcSine_B.itemp = (nmpcSine_B.j_j - 1) * nmpcSine_B.ma;
    nmpcSine_B.ii_k = (nmpcSine_B.itemp + nmpcSine_B.j_j) - 1;
    nmpcSine_B.nmi = n - nmpcSine_B.j_j;
    nmpcSine_B.mmi = m - nmpcSine_B.j_j;
    if (nmpcSine_B.nmi + 1 < 1) {
      nmpcSine_B.idxmax = -2;
    } else {
      nmpcSine_B.idxmax = -1;
      if (nmpcSine_B.nmi + 1 > 1) {
        nmpcSine_B.smax = fabs(nmpcSine_B.vn1_data[nmpcSine_B.j_j - 1]);
        for (nmpcSine_B.pvt = 2; nmpcSine_B.pvt <= nmpcSine_B.nmi + 1;
             nmpcSine_B.pvt++) {
          nmpcSine_B.s_k = fabs(nmpcSine_B.vn1_data[(nmpcSine_B.j_j +
            nmpcSine_B.pvt) - 2]);
          if (nmpcSine_B.s_k > nmpcSine_B.smax) {
            nmpcSine_B.idxmax = nmpcSine_B.pvt - 2;
            nmpcSine_B.smax = nmpcSine_B.s_k;
          }
        }
      }
    }

    nmpcSine_B.pvt = nmpcSine_B.j_j + nmpcSine_B.idxmax;
    if (nmpcSine_B.pvt + 1 != nmpcSine_B.j_j) {
      nmpcSine_xswap(m, A_data, nmpcSine_B.pvt * nmpcSine_B.ma + 1,
                     nmpcSine_B.itemp + 1);
      nmpcSine_B.itemp = jpvt_data[nmpcSine_B.pvt];
      jpvt_data[nmpcSine_B.pvt] = jpvt_data[nmpcSine_B.j_j - 1];
      jpvt_data[nmpcSine_B.j_j - 1] = nmpcSine_B.itemp;
      nmpcSine_B.vn1_data[nmpcSine_B.pvt] = nmpcSine_B.vn1_data[nmpcSine_B.j_j -
        1];
      nmpcSine_B.vn2_data[nmpcSine_B.pvt] = nmpcSine_B.vn2_data[nmpcSine_B.j_j -
        1];
    }

    if (nmpcSine_B.j_j < m) {
      nmpcSine_B.s_k = A_data[nmpcSine_B.ii_k];
      nmpcSine_B.smax = nmpcSine_xzlarfg(nmpcSine_B.mmi + 1, &nmpcSine_B.s_k,
        A_data, nmpcSine_B.ii_k + 2);
      tau_data[nmpcSine_B.j_j - 1] = nmpcSine_B.smax;
      A_data[nmpcSine_B.ii_k] = nmpcSine_B.s_k;
    } else {
      nmpcSine_B.smax = 0.0;
      tau_data[nmpcSine_B.j_j - 1] = 0.0;
    }

    if (nmpcSine_B.j_j < n) {
      nmpcSine_B.s_k = A_data[nmpcSine_B.ii_k];
      A_data[nmpcSine_B.ii_k] = 1.0;
      nmpcSine_xzlarf(nmpcSine_B.mmi + 1, nmpcSine_B.nmi, nmpcSine_B.ii_k + 1,
                      nmpcSine_B.smax, A_data, (nmpcSine_B.ii_k + nmpcSine_B.ma)
                      + 1, nmpcSine_B.ma, nmpcSine_B.work_data_b);
      A_data[nmpcSine_B.ii_k] = nmpcSine_B.s_k;
    }

    for (nmpcSine_B.ii_k = nmpcSine_B.j_j + 1; nmpcSine_B.ii_k <= n;
         nmpcSine_B.ii_k++) {
      nmpcSine_B.nmi = (nmpcSine_B.ii_k - 1) * nmpcSine_B.ma + nmpcSine_B.j_j;
      nmpcSine_B.smax = nmpcSine_B.vn1_data[nmpcSine_B.ii_k - 1];
      if (nmpcSine_B.smax != 0.0) {
        nmpcSine_B.s_k = fabs(A_data[nmpcSine_B.nmi - 1]) / nmpcSine_B.smax;
        nmpcSine_B.s_k = 1.0 - nmpcSine_B.s_k * nmpcSine_B.s_k;
        if (nmpcSine_B.s_k < 0.0) {
          nmpcSine_B.s_k = 0.0;
        }

        nmpcSine_B.temp2 = nmpcSine_B.smax / nmpcSine_B.vn2_data[nmpcSine_B.ii_k
          - 1];
        nmpcSine_B.temp2 = nmpcSine_B.temp2 * nmpcSine_B.temp2 * nmpcSine_B.s_k;
        if (nmpcSine_B.temp2 <= 1.4901161193847656E-8) {
          if (nmpcSine_B.j_j < m) {
            nmpcSine_B.smax = nmpcSine_xnrm2_n(nmpcSine_B.mmi, A_data,
              nmpcSine_B.nmi + 1);
            nmpcSine_B.vn1_data[nmpcSine_B.ii_k - 1] = nmpcSine_B.smax;
            nmpcSine_B.vn2_data[nmpcSine_B.ii_k - 1] = nmpcSine_B.smax;
          } else {
            nmpcSine_B.vn1_data[nmpcSine_B.ii_k - 1] = 0.0;
            nmpcSine_B.vn2_data[nmpcSine_B.ii_k - 1] = 0.0;
          }
        } else {
          nmpcSine_B.vn1_data[nmpcSine_B.ii_k - 1] = nmpcSine_B.smax * sqrt
            (nmpcSine_B.s_k);
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgeqp3(real_T A_data[], const int32_T A_size[2], int32_T
  m, int32_T n, int32_T jpvt_data[], real_T tau_data[], int32_T tau_size[1])
{
  static const int32_T offsets[4] = { 0, 1, 2, 3 };

  nmpcSine_B.ma_tmp = A_size[0];
  if (m <= n) {
    nmpcSine_B.minmn = m;
  } else {
    nmpcSine_B.minmn = n;
  }

  if (A_size[0] <= A_size[1]) {
    nmpcSine_B.nfxd = A_size[0];
  } else {
    nmpcSine_B.nfxd = A_size[1];
  }

  tau_size[0] = nmpcSine_B.nfxd;
  if (nmpcSine_B.nfxd - 1 >= 0) {
    memset(&tau_data[0], 0, static_cast<uint32_T>(nmpcSine_B.nfxd) * sizeof
           (real_T));
  }

  if (nmpcSine_B.minmn < 1) {
    nmpcSine_B.minmn = (n / 4) << 2;
    nmpcSine_B.b_j_e = nmpcSine_B.minmn - 4;
    for (nmpcSine_B.ma_tmp = 0; nmpcSine_B.ma_tmp <= nmpcSine_B.b_j_e;
         nmpcSine_B.ma_tmp += 4) {
      _mm_storeu_si128((__m128i *)&jpvt_data[nmpcSine_B.ma_tmp], _mm_add_epi32
                       (_mm_add_epi32(_mm_set1_epi32(nmpcSine_B.ma_tmp),
        _mm_loadu_si128((const __m128i *)&offsets[0])), _mm_set1_epi32(1)));
    }

    for (nmpcSine_B.ma_tmp = nmpcSine_B.minmn; nmpcSine_B.ma_tmp < n;
         nmpcSine_B.ma_tmp++) {
      jpvt_data[nmpcSine_B.ma_tmp] = nmpcSine_B.ma_tmp + 1;
    }
  } else {
    nmpcSine_B.nfxd = -1;
    for (nmpcSine_B.b_j_e = 0; nmpcSine_B.b_j_e < n; nmpcSine_B.b_j_e++) {
      if (jpvt_data[nmpcSine_B.b_j_e] != 0) {
        nmpcSine_B.nfxd++;
        if (nmpcSine_B.b_j_e + 1 != nmpcSine_B.nfxd + 1) {
          nmpcSine_xswap(m, A_data, nmpcSine_B.b_j_e * nmpcSine_B.ma_tmp + 1,
                         nmpcSine_B.nfxd * nmpcSine_B.ma_tmp + 1);
          jpvt_data[nmpcSine_B.b_j_e] = jpvt_data[nmpcSine_B.nfxd];
          jpvt_data[nmpcSine_B.nfxd] = nmpcSine_B.b_j_e + 1;
        } else {
          jpvt_data[nmpcSine_B.b_j_e] = nmpcSine_B.b_j_e + 1;
        }
      } else {
        jpvt_data[nmpcSine_B.b_j_e] = nmpcSine_B.b_j_e + 1;
      }
    }

    if (nmpcSine_B.nfxd + 1 <= nmpcSine_B.minmn) {
      nmpcSine_B.nfxd++;
    } else {
      nmpcSine_B.nfxd = nmpcSine_B.minmn;
    }

    nmpcSine_qrf(A_data, A_size, m, n, nmpcSine_B.nfxd, tau_data);
    if (nmpcSine_B.nfxd < nmpcSine_B.minmn) {
      nmpcSine_qrpf(A_data, A_size, m, n, nmpcSine_B.nfxd, tau_data, jpvt_data);
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factorQRE(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj, const
  real_T A_data[], int32_T mrows, int32_T ncols, int32_T ldA)
{
  boolean_T guard1;
  nmpcSine_B.idx_k = mrows * ncols;
  guard1 = false;
  if (nmpcSine_B.idx_k > 0) {
    for (nmpcSine_B.idx_k = 0; nmpcSine_B.idx_k < ncols; nmpcSine_B.idx_k++) {
      nmpcSine_B.ix0 = ldA * nmpcSine_B.idx_k;
      nmpcSine_B.iy0 = obj->ldq * nmpcSine_B.idx_k;
      nmpcSine_B.b_b = static_cast<uint16_T>(mrows);
      for (nmpcSine_B.k_hm = 0; nmpcSine_B.k_hm < nmpcSine_B.b_b;
           nmpcSine_B.k_hm++) {
        obj->QR.data[nmpcSine_B.iy0 + nmpcSine_B.k_hm] = A_data[nmpcSine_B.ix0 +
          nmpcSine_B.k_hm];
      }
    }

    guard1 = true;
  } else if (nmpcSine_B.idx_k == 0) {
    obj->mrows = mrows;
    obj->ncols = ncols;
    obj->minRowCol = 0;
  } else {
    guard1 = true;
  }

  if (guard1) {
    obj->usedPivoting = true;
    obj->mrows = mrows;
    obj->ncols = ncols;
    if (mrows <= ncols) {
      obj->minRowCol = mrows;
    } else {
      obj->minRowCol = ncols;
    }

    nmpcSine_xgeqp3(obj->QR.data, obj->QR.size, mrows, ncols, obj->jpvt.data,
                    obj->tau.data, obj->tau.size);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xorgqr(int32_T m, int32_T n, int32_T k, real_T A_data[],
  const int32_T A_size[2], int32_T lda, const real_T tau_data[])
{
  __m128d tmp;
  if (n >= 1) {
    for (nmpcSine_B.itau = k; nmpcSine_B.itau < n; nmpcSine_B.itau++) {
      nmpcSine_B.ia_c = nmpcSine_B.itau * lda;
      memset(&A_data[nmpcSine_B.ia_c], 0, static_cast<uint32_T>(m) * sizeof
             (real_T));
      A_data[nmpcSine_B.ia_c + nmpcSine_B.itau] = 1.0;
    }

    nmpcSine_B.itau = k - 1;
    nmpcSine_B.ia_c = A_size[1];
    if (nmpcSine_B.ia_c - 1 >= 0) {
      memset(&nmpcSine_B.work_data_p[0], 0, static_cast<uint32_T>
             (nmpcSine_B.ia_c) * sizeof(real_T));
    }

    nmpcSine_B.i_e = k;
    while (nmpcSine_B.i_e >= 1) {
      nmpcSine_B.ia_c = (nmpcSine_B.i_e - 1) * lda + nmpcSine_B.i_e;
      if (nmpcSine_B.i_e < n) {
        A_data[nmpcSine_B.ia_c - 1] = 1.0;
        nmpcSine_xzlarf((m - nmpcSine_B.i_e) + 1, n - nmpcSine_B.i_e,
                        nmpcSine_B.ia_c, tau_data[nmpcSine_B.itau], A_data,
                        nmpcSine_B.ia_c + lda, lda, nmpcSine_B.work_data_p);
      }

      if (nmpcSine_B.i_e < m) {
        nmpcSine_B.c_m = (nmpcSine_B.ia_c + m) - nmpcSine_B.i_e;
        nmpcSine_B.scalarLB_a = ((((nmpcSine_B.c_m - nmpcSine_B.ia_c) / 2) << 1)
          + nmpcSine_B.ia_c) + 1;
        nmpcSine_B.vectorUB_j = nmpcSine_B.scalarLB_a - 2;
        for (nmpcSine_B.b_k_o = nmpcSine_B.ia_c + 1; nmpcSine_B.b_k_o <=
             nmpcSine_B.vectorUB_j; nmpcSine_B.b_k_o += 2) {
          tmp = _mm_loadu_pd(&A_data[nmpcSine_B.b_k_o - 1]);
          _mm_storeu_pd(&A_data[nmpcSine_B.b_k_o - 1], _mm_mul_pd(tmp,
            _mm_set1_pd(-tau_data[nmpcSine_B.itau])));
        }

        for (nmpcSine_B.b_k_o = nmpcSine_B.scalarLB_a; nmpcSine_B.b_k_o <=
             nmpcSine_B.c_m; nmpcSine_B.b_k_o++) {
          A_data[nmpcSine_B.b_k_o - 1] *= -tau_data[nmpcSine_B.itau];
        }
      }

      A_data[nmpcSine_B.ia_c - 1] = 1.0 - tau_data[nmpcSine_B.itau];
      nmpcSine_B.c_m = static_cast<uint16_T>(nmpcSine_B.i_e - 1);
      for (nmpcSine_B.b_k_o = 0; nmpcSine_B.b_k_o < nmpcSine_B.c_m;
           nmpcSine_B.b_k_o++) {
        A_data[(nmpcSine_B.ia_c - nmpcSine_B.b_k_o) - 2] = 0.0;
      }

      nmpcSine_B.itau--;
      nmpcSine_B.i_e--;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_sortLambdaQP(real_T lambda_data[], int32_T
  WorkingSet_nActiveConstr, const int32_T WorkingSet_sizes[5], const int32_T
  WorkingSet_isActiveIdx[6], const int32_T WorkingSet_Wid_data[], const int32_T
  WorkingSet_Wlocalidx_data[], real_T workspace_data[])
{
  if (WorkingSet_nActiveConstr != 0) {
    int32_T currentMplier;
    int32_T idxOffset;
    int32_T mAll;
    mAll = (((WorkingSet_sizes[0] + WorkingSet_sizes[3]) + WorkingSet_sizes[4])
            + WorkingSet_sizes[2]) + 119;
    if (static_cast<uint16_T>(mAll + 1) - 1 >= 0) {
      memcpy(&workspace_data[0], &lambda_data[0], static_cast<uint16_T>(mAll + 1)
             * sizeof(real_T));
    }

    if (mAll >= 0) {
      memset(&lambda_data[0], 0, static_cast<uint32_T>(mAll + 1) * sizeof(real_T));
    }

    currentMplier = 0;
    mAll = 0;
    while ((mAll + 1 <= WorkingSet_nActiveConstr) && (WorkingSet_Wid_data[mAll] <=
            2)) {
      if (WorkingSet_Wid_data[mAll] == 1) {
        idxOffset = 1;
      } else {
        idxOffset = WorkingSet_isActiveIdx[1];
      }

      lambda_data[(idxOffset + WorkingSet_Wlocalidx_data[mAll]) - 2] =
        workspace_data[currentMplier];
      currentMplier++;
      mAll++;
    }

    while (mAll + 1 <= WorkingSet_nActiveConstr) {
      switch (WorkingSet_Wid_data[mAll]) {
       case 3:
        idxOffset = WorkingSet_isActiveIdx[2];
        break;

       case 4:
        idxOffset = WorkingSet_isActiveIdx[3];
        break;

       default:
        idxOffset = WorkingSet_isActiveIdx[4];
        break;
      }

      lambda_data[(idxOffset + WorkingSet_Wlocalidx_data[mAll]) - 2] =
        workspace_data[currentMplier];
      currentMplier++;
      mAll++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_test_exit_n(s7RdrPWkr8UPAUyTdDJkLaG_nmpcS_T *Flags,
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, sG8JZ69axY52WWR6RKyApQC_nmpcS_T
  *MeritFunction, const int32_T fscales_lineq_constraint_size[1], const int32_T
  fscales_cineq_constraint_size[1], s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
  *QRManager, const real_T lb[125], const real_T ub[125])
{
  real_T nlpComplErrorTmp;
  real_T s;
  real_T smax;
  real_T tmp;
  real_T tmp_0;
  int32_T c_ix;
  int32_T fullRank_R;
  int32_T iQR0;
  int32_T idx_max;
  int32_T mLambda;
  int32_T nVar;
  int32_T nVar_tmp;
  int32_T rankR;
  boolean_T dxTooSmall;
  boolean_T exitg1;
  boolean_T guard1;
  boolean_T isFeasible;
  nVar_tmp = WorkingSet->nVar;
  mLambda = (((WorkingSet->sizes[0] + WorkingSet->sizes[2]) + WorkingSet->sizes
              [3]) + WorkingSet->sizes[4]) + 119;
  if (mLambda >= 0) {
    memcpy(&TrialState->lambdaStopTest.data[0], &TrialState->lambdasqp.data[0],
           static_cast<uint32_T>(mLambda + 1) * sizeof(real_T));
  }

  nmpcSine_computeGradLag(TrialState->gradLag.data, WorkingSet->ldA,
    WorkingSet->nVar, TrialState->grad.data, WorkingSet->sizes[2],
    WorkingSet->Aineq.data, WorkingSet->Aeq.data, WorkingSet->indexFixed.data,
    WorkingSet->sizes[0], WorkingSet->indexLB.data, WorkingSet->sizes[3],
    WorkingSet->indexUB.data, WorkingSet->sizes[4],
    TrialState->lambdaStopTest.data);
  if (WorkingSet->nVar < 1) {
    idx_max = 0;
  } else {
    idx_max = 1;
    if (WorkingSet->nVar > 1) {
      smax = fabs(TrialState->grad.data[0]);
      for (nVar = 2; nVar <= nVar_tmp; nVar++) {
        s = fabs(TrialState->grad.data[nVar - 1]);
        if (s > smax) {
          idx_max = nVar;
          smax = s;
        }
      }
    }
  }

  smax = fabs(TrialState->grad.data[idx_max - 1]);
  if ((smax <= 1.0) || rtIsNaN(smax)) {
    smax = 1.0;
  }

  if (rtIsInf(smax)) {
    smax = 1.0;
  }

  MeritFunction->nlpPrimalFeasError = nmpcSine_computePrimalFeasError
    (TrialState->xstarsqp, WorkingSet->sizes[2] - TrialState->mNonlinIneq,
     TrialState->mNonlinIneq, TrialState->cIneq.data, TrialState->cEq,
     WorkingSet->indexLB.data, WorkingSet->sizes[3], lb,
     WorkingSet->indexUB.data, WorkingSet->sizes[4], ub);
  if (TrialState->sqpIterations == 0) {
    if ((MeritFunction->nlpPrimalFeasError <= 1.0) || rtIsNaN
        (MeritFunction->nlpPrimalFeasError)) {
      MeritFunction->feasRelativeFactor = 1.0;
    } else {
      MeritFunction->feasRelativeFactor = MeritFunction->nlpPrimalFeasError;
    }
  }

  isFeasible = (MeritFunction->nlpPrimalFeasError <= 0.001 *
                MeritFunction->feasRelativeFactor);
  nmpcSine_computeDualFeasError(WorkingSet->nVar, TrialState->gradLag.data,
    &Flags->gradOK, &MeritFunction->nlpDualFeasError);
  if (!Flags->gradOK) {
    Flags->done = true;
    if (isFeasible) {
      TrialState->sqpExitFlag = 2;
    } else {
      TrialState->sqpExitFlag = -2;
    }
  } else {
    MeritFunction->nlpComplError = nmpcSine_computeComplError
      (fscales_lineq_constraint_size, fscales_cineq_constraint_size,
       TrialState->xstarsqp, WorkingSet->sizes[2], TrialState->cIneq.data,
       WorkingSet->indexLB.data, WorkingSet->sizes[3], lb,
       WorkingSet->indexUB.data, WorkingSet->sizes[4], ub,
       TrialState->lambdaStopTest.data, WorkingSet->sizes[0] + 121);
    if ((MeritFunction->nlpDualFeasError >= MeritFunction->nlpComplError) ||
        rtIsNaN(MeritFunction->nlpComplError)) {
      MeritFunction->firstOrderOpt = MeritFunction->nlpDualFeasError;
    } else {
      MeritFunction->firstOrderOpt = MeritFunction->nlpComplError;
    }

    if (TrialState->sqpIterations > 1) {
      nmpcSine_computeGradLag_e(memspace->workspace_float.data, WorkingSet->ldA,
        WorkingSet->nVar, TrialState->grad.data, WorkingSet->sizes[2],
        WorkingSet->Aineq.data, WorkingSet->Aeq.data,
        WorkingSet->indexFixed.data, WorkingSet->sizes[0],
        WorkingSet->indexLB.data, WorkingSet->sizes[3], WorkingSet->indexUB.data,
        WorkingSet->sizes[4], TrialState->lambdaStopTestPrev.data);
      nmpcSine_computeDualFeasError_k(WorkingSet->nVar,
        memspace->workspace_float.data, &dxTooSmall, &s);
      nlpComplErrorTmp = nmpcSine_computeComplError
        (fscales_lineq_constraint_size, fscales_cineq_constraint_size,
         TrialState->xstarsqp, WorkingSet->sizes[2], TrialState->cIneq.data,
         WorkingSet->indexLB.data, WorkingSet->sizes[3], lb,
         WorkingSet->indexUB.data, WorkingSet->sizes[4], ub,
         TrialState->lambdaStopTestPrev.data, WorkingSet->sizes[0] + 121);
      if ((s < MeritFunction->nlpDualFeasError) && (nlpComplErrorTmp <
           MeritFunction->nlpComplError)) {
        MeritFunction->nlpDualFeasError = s;
        MeritFunction->nlpComplError = nlpComplErrorTmp;
        if (s >= nlpComplErrorTmp) {
          MeritFunction->firstOrderOpt = s;
        } else {
          MeritFunction->firstOrderOpt = nlpComplErrorTmp;
        }

        if (mLambda >= 0) {
          memcpy(&TrialState->lambdaStopTest.data[0],
                 &TrialState->lambdaStopTestPrev.data[0], static_cast<uint32_T>
                 (mLambda + 1) * sizeof(real_T));
        }
      } else if (mLambda >= 0) {
        memcpy(&TrialState->lambdaStopTestPrev.data[0],
               &TrialState->lambdaStopTest.data[0], static_cast<uint32_T>
               (mLambda + 1) * sizeof(real_T));
      }
    } else if (mLambda >= 0) {
      memcpy(&TrialState->lambdaStopTestPrev.data[0],
             &TrialState->lambdaStopTest.data[0], static_cast<uint32_T>(mLambda
              + 1) * sizeof(real_T));
    }

    if (isFeasible && (MeritFunction->nlpDualFeasError <= 0.001 * smax) &&
        (MeritFunction->nlpComplError <= 0.001 * smax)) {
      Flags->done = true;
      TrialState->sqpExitFlag = 1;
    } else {
      Flags->done = false;
      if (isFeasible && (TrialState->sqpFval < -1.0E+20)) {
        Flags->done = true;
        TrialState->sqpExitFlag = -3;
      } else {
        guard1 = false;
        if (TrialState->sqpIterations > 0) {
          dxTooSmall = true;
          nVar = 0;
          exitg1 = false;
          while ((!exitg1) && (nVar <= static_cast<uint16_T>(WorkingSet->nVar) -
                               1)) {
            s = fabs(TrialState->xstarsqp[nVar]);
            if ((s <= 1.0) || rtIsNaN(s)) {
              s = 1.0;
            }

            if (0.001 * s <= fabs(TrialState->delta_x.data[nVar])) {
              dxTooSmall = false;
              exitg1 = true;
            } else {
              nVar++;
            }
          }

          if (dxTooSmall) {
            if (!isFeasible) {
              if (Flags->stepType == 2) {
                Flags->done = true;
                TrialState->sqpExitFlag = -2;
              } else {
                Flags->stepType = 2;
                Flags->failedLineSearch = false;
                Flags->stepAccepted = false;
                guard1 = true;
              }
            } else if (WorkingSet->nActiveConstr == 0) {
              Flags->done = true;
              TrialState->sqpExitFlag = 2;
            } else {
              nmpcSi_updateWorkingSetForNewQP(TrialState->xstarsqp, WorkingSet,
                WorkingSet->sizes[2], TrialState->mNonlinIneq,
                TrialState->cIneq.data, TrialState->cEq, WorkingSet->sizes[3],
                lb, WorkingSet->sizes[4], ub, WorkingSet->sizes[0]);
              if (WorkingSet->nActiveConstr - 1 >= 0) {
                memset(&TrialState->lambda.data[0], 0, static_cast<uint32_T>
                       (WorkingSet->nActiveConstr) * sizeof(real_T));
              }

              nmpcSine_factorQRE(QRManager, WorkingSet->ATwset.data,
                                 WorkingSet->nVar, WorkingSet->nActiveConstr,
                                 WorkingSet->ldA);
              rankR = QRManager->minRowCol;
              for (idx_max = 0; idx_max < rankR; idx_max++) {
                iQR0 = QRManager->ldq * idx_max + idx_max;
                c_ix = QRManager->mrows - idx_max;
                if (c_ix - 2 >= 0) {
                  memcpy(&QRManager->Q.data[iQR0 + 1], &QRManager->QR.data[iQR0
                         + 1], static_cast<uint32_T>(c_ix - 1) * sizeof(real_T));
                }
              }

              nmpcSine_xorgqr(QRManager->mrows, QRManager->mrows,
                              QRManager->minRowCol, QRManager->Q.data,
                              QRManager->Q.size, QRManager->ldq,
                              QRManager->tau.data);
              fullRank_R = QRManager->ldq;
              memset(&memspace->workspace_float.data[0], 0, static_cast<uint16_T>
                     (WorkingSet->nVar) * sizeof(real_T));
              rankR = 0;
              iQR0 = (WorkingSet->nVar - 1) * QRManager->ldq + 1;
              for (idx_max = 1; fullRank_R < 0 ? idx_max >= iQR0 : idx_max <=
                   iQR0; idx_max += fullRank_R) {
                s = 0.0;
                c_ix = (idx_max + nVar_tmp) - 1;
                for (nVar = idx_max; nVar <= c_ix; nVar++) {
                  s += QRManager->Q.data[nVar - 1] * TrialState->grad.data[nVar
                    - idx_max];
                }

                memspace->workspace_float.data[rankR] -= s;
                rankR++;
              }

              if (WorkingSet->nVar >= WorkingSet->nActiveConstr) {
                nVar = WorkingSet->nVar;
              } else {
                nVar = WorkingSet->nActiveConstr;
              }

              s = static_cast<real_T>(nVar) * 2.2204460492503131E-16;
              if (s >= 1.4901161193847656E-8) {
                s = 1.4901161193847656E-8;
              }

              s *= fabs(QRManager->QR.data[0]);
              if (WorkingSet->nVar <= WorkingSet->nActiveConstr) {
                fullRank_R = WorkingSet->nVar;
              } else {
                fullRank_R = WorkingSet->nActiveConstr;
              }

              rankR = 0;
              nVar = 0;
              while ((rankR < fullRank_R) && (fabs(QRManager->QR.data[nVar]) > s))
              {
                rankR++;
                nVar = (nVar + QRManager->ldq) + 1;
              }

              if (rankR != 0) {
                for (nVar = rankR; nVar >= 1; nVar--) {
                  iQR0 = ((nVar - 1) * QRManager->ldq + nVar) - 2;
                  memspace->workspace_float.data[nVar - 1] /= QRManager->
                    QR.data[iQR0 + 1];
                  for (idx_max = 0; idx_max <= nVar - 2; idx_max++) {
                    c_ix = (nVar - idx_max) - 2;
                    memspace->workspace_float.data[c_ix] -=
                      memspace->workspace_float.data[nVar - 1] *
                      QRManager->QR.data[iQR0 - idx_max];
                  }
                }
              }

              if (WorkingSet->nActiveConstr <= fullRank_R) {
                fullRank_R = WorkingSet->nActiveConstr;
              }

              for (nVar = 0; nVar < fullRank_R; nVar++) {
                TrialState->lambda.data[QRManager->jpvt.data[nVar] - 1] =
                  memspace->workspace_float.data[nVar];
              }

              nmpcSine_sortLambdaQP(TrialState->lambda.data,
                                    WorkingSet->nActiveConstr, WorkingSet->sizes,
                                    WorkingSet->isActiveIdx,
                                    WorkingSet->Wid.data,
                                    WorkingSet->Wlocalidx.data,
                                    memspace->workspace_float.data);
              nmpcSine_computeGradLag_e(memspace->workspace_float.data,
                WorkingSet->ldA, WorkingSet->nVar, TrialState->grad.data,
                WorkingSet->sizes[2], WorkingSet->Aineq.data,
                WorkingSet->Aeq.data, WorkingSet->indexFixed.data,
                WorkingSet->sizes[0], WorkingSet->indexLB.data,
                WorkingSet->sizes[3], WorkingSet->indexUB.data,
                WorkingSet->sizes[4], TrialState->lambda.data);
              nmpcSine_computeDualFeasError_k(WorkingSet->nVar,
                memspace->workspace_float.data, &isFeasible, &s);
              nlpComplErrorTmp = nmpcSine_computeComplError
                (fscales_lineq_constraint_size, fscales_cineq_constraint_size,
                 TrialState->xstarsqp, WorkingSet->sizes[2],
                 TrialState->cIneq.data, WorkingSet->indexLB.data,
                 WorkingSet->sizes[3], lb, WorkingSet->indexUB.data,
                 WorkingSet->sizes[4], ub, TrialState->lambda.data,
                 WorkingSet->sizes[0] + 121);
              if ((s >= nlpComplErrorTmp) || rtIsNaN(nlpComplErrorTmp)) {
                tmp = s;
              } else {
                tmp = nlpComplErrorTmp;
              }

              if ((MeritFunction->nlpDualFeasError >=
                   MeritFunction->nlpComplError) || rtIsNaN
                  (MeritFunction->nlpComplError)) {
                tmp_0 = MeritFunction->nlpDualFeasError;
              } else {
                tmp_0 = MeritFunction->nlpComplError;
              }

              if (tmp <= tmp_0) {
                MeritFunction->nlpDualFeasError = s;
                MeritFunction->nlpComplError = nlpComplErrorTmp;
                MeritFunction->firstOrderOpt = tmp;
                if (mLambda >= 0) {
                  memcpy(&TrialState->lambdaStopTest.data[0],
                         &TrialState->lambda.data[0], static_cast<uint32_T>
                         (mLambda + 1) * sizeof(real_T));
                }
              }

              if ((MeritFunction->nlpDualFeasError <= 0.001 * smax) &&
                  (MeritFunction->nlpComplError <= 0.001 * smax)) {
                TrialState->sqpExitFlag = 1;
              } else {
                TrialState->sqpExitFlag = 2;
              }

              Flags->done = true;
              guard1 = true;
            }
          } else {
            guard1 = true;
          }
        } else {
          guard1 = true;
        }

        if (guard1) {
          if (TrialState->sqpIterations >= 18) {
            Flags->done = true;
            TrialState->sqpExitFlag = 0;
          } else if (TrialState->FunctionEvaluations >= 12500) {
            Flags->done = true;
            TrialState->sqpExitFlag = 0;
          }
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
boolean_T nmpcSine::nmpcSine_BFGSUpdate(int32_T nvar, real_T Bk[15625], const
  real_T sk_data[], real_T yk_data[], real_T workspace_data[])
{
  __m128d tmp_0;
  real_T curvatureS;
  real_T dotSY;
  real_T theta;
  int32_T b_iy;
  int32_T b_jA;
  int32_T ia;
  int32_T iac;
  int32_T ix;
  int32_T jA;
  int32_T k;
  boolean_T success;
  dotSY = 0.0;
  if (nvar >= 1) {
    for (k = 0; k < nvar; k++) {
      dotSY += sk_data[k] * yk_data[k];
    }
  }

  k = static_cast<uint16_T>(nvar);
  memset(&workspace_data[0], 0, static_cast<uint16_T>(nvar) * sizeof(real_T));
  ix = 0;
  b_jA = (nvar - 1) * 125 + 1;
  for (iac = 1; iac <= b_jA; iac += 125) {
    jA = (iac + nvar) - 1;
    for (ia = iac; ia <= jA; ia++) {
      b_iy = ia - iac;
      workspace_data[b_iy] += Bk[ia - 1] * sk_data[ix];
    }

    ix++;
  }

  curvatureS = 0.0;
  if (nvar >= 1) {
    for (b_iy = 0; b_iy < nvar; b_iy++) {
      curvatureS += workspace_data[b_iy] * sk_data[b_iy];
    }
  }

  if (dotSY < 0.2 * curvatureS) {
    theta = 0.8 * curvatureS / (curvatureS - dotSY);
    iac = (static_cast<uint16_T>(nvar) / 2) << 1;
    ia = iac - 2;
    for (b_iy = 0; b_iy <= ia; b_iy += 2) {
      tmp_0 = _mm_loadu_pd(&yk_data[b_iy]);
      _mm_storeu_pd(&yk_data[b_iy], _mm_mul_pd(_mm_set1_pd(theta), tmp_0));
    }

    for (b_iy = iac; b_iy < k; b_iy++) {
      yk_data[b_iy] *= theta;
    }

    dotSY = 0.0;
    for (b_iy = 0; b_iy < nvar; b_iy++) {
      if (!(1.0 - theta == 0.0)) {
        yk_data[b_iy] += (1.0 - theta) * workspace_data[b_iy];
      }

      if (nvar >= 1) {
        dotSY += sk_data[b_iy] * yk_data[b_iy];
      }
    }
  }

  success = ((curvatureS > 2.2204460492503131E-16) && (dotSY >
              2.2204460492503131E-16));
  if (success) {
    __m128d tmp;
    curvatureS = -1.0 / curvatureS;
    if (!(curvatureS == 0.0)) {
      jA = 1;
      for (ix = 0; ix < k; ix++) {
        if (workspace_data[ix] != 0.0) {
          theta = workspace_data[ix] * curvatureS;
          b_iy = (nvar + jA) - 1;
          iac = ((((b_iy - jA) + 1) / 2) << 1) + jA;
          ia = iac - 2;
          for (b_jA = jA; b_jA <= ia; b_jA += 2) {
            tmp_0 = _mm_loadu_pd(&workspace_data[b_jA - jA]);
            tmp = _mm_loadu_pd(&Bk[b_jA - 1]);
            _mm_storeu_pd(&Bk[b_jA - 1], _mm_add_pd(_mm_mul_pd(tmp_0,
              _mm_set1_pd(theta)), tmp));
          }

          for (b_jA = iac; b_jA <= b_iy; b_jA++) {
            Bk[b_jA - 1] += workspace_data[b_jA - jA] * theta;
          }
        }

        jA += 125;
      }
    }

    dotSY = 1.0 / dotSY;
    if (!(dotSY == 0.0)) {
      b_jA = 1;
      for (b_iy = 0; b_iy < k; b_iy++) {
        curvatureS = yk_data[b_iy];
        if (curvatureS != 0.0) {
          curvatureS *= dotSY;
          jA = (nvar + b_jA) - 1;
          iac = ((((jA - b_jA) + 1) / 2) << 1) + b_jA;
          ia = iac - 2;
          for (ix = b_jA; ix <= ia; ix += 2) {
            tmp_0 = _mm_loadu_pd(&yk_data[ix - b_jA]);
            tmp = _mm_loadu_pd(&Bk[ix - 1]);
            _mm_storeu_pd(&Bk[ix - 1], _mm_add_pd(_mm_mul_pd(tmp_0, _mm_set1_pd
              (curvatureS)), tmp));
          }

          for (ix = iac; ix <= jA; ix++) {
            Bk[ix - 1] += yk_data[ix - b_jA] * curvatureS;
          }
        }

        b_jA += 125;
      }
    }
  }

  return success;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factorQRE_o(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj,
  int32_T mrows, int32_T ncols)
{
  if (mrows * ncols == 0) {
    obj->mrows = mrows;
    obj->ncols = ncols;
    obj->minRowCol = 0;
  } else {
    obj->usedPivoting = true;
    obj->mrows = mrows;
    obj->ncols = ncols;
    if (mrows <= ncols) {
      obj->minRowCol = mrows;
    } else {
      obj->minRowCol = ncols;
    }

    nmpcSine_xgeqp3(obj->QR.data, obj->QR.size, mrows, ncols, obj->jpvt.data,
                    obj->tau.data, obj->tau.size);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_countsort(int32_T x_data[], int32_T xLen, int32_T
  workspace_data[], int32_T xMin, int32_T xMax)
{
  if ((xLen > 1) && (xMax > xMin)) {
    int32_T b_tmp;
    int32_T idxEnd;
    int32_T idxFill;
    int32_T idxStart;
    int32_T maxOffset;
    b_tmp = xMax - xMin;
    if (b_tmp >= 0) {
      memset(&workspace_data[0], 0, static_cast<uint32_T>(b_tmp + 1) * sizeof
             (int32_T));
    }

    maxOffset = b_tmp - 1;
    for (b_tmp = 0; b_tmp < xLen; b_tmp++) {
      idxFill = x_data[b_tmp] - xMin;
      workspace_data[idxFill]++;
    }

    for (b_tmp = 2; b_tmp <= maxOffset + 2; b_tmp++) {
      workspace_data[b_tmp - 1] += workspace_data[b_tmp - 2];
    }

    idxStart = 1;
    idxEnd = workspace_data[0];
    for (b_tmp = 0; b_tmp <= maxOffset; b_tmp++) {
      for (idxFill = idxStart; idxFill <= idxEnd; idxFill++) {
        x_data[idxFill - 1] = b_tmp + xMin;
      }

      idxStart = workspace_data[b_tmp] + 1;
      idxEnd = workspace_data[b_tmp + 1];
    }

    for (maxOffset = idxStart; maxOffset <= idxEnd; maxOffset++) {
      x_data[maxOffset - 1] = xMax;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_removeConstr(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
  int32_T idx_global)
{
  int32_T TYPE_tmp;
  TYPE_tmp = obj->Wid.data[idx_global - 1] - 1;
  obj->isActiveConstr.data[(obj->isActiveIdx[TYPE_tmp] + obj->
    Wlocalidx.data[idx_global - 1]) - 2] = false;
  if (idx_global < obj->nActiveConstr) {
    int32_T b;
    obj->Wid.data[idx_global - 1] = obj->Wid.data[obj->nActiveConstr - 1];
    obj->Wlocalidx.data[idx_global - 1] = obj->Wlocalidx.data[obj->nActiveConstr
      - 1];
    b = static_cast<uint16_T>(obj->nVar);
    for (int32_T idx = 0; idx < b; idx++) {
      obj->ATwset.data[idx + obj->ldA * (idx_global - 1)] = obj->ATwset.data
        [(obj->nActiveConstr - 1) * obj->ldA + idx];
    }

    obj->bwset.data[idx_global - 1] = obj->bwset.data[obj->nActiveConstr - 1];
  }

  obj->nActiveConstr--;
  obj->nWConstr[TYPE_tmp]--;
}

// Function for MATLAB Function: '<S24>/NLMPC'
int32_T nmpcSine::nmpcSine_RemoveDependentEq_(s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager)
{
  real_T qtb;
  real_T tol;
  int32_T i;
  int32_T iQR0;
  int32_T ix;
  int32_T mTotalWorkingEq_tmp;
  int32_T mWorkingFixed;
  int32_T nDepInd;
  int32_T totalEq_tmp;
  int32_T totalRank;
  boolean_T exitg1;
  mWorkingFixed = workingset->nWConstr[0];
  mTotalWorkingEq_tmp = workingset->nWConstr[0] + workingset->nWConstr[1];
  nDepInd = 0;
  if (mTotalWorkingEq_tmp > 0) {
    totalEq_tmp = static_cast<uint16_T>(workingset->nVar);
    for (totalRank = 0; totalRank < mTotalWorkingEq_tmp; totalRank++) {
      for (iQR0 = 0; iQR0 < totalEq_tmp; iQR0++) {
        qrmanager->QR.data[totalRank + qrmanager->ldq * iQR0] =
          workingset->ATwset.data[workingset->ldA * totalRank + iQR0];
      }
    }

    nDepInd = mTotalWorkingEq_tmp - workingset->nVar;
    if (nDepInd <= 0) {
      nDepInd = 0;
    }

    memset(&qrmanager->jpvt.data[0], 0, static_cast<uint16_T>(workingset->nVar) *
           sizeof(int32_T));
    nmpcSine_factorQRE_o(qrmanager, mTotalWorkingEq_tmp, workingset->nVar);
    if (mTotalWorkingEq_tmp >= workingset->nVar) {
      ix = mTotalWorkingEq_tmp;
    } else {
      ix = workingset->nVar;
    }

    tol = 2.2204460492503131E-15 * static_cast<real_T>(ix);
    if (tol >= 1.4901161193847656E-8) {
      tol = 1.4901161193847656E-8;
    }

    if (workingset->nVar <= mTotalWorkingEq_tmp) {
      totalRank = workingset->nVar;
    } else {
      totalRank = mTotalWorkingEq_tmp;
    }

    totalRank += (totalRank - 1) * qrmanager->ldq;
    while ((totalRank > 0) && (fabs(qrmanager->QR.data[totalRank - 1]) < fabs
            (qrmanager->QR.data[0]) * tol)) {
      totalRank = (totalRank - qrmanager->ldq) - 1;
      nDepInd++;
    }

    if (nDepInd > 0) {
      i = qrmanager->minRowCol;
      for (totalRank = 0; totalRank < i; totalRank++) {
        iQR0 = qrmanager->ldq * totalRank + totalRank;
        ix = qrmanager->mrows - totalRank;
        if (ix - 2 >= 0) {
          memcpy(&qrmanager->Q.data[iQR0 + 1], &qrmanager->QR.data[iQR0 + 1],
                 static_cast<uint32_T>(ix - 1) * sizeof(real_T));
        }
      }

      nmpcSine_xorgqr(qrmanager->mrows, qrmanager->mrows, qrmanager->minRowCol,
                      qrmanager->Q.data, qrmanager->Q.size, qrmanager->ldq,
                      qrmanager->tau.data);
      iQR0 = 0;
      exitg1 = false;
      while ((!exitg1) && (iQR0 <= nDepInd - 1)) {
        ix = ((mTotalWorkingEq_tmp - iQR0) - 1) * qrmanager->ldq;
        qtb = 0.0;
        for (totalRank = 0; totalRank < mTotalWorkingEq_tmp; totalRank++) {
          qtb += qrmanager->Q.data[ix + totalRank] * workingset->
            bwset.data[totalRank];
        }

        if (fabs(qtb) >= tol) {
          nDepInd = -1;
          exitg1 = true;
        } else {
          iQR0++;
        }
      }
    }

    if (nDepInd > 0) {
      for (totalRank = 0; totalRank < mTotalWorkingEq_tmp; totalRank++) {
        ix = qrmanager->ldq * totalRank;
        i = workingset->ldA * totalRank;
        for (iQR0 = 0; iQR0 < totalEq_tmp; iQR0++) {
          qrmanager->QR.data[ix + iQR0] = workingset->ATwset.data[i + iQR0];
        }
      }

      for (iQR0 = 0; iQR0 < mWorkingFixed; iQR0++) {
        qrmanager->jpvt.data[iQR0] = 1;
      }

      iQR0 = workingset->nWConstr[0] + 1;
      if (iQR0 <= mTotalWorkingEq_tmp) {
        memset(&qrmanager->jpvt.data[iQR0 + -1], 0, static_cast<uint32_T>
               ((mTotalWorkingEq_tmp - iQR0) + 1) * sizeof(int32_T));
      }

      nmpcSine_factorQRE_o(qrmanager, workingset->nVar, mTotalWorkingEq_tmp);
      for (mWorkingFixed = 0; mWorkingFixed < nDepInd; mWorkingFixed++) {
        memspace->workspace_int.data[mWorkingFixed] = qrmanager->jpvt.data
          [(mTotalWorkingEq_tmp - nDepInd) + mWorkingFixed];
      }

      nmpcSine_countsort(memspace->workspace_int.data, nDepInd,
                         memspace->workspace_sort.data, 1, mTotalWorkingEq_tmp);
      for (totalRank = nDepInd; totalRank >= 1; totalRank--) {
        mTotalWorkingEq_tmp = workingset->nWConstr[0] + workingset->nWConstr[1];
        if (mTotalWorkingEq_tmp != 0) {
          iQR0 = memspace->workspace_int.data[totalRank - 1];
          if (iQR0 <= mTotalWorkingEq_tmp) {
            if ((mTotalWorkingEq_tmp == workingset->nActiveConstr) ||
                (mTotalWorkingEq_tmp == iQR0)) {
              workingset->mEqRemoved++;
              workingset->indexEqRemoved[workingset->mEqRemoved - 1] =
                workingset->Wlocalidx.data[iQR0 - 1];
              nmpcSine_removeConstr(workingset, memspace->
                                    workspace_int.data[totalRank - 1]);
            } else {
              workingset->mEqRemoved++;
              i = workingset->Wid.data[iQR0 - 1] - 1;
              mWorkingFixed = workingset->Wlocalidx.data[iQR0 - 1];
              workingset->indexEqRemoved[workingset->mEqRemoved - 1] =
                mWorkingFixed;
              workingset->isActiveConstr.data[(workingset->isActiveIdx[i] +
                mWorkingFixed) - 2] = false;
              workingset->Wid.data[iQR0 - 1] = workingset->
                Wid.data[mTotalWorkingEq_tmp - 1];
              workingset->Wlocalidx.data[iQR0 - 1] = workingset->
                Wlocalidx.data[mTotalWorkingEq_tmp - 1];
              for (mWorkingFixed = 0; mWorkingFixed < totalEq_tmp; mWorkingFixed
                   ++) {
                workingset->ATwset.data[mWorkingFixed + workingset->ldA * (iQR0
                  - 1)] = workingset->ATwset.data[(mTotalWorkingEq_tmp - 1) *
                  workingset->ldA + mWorkingFixed];
              }

              workingset->bwset.data[iQR0 - 1] = workingset->
                bwset.data[mTotalWorkingEq_tmp - 1];
              workingset->Wid.data[mTotalWorkingEq_tmp - 1] =
                workingset->Wid.data[workingset->nActiveConstr - 1];
              workingset->Wlocalidx.data[mTotalWorkingEq_tmp - 1] =
                workingset->Wlocalidx.data[workingset->nActiveConstr - 1];
              for (mWorkingFixed = 0; mWorkingFixed < totalEq_tmp; mWorkingFixed
                   ++) {
                workingset->ATwset.data[mWorkingFixed + workingset->ldA *
                  (mTotalWorkingEq_tmp - 1)] = workingset->ATwset.data
                  [(workingset->nActiveConstr - 1) * workingset->ldA +
                  mWorkingFixed];
              }

              workingset->bwset.data[mTotalWorkingEq_tmp - 1] =
                workingset->bwset.data[workingset->nActiveConstr - 1];
              workingset->nActiveConstr--;
              workingset->nWConstr[i]--;
            }
          }
        }
      }
    }
  }

  return nDepInd;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_RemoveDependentIneq_(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager,
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace)
{
  real_T maxDiag;
  real_T tol;
  real_T u1;
  int32_T b_idx;
  int32_T d;
  int32_T idxDiag;
  int32_T idxDiag_tmp;
  int32_T ix0;
  int32_T iy0;
  int32_T nDepIneq;
  int32_T nFixedConstr;
  nDepIneq = workingset->nActiveConstr;
  nFixedConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  idxDiag_tmp = workingset->nVar;
  if ((workingset->nWConstr[2] + workingset->nWConstr[3]) + workingset->
      nWConstr[4] > 0) {
    if (workingset->nVar >= workingset->nActiveConstr) {
      b_idx = workingset->nVar;
    } else {
      b_idx = workingset->nActiveConstr;
    }

    tol = 2.2204460492503131E-15 * static_cast<real_T>(b_idx);
    if (tol >= 1.4901161193847656E-8) {
      tol = 1.4901161193847656E-8;
    }

    for (b_idx = 0; b_idx < nFixedConstr; b_idx++) {
      qrmanager->jpvt.data[b_idx] = 1;
    }

    if (nFixedConstr + 1 <= workingset->nActiveConstr) {
      memset(&qrmanager->jpvt.data[nFixedConstr], 0, static_cast<uint32_T>
             (workingset->nActiveConstr - nFixedConstr) * sizeof(int32_T));
    }

    for (b_idx = 0; b_idx < nDepIneq; b_idx++) {
      iy0 = qrmanager->ldq * b_idx;
      ix0 = workingset->ldA * b_idx;
      d = static_cast<uint16_T>(idxDiag_tmp);
      for (idxDiag = 0; idxDiag < d; idxDiag++) {
        qrmanager->QR.data[iy0 + idxDiag] = workingset->ATwset.data[ix0 +
          idxDiag];
      }
    }

    nmpcSine_factorQRE_o(qrmanager, workingset->nVar, workingset->nActiveConstr);
    nDepIneq = 0;
    for (b_idx = workingset->nActiveConstr - 1; b_idx + 1 > idxDiag_tmp; b_idx--)
    {
      nDepIneq++;
      memspace->workspace_int.data[nDepIneq - 1] = qrmanager->jpvt.data[b_idx];
    }

    maxDiag = fabs(qrmanager->QR.data[0]);
    for (idxDiag = 0; idxDiag < b_idx; idxDiag++) {
      u1 = fabs(qrmanager->QR.data[((idxDiag + 1) * qrmanager->ldq + idxDiag) +
                1]);
      if ((!(maxDiag >= u1)) && (!rtIsNaN(u1))) {
        maxDiag = u1;
      }
    }

    if (b_idx + 1 <= workingset->nVar) {
      idxDiag = qrmanager->ldq * b_idx + b_idx;
      while ((b_idx + 1 > nFixedConstr) && (fabs(qrmanager->QR.data[idxDiag]) <
              tol * maxDiag)) {
        nDepIneq++;
        memspace->workspace_int.data[nDepIneq - 1] = qrmanager->jpvt.data[b_idx];
        b_idx--;
        idxDiag = (idxDiag - qrmanager->ldq) - 1;
      }
    }

    nmpcSine_countsort(memspace->workspace_int.data, nDepIneq,
                       memspace->workspace_sort.data, nFixedConstr + 1,
                       workingset->nActiveConstr);
    for (nFixedConstr = nDepIneq; nFixedConstr >= 1; nFixedConstr--) {
      nmpcSine_removeConstr(workingset, memspace->
                            workspace_int.data[nFixedConstr - 1]);
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
int32_T nmpcSine::nmpcSine_rank(const real_T qrmanager_QR_data[], const int32_T
  qrmanager_QR_size[2], int32_T qrmanager_mrows, int32_T qrmanager_ncols)
{
  int32_T minmn;
  int32_T r;
  r = 0;
  if (qrmanager_mrows <= qrmanager_ncols) {
    minmn = qrmanager_mrows;
  } else {
    minmn = qrmanager_ncols;
  }

  if (minmn > 0) {
    real_T tol;
    int32_T tmp;
    if (qrmanager_mrows >= qrmanager_ncols) {
      tmp = qrmanager_mrows;
    } else {
      tmp = qrmanager_ncols;
    }

    tol = 2.2204460492503131E-15 * static_cast<real_T>(tmp);
    if (tol >= 1.4901161193847656E-8) {
      tol = 1.4901161193847656E-8;
    }

    tol *= fabs(qrmanager_QR_data[0]);
    while ((r < minmn) && (!(fabs(qrmanager_QR_data[qrmanager_QR_size[0] * r + r])
             <= tol))) {
      r++;
    }
  }

  return r;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemv_d(int32_T m, int32_T n, const real_T A_data[],
  int32_T lda, const real_T x_data[], real_T y_data[])
{
  if (n != 0) {
    nmpcSine_B.b_np = static_cast<uint16_T>(n);
    nmpcSine_B.scalarLB_h = (static_cast<uint16_T>(n) / 2) << 1;
    nmpcSine_B.vectorUB_i = nmpcSine_B.scalarLB_h - 2;
    for (nmpcSine_B.b_iy_j5 = 0; nmpcSine_B.b_iy_j5 <= nmpcSine_B.vectorUB_i;
         nmpcSine_B.b_iy_j5 += 2) {
      __m128d tmp;
      tmp = _mm_loadu_pd(&y_data[nmpcSine_B.b_iy_j5]);
      _mm_storeu_pd(&y_data[nmpcSine_B.b_iy_j5], _mm_mul_pd(tmp, _mm_set1_pd
        (-1.0)));
    }

    for (nmpcSine_B.b_iy_j5 = nmpcSine_B.scalarLB_h; nmpcSine_B.b_iy_j5 <
         nmpcSine_B.b_np; nmpcSine_B.b_iy_j5++) {
      y_data[nmpcSine_B.b_iy_j5] = -y_data[nmpcSine_B.b_iy_j5];
    }

    nmpcSine_B.scalarLB_h = 0;
    nmpcSine_B.vectorUB_i = (n - 1) * lda + 1;
    for (nmpcSine_B.b_iy_j5 = 1; lda < 0 ? nmpcSine_B.b_iy_j5 >=
         nmpcSine_B.vectorUB_i : nmpcSine_B.b_iy_j5 <= nmpcSine_B.vectorUB_i;
         nmpcSine_B.b_iy_j5 += lda) {
      nmpcSine_B.c_l = 0.0;
      nmpcSine_B.e_a3 = (nmpcSine_B.b_iy_j5 + m) - 1;
      for (nmpcSine_B.b_np = nmpcSine_B.b_iy_j5; nmpcSine_B.b_np <=
           nmpcSine_B.e_a3; nmpcSine_B.b_np++) {
        nmpcSine_B.c_l += x_data[nmpcSine_B.b_np - nmpcSine_B.b_iy_j5] *
          A_data[nmpcSine_B.b_np - 1];
      }

      y_data[nmpcSine_B.scalarLB_h] += nmpcSine_B.c_l;
      nmpcSine_B.scalarLB_h++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSi_maxConstraintViolation_i(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *obj, const real_T x_data[])
{
  real_T v;
  if (obj->probType == 2) {
    v = 0.0;
    nmpcSine_B.b_mIneq_k = obj->sizes[2];
    if (obj->Aineq.size[0] != 0) {
      if (nmpcSine_B.b_mIneq_k - 1 >= 0) {
        memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0],
               static_cast<uint32_T>(nmpcSine_B.b_mIneq_k) * sizeof(real_T));
      }

      nmpcSine_xgemv_d(125, obj->sizes[2], obj->Aineq.data, obj->ldA, x_data,
                       obj->maxConstrWorkspace.data);
      nmpcSine_B.f_d = static_cast<uint16_T>(obj->sizes[2]);
      for (nmpcSine_B.mIneq_l = 0; nmpcSine_B.mIneq_l < nmpcSine_B.f_d;
           nmpcSine_B.mIneq_l++) {
        nmpcSine_B.u1_g = obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_l] -
          x_data[nmpcSine_B.mIneq_l + 125];
        obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_l] = nmpcSine_B.u1_g;
        if ((!(v >= nmpcSine_B.u1_g)) && (!rtIsNaN(nmpcSine_B.u1_g))) {
          v = nmpcSine_B.u1_g;
        }
      }
    }

    memcpy(&obj->maxConstrWorkspace.data[0], &obj->beq[0], 120U * sizeof(real_T));
    nmpcSine_xgemv_d(125, 120, obj->Aeq.data, obj->ldA, x_data,
                     obj->maxConstrWorkspace.data);
    for (nmpcSine_B.mIneq_l = 0; nmpcSine_B.mIneq_l < 120; nmpcSine_B.mIneq_l++)
    {
      obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_l] =
        (obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_l] - x_data
         [(nmpcSine_B.b_mIneq_k + nmpcSine_B.mIneq_l) + 125]) + x_data
        [(obj->sizes[2] + nmpcSine_B.mIneq_l) + 245];
      nmpcSine_B.u1_g = fabs(obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_l]);
      if ((!(v >= nmpcSine_B.u1_g)) && (!rtIsNaN(nmpcSine_B.u1_g))) {
        v = nmpcSine_B.u1_g;
      }
    }
  } else {
    v = 0.0;
    nmpcSine_B.mIneq_l = obj->sizes[2];
    if (obj->Aineq.size[0] != 0) {
      if (nmpcSine_B.mIneq_l - 1 >= 0) {
        memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0],
               static_cast<uint32_T>(nmpcSine_B.mIneq_l) * sizeof(real_T));
      }

      nmpcSine_xgemv_d(obj->nVar, obj->sizes[2], obj->Aineq.data, obj->ldA,
                       x_data, obj->maxConstrWorkspace.data);
      nmpcSine_B.mIneq_l = static_cast<uint16_T>(obj->sizes[2]);
      for (nmpcSine_B.b_mIneq_k = 0; nmpcSine_B.b_mIneq_k < nmpcSine_B.mIneq_l;
           nmpcSine_B.b_mIneq_k++) {
        nmpcSine_B.u1_g = obj->maxConstrWorkspace.data[nmpcSine_B.b_mIneq_k];
        if ((!(v >= nmpcSine_B.u1_g)) && (!rtIsNaN(nmpcSine_B.u1_g))) {
          v = nmpcSine_B.u1_g;
        }
      }
    }

    memcpy(&obj->maxConstrWorkspace.data[0], &obj->beq[0], 120U * sizeof(real_T));
    nmpcSine_xgemv_d(obj->nVar, 120, obj->Aeq.data, obj->ldA, x_data,
                     obj->maxConstrWorkspace.data);
    for (nmpcSine_B.b_mIneq_k = 0; nmpcSine_B.b_mIneq_k < 120;
         nmpcSine_B.b_mIneq_k++) {
      nmpcSine_B.u1_g = fabs(obj->maxConstrWorkspace.data[nmpcSine_B.b_mIneq_k]);
      if ((!(v >= nmpcSine_B.u1_g)) && (!rtIsNaN(nmpcSine_B.u1_g))) {
        v = nmpcSine_B.u1_g;
      }
    }
  }

  if (obj->sizes[3] > 0) {
    nmpcSine_B.mIneq_l = static_cast<uint16_T>(obj->sizes[3]);
    for (nmpcSine_B.b_mIneq_k = 0; nmpcSine_B.b_mIneq_k < nmpcSine_B.mIneq_l;
         nmpcSine_B.b_mIneq_k++) {
      nmpcSine_B.u1_g = -x_data[obj->indexLB.data[nmpcSine_B.b_mIneq_k] - 1] -
        obj->lb.data[obj->indexLB.data[nmpcSine_B.b_mIneq_k] - 1];
      if ((!(v >= nmpcSine_B.u1_g)) && (!rtIsNaN(nmpcSine_B.u1_g))) {
        v = nmpcSine_B.u1_g;
      }
    }
  }

  if (obj->sizes[4] > 0) {
    nmpcSine_B.mIneq_l = static_cast<uint16_T>(obj->sizes[4]);
    for (nmpcSine_B.b_mIneq_k = 0; nmpcSine_B.b_mIneq_k < nmpcSine_B.mIneq_l;
         nmpcSine_B.b_mIneq_k++) {
      nmpcSine_B.u1_g = x_data[obj->indexUB.data[nmpcSine_B.b_mIneq_k] - 1] -
        obj->ub.data[obj->indexUB.data[nmpcSine_B.b_mIneq_k] - 1];
      if ((!(v >= nmpcSine_B.u1_g)) && (!rtIsNaN(nmpcSine_B.u1_g))) {
        v = nmpcSine_B.u1_g;
      }
    }
  }

  if (obj->sizes[0] > 0) {
    nmpcSine_B.mIneq_l = static_cast<uint16_T>(obj->sizes[0]);
    for (nmpcSine_B.b_mIneq_k = 0; nmpcSine_B.b_mIneq_k < nmpcSine_B.mIneq_l;
         nmpcSine_B.b_mIneq_k++) {
      nmpcSine_B.u1_g = fabs(x_data[obj->indexFixed.data[nmpcSine_B.b_mIneq_k] -
        1] - obj->ub.data[obj->indexFixed.data[nmpcSine_B.b_mIneq_k] - 1]);
      if ((!(v >= nmpcSine_B.u1_g)) && (!rtIsNaN(nmpcSine_B.u1_g))) {
        v = nmpcSine_B.u1_g;
      }
    }
  }

  return v;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemv_df(int32_T m, int32_T n, const real_T A_data[],
  int32_T lda, const real_T x_data[], int32_T ix0, real_T y_data[])
{
  if (n != 0) {
    nmpcSine_B.b_o = static_cast<uint16_T>(n);
    nmpcSine_B.scalarLB_g = (static_cast<uint16_T>(n) / 2) << 1;
    nmpcSine_B.vectorUB_jz = nmpcSine_B.scalarLB_g - 2;
    for (nmpcSine_B.b_iy_h = 0; nmpcSine_B.b_iy_h <= nmpcSine_B.vectorUB_jz;
         nmpcSine_B.b_iy_h += 2) {
      __m128d tmp;
      tmp = _mm_loadu_pd(&y_data[nmpcSine_B.b_iy_h]);
      _mm_storeu_pd(&y_data[nmpcSine_B.b_iy_h], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
    }

    for (nmpcSine_B.b_iy_h = nmpcSine_B.scalarLB_g; nmpcSine_B.b_iy_h <
         nmpcSine_B.b_o; nmpcSine_B.b_iy_h++) {
      y_data[nmpcSine_B.b_iy_h] = -y_data[nmpcSine_B.b_iy_h];
    }

    nmpcSine_B.scalarLB_g = 0;
    nmpcSine_B.vectorUB_jz = (n - 1) * lda + 1;
    for (nmpcSine_B.b_iy_h = 1; lda < 0 ? nmpcSine_B.b_iy_h >=
         nmpcSine_B.vectorUB_jz : nmpcSine_B.b_iy_h <= nmpcSine_B.vectorUB_jz;
         nmpcSine_B.b_iy_h += lda) {
      nmpcSine_B.c_d = 0.0;
      nmpcSine_B.e_j = (nmpcSine_B.b_iy_h + m) - 1;
      for (nmpcSine_B.b_o = nmpcSine_B.b_iy_h; nmpcSine_B.b_o <= nmpcSine_B.e_j;
           nmpcSine_B.b_o++) {
        nmpcSine_B.c_d += x_data[((ix0 + nmpcSine_B.b_o) - nmpcSine_B.b_iy_h) -
          1] * A_data[nmpcSine_B.b_o - 1];
      }

      y_data[nmpcSine_B.scalarLB_g] += nmpcSine_B.c_d;
      nmpcSine_B.scalarLB_g++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcS_maxConstraintViolation_in(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *obj, const real_T x_data[], int32_T ix0)
{
  real_T v;
  if (obj->probType == 2) {
    v = 0.0;
    nmpcSine_B.b_mIneq = obj->sizes[2];
    if (obj->Aineq.size[0] != 0) {
      if (nmpcSine_B.b_mIneq - 1 >= 0) {
        memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0],
               static_cast<uint32_T>(nmpcSine_B.b_mIneq) * sizeof(real_T));
      }

      nmpcSine_xgemv_df(125, obj->sizes[2], obj->Aineq.data, obj->ldA, x_data,
                        ix0, obj->maxConstrWorkspace.data);
      nmpcSine_B.f_h = static_cast<uint16_T>(obj->sizes[2]);
      for (nmpcSine_B.mIneq_a = 0; nmpcSine_B.mIneq_a < nmpcSine_B.f_h;
           nmpcSine_B.mIneq_a++) {
        nmpcSine_B.u1_i = obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_a] -
          x_data[(ix0 + nmpcSine_B.mIneq_a) + 124];
        obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_a] = nmpcSine_B.u1_i;
        if ((!(v >= nmpcSine_B.u1_i)) && (!rtIsNaN(nmpcSine_B.u1_i))) {
          v = nmpcSine_B.u1_i;
        }
      }
    }

    memcpy(&obj->maxConstrWorkspace.data[0], &obj->beq[0], 120U * sizeof(real_T));
    nmpcSine_xgemv_df(125, 120, obj->Aeq.data, obj->ldA, x_data, ix0,
                      obj->maxConstrWorkspace.data);
    for (nmpcSine_B.mIneq_a = 0; nmpcSine_B.mIneq_a < 120; nmpcSine_B.mIneq_a++)
    {
      obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_a] =
        (obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_a] - x_data[((ix0 +
           nmpcSine_B.b_mIneq) + nmpcSine_B.mIneq_a) + 124]) + x_data[((ix0 +
        obj->sizes[2]) + nmpcSine_B.mIneq_a) + 244];
      nmpcSine_B.u1_i = fabs(obj->maxConstrWorkspace.data[nmpcSine_B.mIneq_a]);
      if ((!(v >= nmpcSine_B.u1_i)) && (!rtIsNaN(nmpcSine_B.u1_i))) {
        v = nmpcSine_B.u1_i;
      }
    }
  } else {
    v = 0.0;
    nmpcSine_B.mIneq_a = obj->sizes[2];
    if (obj->Aineq.size[0] != 0) {
      if (nmpcSine_B.mIneq_a - 1 >= 0) {
        memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0],
               static_cast<uint32_T>(nmpcSine_B.mIneq_a) * sizeof(real_T));
      }

      nmpcSine_xgemv_df(obj->nVar, obj->sizes[2], obj->Aineq.data, obj->ldA,
                        x_data, ix0, obj->maxConstrWorkspace.data);
      nmpcSine_B.mIneq_a = static_cast<uint16_T>(obj->sizes[2]);
      for (nmpcSine_B.b_mIneq = 0; nmpcSine_B.b_mIneq < nmpcSine_B.mIneq_a;
           nmpcSine_B.b_mIneq++) {
        nmpcSine_B.u1_i = obj->maxConstrWorkspace.data[nmpcSine_B.b_mIneq];
        if ((!(v >= nmpcSine_B.u1_i)) && (!rtIsNaN(nmpcSine_B.u1_i))) {
          v = nmpcSine_B.u1_i;
        }
      }
    }

    memcpy(&obj->maxConstrWorkspace.data[0], &obj->beq[0], 120U * sizeof(real_T));
    nmpcSine_xgemv_df(obj->nVar, 120, obj->Aeq.data, obj->ldA, x_data, ix0,
                      obj->maxConstrWorkspace.data);
    for (nmpcSine_B.b_mIneq = 0; nmpcSine_B.b_mIneq < 120; nmpcSine_B.b_mIneq++)
    {
      nmpcSine_B.u1_i = fabs(obj->maxConstrWorkspace.data[nmpcSine_B.b_mIneq]);
      if ((!(v >= nmpcSine_B.u1_i)) && (!rtIsNaN(nmpcSine_B.u1_i))) {
        v = nmpcSine_B.u1_i;
      }
    }
  }

  if (obj->sizes[3] > 0) {
    nmpcSine_B.mIneq_a = static_cast<uint16_T>(obj->sizes[3]);
    for (nmpcSine_B.b_mIneq = 0; nmpcSine_B.b_mIneq < nmpcSine_B.mIneq_a;
         nmpcSine_B.b_mIneq++) {
      nmpcSine_B.u1_i = -x_data[(ix0 + obj->indexLB.data[nmpcSine_B.b_mIneq]) -
        2] - obj->lb.data[obj->indexLB.data[nmpcSine_B.b_mIneq] - 1];
      if ((!(v >= nmpcSine_B.u1_i)) && (!rtIsNaN(nmpcSine_B.u1_i))) {
        v = nmpcSine_B.u1_i;
      }
    }
  }

  if (obj->sizes[4] > 0) {
    nmpcSine_B.mIneq_a = static_cast<uint16_T>(obj->sizes[4]);
    for (nmpcSine_B.b_mIneq = 0; nmpcSine_B.b_mIneq < nmpcSine_B.mIneq_a;
         nmpcSine_B.b_mIneq++) {
      nmpcSine_B.u1_i = x_data[(ix0 + obj->indexUB.data[nmpcSine_B.b_mIneq]) - 2]
        - obj->ub.data[obj->indexUB.data[nmpcSine_B.b_mIneq] - 1];
      if ((!(v >= nmpcSine_B.u1_i)) && (!rtIsNaN(nmpcSine_B.u1_i))) {
        v = nmpcSine_B.u1_i;
      }
    }
  }

  if (obj->sizes[0] > 0) {
    nmpcSine_B.mIneq_a = static_cast<uint16_T>(obj->sizes[0]);
    for (nmpcSine_B.b_mIneq = 0; nmpcSine_B.b_mIneq < nmpcSine_B.mIneq_a;
         nmpcSine_B.b_mIneq++) {
      nmpcSine_B.u1_i = fabs(x_data[(ix0 + obj->
        indexFixed.data[nmpcSine_B.b_mIneq]) - 2] - obj->ub.data
        [obj->indexFixed.data[nmpcSine_B.b_mIneq] - 1]);
      if ((!(v >= nmpcSine_B.u1_i)) && (!rtIsNaN(nmpcSine_B.u1_i))) {
        v = nmpcSine_B.u1_i;
      }
    }
  }

  return v;
}

// Function for MATLAB Function: '<S24>/NLMPC'
boolean_T nmpcSine::nmpcSin_feasibleX0ForWorkingSet(real_T workspace_data[],
  const int32_T workspace_size[2], real_T xCurrent_data[],
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
  *qrmanager)
{
  __m128d tmp;
  __m128d tmp_0;
  int32_T exitg1;
  boolean_T nonDegenerateWset;
  nmpcSine_B.mWConstr = workingset->nActiveConstr;
  nmpcSine_B.nVar_n = workingset->nVar;
  nonDegenerateWset = true;
  if (workingset->nActiveConstr != 0) {
    if (workingset->nActiveConstr >= workingset->nVar) {
      nmpcSine_B.i_k_j = static_cast<uint16_T>(workingset->nVar);
      for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR < nmpcSine_B.i_k_j;
           nmpcSine_B.rankQR++) {
        nmpcSine_B.iQR0_f = qrmanager->ldq * nmpcSine_B.rankQR;
        for (nmpcSine_B.ldq = 0; nmpcSine_B.ldq < nmpcSine_B.mWConstr;
             nmpcSine_B.ldq++) {
          qrmanager->QR.data[nmpcSine_B.ldq + nmpcSine_B.iQR0_f] =
            workingset->ATwset.data[workingset->ldA * nmpcSine_B.ldq +
            nmpcSine_B.rankQR];
        }
      }

      memset(&qrmanager->jpvt.data[0], 0, static_cast<uint32_T>(nmpcSine_B.i_k_j)
             * sizeof(int32_T));
      nmpcSine_factorQRE_o(qrmanager, workingset->nActiveConstr,
                           workingset->nVar);
      nmpcSine_B.iQR0_f = qrmanager->minRowCol;
      for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR < nmpcSine_B.iQR0_f;
           nmpcSine_B.rankQR++) {
        nmpcSine_B.b_iQR0 = qrmanager->ldq * nmpcSine_B.rankQR +
          nmpcSine_B.rankQR;
        nmpcSine_B.jBcol = qrmanager->mrows - nmpcSine_B.rankQR;
        if (nmpcSine_B.jBcol - 2 >= 0) {
          memcpy(&qrmanager->Q.data[nmpcSine_B.b_iQR0 + 1], &qrmanager->
                 QR.data[nmpcSine_B.b_iQR0 + 1], static_cast<uint32_T>
                 (nmpcSine_B.jBcol - 1) * sizeof(real_T));
        }
      }

      nmpcSine_xorgqr(qrmanager->mrows, qrmanager->mrows, qrmanager->minRowCol,
                      qrmanager->Q.data, qrmanager->Q.size, qrmanager->ldq,
                      qrmanager->tau.data);
      nmpcSine_B.rankQR = nmpcSine_rank(qrmanager->QR.data, qrmanager->QR.size,
        qrmanager->mrows, qrmanager->ncols);
      for (nmpcSine_B.ldq = 0; nmpcSine_B.ldq < nmpcSine_B.mWConstr;
           nmpcSine_B.ldq++) {
        workspace_data[nmpcSine_B.ldq] = workingset->bwset.data[nmpcSine_B.ldq];
        workspace_data[nmpcSine_B.ldq + workspace_size[0]] =
          workingset->bwset.data[nmpcSine_B.ldq];
      }

      nmpcSine_B.b_iQR0 = workingset->ldA;
      nmpcSine_B.jBcol = 0;
      nmpcSine_B.iAcol = (workingset->nActiveConstr - 1) * workingset->ldA + 1;
      for (nmpcSine_B.ldq = 1; nmpcSine_B.b_iQR0 < 0 ? nmpcSine_B.ldq >=
           nmpcSine_B.iAcol : nmpcSine_B.ldq <= nmpcSine_B.iAcol; nmpcSine_B.ldq
           += nmpcSine_B.b_iQR0) {
        nmpcSine_B.temp = 0.0;
        nmpcSine_B.h_k = (nmpcSine_B.ldq + nmpcSine_B.nVar_n) - 1;
        for (nmpcSine_B.iQR0_f = nmpcSine_B.ldq; nmpcSine_B.iQR0_f <=
             nmpcSine_B.h_k; nmpcSine_B.iQR0_f++) {
          nmpcSine_B.temp += workingset->ATwset.data[nmpcSine_B.iQR0_f - 1] *
            xCurrent_data[nmpcSine_B.iQR0_f - nmpcSine_B.ldq];
        }

        workspace_data[nmpcSine_B.jBcol] -= nmpcSine_B.temp;
        nmpcSine_B.jBcol++;
      }

      nmpcSine_B.ldq = qrmanager->ldq;
      nmpcSine_B.iQR0_f = workspace_size[0];
      nmpcSine_B.b_iQR0 = workspace_size[0] * workspace_size[1];
      if (nmpcSine_B.b_iQR0 - 1 >= 0) {
        memcpy(&nmpcSine_B.B_data[0], &workspace_data[0], static_cast<uint32_T>
               (nmpcSine_B.b_iQR0) * sizeof(real_T));
      }

      for (nmpcSine_B.h_k = 0; nmpcSine_B.iQR0_f < 0 ? nmpcSine_B.h_k >=
           nmpcSine_B.iQR0_f : nmpcSine_B.h_k <= nmpcSine_B.iQR0_f;
           nmpcSine_B.h_k += nmpcSine_B.iQR0_f) {
        nmpcSine_B.iAcol = nmpcSine_B.h_k + nmpcSine_B.nVar_n;
        for (nmpcSine_B.jBcol = nmpcSine_B.h_k + 1; nmpcSine_B.jBcol <=
             nmpcSine_B.iAcol; nmpcSine_B.jBcol++) {
          workspace_data[nmpcSine_B.jBcol - 1] = 0.0;
        }
      }

      nmpcSine_B.b_br = -1;
      for (nmpcSine_B.h_k = 0; nmpcSine_B.iQR0_f < 0 ? nmpcSine_B.h_k >=
           nmpcSine_B.iQR0_f : nmpcSine_B.h_k <= nmpcSine_B.iQR0_f;
           nmpcSine_B.h_k += nmpcSine_B.iQR0_f) {
        nmpcSine_B.br = -1;
        nmpcSine_B.b_iQR0 = nmpcSine_B.h_k + nmpcSine_B.nVar_n;
        for (nmpcSine_B.jBcol = nmpcSine_B.h_k + 1; nmpcSine_B.jBcol <=
             nmpcSine_B.b_iQR0; nmpcSine_B.jBcol++) {
          nmpcSine_B.temp = 0.0;
          for (nmpcSine_B.iAcol = 0; nmpcSine_B.iAcol < nmpcSine_B.mWConstr;
               nmpcSine_B.iAcol++) {
            nmpcSine_B.temp += qrmanager->Q.data[(nmpcSine_B.iAcol +
              nmpcSine_B.br) + 1] * nmpcSine_B.B_data[(nmpcSine_B.iAcol +
              nmpcSine_B.b_br) + 1];
          }

          workspace_data[nmpcSine_B.jBcol - 1] += nmpcSine_B.temp;
          nmpcSine_B.br += nmpcSine_B.ldq;
        }

        nmpcSine_B.b_br += nmpcSine_B.iQR0_f;
      }

      for (nmpcSine_B.mWConstr = 0; nmpcSine_B.mWConstr < 2; nmpcSine_B.mWConstr
           ++) {
        nmpcSine_B.b_br = nmpcSine_B.iQR0_f * nmpcSine_B.mWConstr - 1;
        for (nmpcSine_B.jBcol = nmpcSine_B.rankQR; nmpcSine_B.jBcol >= 1;
             nmpcSine_B.jBcol--) {
          nmpcSine_B.br = (nmpcSine_B.jBcol - 1) * nmpcSine_B.ldq;
          nmpcSine_B.h_k = nmpcSine_B.jBcol + nmpcSine_B.b_br;
          nmpcSine_B.temp = workspace_data[nmpcSine_B.h_k];
          if (nmpcSine_B.temp != 0.0) {
            workspace_data[nmpcSine_B.h_k] = nmpcSine_B.temp /
              qrmanager->QR.data[(nmpcSine_B.jBcol + nmpcSine_B.br) - 1];
            nmpcSine_B.b_iQR0 = static_cast<uint16_T>(nmpcSine_B.jBcol - 1);
            for (nmpcSine_B.iAcol = 0; nmpcSine_B.iAcol < nmpcSine_B.b_iQR0;
                 nmpcSine_B.iAcol++) {
              nmpcSine_B.n_a = (nmpcSine_B.iAcol + nmpcSine_B.b_br) + 1;
              workspace_data[nmpcSine_B.n_a] -= qrmanager->
                QR.data[nmpcSine_B.iAcol + nmpcSine_B.br] *
                workspace_data[nmpcSine_B.h_k];
            }
          }
        }
      }

      for (nmpcSine_B.ldq = nmpcSine_B.rankQR + 1; nmpcSine_B.ldq <=
           nmpcSine_B.nVar_n; nmpcSine_B.ldq++) {
        workspace_data[nmpcSine_B.ldq - 1] = 0.0;
        workspace_data[(nmpcSine_B.ldq + workspace_size[0]) - 1] = 0.0;
      }

      for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR < nmpcSine_B.i_k_j;
           nmpcSine_B.rankQR++) {
        workspace_data[(qrmanager->jpvt.data[nmpcSine_B.rankQR] +
                        (workspace_size[0] << 1)) - 1] =
          workspace_data[nmpcSine_B.rankQR];
      }

      for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR < nmpcSine_B.i_k_j;
           nmpcSine_B.rankQR++) {
        workspace_data[nmpcSine_B.rankQR] = workspace_data[(workspace_size[0] <<
          1) + nmpcSine_B.rankQR];
      }

      for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR < nmpcSine_B.i_k_j;
           nmpcSine_B.rankQR++) {
        workspace_data[(qrmanager->jpvt.data[nmpcSine_B.rankQR] +
                        (workspace_size[0] << 1)) - 1] =
          workspace_data[nmpcSine_B.rankQR + workspace_size[0]];
      }

      for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR < nmpcSine_B.i_k_j;
           nmpcSine_B.rankQR++) {
        workspace_data[nmpcSine_B.rankQR + workspace_size[0]] = workspace_data
          [(workspace_size[0] << 1) + nmpcSine_B.rankQR];
      }
    } else {
      if (nmpcSine_B.mWConstr - 1 >= 0) {
        memset(&qrmanager->jpvt.data[0], 0, static_cast<uint32_T>
               (nmpcSine_B.mWConstr) * sizeof(int32_T));
      }

      nmpcSine_factorQRE(qrmanager, workingset->ATwset.data, workingset->nVar,
                         workingset->nActiveConstr, workingset->ldA);
      nmpcSine_B.ldq = qrmanager->minRowCol;
      for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR < nmpcSine_B.ldq;
           nmpcSine_B.rankQR++) {
        nmpcSine_B.iQR0_f = qrmanager->ldq * nmpcSine_B.rankQR +
          nmpcSine_B.rankQR;
        nmpcSine_B.b_iQR0 = qrmanager->mrows - nmpcSine_B.rankQR;
        if (nmpcSine_B.b_iQR0 - 2 >= 0) {
          memcpy(&qrmanager->Q.data[nmpcSine_B.iQR0_f + 1], &qrmanager->
                 QR.data[nmpcSine_B.iQR0_f + 1], static_cast<uint32_T>
                 (nmpcSine_B.b_iQR0 - 1) * sizeof(real_T));
        }
      }

      nmpcSine_xorgqr(qrmanager->mrows, qrmanager->minRowCol,
                      qrmanager->minRowCol, qrmanager->Q.data, qrmanager->Q.size,
                      qrmanager->ldq, qrmanager->tau.data);
      nmpcSine_B.rankQR = nmpcSine_rank(qrmanager->QR.data, qrmanager->QR.size,
        qrmanager->mrows, qrmanager->ncols);
      for (nmpcSine_B.i_k_j = 0; nmpcSine_B.i_k_j < nmpcSine_B.mWConstr;
           nmpcSine_B.i_k_j++) {
        nmpcSine_B.iQR0_f = (qrmanager->jpvt.data[nmpcSine_B.i_k_j] - 1) *
          workingset->ldA;
        nmpcSine_B.temp = 0.0;
        nmpcSine_B.b_iQR0 = static_cast<uint16_T>(nmpcSine_B.nVar_n);
        for (nmpcSine_B.ldq = 0; nmpcSine_B.ldq < nmpcSine_B.b_iQR0;
             nmpcSine_B.ldq++) {
          nmpcSine_B.temp += workingset->ATwset.data[nmpcSine_B.iQR0_f +
            nmpcSine_B.ldq] * xCurrent_data[nmpcSine_B.ldq];
        }

        nmpcSine_B.constrViolation_basicX = workingset->bwset.data
          [qrmanager->jpvt.data[nmpcSine_B.i_k_j] - 1];
        workspace_data[nmpcSine_B.i_k_j] = nmpcSine_B.constrViolation_basicX -
          nmpcSine_B.temp;
        workspace_data[nmpcSine_B.i_k_j + workspace_size[0]] =
          nmpcSine_B.constrViolation_basicX;
      }

      nmpcSine_B.ldq = qrmanager->ldq;
      nmpcSine_B.iQR0_f = workspace_size[0];
      nmpcSine_B.h_k = static_cast<uint16_T>(nmpcSine_B.rankQR);
      for (nmpcSine_B.i_k_j = 0; nmpcSine_B.i_k_j < 2; nmpcSine_B.i_k_j++) {
        nmpcSine_B.jBcol = nmpcSine_B.iQR0_f * nmpcSine_B.i_k_j;
        for (nmpcSine_B.mWConstr = 0; nmpcSine_B.mWConstr < nmpcSine_B.h_k;
             nmpcSine_B.mWConstr++) {
          nmpcSine_B.iAcol = nmpcSine_B.ldq * nmpcSine_B.mWConstr;
          nmpcSine_B.br = nmpcSine_B.mWConstr + nmpcSine_B.jBcol;
          nmpcSine_B.temp = workspace_data[nmpcSine_B.br];
          for (nmpcSine_B.b_iQR0 = 0; nmpcSine_B.b_iQR0 < nmpcSine_B.mWConstr;
               nmpcSine_B.b_iQR0++) {
            nmpcSine_B.temp -= qrmanager->QR.data[nmpcSine_B.b_iQR0 +
              nmpcSine_B.iAcol] * workspace_data[nmpcSine_B.b_iQR0 +
              nmpcSine_B.jBcol];
          }

          workspace_data[nmpcSine_B.br] = nmpcSine_B.temp / qrmanager->
            QR.data[nmpcSine_B.mWConstr + nmpcSine_B.iAcol];
        }
      }

      nmpcSine_B.mWConstr = workspace_size[0] * workspace_size[1];
      if (nmpcSine_B.mWConstr - 1 >= 0) {
        memcpy(&nmpcSine_B.B_data[0], &workspace_data[0], static_cast<uint32_T>
               (nmpcSine_B.mWConstr) * sizeof(real_T));
      }

      for (nmpcSine_B.i_k_j = 0; nmpcSine_B.iQR0_f < 0 ? nmpcSine_B.i_k_j >=
           nmpcSine_B.iQR0_f : nmpcSine_B.i_k_j <= nmpcSine_B.iQR0_f;
           nmpcSine_B.i_k_j += nmpcSine_B.iQR0_f) {
        nmpcSine_B.h_k = nmpcSine_B.i_k_j + nmpcSine_B.nVar_n;
        for (nmpcSine_B.mWConstr = nmpcSine_B.i_k_j + 1; nmpcSine_B.mWConstr <=
             nmpcSine_B.h_k; nmpcSine_B.mWConstr++) {
          workspace_data[nmpcSine_B.mWConstr - 1] = 0.0;
        }
      }

      nmpcSine_B.br = 1;
      for (nmpcSine_B.i_k_j = 0; nmpcSine_B.iQR0_f < 0 ? nmpcSine_B.i_k_j >=
           nmpcSine_B.iQR0_f : nmpcSine_B.i_k_j <= nmpcSine_B.iQR0_f;
           nmpcSine_B.i_k_j += nmpcSine_B.iQR0_f) {
        nmpcSine_B.b_iQR0 = -1;
        nmpcSine_B.n_a = nmpcSine_B.br + nmpcSine_B.rankQR;
        for (nmpcSine_B.mWConstr = nmpcSine_B.br; nmpcSine_B.mWConstr <
             nmpcSine_B.n_a; nmpcSine_B.mWConstr++) {
          nmpcSine_B.h_k = nmpcSine_B.i_k_j + nmpcSine_B.nVar_n;
          nmpcSine_B.iAcol = ((((nmpcSine_B.h_k - nmpcSine_B.i_k_j) / 2) << 1) +
                              nmpcSine_B.i_k_j) + 1;
          nmpcSine_B.b_br = nmpcSine_B.iAcol - 2;
          for (nmpcSine_B.jBcol = nmpcSine_B.i_k_j + 1; nmpcSine_B.jBcol <=
               nmpcSine_B.b_br; nmpcSine_B.jBcol += 2) {
            tmp = _mm_loadu_pd(&qrmanager->Q.data[(nmpcSine_B.b_iQR0 +
              nmpcSine_B.jBcol) - nmpcSine_B.i_k_j]);
            tmp_0 = _mm_loadu_pd(&workspace_data[nmpcSine_B.jBcol - 1]);
            _mm_storeu_pd(&workspace_data[nmpcSine_B.jBcol - 1], _mm_add_pd
                          (_mm_mul_pd(_mm_set1_pd
              (nmpcSine_B.B_data[nmpcSine_B.mWConstr - 1]), tmp), tmp_0));
          }

          for (nmpcSine_B.jBcol = nmpcSine_B.iAcol; nmpcSine_B.jBcol <=
               nmpcSine_B.h_k; nmpcSine_B.jBcol++) {
            workspace_data[nmpcSine_B.jBcol - 1] += qrmanager->Q.data
              [(nmpcSine_B.b_iQR0 + nmpcSine_B.jBcol) - nmpcSine_B.i_k_j] *
              nmpcSine_B.B_data[nmpcSine_B.mWConstr - 1];
          }

          nmpcSine_B.b_iQR0 += nmpcSine_B.ldq;
        }

        nmpcSine_B.br += nmpcSine_B.iQR0_f;
      }
    }

    nmpcSine_B.rankQR = 0;
    do {
      exitg1 = 0;
      if (nmpcSine_B.rankQR <= static_cast<uint16_T>(nmpcSine_B.nVar_n) - 1) {
        nmpcSine_B.temp = workspace_data[nmpcSine_B.rankQR];
        if (rtIsInf(nmpcSine_B.temp) || rtIsNaN(nmpcSine_B.temp)) {
          nonDegenerateWset = false;
          exitg1 = 1;
        } else {
          nmpcSine_B.temp = workspace_data[nmpcSine_B.rankQR + workspace_size[0]];
          if (rtIsInf(nmpcSine_B.temp) || rtIsNaN(nmpcSine_B.temp)) {
            nonDegenerateWset = false;
            exitg1 = 1;
          } else {
            nmpcSine_B.rankQR++;
          }
        }
      } else {
        nmpcSine_B.iAcol = (nmpcSine_B.nVar_n / 2) << 1;
        nmpcSine_B.b_br = nmpcSine_B.iAcol - 2;
        for (nmpcSine_B.rankQR = 0; nmpcSine_B.rankQR <= nmpcSine_B.b_br;
             nmpcSine_B.rankQR += 2) {
          tmp = _mm_loadu_pd(&workspace_data[nmpcSine_B.rankQR]);
          tmp_0 = _mm_loadu_pd(&xCurrent_data[nmpcSine_B.rankQR]);
          _mm_storeu_pd(&workspace_data[nmpcSine_B.rankQR], _mm_add_pd(tmp,
            tmp_0));
        }

        for (nmpcSine_B.rankQR = nmpcSine_B.iAcol; nmpcSine_B.rankQR <
             nmpcSine_B.nVar_n; nmpcSine_B.rankQR++) {
          workspace_data[nmpcSine_B.rankQR] += xCurrent_data[nmpcSine_B.rankQR];
        }

        nmpcSine_B.temp = nmpcSi_maxConstraintViolation_i(workingset,
          workspace_data);
        nmpcSine_B.constrViolation_basicX = nmpcS_maxConstraintViolation_in
          (workingset, workspace_data, workspace_size[0] + 1);
        if ((nmpcSine_B.temp <= 2.2204460492503131E-16) || (nmpcSine_B.temp <
             nmpcSine_B.constrViolation_basicX)) {
          nmpcSine_B.rankQR = static_cast<uint16_T>(nmpcSine_B.nVar_n);
          memcpy(&xCurrent_data[0], &workspace_data[0], static_cast<uint32_T>
                 (nmpcSine_B.rankQR) * sizeof(real_T));
        } else {
          nmpcSine_B.rankQR = static_cast<uint16_T>(nmpcSine_B.nVar_n);
          for (nmpcSine_B.nVar_n = 0; nmpcSine_B.nVar_n < nmpcSine_B.rankQR;
               nmpcSine_B.nVar_n++) {
            xCurrent_data[nmpcSine_B.nVar_n] = workspace_data[workspace_size[0]
              + nmpcSine_B.nVar_n];
          }
        }

        exitg1 = 1;
      }
    } while (exitg1 == 0);
  }

  return nonDegenerateWset;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_RemoveDependentIneq__b(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager,
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace)
{
  nmpcSine_B.nDepIneq = workingset->nActiveConstr;
  nmpcSine_B.nFixedConstr = workingset->nWConstr[0] + workingset->nWConstr[1];
  nmpcSine_B.idxDiag_tmp = workingset->nVar;
  if ((workingset->nWConstr[2] + workingset->nWConstr[3]) + workingset->
      nWConstr[4] > 0) {
    if (workingset->nVar >= workingset->nActiveConstr) {
      nmpcSine_B.b_idx_j = workingset->nVar;
    } else {
      nmpcSine_B.b_idx_j = workingset->nActiveConstr;
    }

    nmpcSine_B.u1_p = 2.2204460492503131E-15 * static_cast<real_T>
      (nmpcSine_B.b_idx_j);
    if (nmpcSine_B.u1_p >= 1.4901161193847656E-8) {
      nmpcSine_B.u1_p = 1.4901161193847656E-8;
    }

    nmpcSine_B.tol = 10.0 * nmpcSine_B.u1_p;
    for (nmpcSine_B.b_idx_j = 0; nmpcSine_B.b_idx_j < nmpcSine_B.nFixedConstr;
         nmpcSine_B.b_idx_j++) {
      qrmanager->jpvt.data[nmpcSine_B.b_idx_j] = 1;
    }

    if (nmpcSine_B.nFixedConstr + 1 <= nmpcSine_B.nDepIneq) {
      memset(&qrmanager->jpvt.data[nmpcSine_B.nFixedConstr], 0,
             static_cast<uint32_T>(nmpcSine_B.nDepIneq - nmpcSine_B.nFixedConstr)
             * sizeof(int32_T));
    }

    for (nmpcSine_B.b_idx_j = 0; nmpcSine_B.b_idx_j < nmpcSine_B.nDepIneq;
         nmpcSine_B.b_idx_j++) {
      nmpcSine_B.iy0_mi = qrmanager->ldq * nmpcSine_B.b_idx_j;
      nmpcSine_B.ix0_i = workingset->ldA * nmpcSine_B.b_idx_j;
      nmpcSine_B.d_f = static_cast<uint16_T>(nmpcSine_B.idxDiag_tmp);
      for (nmpcSine_B.idxDiag = 0; nmpcSine_B.idxDiag < nmpcSine_B.d_f;
           nmpcSine_B.idxDiag++) {
        qrmanager->QR.data[nmpcSine_B.iy0_mi + nmpcSine_B.idxDiag] =
          workingset->ATwset.data[nmpcSine_B.ix0_i + nmpcSine_B.idxDiag];
      }
    }

    nmpcSine_factorQRE_o(qrmanager, workingset->nVar, workingset->nActiveConstr);
    nmpcSine_B.nDepIneq = 0;
    nmpcSine_B.b_idx_j = workingset->nActiveConstr - 1;
    while (nmpcSine_B.b_idx_j + 1 > nmpcSine_B.idxDiag_tmp) {
      nmpcSine_B.nDepIneq++;
      memspace->workspace_int.data[nmpcSine_B.nDepIneq - 1] =
        qrmanager->jpvt.data[nmpcSine_B.b_idx_j];
      nmpcSine_B.b_idx_j--;
    }

    nmpcSine_B.maxDiag = fabs(qrmanager->QR.data[0]);
    for (nmpcSine_B.idxDiag = 0; nmpcSine_B.idxDiag < nmpcSine_B.b_idx_j;
         nmpcSine_B.idxDiag++) {
      nmpcSine_B.u1_p = fabs(qrmanager->QR.data[((nmpcSine_B.idxDiag + 1) *
        qrmanager->ldq + nmpcSine_B.idxDiag) + 1]);
      if ((!(nmpcSine_B.maxDiag >= nmpcSine_B.u1_p)) && (!rtIsNaN
           (nmpcSine_B.u1_p))) {
        nmpcSine_B.maxDiag = nmpcSine_B.u1_p;
      }
    }

    if (nmpcSine_B.b_idx_j + 1 <= workingset->nVar) {
      nmpcSine_B.idxDiag = qrmanager->ldq * nmpcSine_B.b_idx_j +
        nmpcSine_B.b_idx_j;
      while ((nmpcSine_B.b_idx_j + 1 > nmpcSine_B.nFixedConstr) && (fabs
              (qrmanager->QR.data[nmpcSine_B.idxDiag]) < nmpcSine_B.tol *
              nmpcSine_B.maxDiag)) {
        nmpcSine_B.nDepIneq++;
        memspace->workspace_int.data[nmpcSine_B.nDepIneq - 1] =
          qrmanager->jpvt.data[nmpcSine_B.b_idx_j];
        nmpcSine_B.b_idx_j--;
        nmpcSine_B.idxDiag = (nmpcSine_B.idxDiag - qrmanager->ldq) - 1;
      }
    }

    nmpcSine_countsort(memspace->workspace_int.data, nmpcSine_B.nDepIneq,
                       memspace->workspace_sort.data, nmpcSine_B.nFixedConstr +
                       1, workingset->nActiveConstr);
    for (nmpcSine_B.nFixedConstr = nmpcSine_B.nDepIneq; nmpcSine_B.nFixedConstr >=
         1; nmpcSine_B.nFixedConstr--) {
      nmpcSine_removeConstr(workingset, memspace->
                            workspace_int.data[nmpcSine_B.nFixedConstr - 1]);
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemv_dff(int32_T m, int32_T n, const real_T A_data[],
  int32_T lda, const real_T x_data[], real_T y_data[])
{
  if (n != 0) {
    nmpcSine_B.b_gbr = static_cast<uint16_T>(n);
    nmpcSine_B.scalarLB_k = (static_cast<uint16_T>(n) / 2) << 1;
    nmpcSine_B.vectorUB_d = nmpcSine_B.scalarLB_k - 2;
    for (nmpcSine_B.b_iy_g = 0; nmpcSine_B.b_iy_g <= nmpcSine_B.vectorUB_d;
         nmpcSine_B.b_iy_g += 2) {
      __m128d tmp;
      tmp = _mm_loadu_pd(&y_data[nmpcSine_B.b_iy_g]);
      _mm_storeu_pd(&y_data[nmpcSine_B.b_iy_g], _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
    }

    for (nmpcSine_B.b_iy_g = nmpcSine_B.scalarLB_k; nmpcSine_B.b_iy_g <
         nmpcSine_B.b_gbr; nmpcSine_B.b_iy_g++) {
      y_data[nmpcSine_B.b_iy_g] = -y_data[nmpcSine_B.b_iy_g];
    }

    nmpcSine_B.scalarLB_k = 0;
    nmpcSine_B.vectorUB_d = (n - 1) * lda + 1;
    for (nmpcSine_B.b_iy_g = 1; lda < 0 ? nmpcSine_B.b_iy_g >=
         nmpcSine_B.vectorUB_d : nmpcSine_B.b_iy_g <= nmpcSine_B.vectorUB_d;
         nmpcSine_B.b_iy_g += lda) {
      nmpcSine_B.c_h = 0.0;
      nmpcSine_B.e_c = (nmpcSine_B.b_iy_g + m) - 1;
      for (nmpcSine_B.b_gbr = nmpcSine_B.b_iy_g; nmpcSine_B.b_gbr <=
           nmpcSine_B.e_c; nmpcSine_B.b_gbr++) {
        nmpcSine_B.c_h += A_data[nmpcSine_B.b_gbr - 1] * x_data[nmpcSine_B.b_gbr
          - nmpcSine_B.b_iy_g];
      }

      y_data[nmpcSine_B.scalarLB_k] += nmpcSine_B.c_h;
      nmpcSine_B.scalarLB_k++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::maxConstraintViolation_AMats_no(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *obj, const real_T x_data[])
{
  real_T u1;
  real_T v;
  int32_T k;
  int32_T mIneq;
  v = 0.0;
  mIneq = obj->sizes[2];
  if (obj->Aineq.size[0] != 0) {
    if (mIneq - 1 >= 0) {
      memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0],
             static_cast<uint32_T>(mIneq) * sizeof(real_T));
    }

    nmpcSine_xgemv_dff(obj->nVar, obj->sizes[2], obj->Aineq.data, obj->ldA,
                       x_data, obj->maxConstrWorkspace.data);
    mIneq = static_cast<uint16_T>(obj->sizes[2]);
    for (k = 0; k < mIneq; k++) {
      u1 = obj->maxConstrWorkspace.data[k];
      if ((!(v >= u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  memcpy(&obj->maxConstrWorkspace.data[0], &obj->beq[0], 120U * sizeof(real_T));
  nmpcSine_xgemv_dff(obj->nVar, 120, obj->Aeq.data, obj->ldA, x_data,
                     obj->maxConstrWorkspace.data);
  for (k = 0; k < 120; k++) {
    u1 = fabs(obj->maxConstrWorkspace.data[k]);
    if ((!(v >= u1)) && (!rtIsNaN(u1))) {
      v = u1;
    }
  }

  return v;
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::maxConstraintViolation_AMats_re(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *obj, const real_T x_data[])
{
  real_T v;
  v = 0.0;
  nmpcSine_B.mIneq_iz = obj->sizes[2];
  if (obj->Aineq.size[0] != 0) {
    if (nmpcSine_B.mIneq_iz - 1 >= 0) {
      memcpy(&obj->maxConstrWorkspace.data[0], &obj->bineq.data[0],
             static_cast<uint32_T>(nmpcSine_B.mIneq_iz) * sizeof(real_T));
    }

    nmpcSine_xgemv_dff(125, obj->sizes[2], obj->Aineq.data, obj->ldA, x_data,
                       obj->maxConstrWorkspace.data);
    nmpcSine_B.b_gb = static_cast<uint16_T>(obj->sizes[2]);
    for (nmpcSine_B.k_ge = 0; nmpcSine_B.k_ge < nmpcSine_B.b_gb; nmpcSine_B.k_ge
         ++) {
      nmpcSine_B.obj_maxConstrWorkspace = obj->
        maxConstrWorkspace.data[nmpcSine_B.k_ge] - x_data[nmpcSine_B.k_ge + 125];
      obj->maxConstrWorkspace.data[nmpcSine_B.k_ge] =
        nmpcSine_B.obj_maxConstrWorkspace;
      if ((!(v >= nmpcSine_B.obj_maxConstrWorkspace)) && (!rtIsNaN
           (nmpcSine_B.obj_maxConstrWorkspace))) {
        v = nmpcSine_B.obj_maxConstrWorkspace;
      }
    }
  }

  memcpy(&obj->maxConstrWorkspace.data[0], &obj->beq[0], 120U * sizeof(real_T));
  nmpcSine_xgemv_dff(125, 120, obj->Aeq.data, obj->ldA, x_data,
                     obj->maxConstrWorkspace.data);
  for (nmpcSine_B.k_ge = 0; nmpcSine_B.k_ge < 120; nmpcSine_B.k_ge++) {
    obj->maxConstrWorkspace.data[nmpcSine_B.k_ge] =
      (obj->maxConstrWorkspace.data[nmpcSine_B.k_ge] - x_data
       [(nmpcSine_B.mIneq_iz + nmpcSine_B.k_ge) + 125]) + x_data[(obj->sizes[2]
      + nmpcSine_B.k_ge) + 245];
    nmpcSine_B.obj_maxConstrWorkspace = fabs(obj->
      maxConstrWorkspace.data[nmpcSine_B.k_ge]);
    if ((!(v >= nmpcSine_B.obj_maxConstrWorkspace)) && (!rtIsNaN
         (nmpcSine_B.obj_maxConstrWorkspace))) {
      v = nmpcSine_B.obj_maxConstrWorkspace;
    }
  }

  return v;
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpc_maxConstraintViolation_ini(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *obj, const real_T x_data[])
{
  real_T v;
  if (obj->probType == 2) {
    v = maxConstraintViolation_AMats_re(obj, x_data);
  } else {
    v = maxConstraintViolation_AMats_no(obj, x_data);
  }

  if (obj->sizes[3] > 0) {
    nmpcSine_B.b_gz = static_cast<uint16_T>(obj->sizes[3]);
    for (nmpcSine_B.idx_e = 0; nmpcSine_B.idx_e < nmpcSine_B.b_gz;
         nmpcSine_B.idx_e++) {
      nmpcSine_B.u1_i3 = -x_data[obj->indexLB.data[nmpcSine_B.idx_e] - 1] -
        obj->lb.data[obj->indexLB.data[nmpcSine_B.idx_e] - 1];
      if ((!(v >= nmpcSine_B.u1_i3)) && (!rtIsNaN(nmpcSine_B.u1_i3))) {
        v = nmpcSine_B.u1_i3;
      }
    }
  }

  if (obj->sizes[4] > 0) {
    nmpcSine_B.b_gz = static_cast<uint16_T>(obj->sizes[4]);
    for (nmpcSine_B.idx_e = 0; nmpcSine_B.idx_e < nmpcSine_B.b_gz;
         nmpcSine_B.idx_e++) {
      nmpcSine_B.u1_i3 = x_data[obj->indexUB.data[nmpcSine_B.idx_e] - 1] -
        obj->ub.data[obj->indexUB.data[nmpcSine_B.idx_e] - 1];
      if ((!(v >= nmpcSine_B.u1_i3)) && (!rtIsNaN(nmpcSine_B.u1_i3))) {
        v = nmpcSine_B.u1_i3;
      }
    }
  }

  if (obj->sizes[0] > 0) {
    nmpcSine_B.b_gz = static_cast<uint16_T>(obj->sizes[0]);
    for (nmpcSine_B.idx_e = 0; nmpcSine_B.idx_e < nmpcSine_B.b_gz;
         nmpcSine_B.idx_e++) {
      nmpcSine_B.u1_i3 = fabs(x_data[obj->indexFixed.data[nmpcSine_B.idx_e] - 1]
        - obj->ub.data[obj->indexFixed.data[nmpcSine_B.idx_e] - 1]);
      if ((!(v >= nmpcSine_B.u1_i3)) && (!rtIsNaN(nmpcSine_B.u1_i3))) {
        v = nmpcSine_B.u1_i3;
      }
    }
  }

  return v;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_PresolveWorkingSet(s_2COE1uYisQtyPYvPjrXP9G_nmpc_T
  *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
  *qrmanager)
{
  boolean_T guard1;
  solution->state = 82;
  nmpcSine_B.b_h = nmpcSine_RemoveDependentEq_(memspace, workingset, qrmanager);
  if ((nmpcSine_B.b_h != -1) && (workingset->nActiveConstr <= qrmanager->ldq)) {
    nmpcSine_RemoveDependentIneq_(workingset, qrmanager, memspace);
    nmpcSine_B.okWorkingSet = nmpcSin_feasibleX0ForWorkingSet
      (memspace->workspace_float.data, memspace->workspace_float.size,
       solution->xstar.data, workingset, qrmanager);
    guard1 = false;
    if (!nmpcSine_B.okWorkingSet) {
      nmpcSine_RemoveDependentIneq__b(workingset, qrmanager, memspace);
      nmpcSine_B.okWorkingSet = nmpcSin_feasibleX0ForWorkingSet
        (memspace->workspace_float.data, memspace->workspace_float.size,
         solution->xstar.data, workingset, qrmanager);
      if (!nmpcSine_B.okWorkingSet) {
        solution->state = -7;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      if (workingset->nWConstr[0] + workingset->nWConstr[1] == workingset->nVar)
      {
        nmpcSine_B.constrViolation = nmpc_maxConstraintViolation_ini(workingset,
          solution->xstar.data);
        if (nmpcSine_B.constrViolation > 0.001) {
          solution->state = -2;
        }
      }
    }
  } else {
    solution->state = -3;
    nmpcSine_B.idxStartIneq_tmp_l = workingset->nWConstr[0] +
      workingset->nWConstr[1];
    nmpcSine_B.idxStartIneq_c = nmpcSine_B.idxStartIneq_tmp_l + 1;
    nmpcSine_B.idxEndIneq_a = workingset->nActiveConstr;
    for (nmpcSine_B.b_h = nmpcSine_B.idxStartIneq_c; nmpcSine_B.b_h <=
         nmpcSine_B.idxEndIneq_a; nmpcSine_B.b_h++) {
      workingset->isActiveConstr.data[(workingset->isActiveIdx
        [workingset->Wid.data[nmpcSine_B.b_h - 1] - 1] +
        workingset->Wlocalidx.data[nmpcSine_B.b_h - 1]) - 2] = false;
    }

    workingset->nWConstr[2] = 0;
    workingset->nWConstr[3] = 0;
    workingset->nWConstr[4] = 0;
    workingset->nActiveConstr = nmpcSine_B.idxStartIneq_tmp_l;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemv_dffo(int32_T m, int32_T n, const real_T A[15625],
  int32_T lda, const real_T x_data[], real_T y_data[])
{
  if ((m != 0) && (n != 0)) {
    if (m - 1 >= 0) {
      memset(&y_data[0], 0, static_cast<uint32_T>(m) * sizeof(real_T));
    }

    nmpcSine_B.ix_p = 0;
    nmpcSine_B.b_l = (n - 1) * lda + 1;
    for (nmpcSine_B.b_iy_l = 1; lda < 0 ? nmpcSine_B.b_iy_l >= nmpcSine_B.b_l :
         nmpcSine_B.b_iy_l <= nmpcSine_B.b_l; nmpcSine_B.b_iy_l += lda) {
      nmpcSine_B.c_hb = (nmpcSine_B.b_iy_l + m) - 1;
      for (nmpcSine_B.ia_cg = nmpcSine_B.b_iy_l; nmpcSine_B.ia_cg <=
           nmpcSine_B.c_hb; nmpcSine_B.ia_cg++) {
        nmpcSine_B.i8 = nmpcSine_B.ia_cg - nmpcSine_B.b_iy_l;
        y_data[nmpcSine_B.i8] += A[nmpcSine_B.ia_cg - 1] *
          x_data[nmpcSine_B.ix_p];
      }

      nmpcSine_B.ix_p++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_computeGrad_StoreHx(s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *obj,
  const real_T H[15625], const real_T f_data[], const real_T x_data[])
{
  __m128d tmp;
  switch (obj->objtype) {
   case 5:
    nmpcSine_B.b_j1 = obj->nvar;
    if (nmpcSine_B.b_j1 - 2 >= 0) {
      memset(&obj->grad.data[0], 0, static_cast<uint32_T>(nmpcSine_B.b_j1 - 1) *
             sizeof(real_T));
    }

    obj->grad.data[obj->nvar - 1] = obj->gammaScalar;
    break;

   case 3:
    nmpcSine_xgemv_dffo(obj->nvar, obj->nvar, H, obj->nvar, x_data, obj->Hx.data);
    nmpcSine_B.b_j1 = obj->nvar;
    if (nmpcSine_B.b_j1 - 1 >= 0) {
      memcpy(&obj->grad.data[0], &obj->Hx.data[0], static_cast<uint32_T>
             (nmpcSine_B.b_j1) * sizeof(real_T));
    }

    if (obj->hasLinear) {
      nmpcSine_B.iy_p = obj->grad.size[0];
      nmpcSine_B.i_n = obj->grad.size[0];
      if (nmpcSine_B.i_n - 1 >= 0) {
        memcpy(&nmpcSine_B.y_data_cv[0], &obj->grad.data[0], static_cast<
               uint32_T>(nmpcSine_B.i_n) * sizeof(real_T));
      }

      if (obj->nvar >= 1) {
        nmpcSine_B.ixlast_j = obj->nvar;
        nmpcSine_B.i_n = (obj->nvar / 2) << 1;
        nmpcSine_B.b_j1 = nmpcSine_B.i_n - 2;
        for (nmpcSine_B.k_lc = 0; nmpcSine_B.k_lc <= nmpcSine_B.b_j1;
             nmpcSine_B.k_lc += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.y_data_cv[nmpcSine_B.k_lc]);
          _mm_storeu_pd(&nmpcSine_B.y_data_cv[nmpcSine_B.k_lc], _mm_add_pd(tmp,
            _mm_loadu_pd(&f_data[nmpcSine_B.k_lc])));
        }

        for (nmpcSine_B.k_lc = nmpcSine_B.i_n; nmpcSine_B.k_lc <
             nmpcSine_B.ixlast_j; nmpcSine_B.k_lc++) {
          nmpcSine_B.y_data_cv[nmpcSine_B.k_lc] += f_data[nmpcSine_B.k_lc];
        }
      }

      if (nmpcSine_B.iy_p - 1 >= 0) {
        memcpy(&obj->grad.data[0], &nmpcSine_B.y_data_cv[0],
               static_cast<uint32_T>(nmpcSine_B.iy_p) * sizeof(real_T));
      }
    }
    break;

   case 4:
    nmpcSine_B.ixlast_j = obj->maxVar;
    nmpcSine_xgemv_dffo(obj->nvar, obj->nvar, H, obj->nvar, x_data, obj->Hx.data);
    nmpcSine_B.iy_p = obj->nvar + 1;
    nmpcSine_B.i_tmp = (obj->maxVar - obj->nvar) - 1;
    nmpcSine_B.i_n = (((nmpcSine_B.i_tmp / 2) << 1) + obj->nvar) + 1;
    nmpcSine_B.b_j1 = nmpcSine_B.i_n - 2;
    for (nmpcSine_B.k_lc = nmpcSine_B.iy_p; nmpcSine_B.k_lc <= nmpcSine_B.b_j1;
         nmpcSine_B.k_lc += 2) {
      _mm_storeu_pd(&obj->Hx.data[nmpcSine_B.k_lc - 1], _mm_mul_pd(_mm_loadu_pd(
        &x_data[nmpcSine_B.k_lc - 1]), _mm_set1_pd(obj->beta)));
    }

    for (nmpcSine_B.k_lc = nmpcSine_B.i_n; nmpcSine_B.k_lc < nmpcSine_B.ixlast_j;
         nmpcSine_B.k_lc++) {
      obj->Hx.data[nmpcSine_B.k_lc - 1] = x_data[nmpcSine_B.k_lc - 1] *
        obj->beta;
    }

    nmpcSine_B.b_j1 = static_cast<uint16_T>(obj->maxVar - 1);
    memcpy(&obj->grad.data[0], &obj->Hx.data[0], static_cast<uint32_T>
           (nmpcSine_B.b_j1) * sizeof(real_T));
    if (obj->hasLinear) {
      nmpcSine_B.iy_p = obj->grad.size[0];
      nmpcSine_B.i_n = obj->grad.size[0];
      if (nmpcSine_B.i_n - 1 >= 0) {
        memcpy(&nmpcSine_B.y_data_cv[0], &obj->grad.data[0],
               static_cast<uint32_T>(nmpcSine_B.i_n) * sizeof(real_T));
      }

      if (obj->nvar >= 1) {
        nmpcSine_B.ixlast_j = obj->nvar;
        nmpcSine_B.i_n = (obj->nvar / 2) << 1;
        nmpcSine_B.b_j1 = nmpcSine_B.i_n - 2;
        for (nmpcSine_B.k_lc = 0; nmpcSine_B.k_lc <= nmpcSine_B.b_j1;
             nmpcSine_B.k_lc += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.y_data_cv[nmpcSine_B.k_lc]);
          _mm_storeu_pd(&nmpcSine_B.y_data_cv[nmpcSine_B.k_lc], _mm_add_pd(tmp,
            _mm_loadu_pd(&f_data[nmpcSine_B.k_lc])));
        }

        for (nmpcSine_B.k_lc = nmpcSine_B.i_n; nmpcSine_B.k_lc <
             nmpcSine_B.ixlast_j; nmpcSine_B.k_lc++) {
          nmpcSine_B.y_data_cv[nmpcSine_B.k_lc] += f_data[nmpcSine_B.k_lc];
        }
      }

      if (nmpcSine_B.iy_p - 1 >= 0) {
        memcpy(&obj->grad.data[0], &nmpcSine_B.y_data_cv[0],
               static_cast<uint32_T>(nmpcSine_B.iy_p) * sizeof(real_T));
      }
    }

    if (nmpcSine_B.i_tmp >= 1) {
      nmpcSine_B.iy_p = obj->nvar;
      nmpcSine_B.i_n = (nmpcSine_B.i_tmp / 2) << 1;
      nmpcSine_B.b_j1 = nmpcSine_B.i_n - 2;
      for (nmpcSine_B.ixlast_j = 0; nmpcSine_B.ixlast_j <= nmpcSine_B.b_j1;
           nmpcSine_B.ixlast_j += 2) {
        nmpcSine_B.k_lc = nmpcSine_B.iy_p + nmpcSine_B.ixlast_j;
        tmp = _mm_loadu_pd(&obj->grad.data[nmpcSine_B.k_lc]);
        _mm_storeu_pd(&obj->grad.data[nmpcSine_B.k_lc], _mm_add_pd(tmp,
          _mm_set1_pd(obj->rho)));
      }

      for (nmpcSine_B.ixlast_j = nmpcSine_B.i_n; nmpcSine_B.ixlast_j <
           nmpcSine_B.i_tmp; nmpcSine_B.ixlast_j++) {
        nmpcSine_B.k_lc = nmpcSine_B.iy_p + nmpcSine_B.ixlast_j;
        obj->grad.data[nmpcSine_B.k_lc] += obj->rho;
      }
    }
    break;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_computeFval_ReuseHx(const
  s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *obj, real_T workspace_data[], const real_T
  f_data[], const real_T x_data[])
{
  real_T val;
  val = 0.0;
  switch (obj->objtype) {
   case 5:
    val = x_data[obj->nvar - 1] * obj->gammaScalar;
    break;

   case 3:
    {
      if (obj->hasLinear) {
        int32_T b;
        int32_T ixlast;
        int32_T maxRegVar;
        b = obj->nvar;
        maxRegVar = (obj->nvar / 2) << 1;
        ixlast = maxRegVar - 2;
        for (int32_T b_ixlast = 0; b_ixlast <= ixlast; b_ixlast += 2) {
          __m128d tmp_0;
          tmp_0 = _mm_loadu_pd(&obj->Hx.data[b_ixlast]);
          _mm_storeu_pd(&workspace_data[b_ixlast], _mm_add_pd(_mm_mul_pd
            (_mm_set1_pd(0.5), tmp_0), _mm_loadu_pd(&f_data[b_ixlast])));
        }

        for (int32_T b_ixlast = maxRegVar; b_ixlast < b; b_ixlast++) {
          workspace_data[b_ixlast] = 0.5 * obj->Hx.data[b_ixlast] +
            f_data[b_ixlast];
        }

        if (obj->nvar >= 1) {
          ixlast = obj->nvar;
          for (maxRegVar = 0; maxRegVar < ixlast; maxRegVar++) {
            val += workspace_data[maxRegVar] * x_data[maxRegVar];
          }
        }
      } else {
        if (obj->nvar >= 1) {
          int32_T ixlast;
          ixlast = obj->nvar;
          for (int32_T maxRegVar = 0; maxRegVar < ixlast; maxRegVar++) {
            val += x_data[maxRegVar] * obj->Hx.data[maxRegVar];
          }
        }

        val *= 0.5;
      }
    }
    break;

   case 4:
    {
      int32_T b;
      b = obj->maxVar;
      if (obj->hasLinear) {
        int32_T b_tmp;
        int32_T ixlast;
        int32_T maxRegVar;
        if (obj->nvar - 1 >= 0) {
          memcpy(&workspace_data[0], &f_data[0], static_cast<uint32_T>(obj->nvar)
                 * sizeof(real_T));
        }

        ixlast = obj->maxVar - obj->nvar;
        for (maxRegVar = 0; maxRegVar <= ixlast - 2; maxRegVar++) {
          workspace_data[obj->nvar + maxRegVar] = obj->rho;
        }

        b_tmp = static_cast<uint16_T>(obj->maxVar - 1);
        maxRegVar = (static_cast<uint16_T>(obj->maxVar - 1) / 2) << 1;
        ixlast = maxRegVar - 2;
        for (int32_T b_ixlast = 0; b_ixlast <= ixlast; b_ixlast += 2) {
          __m128d tmp;
          __m128d tmp_0;
          tmp_0 = _mm_loadu_pd(&obj->Hx.data[b_ixlast]);
          tmp = _mm_loadu_pd(&workspace_data[b_ixlast]);
          _mm_storeu_pd(&workspace_data[b_ixlast], _mm_add_pd(tmp, _mm_mul_pd
            (_mm_set1_pd(0.5), tmp_0)));
        }

        for (int32_T b_ixlast = maxRegVar; b_ixlast < b_tmp; b_ixlast++) {
          workspace_data[b_ixlast] += 0.5 * obj->Hx.data[b_ixlast];
        }

        if (obj->maxVar - 1 >= 1) {
          for (maxRegVar = 0; maxRegVar <= b - 2; maxRegVar++) {
            val += workspace_data[maxRegVar] * x_data[maxRegVar];
          }
        }
      } else {
        int32_T b_ixlast;
        if (obj->maxVar - 1 >= 1) {
          for (int32_T ixlast = 0; ixlast <= b - 2; ixlast++) {
            val += x_data[ixlast] * obj->Hx.data[ixlast];
          }
        }

        val *= 0.5;
        b_ixlast = obj->nvar + 1;
        for (int32_T ixlast = b_ixlast; ixlast < b; ixlast++) {
          val += x_data[ixlast - 1] * obj->rho;
        }
      }
    }
    break;
  }

  return val;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgeqrf(real_T A_data[], const int32_T A_size[2], int32_T
  m, int32_T n, real_T tau_data[], int32_T tau_size[1])
{
  int32_T loop_ub;
  int32_T minmn;
  if (m <= n) {
    minmn = m;
  } else {
    minmn = n;
  }

  if (A_size[0] <= A_size[1]) {
    loop_ub = A_size[0];
  } else {
    loop_ub = A_size[1];
  }

  tau_size[0] = loop_ub;
  if (loop_ub - 1 >= 0) {
    memset(&tau_data[0], 0, static_cast<uint32_T>(loop_ub) * sizeof(real_T));
  }

  if (minmn >= 1) {
    nmpcSine_qrf(A_data, A_size, m, n, minmn, tau_data);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factorQR(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj, const
  real_T A_data[], int32_T mrows, int32_T ncols, int32_T ldA)
{
  int32_T b;
  int32_T idx;
  int32_T ix0;
  int32_T iy0;
  int32_T k;
  static const int32_T offsets[4] = { 0, 1, 2, 3 };

  boolean_T guard1;
  idx = mrows * ncols;
  guard1 = false;
  if (idx > 0) {
    for (idx = 0; idx < ncols; idx++) {
      ix0 = ldA * idx;
      iy0 = obj->ldq * idx;
      b = static_cast<uint16_T>(mrows);
      for (k = 0; k < b; k++) {
        obj->QR.data[iy0 + k] = A_data[ix0 + k];
      }
    }

    guard1 = true;
  } else if (idx == 0) {
    obj->mrows = mrows;
    obj->ncols = ncols;
    obj->minRowCol = 0;
  } else {
    guard1 = true;
  }

  if (guard1) {
    obj->usedPivoting = false;
    obj->mrows = mrows;
    obj->ncols = ncols;
    k = (ncols / 4) << 2;
    ix0 = k - 4;
    for (idx = 0; idx <= ix0; idx += 4) {
      _mm_storeu_si128((__m128i *)&obj->jpvt.data[idx], _mm_add_epi32
                       (_mm_add_epi32(_mm_set1_epi32(idx), _mm_loadu_si128((
        const __m128i *)&offsets[0])), _mm_set1_epi32(1)));
    }

    for (idx = k; idx < ncols; idx++) {
      obj->jpvt.data[idx] = idx + 1;
    }

    if (mrows <= ncols) {
      obj->minRowCol = mrows;
    } else {
      obj->minRowCol = ncols;
    }

    nmpcSine_xgeqrf(obj->QR.data, obj->QR.size, mrows, ncols, obj->tau.data,
                    obj->tau.size);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xrotg(real_T *a, real_T *b, real_T *c, real_T *s)
{
  nmpcSine_B.roe = *b;
  nmpcSine_B.absa = fabs(*a);
  nmpcSine_B.absb = fabs(*b);
  if (nmpcSine_B.absa > nmpcSine_B.absb) {
    nmpcSine_B.roe = *a;
  }

  nmpcSine_B.scale_c = nmpcSine_B.absa + nmpcSine_B.absb;
  if (nmpcSine_B.scale_c == 0.0) {
    *s = 0.0;
    *c = 1.0;
    *a = 0.0;
    *b = 0.0;
  } else {
    nmpcSine_B.ads = nmpcSine_B.absa / nmpcSine_B.scale_c;
    nmpcSine_B.bds = nmpcSine_B.absb / nmpcSine_B.scale_c;
    nmpcSine_B.scale_c *= sqrt(nmpcSine_B.ads * nmpcSine_B.ads + nmpcSine_B.bds *
      nmpcSine_B.bds);
    if (nmpcSine_B.roe < 0.0) {
      nmpcSine_B.scale_c = -nmpcSine_B.scale_c;
    }

    *c = *a / nmpcSine_B.scale_c;
    *s = *b / nmpcSine_B.scale_c;
    if (nmpcSine_B.absa > nmpcSine_B.absb) {
      *b = *s;
    } else if (*c != 0.0) {
      *b = 1.0 / *c;
    } else {
      *b = 1.0;
    }

    *a = nmpcSine_B.scale_c;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_squareQ_appendCol(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj,
  const real_T vec_data[], int32_T iv0)
{
  real_T b_c;
  real_T s;
  real_T temp;
  int32_T Qk0;
  int32_T b_iy;
  int32_T e;
  int32_T idx;
  int32_T iy;
  int32_T iyend;
  int32_T temp_tmp;
  if (obj->mrows <= obj->ncols + 1) {
    obj->minRowCol = obj->mrows;
  } else {
    obj->minRowCol = obj->ncols + 1;
  }

  b_iy = obj->ldq * obj->ncols;
  idx = obj->ldq;
  if (obj->mrows != 0) {
    iyend = b_iy + obj->mrows;
    if (b_iy + 1 <= iyend) {
      memset(&obj->QR.data[b_iy], 0, static_cast<uint32_T>(iyend - b_iy) *
             sizeof(real_T));
    }

    iy = (obj->mrows - 1) * obj->ldq + 1;
    for (Qk0 = 1; idx < 0 ? Qk0 >= iy : Qk0 <= iy; Qk0 += idx) {
      b_c = 0.0;
      e = (Qk0 + obj->mrows) - 1;
      for (iyend = Qk0; iyend <= e; iyend++) {
        b_c += vec_data[((iv0 + iyend) - Qk0) - 1] * obj->Q.data[iyend - 1];
      }

      obj->QR.data[b_iy] += b_c;
      b_iy++;
    }
  }

  obj->ncols++;
  obj->jpvt.data[obj->ncols - 1] = obj->ncols;
  for (idx = obj->mrows - 2; idx + 2 > obj->ncols; idx--) {
    b_iy = (obj->ncols - 1) * obj->ldq + idx;
    temp = obj->QR.data[b_iy + 1];
    nmpcSine_xrotg(&obj->QR.data[b_iy], &temp, &b_c, &s);
    obj->QR.data[b_iy + 1] = temp;
    Qk0 = obj->ldq * idx;
    iyend = obj->mrows;
    if (obj->mrows >= 1) {
      iy = obj->ldq + Qk0;
      for (b_iy = 0; b_iy < iyend; b_iy++) {
        e = iy + b_iy;
        temp_tmp = Qk0 + b_iy;
        temp = obj->Q.data[temp_tmp] * b_c + obj->Q.data[e] * s;
        obj->Q.data[e] = obj->Q.data[e] * b_c - obj->Q.data[temp_tmp] * s;
        obj->Q.data[temp_tmp] = temp;
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_deleteColMoveEnd(s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *obj,
  int32_T idx)
{
  if (obj->usedPivoting) {
    nmpcSine_B.i_o4 = 1;
    while ((nmpcSine_B.i_o4 <= obj->ncols) && (obj->jpvt.data[nmpcSine_B.i_o4 -
            1] != idx)) {
      nmpcSine_B.i_o4++;
    }

    idx = nmpcSine_B.i_o4;
  }

  if (idx >= obj->ncols) {
    obj->ncols--;
  } else {
    obj->jpvt.data[idx - 1] = obj->jpvt.data[obj->ncols - 1];
    nmpcSine_B.QRk0 = obj->minRowCol;
    for (nmpcSine_B.i_o4 = 0; nmpcSine_B.i_o4 < nmpcSine_B.QRk0; nmpcSine_B.i_o4
         ++) {
      obj->QR.data[nmpcSine_B.i_o4 + obj->ldq * (idx - 1)] = obj->QR.data
        [(obj->ncols - 1) * obj->ldq + nmpcSine_B.i_o4];
    }

    obj->ncols--;
    if (obj->mrows <= obj->ncols) {
      obj->minRowCol = obj->mrows;
    } else {
      obj->minRowCol = obj->ncols;
    }

    if (idx < obj->mrows) {
      if (obj->mrows - 1 <= obj->ncols) {
        nmpcSine_B.i_o4 = obj->mrows - 1;
      } else {
        nmpcSine_B.i_o4 = obj->ncols;
      }

      nmpcSine_B.k_fr = nmpcSine_B.i_o4;
      nmpcSine_B.idxRotGCol = (idx - 1) * obj->ldq;
      while (nmpcSine_B.k_fr >= idx) {
        nmpcSine_B.QRk0 = nmpcSine_B.k_fr + nmpcSine_B.idxRotGCol;
        nmpcSine_B.b_temp = obj->QR.data[nmpcSine_B.QRk0];
        nmpcSine_xrotg(&obj->QR.data[nmpcSine_B.QRk0 - 1], &nmpcSine_B.b_temp,
                       &nmpcSine_B.c_c_c, &nmpcSine_B.b_s);
        obj->QR.data[nmpcSine_B.QRk0] = nmpcSine_B.b_temp;
        obj->QR.data[nmpcSine_B.k_fr + obj->ldq * (nmpcSine_B.k_fr - 1)] = 0.0;
        nmpcSine_B.QRk0 = obj->ldq * idx + nmpcSine_B.k_fr;
        nmpcSine_B.b_ix = obj->ncols - idx;
        if (nmpcSine_B.b_ix >= 1) {
          nmpcSine_B.ix_m = nmpcSine_B.QRk0 - 1;
          for (nmpcSine_B.b_n = 0; nmpcSine_B.b_n < nmpcSine_B.b_ix;
               nmpcSine_B.b_n++) {
            nmpcSine_B.b_temp = obj->QR.data[nmpcSine_B.ix_m] * nmpcSine_B.c_c_c
              + obj->QR.data[nmpcSine_B.QRk0] * nmpcSine_B.b_s;
            obj->QR.data[nmpcSine_B.QRk0] = obj->QR.data[nmpcSine_B.QRk0] *
              nmpcSine_B.c_c_c - obj->QR.data[nmpcSine_B.ix_m] * nmpcSine_B.b_s;
            obj->QR.data[nmpcSine_B.ix_m] = nmpcSine_B.b_temp;
            nmpcSine_B.QRk0 += obj->ldq;
            nmpcSine_B.ix_m += obj->ldq;
          }
        }

        nmpcSine_B.QRk0 = (nmpcSine_B.k_fr - 1) * obj->ldq;
        nmpcSine_B.b_ix = obj->mrows;
        if (obj->mrows >= 1) {
          nmpcSine_B.ix_m = obj->ldq + nmpcSine_B.QRk0;
          for (nmpcSine_B.b_n = 0; nmpcSine_B.b_n < nmpcSine_B.b_ix;
               nmpcSine_B.b_n++) {
            nmpcSine_B.d_temp_tmp = nmpcSine_B.ix_m + nmpcSine_B.b_n;
            nmpcSine_B.QRk0_tmp = nmpcSine_B.QRk0 + nmpcSine_B.b_n;
            nmpcSine_B.b_temp = obj->Q.data[nmpcSine_B.QRk0_tmp] *
              nmpcSine_B.c_c_c + obj->Q.data[nmpcSine_B.d_temp_tmp] *
              nmpcSine_B.b_s;
            obj->Q.data[nmpcSine_B.d_temp_tmp] = obj->
              Q.data[nmpcSine_B.d_temp_tmp] * nmpcSine_B.c_c_c - obj->
              Q.data[nmpcSine_B.QRk0_tmp] * nmpcSine_B.b_s;
            obj->Q.data[nmpcSine_B.QRk0_tmp] = nmpcSine_B.b_temp;
          }
        }

        nmpcSine_B.k_fr--;
      }

      for (nmpcSine_B.k_fr = idx + 1; nmpcSine_B.k_fr <= nmpcSine_B.i_o4;
           nmpcSine_B.k_fr++) {
        nmpcSine_B.QRk0_tmp = (nmpcSine_B.k_fr - 1) * obj->ldq;
        nmpcSine_B.QRk0 = nmpcSine_B.k_fr + nmpcSine_B.QRk0_tmp;
        nmpcSine_B.b_temp = obj->QR.data[nmpcSine_B.QRk0];
        nmpcSine_xrotg(&obj->QR.data[nmpcSine_B.QRk0 - 1], &nmpcSine_B.b_temp,
                       &nmpcSine_B.c_c_c, &nmpcSine_B.b_s);
        obj->QR.data[nmpcSine_B.QRk0] = nmpcSine_B.b_temp;
        nmpcSine_B.QRk0 = (obj->ldq + 1) * nmpcSine_B.k_fr;
        nmpcSine_B.b_n = obj->ncols - nmpcSine_B.k_fr;
        if (nmpcSine_B.b_n >= 1) {
          nmpcSine_B.b_ix = nmpcSine_B.QRk0 - 1;
          for (nmpcSine_B.idxRotGCol = 0; nmpcSine_B.idxRotGCol < nmpcSine_B.b_n;
               nmpcSine_B.idxRotGCol++) {
            nmpcSine_B.b_temp = obj->QR.data[nmpcSine_B.b_ix] * nmpcSine_B.c_c_c
              + obj->QR.data[nmpcSine_B.QRk0] * nmpcSine_B.b_s;
            obj->QR.data[nmpcSine_B.QRk0] = obj->QR.data[nmpcSine_B.QRk0] *
              nmpcSine_B.c_c_c - obj->QR.data[nmpcSine_B.b_ix] * nmpcSine_B.b_s;
            obj->QR.data[nmpcSine_B.b_ix] = nmpcSine_B.b_temp;
            nmpcSine_B.QRk0 += obj->ldq;
            nmpcSine_B.b_ix += obj->ldq;
          }
        }

        nmpcSine_B.b_n = obj->mrows;
        if (obj->mrows >= 1) {
          nmpcSine_B.b_ix = obj->ldq + nmpcSine_B.QRk0_tmp;
          for (nmpcSine_B.idxRotGCol = 0; nmpcSine_B.idxRotGCol < nmpcSine_B.b_n;
               nmpcSine_B.idxRotGCol++) {
            nmpcSine_B.ix_m = nmpcSine_B.b_ix + nmpcSine_B.idxRotGCol;
            nmpcSine_B.d_temp_tmp = nmpcSine_B.QRk0_tmp + nmpcSine_B.idxRotGCol;
            nmpcSine_B.b_temp = obj->Q.data[nmpcSine_B.d_temp_tmp] *
              nmpcSine_B.c_c_c + obj->Q.data[nmpcSine_B.ix_m] * nmpcSine_B.b_s;
            obj->Q.data[nmpcSine_B.ix_m] = obj->Q.data[nmpcSine_B.ix_m] *
              nmpcSine_B.c_c_c - obj->Q.data[nmpcSine_B.d_temp_tmp] *
              nmpcSine_B.b_s;
            obj->Q.data[nmpcSine_B.d_temp_tmp] = nmpcSine_B.b_temp;
          }
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
boolean_T nmpcSine::nmpcSine_strcmp(const char_T a[7])
{
  int32_T ret;
  static const char_T b[7] = { 'f', 'm', 'i', 'n', 'c', 'o', 'n' };

  ret = memcmp(&a[0], &b[0], 7);
  return ret == 0;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemm(int32_T m, int32_T n, int32_T k, const real_T A
  [15625], int32_T lda, const real_T B_data[], int32_T ib0, int32_T ldb, real_T
  C_data[], int32_T ldc)
{
  if ((m != 0) && (n != 0)) {
    nmpcSine_B.br_a = ib0;
    nmpcSine_B.lastColC_k = (n - 1) * ldc;
    for (nmpcSine_B.cr = 0; ldc < 0 ? nmpcSine_B.cr >= nmpcSine_B.lastColC_k :
         nmpcSine_B.cr <= nmpcSine_B.lastColC_k; nmpcSine_B.cr += ldc) {
      nmpcSine_B.b_f = nmpcSine_B.cr + m;
      for (nmpcSine_B.ic_c = nmpcSine_B.cr + 1; nmpcSine_B.ic_c <=
           nmpcSine_B.b_f; nmpcSine_B.ic_c++) {
        C_data[nmpcSine_B.ic_c - 1] = 0.0;
      }
    }

    for (nmpcSine_B.cr = 0; ldc < 0 ? nmpcSine_B.cr >= nmpcSine_B.lastColC_k :
         nmpcSine_B.cr <= nmpcSine_B.lastColC_k; nmpcSine_B.cr += ldc) {
      nmpcSine_B.ar = -1;
      nmpcSine_B.c_jk = nmpcSine_B.br_a + k;
      for (nmpcSine_B.ic_c = nmpcSine_B.br_a; nmpcSine_B.ic_c < nmpcSine_B.c_jk;
           nmpcSine_B.ic_c++) {
        nmpcSine_B.d_k = nmpcSine_B.cr + m;
        nmpcSine_B.scalarLB_hj = ((((nmpcSine_B.d_k - nmpcSine_B.cr) / 2) << 1)
          + nmpcSine_B.cr) + 1;
        nmpcSine_B.vectorUB_d1 = nmpcSine_B.scalarLB_hj - 2;
        for (nmpcSine_B.b_f = nmpcSine_B.cr + 1; nmpcSine_B.b_f <=
             nmpcSine_B.vectorUB_d1; nmpcSine_B.b_f += 2) {
          __m128d tmp;
          tmp = _mm_loadu_pd(&C_data[nmpcSine_B.b_f - 1]);
          _mm_storeu_pd(&C_data[nmpcSine_B.b_f - 1], _mm_add_pd(_mm_mul_pd
            (_mm_set1_pd(B_data[nmpcSine_B.ic_c - 1]), _mm_loadu_pd(&A
            [(nmpcSine_B.ar + nmpcSine_B.b_f) - nmpcSine_B.cr])), tmp));
        }

        for (nmpcSine_B.b_f = nmpcSine_B.scalarLB_hj; nmpcSine_B.b_f <=
             nmpcSine_B.d_k; nmpcSine_B.b_f++) {
          C_data[nmpcSine_B.b_f - 1] += A[(nmpcSine_B.ar + nmpcSine_B.b_f) -
            nmpcSine_B.cr] * B_data[nmpcSine_B.ic_c - 1];
        }

        nmpcSine_B.ar += lda;
      }

      nmpcSine_B.br_a += ldb;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemm_k(int32_T m, int32_T n, int32_T k, const real_T
  A_data[], int32_T ia0, int32_T lda, const real_T B_data[], int32_T ldb, real_T
  C_data[], int32_T ldc)
{
  if ((m != 0) && (n != 0)) {
    int32_T br;
    int32_T lastColC;
    lastColC = (n - 1) * ldc;
    for (int32_T cr = 0; ldc < 0 ? cr >= lastColC : cr <= lastColC; cr += ldc) {
      br = cr + m;
      for (int32_T ic = cr + 1; ic <= br; ic++) {
        C_data[ic - 1] = 0.0;
      }
    }

    br = -1;
    for (int32_T cr = 0; ldc < 0 ? cr >= lastColC : cr <= lastColC; cr += ldc) {
      int32_T ar;
      int32_T c;
      ar = ia0;
      c = cr + m;
      for (int32_T ic = cr + 1; ic <= c; ic++) {
        real_T temp;
        temp = 0.0;
        for (int32_T b_w = 0; b_w < k; b_w++) {
          temp += A_data[(b_w + ar) - 1] * B_data[(b_w + br) + 1];
        }

        C_data[ic - 1] += temp;
        ar += lda;
      }

      br += ldb;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_fullColLDL2_(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj,
  int32_T LD_offset, int32_T NColsRemain)
{
  nmpcSine_B.LDimSizeP1_d = obj->ldm;
  for (nmpcSine_B.k_hj = 0; nmpcSine_B.k_hj < NColsRemain; nmpcSine_B.k_hj++) {
    __m128d tmp;
    nmpcSine_B.LD_diagOffset_b = (nmpcSine_B.LDimSizeP1_d + 1) * nmpcSine_B.k_hj
      + LD_offset;
    nmpcSine_B.alpha1 = -1.0 / obj->FMat.data[nmpcSine_B.LD_diagOffset_b - 1];
    nmpcSine_B.subMatrixDim = (NColsRemain - nmpcSine_B.k_hj) - 2;
    for (nmpcSine_B.b_k_p = 0; nmpcSine_B.b_k_p <= nmpcSine_B.subMatrixDim;
         nmpcSine_B.b_k_p++) {
      obj->workspace_ = obj->FMat.data[nmpcSine_B.LD_diagOffset_b +
        nmpcSine_B.b_k_p];
    }

    nmpcSine_B.y_f = obj->workspace_;
    if (!(nmpcSine_B.alpha1 == 0.0)) {
      nmpcSine_B.jA_n = nmpcSine_B.LD_diagOffset_b + nmpcSine_B.LDimSizeP1_d;
      for (nmpcSine_B.j_jz = 0; nmpcSine_B.j_jz <= nmpcSine_B.subMatrixDim;
           nmpcSine_B.j_jz++) {
        if (nmpcSine_B.y_f != 0.0) {
          nmpcSine_B.temp_d = nmpcSine_B.y_f * nmpcSine_B.alpha1;
          nmpcSine_B.b_ot = (nmpcSine_B.subMatrixDim + nmpcSine_B.jA_n) + 1;
          nmpcSine_B.b_k_p = ((((nmpcSine_B.b_ot - nmpcSine_B.jA_n) / 2) << 1) +
                              nmpcSine_B.jA_n) + 1;
          nmpcSine_B.alpha1_tmp = nmpcSine_B.b_k_p - 2;
          for (nmpcSine_B.ijA = nmpcSine_B.jA_n + 1; nmpcSine_B.ijA <=
               nmpcSine_B.alpha1_tmp; nmpcSine_B.ijA += 2) {
            tmp = _mm_loadu_pd(&obj->FMat.data[nmpcSine_B.ijA - 1]);
            _mm_storeu_pd(&obj->FMat.data[nmpcSine_B.ijA - 1], _mm_add_pd(tmp,
              _mm_set1_pd(obj->workspace_ * nmpcSine_B.temp_d)));
          }

          for (nmpcSine_B.ijA = nmpcSine_B.b_k_p; nmpcSine_B.ijA <=
               nmpcSine_B.b_ot; nmpcSine_B.ijA++) {
            obj->FMat.data[nmpcSine_B.ijA - 1] += obj->workspace_ *
              nmpcSine_B.temp_d;
          }
        }

        nmpcSine_B.jA_n += obj->ldm;
      }
    }

    nmpcSine_B.alpha1 = 1.0 / obj->FMat.data[nmpcSine_B.LD_diagOffset_b - 1];
    nmpcSine_B.j_jz = (nmpcSine_B.LD_diagOffset_b + nmpcSine_B.subMatrixDim) + 1;
    nmpcSine_B.b_k_p = ((((nmpcSine_B.j_jz - nmpcSine_B.LD_diagOffset_b) / 2) <<
                         1) + nmpcSine_B.LD_diagOffset_b) + 1;
    nmpcSine_B.alpha1_tmp = nmpcSine_B.b_k_p - 2;
    for (nmpcSine_B.subMatrixDim = nmpcSine_B.LD_diagOffset_b + 1;
         nmpcSine_B.subMatrixDim <= nmpcSine_B.alpha1_tmp;
         nmpcSine_B.subMatrixDim += 2) {
      tmp = _mm_loadu_pd(&obj->FMat.data[nmpcSine_B.subMatrixDim - 1]);
      _mm_storeu_pd(&obj->FMat.data[nmpcSine_B.subMatrixDim - 1], _mm_mul_pd(tmp,
        _mm_set1_pd(nmpcSine_B.alpha1)));
    }

    for (nmpcSine_B.subMatrixDim = nmpcSine_B.b_k_p; nmpcSine_B.subMatrixDim <=
         nmpcSine_B.j_jz; nmpcSine_B.subMatrixDim++) {
      obj->FMat.data[nmpcSine_B.subMatrixDim - 1] *= nmpcSine_B.alpha1;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_partialColLDL3_(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj,
  int32_T LD_offset, int32_T NColsRemain)
{
  __m128d tmp;
  nmpcSine_B.LDimSizeP1_k = obj->ldm + 1;
  for (nmpcSine_B.k_jk = 0; nmpcSine_B.k_jk < 48; nmpcSine_B.k_jk++) {
    nmpcSine_B.subRows = (NColsRemain - nmpcSine_B.k_jk) - 1;
    nmpcSine_B.LD_diagOffset = (nmpcSine_B.LDimSizeP1_k * nmpcSine_B.k_jk +
      LD_offset) - 1;
    for (nmpcSine_B.subBlockSize = 0; nmpcSine_B.subBlockSize <=
         nmpcSine_B.subRows; nmpcSine_B.subBlockSize++) {
      obj->workspace_ = obj->FMat.data[nmpcSine_B.LD_diagOffset +
        nmpcSine_B.subBlockSize];
    }

    for (nmpcSine_B.subBlockSize = 0; nmpcSine_B.subBlockSize < NColsRemain;
         nmpcSine_B.subBlockSize++) {
      obj->workspace2_ = obj->workspace_;
    }

    nmpcSine_B.br_k = obj->ldm;
    nmpcSine_B.y_p = obj->workspace2_;
    if ((NColsRemain != 0) && (nmpcSine_B.k_jk != 0)) {
      nmpcSine_B.ix_h = LD_offset + nmpcSine_B.k_jk;
      nmpcSine_B.c_f = (nmpcSine_B.k_jk - 1) * obj->ldm + 1;
      for (nmpcSine_B.subBlockSize = 1; nmpcSine_B.br_k < 0 ?
           nmpcSine_B.subBlockSize >= nmpcSine_B.c_f : nmpcSine_B.subBlockSize <=
           nmpcSine_B.c_f; nmpcSine_B.subBlockSize += nmpcSine_B.br_k) {
        nmpcSine_B.d_d = nmpcSine_B.subBlockSize + NColsRemain;
        for (nmpcSine_B.ia = nmpcSine_B.subBlockSize; nmpcSine_B.ia <
             nmpcSine_B.d_d; nmpcSine_B.ia++) {
          nmpcSine_B.y_p += -obj->FMat.data[nmpcSine_B.ix_h - 1] *
            obj->workspace_;
        }

        nmpcSine_B.ix_h += obj->ldm;
      }
    }

    obj->workspace2_ = nmpcSine_B.y_p;
    for (nmpcSine_B.subBlockSize = 0; nmpcSine_B.subBlockSize < NColsRemain;
         nmpcSine_B.subBlockSize++) {
      obj->workspace_ = nmpcSine_B.y_p;
    }

    for (nmpcSine_B.subBlockSize = 0; nmpcSine_B.subBlockSize <=
         nmpcSine_B.subRows; nmpcSine_B.subBlockSize++) {
      obj->FMat.data[nmpcSine_B.LD_diagOffset + nmpcSine_B.subBlockSize] =
        obj->workspace_;
    }

    nmpcSine_B.subBlockSize = (nmpcSine_B.subRows / 2) << 1;
    nmpcSine_B.ia = nmpcSine_B.subBlockSize - 2;
    for (nmpcSine_B.ix_h = 0; nmpcSine_B.ix_h <= nmpcSine_B.ia; nmpcSine_B.ix_h +=
         2) {
      nmpcSine_B.c_f = (nmpcSine_B.ix_h + nmpcSine_B.LD_diagOffset) + 1;
      tmp = _mm_loadu_pd(&obj->FMat.data[nmpcSine_B.c_f]);
      _mm_storeu_pd(&obj->FMat.data[nmpcSine_B.c_f], _mm_div_pd(tmp, _mm_set1_pd
        (obj->FMat.data[nmpcSine_B.LD_diagOffset])));
    }

    for (nmpcSine_B.ix_h = nmpcSine_B.subBlockSize; nmpcSine_B.ix_h <
         nmpcSine_B.subRows; nmpcSine_B.ix_h++) {
      nmpcSine_B.c_f = (nmpcSine_B.ix_h + nmpcSine_B.LD_diagOffset) + 1;
      obj->FMat.data[nmpcSine_B.c_f] /= obj->FMat.data[nmpcSine_B.LD_diagOffset];
    }
  }

  for (nmpcSine_B.k_jk = 48; nmpcSine_B.k_jk <= NColsRemain - 1; nmpcSine_B.k_jk
       += 48) {
    nmpcSine_B.br_k = NColsRemain - nmpcSine_B.k_jk;
    if (nmpcSine_B.br_k >= 48) {
      nmpcSine_B.subBlockSize = 48;
    } else {
      nmpcSine_B.subBlockSize = nmpcSine_B.br_k;
    }

    nmpcSine_B.subRows = nmpcSine_B.k_jk + nmpcSine_B.subBlockSize;
    for (nmpcSine_B.ia = nmpcSine_B.k_jk; nmpcSine_B.ia < nmpcSine_B.subRows;
         nmpcSine_B.ia++) {
      nmpcSine_B.LD_diagOffset = nmpcSine_B.subRows - nmpcSine_B.ia;
      for (nmpcSine_B.ix_h = 0; nmpcSine_B.ix_h < 48; nmpcSine_B.ix_h++) {
        obj->workspace2_ = obj->FMat.data[((LD_offset + nmpcSine_B.ia) +
          nmpcSine_B.ix_h * obj->ldm) - 1];
      }

      nmpcSine_B.d_d = obj->ldm;
      if (nmpcSine_B.LD_diagOffset != 0) {
        nmpcSine_B.e_li = (obj->ldm * 47 + nmpcSine_B.ia) + 1;
        for (nmpcSine_B.ix_h = nmpcSine_B.ia + 1; nmpcSine_B.d_d < 0 ?
             nmpcSine_B.ix_h >= nmpcSine_B.e_li : nmpcSine_B.ix_h <=
             nmpcSine_B.e_li; nmpcSine_B.ix_h += nmpcSine_B.d_d) {
          nmpcSine_B.lastColC = nmpcSine_B.ix_h + nmpcSine_B.LD_diagOffset;
          for (nmpcSine_B.c_f = nmpcSine_B.ix_h; nmpcSine_B.c_f <
               nmpcSine_B.lastColC; nmpcSine_B.c_f++) {
            // Check node always fails. would cause program termination and was eliminated 
          }
        }
      }
    }

    if (nmpcSine_B.subRows < NColsRemain) {
      nmpcSine_B.subRows = nmpcSine_B.br_k - nmpcSine_B.subBlockSize;
      nmpcSine_B.LD_diagOffset = ((LD_offset + nmpcSine_B.subBlockSize) +
        nmpcSine_B.LDimSizeP1_k * nmpcSine_B.k_jk) - 1;
      for (nmpcSine_B.ia = 0; nmpcSine_B.ia < 48; nmpcSine_B.ia++) {
        nmpcSine_B.ix_h = (LD_offset + nmpcSine_B.k_jk) + nmpcSine_B.ia *
          obj->ldm;
        for (nmpcSine_B.br_k = 0; nmpcSine_B.br_k < nmpcSine_B.subBlockSize;
             nmpcSine_B.br_k++) {
          obj->workspace2_ = obj->FMat.data[(nmpcSine_B.ix_h + nmpcSine_B.br_k)
            - 1];
        }
      }

      nmpcSine_B.ix_h = obj->ldm;
      if ((nmpcSine_B.subRows != 0) && (nmpcSine_B.subBlockSize != 0)) {
        nmpcSine_B.lastColC = (nmpcSine_B.subBlockSize - 1) * obj->ldm +
          nmpcSine_B.LD_diagOffset;
        nmpcSine_B.br_k = 0;
        for (nmpcSine_B.c_f = nmpcSine_B.LD_diagOffset; nmpcSine_B.ix_h < 0 ?
             nmpcSine_B.c_f >= nmpcSine_B.lastColC : nmpcSine_B.c_f <=
             nmpcSine_B.lastColC; nmpcSine_B.c_f += nmpcSine_B.ix_h) {
          nmpcSine_B.br_k++;
          nmpcSine_B.g_i = nmpcSine_B.ix_h * 47 + nmpcSine_B.br_k;
          for (nmpcSine_B.d_d = nmpcSine_B.br_k; nmpcSine_B.ix_h < 0 ?
               nmpcSine_B.d_d >= nmpcSine_B.g_i : nmpcSine_B.d_d <=
               nmpcSine_B.g_i; nmpcSine_B.d_d += nmpcSine_B.ix_h) {
            nmpcSine_B.h_h = nmpcSine_B.c_f + nmpcSine_B.subRows;
            nmpcSine_B.subBlockSize = ((((nmpcSine_B.h_h - nmpcSine_B.c_f) / 2) <<
              1) + nmpcSine_B.c_f) + 1;
            nmpcSine_B.ia = nmpcSine_B.subBlockSize - 2;
            for (nmpcSine_B.e_li = nmpcSine_B.c_f + 1; nmpcSine_B.e_li <=
                 nmpcSine_B.ia; nmpcSine_B.e_li += 2) {
              tmp = _mm_loadu_pd(&obj->FMat.data[nmpcSine_B.e_li - 1]);
              _mm_storeu_pd(&obj->FMat.data[nmpcSine_B.e_li - 1], _mm_add_pd(tmp,
                _mm_set1_pd(-obj->workspace2_ * obj->workspace_)));
            }

            for (nmpcSine_B.e_li = nmpcSine_B.subBlockSize; nmpcSine_B.e_li <=
                 nmpcSine_B.h_h; nmpcSine_B.e_li++) {
              obj->FMat.data[nmpcSine_B.e_li - 1] += -obj->workspace2_ *
                obj->workspace_;
            }
          }
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
int32_T nmpcSine::nmpcSine_xpotrf(int32_T n, real_T A_data[], int32_T lda)
{
  int32_T info;
  boolean_T exitg1;
  info = 0;
  nmpcSine_B.b_j_n = 0;
  exitg1 = false;
  while ((!exitg1) && (nmpcSine_B.b_j_n <= n - 1)) {
    nmpcSine_B.idxA1j = nmpcSine_B.b_j_n * lda;
    nmpcSine_B.idxAjj = nmpcSine_B.idxA1j + nmpcSine_B.b_j_n;
    nmpcSine_B.ssq = 0.0;
    if (nmpcSine_B.b_j_n >= 1) {
      for (nmpcSine_B.b_k_ot = 0; nmpcSine_B.b_k_ot < nmpcSine_B.b_j_n;
           nmpcSine_B.b_k_ot++) {
        nmpcSine_B.c_j = A_data[nmpcSine_B.idxA1j + nmpcSine_B.b_k_ot];
        nmpcSine_B.ssq += nmpcSine_B.c_j * nmpcSine_B.c_j;
      }
    }

    nmpcSine_B.ssq = A_data[nmpcSine_B.idxAjj] - nmpcSine_B.ssq;
    if (nmpcSine_B.ssq > 0.0) {
      nmpcSine_B.ssq = sqrt(nmpcSine_B.ssq);
      A_data[nmpcSine_B.idxAjj] = nmpcSine_B.ssq;
      if (nmpcSine_B.b_j_n + 1 < n) {
        nmpcSine_B.nmj = (n - nmpcSine_B.b_j_n) - 2;
        nmpcSine_B.ia0 = (nmpcSine_B.idxA1j + lda) + 1;
        nmpcSine_B.idxAjj += lda;
        if ((nmpcSine_B.b_j_n != 0) && (nmpcSine_B.nmj + 1 != 0)) {
          nmpcSine_B.iy_i = nmpcSine_B.idxAjj;
          nmpcSine_B.b_p = lda * nmpcSine_B.nmj + nmpcSine_B.ia0;
          for (nmpcSine_B.b_k_ot = nmpcSine_B.ia0; lda < 0 ? nmpcSine_B.b_k_ot >=
               nmpcSine_B.b_p : nmpcSine_B.b_k_ot <= nmpcSine_B.b_p;
               nmpcSine_B.b_k_ot += lda) {
            nmpcSine_B.c_j = 0.0;
            nmpcSine_B.d_o = (nmpcSine_B.b_k_ot + nmpcSine_B.b_j_n) - 1;
            for (nmpcSine_B.ia_m = nmpcSine_B.b_k_ot; nmpcSine_B.ia_m <=
                 nmpcSine_B.d_o; nmpcSine_B.ia_m++) {
              nmpcSine_B.c_j += A_data[(nmpcSine_B.idxA1j + nmpcSine_B.ia_m) -
                nmpcSine_B.b_k_ot] * A_data[nmpcSine_B.ia_m - 1];
            }

            A_data[nmpcSine_B.iy_i] -= nmpcSine_B.c_j;
            nmpcSine_B.iy_i += lda;
          }
        }

        nmpcSine_B.ssq = 1.0 / nmpcSine_B.ssq;
        nmpcSine_B.nmj = (lda * nmpcSine_B.nmj + nmpcSine_B.idxAjj) + 1;
        for (nmpcSine_B.idxA1j = nmpcSine_B.idxAjj + 1; lda < 0 ?
             nmpcSine_B.idxA1j >= nmpcSine_B.nmj : nmpcSine_B.idxA1j <=
             nmpcSine_B.nmj; nmpcSine_B.idxA1j += lda) {
          A_data[nmpcSine_B.idxA1j - 1] *= nmpcSine_B.ssq;
        }
      }

      nmpcSine_B.b_j_n++;
    } else {
      A_data[nmpcSine_B.idxAjj] = nmpcSine_B.ssq;
      info = nmpcSine_B.b_j_n + 1;
      exitg1 = true;
    }
  }

  return info;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemv_dffoy(int32_T m, int32_T n, const real_T A_data[],
  int32_T ia0, int32_T lda, const real_T x_data[], real_T y_data[])
{
  if (m != 0) {
    int32_T b;
    int32_T ix;
    memset(&y_data[0], 0, static_cast<uint32_T>(m) * sizeof(real_T));
    ix = 0;
    b = (n - 1) * lda + ia0;
    for (int32_T b_iy = ia0; lda < 0 ? b_iy >= b : b_iy <= b; b_iy += lda) {
      int32_T c;
      c = (b_iy + m) - 1;
      for (int32_T ia = b_iy; ia <= c; ia++) {
        int32_T tmp;
        tmp = ia - b_iy;
        y_data[tmp] += A_data[ia - 1] * x_data[ix];
      }

      ix++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factor_g(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, const
  real_T A[15625], int32_T ndims, int32_T ldA)
{
  int32_T exitg2;
  boolean_T exitg1;
  nmpcSine_B.LDimSizeP1_h = obj->ldm + 1;
  obj->ndims = ndims;
  for (nmpcSine_B.k_o = 0; nmpcSine_B.k_o < ndims; nmpcSine_B.k_o++) {
    nmpcSine_B.order_n = ldA * nmpcSine_B.k_o;
    nmpcSine_B.iy0_m = obj->ldm * nmpcSine_B.k_o;
    for (nmpcSine_B.A_maxDiag_idx_b = 0; nmpcSine_B.A_maxDiag_idx_b < ndims;
         nmpcSine_B.A_maxDiag_idx_b++) {
      obj->FMat.data[nmpcSine_B.iy0_m + nmpcSine_B.A_maxDiag_idx_b] =
        A[nmpcSine_B.A_maxDiag_idx_b + nmpcSine_B.order_n];
    }
  }

  if (ndims < 1) {
    nmpcSine_B.A_maxDiag_idx_b = -1;
  } else {
    nmpcSine_B.A_maxDiag_idx_b = 0;
    if (ndims > 1) {
      nmpcSine_B.smax_c = fabs(obj->FMat.data[0]);
      for (nmpcSine_B.k_o = 2; nmpcSine_B.k_o <= ndims; nmpcSine_B.k_o++) {
        nmpcSine_B.s_f = fabs(obj->FMat.data[(nmpcSine_B.k_o - 1) *
                              nmpcSine_B.LDimSizeP1_h]);
        if (nmpcSine_B.s_f > nmpcSine_B.smax_c) {
          nmpcSine_B.A_maxDiag_idx_b = nmpcSine_B.k_o - 1;
          nmpcSine_B.smax_c = nmpcSine_B.s_f;
        }
      }
    }
  }

  nmpcSine_B.smax_c = fabs(obj->FMat.data[obj->ldm * nmpcSine_B.A_maxDiag_idx_b
    + nmpcSine_B.A_maxDiag_idx_b]) * 2.2204460492503131E-16;
  if (nmpcSine_B.smax_c >= 0.0) {
    obj->regTol_ = nmpcSine_B.smax_c;
  } else {
    obj->regTol_ = 0.0;
  }

  if (ndims > 128) {
    nmpcSine_B.k_o = 0;
    exitg1 = false;
    while ((!exitg1) && (nmpcSine_B.k_o < ndims)) {
      nmpcSine_B.A_maxDiag_idx_b = nmpcSine_B.LDimSizeP1_h * nmpcSine_B.k_o + 1;
      nmpcSine_B.order_n = ndims - nmpcSine_B.k_o;
      if (nmpcSine_B.k_o + 48 <= ndims) {
        nmpcSine_partialColLDL3_(obj, nmpcSine_B.A_maxDiag_idx_b,
          nmpcSine_B.order_n);
        nmpcSine_B.k_o += 48;
      } else {
        nmpcSine_fullColLDL2_(obj, nmpcSine_B.A_maxDiag_idx_b,
                              nmpcSine_B.order_n);
        exitg1 = true;
      }
    }
  } else {
    nmpcSine_fullColLDL2_(obj, 1, ndims);
  }

  if (obj->ConvexCheck) {
    nmpcSine_B.LDimSizeP1_h = 0;
    do {
      exitg2 = 0;
      if (nmpcSine_B.LDimSizeP1_h <= ndims - 1) {
        if (obj->FMat.data[obj->ldm * nmpcSine_B.LDimSizeP1_h +
            nmpcSine_B.LDimSizeP1_h] <= 0.0) {
          obj->info = -nmpcSine_B.LDimSizeP1_h - 1;
          exitg2 = 1;
        } else {
          nmpcSine_B.LDimSizeP1_h++;
        }
      } else {
        obj->ConvexCheck = false;
        exitg2 = 1;
      }
    } while (exitg2 == 0);
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_factor(s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, const
  real_T A[15625], int32_T ndims, int32_T ldA)
{
  obj->ndims = ndims;
  for (nmpcSine_B.idx_b = 0; nmpcSine_B.idx_b < ndims; nmpcSine_B.idx_b++) {
    nmpcSine_B.ix0_j = ldA * nmpcSine_B.idx_b;
    nmpcSine_B.iy0_e = obj->ldm * nmpcSine_B.idx_b;
    for (nmpcSine_B.b_k_i = 0; nmpcSine_B.b_k_i < ndims; nmpcSine_B.b_k_i++) {
      obj->FMat.data[nmpcSine_B.iy0_e + nmpcSine_B.b_k_i] = A[nmpcSine_B.b_k_i +
        nmpcSine_B.ix0_j];
    }
  }

  obj->info = nmpcSine_xpotrf(ndims, obj->FMat.data, obj->ldm);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_solve_i(const s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj,
  real_T rhs_data[])
{
  int32_T jjA;
  int32_T n_tmp;
  n_tmp = obj->ndims;
  if (obj->ndims != 0) {
    for (int32_T b_j = 0; b_j < n_tmp; b_j++) {
      int32_T c;
      jjA = b_j * obj->ldm + b_j;
      c = (n_tmp - b_j) - 2;
      for (int32_T b_i = 0; b_i <= c; b_i++) {
        int32_T ix;
        ix = (b_i + b_j) + 1;
        rhs_data[ix] -= obj->FMat.data[(b_i + jjA) + 1] * rhs_data[b_j];
      }
    }
  }

  for (int32_T b_j = 0; b_j < n_tmp; b_j++) {
    rhs_data[b_j] /= obj->FMat.data[obj->ldm * b_j + b_j];
  }

  if (obj->ndims != 0) {
    for (int32_T b_j = n_tmp; b_j >= 1; b_j--) {
      real_T temp;
      jjA = (b_j - 1) * obj->ldm;
      temp = rhs_data[b_j - 1];
      for (int32_T b_i = n_tmp; b_i >= b_j + 1; b_i--) {
        temp -= obj->FMat.data[(jjA + b_i) - 1] * rhs_data[b_i - 1];
      }

      rhs_data[b_j - 1] = temp;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_solve(const s_mDApTYzDBpxuvxemclsuEF_nmpc_T *obj, real_T
  rhs_data[])
{
  int32_T n_tmp;
  n_tmp = obj->ndims;
  if (obj->ndims != 0) {
    int32_T jA;
    for (int32_T j = 0; j < n_tmp; j++) {
      real_T temp;
      jA = j * obj->ldm;
      temp = rhs_data[j];
      for (int32_T i = 0; i < j; i++) {
        temp -= obj->FMat.data[jA + i] * rhs_data[i];
      }

      rhs_data[j] = temp / obj->FMat.data[jA + j];
    }

    for (int32_T j = n_tmp; j >= 1; j--) {
      jA = ((j - 1) * obj->ldm + j) - 2;
      rhs_data[j - 1] /= obj->FMat.data[jA + 1];
      for (int32_T i = 0; i <= j - 2; i++) {
        int32_T ix;
        ix = (j - i) - 2;
        rhs_data[ix] -= obj->FMat.data[jA - i] * rhs_data[j - 1];
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_compute_deltax(const real_T H[15625],
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, const s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager,
  s_mDApTYzDBpxuvxemclsuEF_nmpc_T *cholmanager, const
  s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, boolean_T alwaysPositiveDef)
{
  __m128d tmp;
  int32_T exitg2;
  boolean_T exitg1;
  nmpcSine_B.nVar_he = qrmanager->mrows - 1;
  nmpcSine_B.mNull_tmp = qrmanager->mrows - qrmanager->ncols;
  if (nmpcSine_B.mNull_tmp <= 0) {
    if (nmpcSine_B.nVar_he >= 0) {
      memset(&solution->searchDir.data[0], 0, static_cast<uint32_T>
             (nmpcSine_B.nVar_he + 1) * sizeof(real_T));
    }
  } else {
    nmpcSine_B.A_maxDiag_idx = (qrmanager->mrows / 2) << 1;
    nmpcSine_B.order = nmpcSine_B.A_maxDiag_idx - 2;
    for (nmpcSine_B.b_idx = 0; nmpcSine_B.b_idx <= nmpcSine_B.order;
         nmpcSine_B.b_idx += 2) {
      tmp = _mm_loadu_pd(&objective->grad.data[nmpcSine_B.b_idx]);
      _mm_storeu_pd(&solution->searchDir.data[nmpcSine_B.b_idx], _mm_mul_pd(tmp,
        _mm_set1_pd(-1.0)));
    }

    for (nmpcSine_B.b_idx = nmpcSine_B.A_maxDiag_idx; nmpcSine_B.b_idx <=
         nmpcSine_B.nVar_he; nmpcSine_B.b_idx++) {
      solution->searchDir.data[nmpcSine_B.b_idx] = -objective->
        grad.data[nmpcSine_B.b_idx];
    }

    if (qrmanager->ncols <= 0) {
      switch (objective->objtype) {
       case 5:
        break;

       case 3:
        if (alwaysPositiveDef) {
          nmpcSine_factor(cholmanager, H, qrmanager->mrows, qrmanager->mrows);
        } else {
          nmpcSine_factor_g(cholmanager, H, qrmanager->mrows, qrmanager->mrows);
        }

        if (cholmanager->info != 0) {
          solution->state = -6;
        } else if (alwaysPositiveDef) {
          nmpcSine_solve(cholmanager, solution->searchDir.data);
        } else {
          nmpcSine_solve_i(cholmanager, solution->searchDir.data);
        }
        break;

       case 4:
        if (alwaysPositiveDef) {
          nmpcSine_factor(cholmanager, H, objective->nvar, objective->nvar);
          if (cholmanager->info != 0) {
            solution->state = -6;
          } else {
            nmpcSine_solve(cholmanager, solution->searchDir.data);
            nmpcSine_B.smax_d = 1.0 / objective->beta;
            nmpcSine_B.b_idx = objective->nvar + 1;
            nmpcSine_B.nVar_he = qrmanager->mrows;
            nmpcSine_B.A_maxDiag_idx = ((((qrmanager->mrows - objective->nvar) /
              2) << 1) + objective->nvar) + 1;
            nmpcSine_B.order = nmpcSine_B.A_maxDiag_idx - 2;
            for (nmpcSine_B.mNull_tmp = nmpcSine_B.b_idx; nmpcSine_B.mNull_tmp <=
                 nmpcSine_B.order; nmpcSine_B.mNull_tmp += 2) {
              tmp = _mm_loadu_pd(&solution->searchDir.data[nmpcSine_B.mNull_tmp
                                 - 1]);
              _mm_storeu_pd(&solution->searchDir.data[nmpcSine_B.mNull_tmp - 1],
                            _mm_mul_pd(tmp, _mm_set1_pd(nmpcSine_B.smax_d)));
            }

            for (nmpcSine_B.mNull_tmp = nmpcSine_B.A_maxDiag_idx;
                 nmpcSine_B.mNull_tmp <= nmpcSine_B.nVar_he;
                 nmpcSine_B.mNull_tmp++) {
              solution->searchDir.data[nmpcSine_B.mNull_tmp - 1] *=
                nmpcSine_B.smax_d;
            }
          }
        }
        break;
      }
    } else {
      nmpcSine_B.b_idx = qrmanager->ldq * qrmanager->ncols + 1;
      if (objective->objtype == 5) {
        for (nmpcSine_B.k_k = 0; nmpcSine_B.k_k < nmpcSine_B.mNull_tmp;
             nmpcSine_B.k_k++) {
          memspace->workspace_float.data[nmpcSine_B.k_k] = -qrmanager->Q.data
            [(qrmanager->ncols + nmpcSine_B.k_k) * qrmanager->ldq +
            nmpcSine_B.nVar_he];
        }

        nmpcSine_xgemv_dffoy(qrmanager->mrows, nmpcSine_B.mNull_tmp,
                             qrmanager->Q.data, nmpcSine_B.b_idx, qrmanager->ldq,
                             memspace->workspace_float.data,
                             solution->searchDir.data);
      } else {
        if (objective->objtype == 3) {
          nmpcSine_xgemm(qrmanager->mrows, nmpcSine_B.mNull_tmp,
                         qrmanager->mrows, H, qrmanager->mrows,
                         qrmanager->Q.data, nmpcSine_B.b_idx, qrmanager->ldq,
                         memspace->workspace_float.data,
                         memspace->workspace_float.size[0]);
          nmpcSine_xgemm_k(nmpcSine_B.mNull_tmp, nmpcSine_B.mNull_tmp,
                           qrmanager->mrows, qrmanager->Q.data, nmpcSine_B.b_idx,
                           qrmanager->ldq, memspace->workspace_float.data,
                           memspace->workspace_float.size[0],
                           cholmanager->FMat.data, cholmanager->ldm);
        } else if (alwaysPositiveDef) {
          nmpcSine_B.nVars = qrmanager->mrows;
          nmpcSine_xgemm(objective->nvar, nmpcSine_B.mNull_tmp, objective->nvar,
                         H, objective->nvar, qrmanager->Q.data, nmpcSine_B.b_idx,
                         qrmanager->ldq, memspace->workspace_float.data,
                         memspace->workspace_float.size[0]);
          for (nmpcSine_B.k_k = 0; nmpcSine_B.k_k < nmpcSine_B.mNull_tmp;
               nmpcSine_B.k_k++) {
            nmpcSine_B.c_ix = objective->nvar + 1;
            nmpcSine_B.A_maxDiag_idx = ((((nmpcSine_B.nVars - objective->nvar) /
              2) << 1) + objective->nvar) + 1;
            nmpcSine_B.order = nmpcSine_B.A_maxDiag_idx - 2;
            for (nmpcSine_B.LDimSizeP1 = nmpcSine_B.c_ix; nmpcSine_B.LDimSizeP1 <=
                 nmpcSine_B.order; nmpcSine_B.LDimSizeP1 += 2) {
              tmp = _mm_loadu_pd(&qrmanager->Q.data[((nmpcSine_B.k_k +
                qrmanager->ncols) * qrmanager->Q.size[0] + nmpcSine_B.LDimSizeP1)
                                 - 1]);
              _mm_storeu_pd(&memspace->workspace_float.data
                            [(nmpcSine_B.LDimSizeP1 +
                              memspace->workspace_float.size[0] * nmpcSine_B.k_k)
                            - 1], _mm_mul_pd(tmp, _mm_set1_pd(objective->beta)));
            }

            for (nmpcSine_B.LDimSizeP1 = nmpcSine_B.A_maxDiag_idx;
                 nmpcSine_B.LDimSizeP1 <= nmpcSine_B.nVars;
                 nmpcSine_B.LDimSizeP1++) {
              memspace->workspace_float.data[(nmpcSine_B.LDimSizeP1 +
                memspace->workspace_float.size[0] * nmpcSine_B.k_k) - 1] =
                qrmanager->Q.data[((nmpcSine_B.k_k + qrmanager->ncols) *
                                   qrmanager->Q.size[0] + nmpcSine_B.LDimSizeP1)
                - 1] * objective->beta;
            }
          }

          nmpcSine_xgemm_k(nmpcSine_B.mNull_tmp, nmpcSine_B.mNull_tmp,
                           qrmanager->mrows, qrmanager->Q.data, nmpcSine_B.b_idx,
                           qrmanager->ldq, memspace->workspace_float.data,
                           memspace->workspace_float.size[0],
                           cholmanager->FMat.data, cholmanager->ldm);
        }

        if (alwaysPositiveDef) {
          cholmanager->ndims = nmpcSine_B.mNull_tmp;
          cholmanager->info = nmpcSine_xpotrf(nmpcSine_B.mNull_tmp,
            cholmanager->FMat.data, cholmanager->ldm);
        } else {
          nmpcSine_B.LDimSizeP1 = cholmanager->ldm + 1;
          cholmanager->ndims = nmpcSine_B.mNull_tmp;
          nmpcSine_B.A_maxDiag_idx = 0;
          if (nmpcSine_B.mNull_tmp > 1) {
            nmpcSine_B.smax_d = fabs(cholmanager->FMat.data[0]);
            for (nmpcSine_B.k_k = 2; nmpcSine_B.k_k <= nmpcSine_B.mNull_tmp;
                 nmpcSine_B.k_k++) {
              nmpcSine_B.s_n = fabs(cholmanager->FMat.data[(nmpcSine_B.k_k - 1) *
                                    nmpcSine_B.LDimSizeP1]);
              if (nmpcSine_B.s_n > nmpcSine_B.smax_d) {
                nmpcSine_B.A_maxDiag_idx = nmpcSine_B.k_k - 1;
                nmpcSine_B.smax_d = nmpcSine_B.s_n;
              }
            }
          }

          nmpcSine_B.smax_d = fabs(cholmanager->FMat.data[cholmanager->ldm *
            nmpcSine_B.A_maxDiag_idx + nmpcSine_B.A_maxDiag_idx]) *
            2.2204460492503131E-16;
          if (nmpcSine_B.smax_d >= 0.0) {
            cholmanager->regTol_ = nmpcSine_B.smax_d;
          } else {
            cholmanager->regTol_ = 0.0;
          }

          if (nmpcSine_B.mNull_tmp > 128) {
            nmpcSine_B.k_k = 0;
            exitg1 = false;
            while ((!exitg1) && (nmpcSine_B.k_k < nmpcSine_B.mNull_tmp)) {
              nmpcSine_B.A_maxDiag_idx = nmpcSine_B.LDimSizeP1 * nmpcSine_B.k_k
                + 1;
              nmpcSine_B.order = nmpcSine_B.mNull_tmp - nmpcSine_B.k_k;
              if (nmpcSine_B.k_k + 48 <= nmpcSine_B.mNull_tmp) {
                nmpcSine_partialColLDL3_(cholmanager, nmpcSine_B.A_maxDiag_idx,
                  nmpcSine_B.order);
                nmpcSine_B.k_k += 48;
              } else {
                nmpcSine_fullColLDL2_(cholmanager, nmpcSine_B.A_maxDiag_idx,
                                      nmpcSine_B.order);
                exitg1 = true;
              }
            }
          } else {
            nmpcSine_fullColLDL2_(cholmanager, 1, nmpcSine_B.mNull_tmp);
          }

          if (cholmanager->ConvexCheck) {
            nmpcSine_B.k_k = 0;
            do {
              exitg2 = 0;
              if (nmpcSine_B.k_k <= nmpcSine_B.mNull_tmp - 1) {
                if (cholmanager->FMat.data[cholmanager->ldm * nmpcSine_B.k_k +
                    nmpcSine_B.k_k] <= 0.0) {
                  cholmanager->info = -nmpcSine_B.k_k - 1;
                  exitg2 = 1;
                } else {
                  nmpcSine_B.k_k++;
                }
              } else {
                cholmanager->ConvexCheck = false;
                exitg2 = 1;
              }
            } while (exitg2 == 0);
          }
        }

        if (cholmanager->info != 0) {
          solution->state = -6;
        } else {
          nmpcSine_B.A_maxDiag_idx = qrmanager->ldq;
          if (qrmanager->mrows != 0) {
            memset(&memspace->workspace_float.data[0], 0, static_cast<uint32_T>
                   (nmpcSine_B.mNull_tmp) * sizeof(real_T));
            nmpcSine_B.order = 0;
            nmpcSine_B.nVars = (nmpcSine_B.mNull_tmp - 1) * qrmanager->ldq +
              nmpcSine_B.b_idx;
            for (nmpcSine_B.LDimSizeP1 = nmpcSine_B.b_idx;
                 nmpcSine_B.A_maxDiag_idx < 0 ? nmpcSine_B.LDimSizeP1 >=
                 nmpcSine_B.nVars : nmpcSine_B.LDimSizeP1 <= nmpcSine_B.nVars;
                 nmpcSine_B.LDimSizeP1 += nmpcSine_B.A_maxDiag_idx) {
              nmpcSine_B.smax_d = 0.0;
              nmpcSine_B.c_ix = nmpcSine_B.LDimSizeP1 + nmpcSine_B.nVar_he;
              for (nmpcSine_B.k_k = nmpcSine_B.LDimSizeP1; nmpcSine_B.k_k <=
                   nmpcSine_B.c_ix; nmpcSine_B.k_k++) {
                nmpcSine_B.smax_d += qrmanager->Q.data[nmpcSine_B.k_k - 1] *
                  objective->grad.data[nmpcSine_B.k_k - nmpcSine_B.LDimSizeP1];
              }

              memspace->workspace_float.data[nmpcSine_B.order] -=
                nmpcSine_B.smax_d;
              nmpcSine_B.order++;
            }
          }

          if (alwaysPositiveDef) {
            nmpcSine_B.c_ix = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (nmpcSine_B.nVar_he = 0; nmpcSine_B.nVar_he < nmpcSine_B.c_ix;
                   nmpcSine_B.nVar_he++) {
                nmpcSine_B.A_maxDiag_idx = nmpcSine_B.nVar_he * cholmanager->ldm;
                nmpcSine_B.smax_d = memspace->
                  workspace_float.data[nmpcSine_B.nVar_he];
                for (nmpcSine_B.k_k = 0; nmpcSine_B.k_k < nmpcSine_B.nVar_he;
                     nmpcSine_B.k_k++) {
                  nmpcSine_B.smax_d -= cholmanager->
                    FMat.data[nmpcSine_B.A_maxDiag_idx + nmpcSine_B.k_k] *
                    memspace->workspace_float.data[nmpcSine_B.k_k];
                }

                memspace->workspace_float.data[nmpcSine_B.nVar_he] =
                  nmpcSine_B.smax_d / cholmanager->
                  FMat.data[nmpcSine_B.A_maxDiag_idx + nmpcSine_B.nVar_he];
              }
            }

            if (cholmanager->ndims != 0) {
              for (nmpcSine_B.nVar_he = nmpcSine_B.c_ix; nmpcSine_B.nVar_he >= 1;
                   nmpcSine_B.nVar_he--) {
                nmpcSine_B.order = ((nmpcSine_B.nVar_he - 1) * cholmanager->ldm
                                    + nmpcSine_B.nVar_he) - 2;
                memspace->workspace_float.data[nmpcSine_B.nVar_he - 1] /=
                  cholmanager->FMat.data[nmpcSine_B.order + 1];
                for (nmpcSine_B.LDimSizeP1 = 0; nmpcSine_B.LDimSizeP1 <=
                     nmpcSine_B.nVar_he - 2; nmpcSine_B.LDimSizeP1++) {
                  nmpcSine_B.nVars = (nmpcSine_B.nVar_he - nmpcSine_B.LDimSizeP1)
                    - 2;
                  memspace->workspace_float.data[nmpcSine_B.nVars] -=
                    memspace->workspace_float.data[nmpcSine_B.nVar_he - 1] *
                    cholmanager->FMat.data[nmpcSine_B.order -
                    nmpcSine_B.LDimSizeP1];
                }
              }
            }
          } else {
            nmpcSine_B.A_maxDiag_idx_tmp = cholmanager->ndims;
            if (cholmanager->ndims != 0) {
              for (nmpcSine_B.nVar_he = 0; nmpcSine_B.nVar_he <
                   nmpcSine_B.A_maxDiag_idx_tmp; nmpcSine_B.nVar_he++) {
                nmpcSine_B.order = nmpcSine_B.nVar_he * cholmanager->ldm +
                  nmpcSine_B.nVar_he;
                nmpcSine_B.nVars = (nmpcSine_B.A_maxDiag_idx_tmp -
                                    nmpcSine_B.nVar_he) - 2;
                for (nmpcSine_B.LDimSizeP1 = 0; nmpcSine_B.LDimSizeP1 <=
                     nmpcSine_B.nVars; nmpcSine_B.LDimSizeP1++) {
                  nmpcSine_B.c_ix = (nmpcSine_B.LDimSizeP1 + nmpcSine_B.nVar_he)
                    + 1;
                  memspace->workspace_float.data[nmpcSine_B.c_ix] -=
                    cholmanager->FMat.data[(nmpcSine_B.LDimSizeP1 +
                    nmpcSine_B.order) + 1] * memspace->
                    workspace_float.data[nmpcSine_B.nVar_he];
                }
              }
            }

            for (nmpcSine_B.nVar_he = 0; nmpcSine_B.nVar_he <
                 nmpcSine_B.A_maxDiag_idx_tmp; nmpcSine_B.nVar_he++) {
              memspace->workspace_float.data[nmpcSine_B.nVar_he] /=
                cholmanager->FMat.data[cholmanager->ldm * nmpcSine_B.nVar_he +
                nmpcSine_B.nVar_he];
            }

            if (cholmanager->ndims != 0) {
              for (nmpcSine_B.nVar_he = nmpcSine_B.A_maxDiag_idx_tmp;
                   nmpcSine_B.nVar_he >= 1; nmpcSine_B.nVar_he--) {
                nmpcSine_B.A_maxDiag_idx = (nmpcSine_B.nVar_he - 1) *
                  cholmanager->ldm;
                nmpcSine_B.smax_d = memspace->
                  workspace_float.data[nmpcSine_B.nVar_he - 1];
                for (nmpcSine_B.k_k = nmpcSine_B.A_maxDiag_idx_tmp;
                     nmpcSine_B.k_k >= nmpcSine_B.nVar_he + 1; nmpcSine_B.k_k--)
                {
                  nmpcSine_B.smax_d -= cholmanager->FMat.data
                    [(nmpcSine_B.A_maxDiag_idx + nmpcSine_B.k_k) - 1] *
                    memspace->workspace_float.data[nmpcSine_B.k_k - 1];
                }

                memspace->workspace_float.data[nmpcSine_B.nVar_he - 1] =
                  nmpcSine_B.smax_d;
              }
            }
          }

          nmpcSine_xgemv_dffoy(qrmanager->mrows, nmpcSine_B.mNull_tmp,
                               qrmanager->Q.data, nmpcSine_B.b_idx,
                               qrmanager->ldq, memspace->workspace_float.data,
                               solution->searchDir.data);
        }
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_xnrm2_nl(int32_T n, const real_T x_data[])
{
  real_T y;
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x_data[0]);
    } else {
      real_T scale;
      scale = 3.3121686421112381E-170;
      for (int32_T k = 0; k < n; k++) {
        real_T absxk;
        absxk = fabs(x_data[k]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  return y;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_xgemv_dffoyl(int32_T m, int32_T n, const real_T A_data[],
  int32_T lda, const real_T x_data[], real_T y_data[])
{
  if (n != 0) {
    nmpcSine_B.b_d = static_cast<uint16_T>(n);
    nmpcSine_B.scalarLB_f = (static_cast<uint16_T>(n) / 2) << 1;
    nmpcSine_B.vectorUB_gd = nmpcSine_B.scalarLB_f - 2;
    for (nmpcSine_B.b_iy_mq = 0; nmpcSine_B.b_iy_mq <= nmpcSine_B.vectorUB_gd;
         nmpcSine_B.b_iy_mq += 2) {
      __m128d tmp;
      tmp = _mm_loadu_pd(&y_data[nmpcSine_B.b_iy_mq]);
      _mm_storeu_pd(&y_data[nmpcSine_B.b_iy_mq], _mm_mul_pd(tmp, _mm_set1_pd
        (-1.0)));
    }

    for (nmpcSine_B.b_iy_mq = nmpcSine_B.scalarLB_f; nmpcSine_B.b_iy_mq <
         nmpcSine_B.b_d; nmpcSine_B.b_iy_mq++) {
      y_data[nmpcSine_B.b_iy_mq] = -y_data[nmpcSine_B.b_iy_mq];
    }

    nmpcSine_B.b_d = 0;
    nmpcSine_B.scalarLB_f = (n - 1) * lda + 1;
    for (nmpcSine_B.b_iy_mq = 1; lda < 0 ? nmpcSine_B.b_iy_mq >=
         nmpcSine_B.scalarLB_f : nmpcSine_B.b_iy_mq <= nmpcSine_B.scalarLB_f;
         nmpcSine_B.b_iy_mq += lda) {
      nmpcSine_B.c_n = 0.0;
      nmpcSine_B.vectorUB_gd = (nmpcSine_B.b_iy_mq + m) - 1;
      for (nmpcSine_B.y_tmp = nmpcSine_B.b_iy_mq; nmpcSine_B.y_tmp <=
           nmpcSine_B.vectorUB_gd; nmpcSine_B.y_tmp++) {
        nmpcSine_B.c_n += A_data[nmpcSine_B.y_tmp - 1] * x_data[nmpcSine_B.y_tmp
          - nmpcSine_B.b_iy_mq];
      }

      y_data[nmpcSine_B.b_d] += nmpcSine_B.c_n;
      nmpcSine_B.b_d++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_feasibleratiotest(const real_T solution_xstar_data[],
  const real_T solution_searchDir_data[], real_T workspace_data[], const int32_T
  workspace_size[2], int32_T workingset_nVar, int32_T workingset_ldA, const
  real_T workingset_Aineq_data[], const real_T workingset_bineq_data[], const
  real_T workingset_lb_data[], const real_T workingset_ub_data[], const int32_T
  workingset_indexLB_data[], const int32_T workingset_indexUB_data[], const
  int32_T workingset_sizes[5], const int32_T workingset_isActiveIdx[6], const
  boolean_T workingset_isActiveConstr_data[], const int32_T workingset_nWConstr
  [5], boolean_T isPhaseOne, real_T *alpha, boolean_T *newBlocking, int32_T
  *constrType, int32_T *constrIdx)
{
  *alpha = 1.0E+30;
  *newBlocking = false;
  *constrType = 0;
  *constrIdx = 0;
  nmpcSine_B.denomTol = 2.2204460492503131E-13 * nmpcSine_xnrm2_nl
    (workingset_nVar, solution_searchDir_data);
  if (workingset_nWConstr[2] < workingset_sizes[2]) {
    nmpcSine_B.e_tmp = static_cast<uint16_T>(workingset_sizes[2]);
    if (nmpcSine_B.e_tmp - 1 >= 0) {
      memcpy(&workspace_data[0], &workingset_bineq_data[0], static_cast<uint32_T>
             (nmpcSine_B.e_tmp) * sizeof(real_T));
    }

    nmpcSine_xgemv_dffoyl(workingset_nVar, workingset_sizes[2],
                          workingset_Aineq_data, workingset_ldA,
                          solution_xstar_data, workspace_data);
    nmpcSine_B.k_l = workspace_size[0];
    if (workingset_sizes[2] != 0) {
      nmpcSine_B.iyend = workspace_size[0] + workingset_sizes[2];
      if (nmpcSine_B.k_l + 1 <= nmpcSine_B.iyend) {
        memset(&workspace_data[nmpcSine_B.k_l], 0, static_cast<uint32_T>
               (nmpcSine_B.iyend - nmpcSine_B.k_l) * sizeof(real_T));
      }

      nmpcSine_B.iyend = workspace_size[0];
      nmpcSine_B.f = (workingset_sizes[2] - 1) * workingset_ldA + 1;
      for (nmpcSine_B.b_iy_m = 1; workingset_ldA < 0 ? nmpcSine_B.b_iy_m >=
           nmpcSine_B.f : nmpcSine_B.b_iy_m <= nmpcSine_B.f; nmpcSine_B.b_iy_m +=
           workingset_ldA) {
        nmpcSine_B.b_c_p = 0.0;
        nmpcSine_B.g_nt = (nmpcSine_B.b_iy_m + workingset_nVar) - 1;
        for (nmpcSine_B.ia_g = nmpcSine_B.b_iy_m; nmpcSine_B.ia_g <=
             nmpcSine_B.g_nt; nmpcSine_B.ia_g++) {
          nmpcSine_B.b_c_p += workingset_Aineq_data[nmpcSine_B.ia_g - 1] *
            solution_searchDir_data[nmpcSine_B.ia_g - nmpcSine_B.b_iy_m];
        }

        workspace_data[nmpcSine_B.iyend] += nmpcSine_B.b_c_p;
        nmpcSine_B.iyend++;
      }
    }

    for (nmpcSine_B.b_iy_m = 0; nmpcSine_B.b_iy_m < nmpcSine_B.e_tmp;
         nmpcSine_B.b_iy_m++) {
      nmpcSine_B.phaseOneCorrectionX = workspace_data[nmpcSine_B.k_l +
        nmpcSine_B.b_iy_m];
      if ((nmpcSine_B.phaseOneCorrectionX > nmpcSine_B.denomTol) &&
          (!workingset_isActiveConstr_data[(workingset_isActiveIdx[2] +
            nmpcSine_B.b_iy_m) - 1])) {
        nmpcSine_B.b_c_p = fabs(workspace_data[nmpcSine_B.b_iy_m]);
        nmpcSine_B.phaseOneCorrectionP = 0.001 -
          workspace_data[nmpcSine_B.b_iy_m];
        if ((nmpcSine_B.b_c_p <= nmpcSine_B.phaseOneCorrectionP) || rtIsNaN
            (nmpcSine_B.phaseOneCorrectionP)) {
          nmpcSine_B.phaseOneCorrectionP = nmpcSine_B.b_c_p;
        }

        nmpcSine_B.b_c_p = nmpcSine_B.phaseOneCorrectionP /
          nmpcSine_B.phaseOneCorrectionX;
        if (nmpcSine_B.b_c_p < *alpha) {
          *alpha = nmpcSine_B.b_c_p;
          *constrType = 3;
          *constrIdx = nmpcSine_B.b_iy_m + 1;
          *newBlocking = true;
        }
      }
    }
  }

  if (workingset_nWConstr[3] < workingset_sizes[3]) {
    _mm_storeu_pd(&nmpcSine_B.dv10[0], _mm_mul_pd(_mm_set_pd
      (solution_searchDir_data[workingset_nVar - 1],
       solution_xstar_data[workingset_nVar - 1]), _mm_set1_pd(static_cast<real_T>
      (isPhaseOne))));
    nmpcSine_B.phaseOneCorrectionX = nmpcSine_B.dv10[0];
    nmpcSine_B.phaseOneCorrectionP = nmpcSine_B.dv10[1];
    nmpcSine_B.k_l = workingset_sizes[3];
    for (nmpcSine_B.e_tmp = 0; nmpcSine_B.e_tmp <= nmpcSine_B.k_l - 2;
         nmpcSine_B.e_tmp++) {
      nmpcSine_B.b_iy_m = workingset_indexLB_data[nmpcSine_B.e_tmp];
      nmpcSine_B.pk_corrected = -solution_searchDir_data[nmpcSine_B.b_iy_m - 1]
        - nmpcSine_B.phaseOneCorrectionP;
      if ((nmpcSine_B.pk_corrected > nmpcSine_B.denomTol) &&
          (!workingset_isActiveConstr_data[(workingset_isActiveIdx[3] +
            nmpcSine_B.e_tmp) - 1])) {
        nmpcSine_B.ratio = (-solution_xstar_data[nmpcSine_B.b_iy_m - 1] -
                            workingset_lb_data[nmpcSine_B.b_iy_m - 1]) -
          nmpcSine_B.phaseOneCorrectionX;
        nmpcSine_B.b_c_p = fabs(nmpcSine_B.ratio);
        if ((!(nmpcSine_B.b_c_p <= 0.001 - nmpcSine_B.ratio)) && (!rtIsNaN(0.001
              - nmpcSine_B.ratio))) {
          nmpcSine_B.b_c_p = 0.001 - nmpcSine_B.ratio;
        }

        nmpcSine_B.b_c_p /= nmpcSine_B.pk_corrected;
        if (nmpcSine_B.b_c_p < *alpha) {
          *alpha = nmpcSine_B.b_c_p;
          *constrType = 4;
          *constrIdx = nmpcSine_B.e_tmp + 1;
          *newBlocking = true;
        }
      }
    }

    nmpcSine_B.iyend = workingset_indexLB_data[workingset_sizes[3] - 1] - 1;
    nmpcSine_B.phaseOneCorrectionX = -solution_searchDir_data[nmpcSine_B.iyend];
    if ((nmpcSine_B.phaseOneCorrectionX > nmpcSine_B.denomTol) &&
        (!workingset_isActiveConstr_data[(workingset_isActiveIdx[3] +
          workingset_sizes[3]) - 2])) {
      nmpcSine_B.ratio = -solution_xstar_data[nmpcSine_B.iyend] -
        workingset_lb_data[nmpcSine_B.iyend];
      nmpcSine_B.b_c_p = fabs(nmpcSine_B.ratio);
      if ((!(nmpcSine_B.b_c_p <= 0.001 - nmpcSine_B.ratio)) && (!rtIsNaN(0.001 -
            nmpcSine_B.ratio))) {
        nmpcSine_B.b_c_p = 0.001 - nmpcSine_B.ratio;
      }

      nmpcSine_B.b_c_p /= nmpcSine_B.phaseOneCorrectionX;
      if (nmpcSine_B.b_c_p < *alpha) {
        *alpha = nmpcSine_B.b_c_p;
        *constrType = 4;
        *constrIdx = workingset_sizes[3];
        *newBlocking = true;
      }
    }
  }

  if (workingset_nWConstr[4] < workingset_sizes[4]) {
    _mm_storeu_pd(&nmpcSine_B.dv10[0], _mm_mul_pd(_mm_set_pd
      (solution_searchDir_data[workingset_nVar - 1],
       solution_xstar_data[workingset_nVar - 1]), _mm_set1_pd(static_cast<real_T>
      (isPhaseOne))));
    nmpcSine_B.phaseOneCorrectionX = nmpcSine_B.dv10[0];
    nmpcSine_B.phaseOneCorrectionP = nmpcSine_B.dv10[1];
    nmpcSine_B.k_l = static_cast<uint16_T>(workingset_sizes[4]);
    for (nmpcSine_B.e_tmp = 0; nmpcSine_B.e_tmp < nmpcSine_B.k_l;
         nmpcSine_B.e_tmp++) {
      nmpcSine_B.b_iy_m = workingset_indexUB_data[nmpcSine_B.e_tmp];
      nmpcSine_B.pk_corrected = solution_searchDir_data[nmpcSine_B.b_iy_m - 1] -
        nmpcSine_B.phaseOneCorrectionP;
      if ((nmpcSine_B.pk_corrected > nmpcSine_B.denomTol) &&
          (!workingset_isActiveConstr_data[(workingset_isActiveIdx[4] +
            nmpcSine_B.e_tmp) - 1])) {
        nmpcSine_B.ratio = (solution_xstar_data[nmpcSine_B.b_iy_m - 1] -
                            workingset_ub_data[nmpcSine_B.b_iy_m - 1]) -
          nmpcSine_B.phaseOneCorrectionX;
        nmpcSine_B.b_c_p = fabs(nmpcSine_B.ratio);
        if ((!(nmpcSine_B.b_c_p <= 0.001 - nmpcSine_B.ratio)) && (!rtIsNaN(0.001
              - nmpcSine_B.ratio))) {
          nmpcSine_B.b_c_p = 0.001 - nmpcSine_B.ratio;
        }

        nmpcSine_B.b_c_p /= nmpcSine_B.pk_corrected;
        if (nmpcSine_B.b_c_p < *alpha) {
          *alpha = nmpcSine_B.b_c_p;
          *constrType = 5;
          *constrIdx = nmpcSine_B.e_tmp + 1;
          *newBlocking = true;
        }
      }
    }
  }

  if (!isPhaseOne) {
    *newBlocking = (((!*newBlocking) || (!(*alpha > 1.0))) && (*newBlocking));
    if (!(*alpha <= 1.0)) {
      *alpha = 1.0;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSi_checkUnboundedOrIllPosed(s_2COE1uYisQtyPYvPjrXP9G_nmpc_T
  *solution, const s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective)
{
  if (objective->objtype == 5) {
    if (nmpcSine_xnrm2_nl(objective->nvar, solution->searchDir.data) > 100.0 *
        static_cast<real_T>(objective->nvar) * 1.4901161193847656E-8) {
      solution->state = 3;
    } else {
      solution->state = 4;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpc_addBoundToActiveSetMatrix_(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *obj, int32_T TYPE, int32_T idx_local)
{
  int32_T colOffset;
  int32_T idx_bnd_local;
  obj->nWConstr[TYPE - 1]++;
  obj->isActiveConstr.data[(obj->isActiveIdx[TYPE - 1] + idx_local) - 2] = true;
  obj->nActiveConstr++;
  obj->Wid.data[obj->nActiveConstr - 1] = TYPE;
  obj->Wlocalidx.data[obj->nActiveConstr - 1] = idx_local;
  colOffset = (obj->nActiveConstr - 1) * obj->ldA - 1;
  if (TYPE == 5) {
    idx_bnd_local = obj->indexUB.data[idx_local - 1];
    obj->bwset.data[obj->nActiveConstr - 1] = obj->ub.data[idx_bnd_local - 1];
  } else {
    idx_bnd_local = obj->indexLB.data[idx_local - 1];
    obj->bwset.data[obj->nActiveConstr - 1] = obj->lb.data[idx_bnd_local - 1];
  }

  if (static_cast<uint16_T>(idx_bnd_local - 1) - 1 >= 0) {
    memset(&obj->ATwset.data[colOffset + 1], 0, static_cast<uint16_T>
           (idx_bnd_local - 1) * sizeof(real_T));
  }

  obj->ATwset.data[idx_bnd_local + colOffset] = static_cast<real_T>(TYPE == 5) *
    2.0 - 1.0;
  if (idx_bnd_local + 1 <= obj->nVar) {
    memset(&obj->ATwset.data[(idx_bnd_local + colOffset) + 1], 0,
           static_cast<uint32_T>(((obj->nVar + colOffset) - idx_bnd_local) -
            colOffset) * sizeof(real_T));
  }

  switch (obj->probType) {
   case 3:
   case 2:
    break;

   default:
    obj->ATwset.data[obj->nVar + colOffset] = -1.0;
    break;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_addAineqConstr(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
  int32_T idx_local)
{
  int32_T b;
  int32_T iAineq0;
  int32_T iAw0;
  obj->nWConstr[2]++;
  obj->isActiveConstr.data[(obj->isActiveIdx[2] + idx_local) - 2] = true;
  obj->nActiveConstr++;
  obj->Wid.data[obj->nActiveConstr - 1] = 3;
  obj->Wlocalidx.data[obj->nActiveConstr - 1] = idx_local;
  iAineq0 = (idx_local - 1) * obj->ldA;
  iAw0 = (obj->nActiveConstr - 1) * obj->ldA;
  b = obj->nVar;
  for (int32_T idx = 0; idx < b; idx++) {
    obj->ATwset.data[iAw0 + idx] = obj->Aineq.data[iAineq0 + idx];
  }

  obj->bwset.data[obj->nActiveConstr - 1] = obj->bineq.data[idx_local - 1];
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_compute_lambda(real_T workspace_data[],
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, const
  s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager)
{
  int32_T nActiveConstr_tmp_tmp;
  nActiveConstr_tmp_tmp = qrmanager->ncols;
  if (qrmanager->ncols > 0) {
    real_T c;
    int32_T b_idx;
    int32_T b_ix;
    int32_T idxQR;
    boolean_T guard1;
    guard1 = false;
    if (objective->objtype != 4) {
      boolean_T nonDegenerate;
      if (qrmanager->mrows >= qrmanager->ncols) {
        b_ix = qrmanager->mrows;
      } else {
        b_ix = qrmanager->ncols;
      }

      c = 2.2204460492503131E-15 * static_cast<real_T>(b_ix);
      if (c >= 1.4901161193847656E-8) {
        c = 1.4901161193847656E-8;
      }

      nonDegenerate = ((qrmanager->mrows > 0) && (qrmanager->ncols > 0));
      if (nonDegenerate) {
        boolean_T guard2;
        b_idx = qrmanager->ncols;
        guard2 = false;
        if (qrmanager->mrows < qrmanager->ncols) {
          idxQR = (qrmanager->ncols - 1) * qrmanager->ldq + qrmanager->mrows;
          while ((b_idx > qrmanager->mrows) && (fabs(qrmanager->QR.data[idxQR -
                   1]) >= c)) {
            b_idx--;
            idxQR -= qrmanager->ldq;
          }

          nonDegenerate = (b_idx == qrmanager->mrows);
          if (!nonDegenerate) {
          } else {
            guard2 = true;
          }
        } else {
          guard2 = true;
        }

        if (guard2) {
          idxQR = (b_idx - 1) * qrmanager->ldq + b_idx;
          while ((b_idx >= 1) && (fabs(qrmanager->QR.data[idxQR - 1]) >= c)) {
            b_idx--;
            idxQR = (idxQR - qrmanager->ldq) - 1;
          }

          nonDegenerate = (b_idx == 0);
        }
      }

      if (!nonDegenerate) {
        solution->state = -7;
      } else {
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      int32_T ia;
      int32_T jjA;
      b_idx = qrmanager->ldq;
      if (qrmanager->mrows != 0) {
        memset(&workspace_data[0], 0, static_cast<uint32_T>(qrmanager->ncols) *
               sizeof(real_T));
        b_ix = 0;
        jjA = (qrmanager->ncols - 1) * qrmanager->ldq + 1;
        for (idxQR = 1; b_idx < 0 ? idxQR >= jjA : idxQR <= jjA; idxQR += b_idx)
        {
          int32_T d;
          c = 0.0;
          d = (idxQR + qrmanager->mrows) - 1;
          for (ia = idxQR; ia <= d; ia++) {
            c += qrmanager->Q.data[ia - 1] * objective->grad.data[ia - idxQR];
          }

          workspace_data[b_ix] += c;
          b_ix++;
        }
      }

      if (qrmanager->ncols != 0) {
        for (idxQR = nActiveConstr_tmp_tmp; idxQR >= 1; idxQR--) {
          jjA = ((idxQR - 1) * b_idx + idxQR) - 2;
          workspace_data[idxQR - 1] /= qrmanager->QR.data[jjA + 1];
          for (ia = 0; ia <= idxQR - 2; ia++) {
            b_ix = (idxQR - ia) - 2;
            workspace_data[b_ix] -= workspace_data[idxQR - 1] *
              qrmanager->QR.data[jjA - ia];
          }
        }
      }

      idxQR = (qrmanager->ncols / 2) << 1;
      ia = idxQR - 2;
      for (b_idx = 0; b_idx <= ia; b_idx += 2) {
        __m128d tmp;
        tmp = _mm_loadu_pd(&workspace_data[b_idx]);
        _mm_storeu_pd(&solution->lambda.data[b_idx], _mm_mul_pd(tmp, _mm_set1_pd
          (-1.0)));
      }

      for (b_idx = idxQR; b_idx < nActiveConstr_tmp_tmp; b_idx++) {
        solution->lambda.data[b_idx] = -workspace_data[b_idx];
      }
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nm_checkStoppingAndUpdateFval_f(int32_T *activeSetChangeID,
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, const s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective,
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
  *qrmanager, int32_T runTimeOptions_MaxIterations, boolean_T *updateFval)
{
  solution->iterations++;
  nmpcSine_B.nVar_jc = objective->nvar;
  if ((solution->iterations >= runTimeOptions_MaxIterations) &&
      ((solution->state != 1) || (objective->objtype == 5))) {
    solution->state = 0;
  }

  if (solution->iterations - solution->iterations / 50 * 50 == 0) {
    nmpcSine_B.tempMaxConstr_k = nmpc_maxConstraintViolation_ini(workingset,
      solution->xstar.data);
    solution->maxConstr = nmpcSine_B.tempMaxConstr_k;
    if (objective->objtype == 5) {
      nmpcSine_B.tempMaxConstr_k = solution->maxConstr - solution->
        xstar.data[objective->nvar - 1];
    }

    if (nmpcSine_B.tempMaxConstr_k > 0.001) {
      if (nmpcSine_B.nVar_jc - 1 >= 0) {
        memcpy(&solution->searchDir.data[0], &solution->xstar.data[0],
               static_cast<uint32_T>(nmpcSine_B.nVar_jc) * sizeof(real_T));
      }

      nmpcSine_B.nonDegenerateWset_i = nmpcSin_feasibleX0ForWorkingSet
        (memspace->workspace_float.data, memspace->workspace_float.size,
         solution->searchDir.data, workingset, qrmanager);
      if ((!nmpcSine_B.nonDegenerateWset_i) && (solution->state != 0)) {
        solution->state = -2;
      }

      *activeSetChangeID = 0;
      nmpcSine_B.tempMaxConstr_k = nmpc_maxConstraintViolation_ini(workingset,
        solution->searchDir.data);
      if (nmpcSine_B.tempMaxConstr_k < solution->maxConstr) {
        if (nmpcSine_B.nVar_jc - 1 >= 0) {
          memcpy(&solution->xstar.data[0], &solution->searchDir.data[0],
                 static_cast<uint32_T>(nmpcSine_B.nVar_jc) * sizeof(real_T));
        }

        solution->maxConstr = nmpcSine_B.tempMaxConstr_k;
      }
    }
  }

  if (*updateFval) {
    *updateFval = false;
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_iterate_d(const real_T H[15625], const real_T f_data[],
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T
  *cholmanager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const char_T
  options_SolverName[7], int32_T runTimeOptions_MaxIterations)
{
  __m128d tmp;
  __m128d tmp_0;
  int32_T exitg1;
  boolean_T guard1;
  nmpcSine_B.subProblemChanged_b = true;
  nmpcSine_B.updateFval_a = true;
  nmpcSine_B.activeSetChangeID_h = 0;
  nmpcSine_B.nVar_h = workingset->nVar;
  nmpcSine_B.globalActiveConstrIdx_f = 0;
  nmpcSine_computeGrad_StoreHx(objective, H, f_data, solution->xstar.data);
  solution->fstar = nmpcSine_computeFval_ReuseHx(objective,
    memspace->workspace_float.data, f_data, solution->xstar.data);
  if (solution->iterations < runTimeOptions_MaxIterations) {
    solution->state = -5;
  } else {
    solution->state = 0;
  }

  nmpcSine_B.idxMinLambda_i = workingset->mConstrMax;
  if (nmpcSine_B.idxMinLambda_i - 1 >= 0) {
    memset(&solution->lambda.data[0], 0, static_cast<uint32_T>
           (nmpcSine_B.idxMinLambda_i) * sizeof(real_T));
  }

  do {
    exitg1 = 0;
    if (solution->state == -5) {
      guard1 = false;
      if (nmpcSine_B.subProblemChanged_b) {
        switch (nmpcSine_B.activeSetChangeID_h) {
         case 1:
          nmpcSine_squareQ_appendCol(qrmanager, workingset->ATwset.data,
            workingset->ldA * (workingset->nActiveConstr - 1) + 1);
          break;

         case -1:
          nmpcSine_deleteColMoveEnd(qrmanager,
            nmpcSine_B.globalActiveConstrIdx_f);
          break;

         default:
          nmpcSine_factorQR(qrmanager, workingset->ATwset.data,
                            nmpcSine_B.nVar_h, workingset->nActiveConstr,
                            workingset->ldA);
          nmpcSine_B.g_n = qrmanager->minRowCol;
          for (nmpcSine_B.k_f4 = 0; nmpcSine_B.k_f4 < nmpcSine_B.g_n;
               nmpcSine_B.k_f4++) {
            nmpcSine_B.iQR0_c = qrmanager->ldq * nmpcSine_B.k_f4 +
              nmpcSine_B.k_f4;
            nmpcSine_B.idxMinLambda_i = qrmanager->mrows - nmpcSine_B.k_f4;
            if (nmpcSine_B.idxMinLambda_i - 2 >= 0) {
              memcpy(&qrmanager->Q.data[nmpcSine_B.iQR0_c + 1],
                     &qrmanager->QR.data[nmpcSine_B.iQR0_c + 1],
                     static_cast<uint32_T>(nmpcSine_B.idxMinLambda_i - 1) *
                     sizeof(real_T));
            }
          }

          nmpcSine_xorgqr(qrmanager->mrows, qrmanager->mrows,
                          qrmanager->minRowCol, qrmanager->Q.data,
                          qrmanager->Q.size, qrmanager->ldq, qrmanager->tau.data);
          break;
        }

        nmpcSine_compute_deltax(H, solution, memspace, qrmanager, cholmanager,
          objective, nmpcSine_strcmp(options_SolverName));
        if (solution->state != -5) {
          exitg1 = 1;
        } else {
          nmpcSine_B.normDelta_g = nmpcSine_xnrm2_nl(nmpcSine_B.nVar_h,
            solution->searchDir.data);
          guard1 = true;
        }
      } else {
        if (nmpcSine_B.nVar_h - 1 >= 0) {
          memset(&solution->searchDir.data[0], 0, static_cast<uint32_T>
                 (nmpcSine_B.nVar_h) * sizeof(real_T));
        }

        nmpcSine_B.normDelta_g = 0.0;
        guard1 = true;
      }

      if (guard1) {
        if ((!nmpcSine_B.subProblemChanged_b) || (nmpcSine_B.normDelta_g <
             1.0E-6) || (workingset->nActiveConstr >= nmpcSine_B.nVar_h)) {
          nmpcSine_compute_lambda(memspace->workspace_float.data, solution,
            objective, qrmanager);
          if ((solution->state != -7) || (workingset->nActiveConstr >
               nmpcSine_B.nVar_h)) {
            nmpcSine_B.idxMinLambda_i = 0;
            nmpcSine_B.normDelta_g = 0.0;
            nmpcSine_B.g_n = (workingset->nWConstr[0] + workingset->nWConstr[1])
              + 1;
            nmpcSine_B.iQR0_c = workingset->nActiveConstr;
            for (nmpcSine_B.k_f4 = nmpcSine_B.g_n; nmpcSine_B.k_f4 <=
                 nmpcSine_B.iQR0_c; nmpcSine_B.k_f4++) {
              nmpcSine_B.solution_lambda_n = solution->
                lambda.data[nmpcSine_B.k_f4 - 1];
              if (nmpcSine_B.solution_lambda_n < nmpcSine_B.normDelta_g) {
                nmpcSine_B.normDelta_g = nmpcSine_B.solution_lambda_n;
                nmpcSine_B.idxMinLambda_i = nmpcSine_B.k_f4;
              }
            }

            if (nmpcSine_B.idxMinLambda_i == 0) {
              solution->state = 1;
            } else {
              nmpcSine_B.activeSetChangeID_h = -1;
              nmpcSine_B.globalActiveConstrIdx_f = nmpcSine_B.idxMinLambda_i;
              nmpcSine_B.subProblemChanged_b = true;
              nmpcSine_removeConstr(workingset, nmpcSine_B.idxMinLambda_i);
              if (nmpcSine_B.idxMinLambda_i < workingset->nActiveConstr + 1) {
                solution->lambda.data[nmpcSine_B.idxMinLambda_i - 1] =
                  solution->lambda.data[workingset->nActiveConstr];
              }

              solution->lambda.data[workingset->nActiveConstr] = 0.0;
            }
          } else {
            nmpcSine_B.idxMinLambda_i = workingset->nActiveConstr;
            nmpcSine_B.activeSetChangeID_h = 0;
            nmpcSine_B.globalActiveConstrIdx_f = workingset->nActiveConstr;
            nmpcSine_B.subProblemChanged_b = true;
            nmpcSine_removeConstr(workingset, workingset->nActiveConstr);
            solution->lambda.data[nmpcSine_B.idxMinLambda_i - 1] = 0.0;
          }

          nmpcSine_B.updateFval_a = false;
        } else {
          nmpcSine_feasibleratiotest(solution->xstar.data,
            solution->searchDir.data, memspace->workspace_float.data,
            memspace->workspace_float.size, workingset->nVar, workingset->ldA,
            workingset->Aineq.data, workingset->bineq.data, workingset->lb.data,
            workingset->ub.data, workingset->indexLB.data,
            workingset->indexUB.data, workingset->sizes, workingset->isActiveIdx,
            workingset->isActiveConstr.data, workingset->nWConstr,
            (objective->objtype == 5), &nmpcSine_B.normDelta_g,
            &nmpcSine_B.updateFval_a, &nmpcSine_B.k_f4,
            &nmpcSine_B.idxMinLambda_i);
          if (nmpcSine_B.updateFval_a) {
            switch (nmpcSine_B.k_f4) {
             case 3:
              nmpcSine_addAineqConstr(workingset, nmpcSine_B.idxMinLambda_i);
              break;

             case 4:
              nmpc_addBoundToActiveSetMatrix_(workingset, 4,
                nmpcSine_B.idxMinLambda_i);
              break;

             default:
              nmpc_addBoundToActiveSetMatrix_(workingset, 5,
                nmpcSine_B.idxMinLambda_i);
              break;
            }

            nmpcSine_B.activeSetChangeID_h = 1;
          } else {
            nmpcSi_checkUnboundedOrIllPosed(solution, objective);
            nmpcSine_B.subProblemChanged_b = false;
            if (workingset->nActiveConstr == 0) {
              solution->state = 1;
            }
          }

          if ((nmpcSine_B.nVar_h >= 1) && (!(nmpcSine_B.normDelta_g == 0.0))) {
            nmpcSine_B.idxMinLambda_i = (nmpcSine_B.nVar_h / 2) << 1;
            nmpcSine_B.g_n = nmpcSine_B.idxMinLambda_i - 2;
            for (nmpcSine_B.k_f4 = 0; nmpcSine_B.k_f4 <= nmpcSine_B.g_n;
                 nmpcSine_B.k_f4 += 2) {
              tmp = _mm_loadu_pd(&solution->searchDir.data[nmpcSine_B.k_f4]);
              tmp_0 = _mm_loadu_pd(&solution->xstar.data[nmpcSine_B.k_f4]);
              _mm_storeu_pd(&solution->xstar.data[nmpcSine_B.k_f4], _mm_add_pd
                            (_mm_mul_pd(_mm_set1_pd(nmpcSine_B.normDelta_g), tmp),
                             tmp_0));
            }

            for (nmpcSine_B.k_f4 = nmpcSine_B.idxMinLambda_i; nmpcSine_B.k_f4 <
                 nmpcSine_B.nVar_h; nmpcSine_B.k_f4++) {
              solution->xstar.data[nmpcSine_B.k_f4] += nmpcSine_B.normDelta_g *
                solution->searchDir.data[nmpcSine_B.k_f4];
            }
          }

          nmpcSine_computeGrad_StoreHx(objective, H, f_data,
            solution->xstar.data);
          nmpcSine_B.updateFval_a = true;
        }

        nm_checkStoppingAndUpdateFval_f(&nmpcSine_B.activeSetChangeID_h,
          solution, memspace, objective, workingset, qrmanager,
          runTimeOptions_MaxIterations, &nmpcSine_B.updateFval_a);
      }
    } else {
      if (!nmpcSine_B.updateFval_a) {
        solution->fstar = nmpcSine_computeFval_ReuseHx(objective,
          memspace->workspace_float.data, f_data, solution->xstar.data);
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpc_checkStoppingAndUpdateFval(int32_T *activeSetChangeID, const
  real_T f_data[], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution,
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, const
  s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *workingset, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, int32_T
  runTimeOptions_MaxIterations, const boolean_T *updateFval)
{
  solution->iterations++;
  nmpcSine_B.nVar_c = objective->nvar;
  if ((solution->iterations >= runTimeOptions_MaxIterations) &&
      ((solution->state != 1) || (objective->objtype == 5))) {
    solution->state = 0;
  }

  if (solution->iterations - solution->iterations / 50 * 50 == 0) {
    solution->maxConstr = nmpc_maxConstraintViolation_ini(workingset,
      solution->xstar.data);
    nmpcSine_B.tempMaxConstr = solution->maxConstr;
    if (objective->objtype == 5) {
      nmpcSine_B.tempMaxConstr = solution->maxConstr - solution->
        xstar.data[objective->nvar - 1];
    }

    if (nmpcSine_B.tempMaxConstr > 0.001) {
      if (nmpcSine_B.nVar_c - 1 >= 0) {
        memcpy(&solution->searchDir.data[0], &solution->xstar.data[0],
               static_cast<uint32_T>(nmpcSine_B.nVar_c) * sizeof(real_T));
      }

      nmpcSine_B.nonDegenerateWset = nmpcSin_feasibleX0ForWorkingSet
        (memspace->workspace_float.data, memspace->workspace_float.size,
         solution->searchDir.data, workingset, qrmanager);
      if ((!nmpcSine_B.nonDegenerateWset) && (solution->state != 0)) {
        solution->state = -2;
      }

      *activeSetChangeID = 0;
      nmpcSine_B.tempMaxConstr = nmpc_maxConstraintViolation_ini(workingset,
        solution->searchDir.data);
      if (nmpcSine_B.tempMaxConstr < solution->maxConstr) {
        if (nmpcSine_B.nVar_c - 1 >= 0) {
          memcpy(&solution->xstar.data[0], &solution->searchDir.data[0],
                 static_cast<uint32_T>(nmpcSine_B.nVar_c) * sizeof(real_T));
        }

        solution->maxConstr = nmpcSine_B.tempMaxConstr;
      }
    }
  }

  if (*updateFval) {
    solution->fstar = nmpcSine_computeFval_ReuseHx(objective,
      memspace->workspace_float.data, f_data, solution->xstar.data);
    if ((solution->fstar < 0.001) && ((solution->state != 0) ||
         (objective->objtype != 5))) {
      solution->state = 2;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_iterate(const real_T H[15625], const real_T f_data[],
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T
  *cholmanager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const char_T
  options_SolverName[7], int32_T runTimeOptions_MaxIterations)
{
  __m128d tmp;
  __m128d tmp_0;
  int32_T exitg1;
  boolean_T guard1;
  nmpcSine_B.subProblemChanged = true;
  nmpcSine_B.updateFval = true;
  nmpcSine_B.activeSetChangeID = 0;
  nmpcSine_B.nVar_js = workingset->nVar;
  nmpcSine_B.globalActiveConstrIdx = 0;
  nmpcSine_computeGrad_StoreHx(objective, H, f_data, solution->xstar.data);
  solution->fstar = nmpcSine_computeFval_ReuseHx(objective,
    memspace->workspace_float.data, f_data, solution->xstar.data);
  solution->state = -5;
  nmpcSine_B.idxMinLambda = workingset->mConstrMax;
  if (nmpcSine_B.idxMinLambda - 1 >= 0) {
    memset(&solution->lambda.data[0], 0, static_cast<uint32_T>
           (nmpcSine_B.idxMinLambda) * sizeof(real_T));
  }

  do {
    exitg1 = 0;
    if (solution->state == -5) {
      guard1 = false;
      if (nmpcSine_B.subProblemChanged) {
        switch (nmpcSine_B.activeSetChangeID) {
         case 1:
          nmpcSine_squareQ_appendCol(qrmanager, workingset->ATwset.data,
            workingset->ldA * (workingset->nActiveConstr - 1) + 1);
          break;

         case -1:
          nmpcSine_deleteColMoveEnd(qrmanager, nmpcSine_B.globalActiveConstrIdx);
          break;

         default:
          nmpcSine_factorQR(qrmanager, workingset->ATwset.data,
                            nmpcSine_B.nVar_js, workingset->nActiveConstr,
                            workingset->ldA);
          nmpcSine_B.g = qrmanager->minRowCol;
          for (nmpcSine_B.k_ho = 0; nmpcSine_B.k_ho < nmpcSine_B.g;
               nmpcSine_B.k_ho++) {
            nmpcSine_B.iQR0 = qrmanager->ldq * nmpcSine_B.k_ho + nmpcSine_B.k_ho;
            nmpcSine_B.idxMinLambda = qrmanager->mrows - nmpcSine_B.k_ho;
            if (nmpcSine_B.idxMinLambda - 2 >= 0) {
              memcpy(&qrmanager->Q.data[nmpcSine_B.iQR0 + 1],
                     &qrmanager->QR.data[nmpcSine_B.iQR0 + 1],
                     static_cast<uint32_T>(nmpcSine_B.idxMinLambda - 1) * sizeof
                     (real_T));
            }
          }

          nmpcSine_xorgqr(qrmanager->mrows, qrmanager->mrows,
                          qrmanager->minRowCol, qrmanager->Q.data,
                          qrmanager->Q.size, qrmanager->ldq, qrmanager->tau.data);
          break;
        }

        nmpcSine_compute_deltax(H, solution, memspace, qrmanager, cholmanager,
          objective, nmpcSine_strcmp(options_SolverName));
        if (solution->state != -5) {
          exitg1 = 1;
        } else {
          nmpcSine_B.normDelta = nmpcSine_xnrm2_nl(nmpcSine_B.nVar_js,
            solution->searchDir.data);
          guard1 = true;
        }
      } else {
        if (nmpcSine_B.nVar_js - 1 >= 0) {
          memset(&solution->searchDir.data[0], 0, static_cast<uint32_T>
                 (nmpcSine_B.nVar_js) * sizeof(real_T));
        }

        nmpcSine_B.normDelta = 0.0;
        guard1 = true;
      }

      if (guard1) {
        if ((!nmpcSine_B.subProblemChanged) || (nmpcSine_B.normDelta <
             1.4901161193847657E-10) || (workingset->nActiveConstr >=
             nmpcSine_B.nVar_js)) {
          nmpcSine_compute_lambda(memspace->workspace_float.data, solution,
            objective, qrmanager);
          if ((solution->state != -7) || (workingset->nActiveConstr >
               nmpcSine_B.nVar_js)) {
            nmpcSine_B.idxMinLambda = 0;
            nmpcSine_B.normDelta = 0.0;
            nmpcSine_B.g = (workingset->nWConstr[0] + workingset->nWConstr[1]) +
              1;
            nmpcSine_B.iQR0 = workingset->nActiveConstr;
            for (nmpcSine_B.k_ho = nmpcSine_B.g; nmpcSine_B.k_ho <=
                 nmpcSine_B.iQR0; nmpcSine_B.k_ho++) {
              nmpcSine_B.solution_lambda = solution->lambda.data[nmpcSine_B.k_ho
                - 1];
              if (nmpcSine_B.solution_lambda < nmpcSine_B.normDelta) {
                nmpcSine_B.normDelta = nmpcSine_B.solution_lambda;
                nmpcSine_B.idxMinLambda = nmpcSine_B.k_ho;
              }
            }

            if (nmpcSine_B.idxMinLambda == 0) {
              solution->state = 1;
            } else {
              nmpcSine_B.activeSetChangeID = -1;
              nmpcSine_B.globalActiveConstrIdx = nmpcSine_B.idxMinLambda;
              nmpcSine_B.subProblemChanged = true;
              nmpcSine_removeConstr(workingset, nmpcSine_B.idxMinLambda);
              if (nmpcSine_B.idxMinLambda < workingset->nActiveConstr + 1) {
                solution->lambda.data[nmpcSine_B.idxMinLambda - 1] =
                  solution->lambda.data[workingset->nActiveConstr];
              }

              solution->lambda.data[workingset->nActiveConstr] = 0.0;
            }
          } else {
            nmpcSine_B.idxMinLambda = workingset->nActiveConstr;
            nmpcSine_B.activeSetChangeID = 0;
            nmpcSine_B.globalActiveConstrIdx = workingset->nActiveConstr;
            nmpcSine_B.subProblemChanged = true;
            nmpcSine_removeConstr(workingset, workingset->nActiveConstr);
            solution->lambda.data[nmpcSine_B.idxMinLambda - 1] = 0.0;
          }

          nmpcSine_B.updateFval = false;
        } else {
          nmpcSine_feasibleratiotest(solution->xstar.data,
            solution->searchDir.data, memspace->workspace_float.data,
            memspace->workspace_float.size, workingset->nVar, workingset->ldA,
            workingset->Aineq.data, workingset->bineq.data, workingset->lb.data,
            workingset->ub.data, workingset->indexLB.data,
            workingset->indexUB.data, workingset->sizes, workingset->isActiveIdx,
            workingset->isActiveConstr.data, workingset->nWConstr, true,
            &nmpcSine_B.normDelta, &nmpcSine_B.updateFval, &nmpcSine_B.k_ho,
            &nmpcSine_B.idxMinLambda);
          if (nmpcSine_B.updateFval) {
            switch (nmpcSine_B.k_ho) {
             case 3:
              nmpcSine_addAineqConstr(workingset, nmpcSine_B.idxMinLambda);
              break;

             case 4:
              nmpc_addBoundToActiveSetMatrix_(workingset, 4,
                nmpcSine_B.idxMinLambda);
              break;

             default:
              nmpc_addBoundToActiveSetMatrix_(workingset, 5,
                nmpcSine_B.idxMinLambda);
              break;
            }

            nmpcSine_B.activeSetChangeID = 1;
          } else {
            nmpcSi_checkUnboundedOrIllPosed(solution, objective);
            nmpcSine_B.subProblemChanged = false;
            if (workingset->nActiveConstr == 0) {
              solution->state = 1;
            }
          }

          if ((nmpcSine_B.nVar_js >= 1) && (!(nmpcSine_B.normDelta == 0.0))) {
            nmpcSine_B.idxMinLambda = (nmpcSine_B.nVar_js / 2) << 1;
            nmpcSine_B.g = nmpcSine_B.idxMinLambda - 2;
            for (nmpcSine_B.k_ho = 0; nmpcSine_B.k_ho <= nmpcSine_B.g;
                 nmpcSine_B.k_ho += 2) {
              tmp = _mm_loadu_pd(&solution->searchDir.data[nmpcSine_B.k_ho]);
              tmp_0 = _mm_loadu_pd(&solution->xstar.data[nmpcSine_B.k_ho]);
              _mm_storeu_pd(&solution->xstar.data[nmpcSine_B.k_ho], _mm_add_pd
                            (_mm_mul_pd(_mm_set1_pd(nmpcSine_B.normDelta), tmp),
                             tmp_0));
            }

            for (nmpcSine_B.k_ho = nmpcSine_B.idxMinLambda; nmpcSine_B.k_ho <
                 nmpcSine_B.nVar_js; nmpcSine_B.k_ho++) {
              solution->xstar.data[nmpcSine_B.k_ho] += nmpcSine_B.normDelta *
                solution->searchDir.data[nmpcSine_B.k_ho];
            }
          }

          nmpcSine_computeGrad_StoreHx(objective, H, f_data,
            solution->xstar.data);
          nmpcSine_B.updateFval = true;
        }

        nmpc_checkStoppingAndUpdateFval(&nmpcSine_B.activeSetChangeID, f_data,
          solution, memspace, objective, workingset, qrmanager,
          runTimeOptions_MaxIterations, &nmpcSine_B.updateFval);
      }
    } else {
      if (!nmpcSine_B.updateFval) {
        solution->fstar = nmpcSine_computeFval_ReuseHx(objective,
          memspace->workspace_float.data, f_data, solution->xstar.data);
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_phaseone(const real_T H[15625], const real_T f_data[],
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T
  *cholmanager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const char_T
  options_SolverName[7], const somzaGboVhDG7PNQS6E98jD_nmpcS_T *runTimeOptions)
{
  boolean_T exitg1;
  nmpcSine_B.PROBTYPE_ORIG = workingset->probType;
  nmpcSine_B.nVar_tmp = workingset->nVar;
  solution->xstar.data[workingset->nVar] = solution->maxConstr + 1.0;
  if (workingset->probType == 3) {
    nmpcSine_B.mConstr_d = 1;
  } else {
    nmpcSine_B.mConstr_d = 4;
  }

  nmpcSine_setProblemType(workingset, nmpcSine_B.mConstr_d);
  nmpcSine_B.idxStartIneq_tmp_f = workingset->nWConstr[0] + workingset->
    nWConstr[1];
  nmpcSine_B.idxStartIneq_d = nmpcSine_B.idxStartIneq_tmp_f + 1;
  nmpcSine_B.idxEndIneq_j = workingset->nActiveConstr;
  for (nmpcSine_B.mConstr_d = nmpcSine_B.idxStartIneq_d; nmpcSine_B.mConstr_d <=
       nmpcSine_B.idxEndIneq_j; nmpcSine_B.mConstr_d++) {
    workingset->isActiveConstr.data[(workingset->isActiveIdx
      [workingset->Wid.data[nmpcSine_B.mConstr_d - 1] - 1] +
      workingset->Wlocalidx.data[nmpcSine_B.mConstr_d - 1]) - 2] = false;
  }

  workingset->nWConstr[2] = 0;
  workingset->nWConstr[3] = 0;
  workingset->nWConstr[4] = 0;
  workingset->nActiveConstr = nmpcSine_B.idxStartIneq_tmp_f;
  objective->prev_objtype = objective->objtype;
  objective->prev_nvar = objective->nvar;
  objective->prev_hasLinear = objective->hasLinear;
  objective->objtype = 5;
  objective->nvar = nmpcSine_B.nVar_tmp + 1;
  objective->gammaScalar = 1.0;
  objective->hasLinear = true;
  solution->fstar = solution->xstar.data[nmpcSine_B.nVar_tmp];
  solution->state = 5;
  nmpcSine_iterate(H, f_data, solution, memspace, workingset, qrmanager,
                   cholmanager, objective, options_SolverName,
                   runTimeOptions->MaxIterations);
  if (workingset->isActiveConstr.data[(workingset->isActiveIdx[3] +
       workingset->sizes[3]) - 2]) {
    nmpcSine_B.mConstr_d = workingset->sizes[0] + 121;
    exitg1 = false;
    while ((!exitg1) && (nmpcSine_B.mConstr_d <= workingset->nActiveConstr)) {
      if ((workingset->Wid.data[nmpcSine_B.mConstr_d - 1] == 4) &&
          (workingset->Wlocalidx.data[nmpcSine_B.mConstr_d - 1] ==
           workingset->sizes[3])) {
        nmpcSine_removeConstr(workingset, nmpcSine_B.mConstr_d);
        exitg1 = true;
      } else {
        nmpcSine_B.mConstr_d++;
      }
    }
  }

  nmpcSine_B.mConstr_d = workingset->nActiveConstr;
  while ((nmpcSine_B.mConstr_d > workingset->sizes[0] + 120) &&
         (nmpcSine_B.mConstr_d > nmpcSine_B.nVar_tmp)) {
    nmpcSine_removeConstr(workingset, nmpcSine_B.mConstr_d);
    nmpcSine_B.mConstr_d--;
  }

  solution->maxConstr = solution->xstar.data[nmpcSine_B.nVar_tmp];
  nmpcSine_setProblemType(workingset, nmpcSine_B.PROBTYPE_ORIG);
  objective->objtype = objective->prev_objtype;
  objective->nvar = objective->prev_nvar;
  objective->hasLinear = objective->prev_hasLinear;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_linearForm_(boolean_T obj_hasLinear, int32_T obj_nvar,
  real_T workspace_data[], const real_T H[15625], const real_T f_data[], const
  real_T x_data[])
{
  int32_T beta1;
  beta1 = 0;
  if (obj_hasLinear) {
    if (obj_nvar - 1 >= 0) {
      memcpy(&workspace_data[0], &f_data[0], static_cast<uint32_T>(obj_nvar) *
             sizeof(real_T));
    }

    beta1 = 1;
  }

  if (obj_nvar != 0) {
    int32_T d;
    int32_T ix;
    if ((beta1 != 1) && (obj_nvar - 1 >= 0)) {
      memset(&workspace_data[0], 0, static_cast<uint32_T>(obj_nvar) * sizeof
             (real_T));
    }

    ix = 0;
    d = (obj_nvar - 1) * obj_nvar + 1;
    for (beta1 = 1; obj_nvar < 0 ? beta1 >= d : beta1 <= d; beta1 += obj_nvar) {
      real_T c;
      int32_T e;
      c = 0.5 * x_data[ix];
      e = (beta1 + obj_nvar) - 1;
      for (int32_T ia = beta1; ia <= e; ia++) {
        int32_T tmp;
        tmp = ia - beta1;
        workspace_data[tmp] += H[ia - 1] * c;
      }

      ix++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_driver_o(const real_T H[15625], const real_T f_data[],
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *solution, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *workingset,
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *qrmanager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T
  *cholmanager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *objective, const
  somzaGboVhDG7PNQS6E98jD_nmpcS_T *options, const
  somzaGboVhDG7PNQS6E98jD_nmpcS_T *runTimeOptions)
{
  __m128d tmp;
  boolean_T guard1;
  solution->iterations = 0;
  nmpcSine_B.nVar_o = workingset->nVar;
  guard1 = false;
  if (workingset->probType == 3) {
    nmpcSine_B.c_ln = static_cast<uint16_T>(workingset->sizes[0]);
    for (nmpcSine_B.ixlast = 0; nmpcSine_B.ixlast < nmpcSine_B.c_ln;
         nmpcSine_B.ixlast++) {
      solution->xstar.data[workingset->indexFixed.data[nmpcSine_B.ixlast] - 1] =
        workingset->ub.data[workingset->indexFixed.data[nmpcSine_B.ixlast] - 1];
    }

    nmpcSine_B.c_ln = static_cast<uint16_T>(workingset->sizes[3]);
    for (nmpcSine_B.ixlast = 0; nmpcSine_B.ixlast < nmpcSine_B.c_ln;
         nmpcSine_B.ixlast++) {
      if (workingset->isActiveConstr.data[(workingset->isActiveIdx[3] +
           nmpcSine_B.ixlast) - 1]) {
        solution->xstar.data[workingset->indexLB.data[nmpcSine_B.ixlast] - 1] =
          -workingset->lb.data[workingset->indexLB.data[nmpcSine_B.ixlast] - 1];
      }
    }

    nmpcSine_B.c_ln = static_cast<uint16_T>(workingset->sizes[4]);
    for (nmpcSine_B.ixlast = 0; nmpcSine_B.ixlast < nmpcSine_B.c_ln;
         nmpcSine_B.ixlast++) {
      if (workingset->isActiveConstr.data[(workingset->isActiveIdx[4] +
           nmpcSine_B.ixlast) - 1]) {
        solution->xstar.data[workingset->indexUB.data[nmpcSine_B.ixlast] - 1] =
          workingset->ub.data[workingset->indexUB.data[nmpcSine_B.ixlast] - 1];
      }
    }

    nmpcSine_PresolveWorkingSet(solution, memspace, workingset, qrmanager);
    if (solution->state < 0) {
    } else {
      guard1 = true;
    }
  } else {
    solution->state = 82;
    guard1 = true;
  }

  if (guard1) {
    solution->iterations = 0;
    solution->maxConstr = nmpc_maxConstraintViolation_ini(workingset,
      solution->xstar.data);
    if (solution->maxConstr > 0.001) {
      nmpcSine_phaseone(H, f_data, solution, memspace, workingset, qrmanager,
                        cholmanager, objective, options->SolverName,
                        runTimeOptions);
      if (solution->state != 0) {
        solution->maxConstr = nmpc_maxConstraintViolation_ini(workingset,
          solution->xstar.data);
        if (solution->maxConstr > 0.001) {
          nmpcSine_B.ixlast = workingset->mConstrMax;
          if (nmpcSine_B.ixlast - 1 >= 0) {
            memset(&solution->lambda.data[0], 0, static_cast<uint32_T>
                   (nmpcSine_B.ixlast) * sizeof(real_T));
          }

          nmpcSine_B.maxConstr_new = 0.0;
          switch (objective->objtype) {
           case 5:
            nmpcSine_B.maxConstr_new = solution->xstar.data[objective->nvar - 1]
              * objective->gammaScalar;
            break;

           case 3:
            nmpcSine_linearForm_(objective->hasLinear, objective->nvar,
                                 memspace->workspace_float.data, H, f_data,
                                 solution->xstar.data);
            if (objective->nvar >= 1) {
              nmpcSine_B.ixlast = objective->nvar;
              for (nmpcSine_B.nVar_o = 0; nmpcSine_B.nVar_o < nmpcSine_B.ixlast;
                   nmpcSine_B.nVar_o++) {
                nmpcSine_B.maxConstr_new += memspace->
                  workspace_float.data[nmpcSine_B.nVar_o] * solution->
                  xstar.data[nmpcSine_B.nVar_o];
              }
            }
            break;

           case 4:
            nmpcSine_linearForm_(objective->hasLinear, objective->nvar,
                                 memspace->workspace_float.data, H, f_data,
                                 solution->xstar.data);
            nmpcSine_B.ixlast = objective->nvar + 1;
            nmpcSine_B.c_ln = objective->maxVar;
            nmpcSine_B.scalarLB_l = (((((objective->maxVar - objective->nvar) -
              1) / 2) << 1) + objective->nvar) + 1;
            nmpcSine_B.vectorUB_g = nmpcSine_B.scalarLB_l - 2;
            for (nmpcSine_B.nVar_o = nmpcSine_B.ixlast; nmpcSine_B.nVar_o <=
                 nmpcSine_B.vectorUB_g; nmpcSine_B.nVar_o += 2) {
              tmp = _mm_loadu_pd(&solution->xstar.data[nmpcSine_B.nVar_o - 1]);
              _mm_storeu_pd(&memspace->workspace_float.data[nmpcSine_B.nVar_o -
                            1], _mm_add_pd(_mm_mul_pd(_mm_set1_pd(0.5 *
                objective->beta), tmp), _mm_set1_pd(objective->rho)));
            }

            for (nmpcSine_B.nVar_o = nmpcSine_B.scalarLB_l; nmpcSine_B.nVar_o <
                 nmpcSine_B.c_ln; nmpcSine_B.nVar_o++) {
              memspace->workspace_float.data[nmpcSine_B.nVar_o - 1] = 0.5 *
                objective->beta * solution->xstar.data[nmpcSine_B.nVar_o - 1] +
                objective->rho;
            }

            if (objective->maxVar - 1 >= 1) {
              for (nmpcSine_B.nVar_o = 0; nmpcSine_B.nVar_o <= nmpcSine_B.c_ln -
                   2; nmpcSine_B.nVar_o++) {
                nmpcSine_B.maxConstr_new += memspace->
                  workspace_float.data[nmpcSine_B.nVar_o] * solution->
                  xstar.data[nmpcSine_B.nVar_o];
              }
            }
            break;
          }

          solution->fstar = nmpcSine_B.maxConstr_new;
          solution->state = -2;
        } else {
          if (solution->maxConstr > 0.0) {
            if (nmpcSine_B.nVar_o - 1 >= 0) {
              memcpy(&solution->searchDir.data[0], &solution->xstar.data[0],
                     static_cast<uint32_T>(nmpcSine_B.nVar_o) * sizeof(real_T));
            }

            nmpcSine_PresolveWorkingSet(solution, memspace, workingset,
              qrmanager);
            nmpcSine_B.maxConstr_new = nmpc_maxConstraintViolation_ini
              (workingset, solution->xstar.data);
            if (nmpcSine_B.maxConstr_new >= solution->maxConstr) {
              solution->maxConstr = nmpcSine_B.maxConstr_new;
              if (nmpcSine_B.nVar_o - 1 >= 0) {
                memcpy(&solution->xstar.data[0], &solution->searchDir.data[0],
                       static_cast<uint32_T>(nmpcSine_B.nVar_o) * sizeof(real_T));
              }
            }
          }

          nmpcSine_iterate_d(H, f_data, solution, memspace, workingset,
                             qrmanager, cholmanager, objective,
                             options->SolverName, runTimeOptions->MaxIterations);
        }
      }
    } else {
      nmpcSine_iterate_d(H, f_data, solution, memspace, workingset, qrmanager,
                         cholmanager, objective, options->SolverName,
                         runTimeOptions->MaxIterations);
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_addAeqConstr(s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj,
  int32_T idx_local)
{
  int32_T totalEq;
  totalEq = obj->nWConstr[0] + obj->nWConstr[1];
  if ((obj->nActiveConstr == totalEq) && (idx_local > obj->nWConstr[1])) {
    int32_T b_idx;
    int32_T iAeq0;
    int32_T iAw0;
    obj->nWConstr[1]++;
    obj->isActiveConstr.data[(obj->isActiveIdx[1] + idx_local) - 2] = true;
    obj->nActiveConstr++;
    obj->Wid.data[obj->nActiveConstr - 1] = 2;
    obj->Wlocalidx.data[obj->nActiveConstr - 1] = idx_local;
    iAeq0 = (idx_local - 1) * obj->ldA;
    iAw0 = (obj->nActiveConstr - 1) * obj->ldA;
    b_idx = static_cast<uint16_T>(obj->nVar);
    for (totalEq = 0; totalEq < b_idx; totalEq++) {
      obj->ATwset.data[iAw0 + totalEq] = obj->Aeq.data[iAeq0 + totalEq];
    }

    obj->bwset.data[obj->nActiveConstr - 1] = obj->beq[idx_local - 1];
  } else {
    int32_T iAeq0;
    int32_T iAw0;
    int32_T iAw0_tmp;
    obj->nActiveConstr++;
    obj->Wid.data[obj->nActiveConstr - 1] = obj->Wid.data[totalEq];
    obj->Wlocalidx.data[obj->nActiveConstr - 1] = obj->Wlocalidx.data[totalEq];
    iAw0_tmp = static_cast<uint16_T>(obj->nVar);
    for (iAeq0 = 0; iAeq0 < iAw0_tmp; iAeq0++) {
      obj->ATwset.data[iAeq0 + obj->ldA * (obj->nActiveConstr - 1)] =
        obj->ATwset.data[obj->ldA * totalEq + iAeq0];
    }

    obj->bwset.data[obj->nActiveConstr - 1] = obj->bwset.data[totalEq];
    obj->nWConstr[1]++;
    obj->isActiveConstr.data[(obj->isActiveIdx[1] + idx_local) - 2] = true;
    obj->Wid.data[totalEq] = 2;
    obj->Wlocalidx.data[totalEq] = idx_local;
    iAeq0 = (idx_local - 1) * obj->ldA;
    iAw0 = obj->ldA * totalEq;
    for (int32_T b_idx = 0; b_idx < iAw0_tmp; b_idx++) {
      obj->ATwset.data[iAw0 + b_idx] = obj->Aeq.data[iAeq0 + b_idx];
    }

    obj->bwset.data[totalEq] = obj->beq[idx_local - 1];
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
boolean_T nmpcSine::nmpcSine_soc(const real_T Hessian[15625], const real_T
  grad_data[], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *WorkingSet, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager,
  s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *
  QPObjective, const somzaGboVhDG7PNQS6E98jD_nmpcS_T *qpoptions)
{
  __m128d tmp;
  __m128d tmp_0;
  boolean_T success;
  nmpcSine_B.nWIneq_old = WorkingSet->nWConstr[2];
  nmpcSine_B.nWLower_old = WorkingSet->nWConstr[3];
  nmpcSine_B.nWUpper_old = WorkingSet->nWConstr[4];
  nmpcSine_B.nVar_j = WorkingSet->nVar;
  nmpcSine_B.mConstrMax_o = WorkingSet->mConstrMax;
  nmpcSine_B.idxStartIneq_tmp = static_cast<uint16_T>(WorkingSet->nVar);
  memcpy(&TrialState->xstarsqp[0], &TrialState->xstarsqp_old[0],
         static_cast<uint32_T>(nmpcSine_B.idxStartIneq_tmp) * sizeof(real_T));
  memcpy(&TrialState->socDirection.data[0], &TrialState->xstar.data[0],
         static_cast<uint32_T>(nmpcSine_B.idxStartIneq_tmp) * sizeof(real_T));
  if (nmpcSine_B.mConstrMax_o - 1 >= 0) {
    memcpy(&TrialState->lambdaStopTest.data[0], &TrialState->lambda.data[0],
           static_cast<uint32_T>(nmpcSine_B.mConstrMax_o) * sizeof(real_T));
  }

  nmpcSine_B.idxIneqOffset = WorkingSet->isActiveIdx[2];
  for (nmpcSine_B.idxStartIneq_f = 0; nmpcSine_B.idxStartIneq_f <= 118;
       nmpcSine_B.idxStartIneq_f += 2) {
    tmp_0 = _mm_loadu_pd(&TrialState->cEq[nmpcSine_B.idxStartIneq_f]);
    _mm_storeu_pd(&WorkingSet->beq[nmpcSine_B.idxStartIneq_f], _mm_mul_pd(tmp_0,
      _mm_set1_pd(-1.0)));
  }

  nmpcSine_B.idx_lower_tmp = WorkingSet->ldA;
  nmpcSine_B.iy = 0;
  nmpcSine_B.l = WorkingSet->ldA * 119 + 1;
  for (nmpcSine_B.idxStartIneq_f = 1; nmpcSine_B.idx_lower_tmp < 0 ?
       nmpcSine_B.idxStartIneq_f >= nmpcSine_B.l : nmpcSine_B.idxStartIneq_f <=
       nmpcSine_B.l; nmpcSine_B.idxStartIneq_f += nmpcSine_B.idx_lower_tmp) {
    nmpcSine_B.b_c_a = 0.0;
    nmpcSine_B.idx_Partition = (nmpcSine_B.idxStartIneq_f + WorkingSet->nVar) -
      1;
    for (nmpcSine_B.idx_Aineq = nmpcSine_B.idxStartIneq_f; nmpcSine_B.idx_Aineq <=
         nmpcSine_B.idx_Partition; nmpcSine_B.idx_Aineq++) {
      nmpcSine_B.b_c_a += WorkingSet->Aeq.data[nmpcSine_B.idx_Aineq - 1] *
        TrialState->searchDir.data[nmpcSine_B.idx_Aineq -
        nmpcSine_B.idxStartIneq_f];
    }

    WorkingSet->beq[nmpcSine_B.iy] += nmpcSine_B.b_c_a;
    nmpcSine_B.iy++;
  }

  for (nmpcSine_B.idxStartIneq_f = 0; nmpcSine_B.idxStartIneq_f < 120;
       nmpcSine_B.idxStartIneq_f++) {
    WorkingSet->bwset.data[WorkingSet->sizes[0] + nmpcSine_B.idxStartIneq_f] =
      WorkingSet->beq[nmpcSine_B.idxStartIneq_f];
  }

  if (WorkingSet->sizes[2] > 0) {
    nmpcSine_B.iy = static_cast<uint16_T>(WorkingSet->sizes[2]);
    nmpcSine_B.idx_Aineq = (static_cast<uint16_T>(WorkingSet->sizes[2]) / 2) <<
      1;
    nmpcSine_B.idx_lower = nmpcSine_B.idx_Aineq - 2;
    for (nmpcSine_B.idxStartIneq_f = 0; nmpcSine_B.idxStartIneq_f <=
         nmpcSine_B.idx_lower; nmpcSine_B.idxStartIneq_f += 2) {
      tmp_0 = _mm_loadu_pd(&TrialState->cIneq.data[nmpcSine_B.idxStartIneq_f]);
      _mm_storeu_pd(&WorkingSet->bineq.data[nmpcSine_B.idxStartIneq_f],
                    _mm_mul_pd(tmp_0, _mm_set1_pd(-1.0)));
    }

    for (nmpcSine_B.idxStartIneq_f = nmpcSine_B.idx_Aineq;
         nmpcSine_B.idxStartIneq_f < nmpcSine_B.iy; nmpcSine_B.idxStartIneq_f++)
    {
      WorkingSet->bineq.data[nmpcSine_B.idxStartIneq_f] =
        -TrialState->cIneq.data[nmpcSine_B.idxStartIneq_f];
    }

    nmpcSine_B.iy = 0;
    nmpcSine_B.l = (WorkingSet->sizes[2] - 1) * WorkingSet->ldA + 1;
    for (nmpcSine_B.idxStartIneq_f = 1; nmpcSine_B.idx_lower_tmp < 0 ?
         nmpcSine_B.idxStartIneq_f >= nmpcSine_B.l : nmpcSine_B.idxStartIneq_f <=
         nmpcSine_B.l; nmpcSine_B.idxStartIneq_f += nmpcSine_B.idx_lower_tmp) {
      nmpcSine_B.b_c_a = 0.0;
      nmpcSine_B.idx_Partition = (nmpcSine_B.idxStartIneq_f + WorkingSet->nVar)
        - 1;
      for (nmpcSine_B.idx_Aineq = nmpcSine_B.idxStartIneq_f;
           nmpcSine_B.idx_Aineq <= nmpcSine_B.idx_Partition;
           nmpcSine_B.idx_Aineq++) {
        nmpcSine_B.b_c_a += WorkingSet->Aineq.data[nmpcSine_B.idx_Aineq - 1] *
          TrialState->searchDir.data[nmpcSine_B.idx_Aineq -
          nmpcSine_B.idxStartIneq_f];
      }

      WorkingSet->bineq.data[nmpcSine_B.iy] += nmpcSine_B.b_c_a;
      nmpcSine_B.iy++;
    }

    nmpcSine_B.idx_Aineq = 1;
    nmpcSine_B.idx_lower = WorkingSet->sizes[2] + 1;
    nmpcSine_B.iy = (WorkingSet->sizes[2] + WorkingSet->sizes[3]) + 1;
    nmpcSine_B.l = WorkingSet->nActiveConstr;
    for (nmpcSine_B.idxStartIneq_f = nmpcSine_B.idxIneqOffset;
         nmpcSine_B.idxStartIneq_f <= nmpcSine_B.l; nmpcSine_B.idxStartIneq_f++)
    {
      switch (WorkingSet->Wid.data[nmpcSine_B.idxStartIneq_f - 1]) {
       case 3:
        nmpcSine_B.idx_Partition = nmpcSine_B.idx_Aineq;
        nmpcSine_B.idx_Aineq++;
        WorkingSet->bwset.data[nmpcSine_B.idxStartIneq_f - 1] =
          WorkingSet->bineq.data[WorkingSet->
          Wlocalidx.data[nmpcSine_B.idxStartIneq_f - 1] - 1];
        break;

       case 4:
        nmpcSine_B.idx_Partition = nmpcSine_B.idx_lower;
        nmpcSine_B.idx_lower++;
        break;

       default:
        nmpcSine_B.idx_Partition = nmpcSine_B.iy;
        nmpcSine_B.iy++;
        break;
      }

      TrialState->workingset_old.data[nmpcSine_B.idx_Partition - 1] =
        WorkingSet->Wlocalidx.data[nmpcSine_B.idxStartIneq_f - 1];
    }
  }

  memcpy(&TrialState->xstar.data[0], &TrialState->xstarsqp[0],
         static_cast<uint32_T>(nmpcSine_B.idxStartIneq_tmp) * sizeof(real_T));
  nmpcSine_driver_o(Hessian, grad_data, TrialState, memspace, WorkingSet,
                    QRManager, CholManager, QPObjective, qpoptions, qpoptions);
  while ((WorkingSet->mEqRemoved > 0) && (WorkingSet->indexEqRemoved
          [WorkingSet->mEqRemoved - 1] >= 1)) {
    nmpcSine_addAeqConstr(WorkingSet, WorkingSet->indexEqRemoved
                          [WorkingSet->mEqRemoved - 1]);
    WorkingSet->mEqRemoved--;
  }

  nmpcSine_B.idxStartIneq_f = static_cast<uint16_T>(nmpcSine_B.nVar_j);
  nmpcSine_B.idx_Aineq = (static_cast<uint16_T>(nmpcSine_B.nVar_j) / 2) << 1;
  nmpcSine_B.idx_lower = nmpcSine_B.idx_Aineq - 2;
  for (nmpcSine_B.idxIneqOffset = 0; nmpcSine_B.idxIneqOffset <=
       nmpcSine_B.idx_lower; nmpcSine_B.idxIneqOffset += 2) {
    tmp_0 = _mm_loadu_pd(&TrialState->socDirection.data[nmpcSine_B.idxIneqOffset]);
    tmp = _mm_loadu_pd(&TrialState->xstar.data[nmpcSine_B.idxIneqOffset]);
    _mm_storeu_pd(&TrialState->socDirection.data[nmpcSine_B.idxIneqOffset],
                  _mm_sub_pd(tmp, tmp_0));
    _mm_storeu_pd(&TrialState->xstar.data[nmpcSine_B.idxIneqOffset], tmp_0);
  }

  for (nmpcSine_B.idxIneqOffset = nmpcSine_B.idx_Aineq; nmpcSine_B.idxIneqOffset
       < nmpcSine_B.idxStartIneq_f; nmpcSine_B.idxIneqOffset++) {
    nmpcSine_B.b_c_a = TrialState->socDirection.data[nmpcSine_B.idxIneqOffset];
    TrialState->socDirection.data[nmpcSine_B.idxIneqOffset] =
      TrialState->xstar.data[nmpcSine_B.idxIneqOffset] - nmpcSine_B.b_c_a;
    TrialState->xstar.data[nmpcSine_B.idxIneqOffset] = nmpcSine_B.b_c_a;
  }

  success = (nmpcSine_xnrm2_nl(nmpcSine_B.nVar_j, TrialState->socDirection.data)
             <= 2.0 * nmpcSine_xnrm2_nl(nmpcSine_B.nVar_j,
              TrialState->xstar.data));
  nmpcSine_B.nVar_j = WorkingSet->sizes[2];
  for (nmpcSine_B.idxIneqOffset = 0; nmpcSine_B.idxIneqOffset <= 118;
       nmpcSine_B.idxIneqOffset += 2) {
    tmp_0 = _mm_loadu_pd(&TrialState->cEq[nmpcSine_B.idxIneqOffset]);
    _mm_storeu_pd(&WorkingSet->beq[nmpcSine_B.idxIneqOffset], _mm_mul_pd(tmp_0,
      _mm_set1_pd(-1.0)));
  }

  for (nmpcSine_B.idxIneqOffset = 0; nmpcSine_B.idxIneqOffset < 120;
       nmpcSine_B.idxIneqOffset++) {
    WorkingSet->bwset.data[WorkingSet->sizes[0] + nmpcSine_B.idxIneqOffset] =
      WorkingSet->beq[nmpcSine_B.idxIneqOffset];
  }

  if (WorkingSet->sizes[2] > 0) {
    nmpcSine_B.idxStartIneq_f = static_cast<uint16_T>(WorkingSet->sizes[2]);
    nmpcSine_B.idx_Aineq = (static_cast<uint16_T>(WorkingSet->sizes[2]) / 2) <<
      1;
    nmpcSine_B.idx_lower = nmpcSine_B.idx_Aineq - 2;
    for (nmpcSine_B.idxIneqOffset = 0; nmpcSine_B.idxIneqOffset <=
         nmpcSine_B.idx_lower; nmpcSine_B.idxIneqOffset += 2) {
      tmp_0 = _mm_loadu_pd(&TrialState->cIneq.data[nmpcSine_B.idxIneqOffset]);
      _mm_storeu_pd(&WorkingSet->bineq.data[nmpcSine_B.idxIneqOffset],
                    _mm_mul_pd(tmp_0, _mm_set1_pd(-1.0)));
    }

    for (nmpcSine_B.idxIneqOffset = nmpcSine_B.idx_Aineq;
         nmpcSine_B.idxIneqOffset < nmpcSine_B.idxStartIneq_f;
         nmpcSine_B.idxIneqOffset++) {
      WorkingSet->bineq.data[nmpcSine_B.idxIneqOffset] = -TrialState->
        cIneq.data[nmpcSine_B.idxIneqOffset];
    }

    if (!success) {
      nmpcSine_B.idxStartIneq_tmp = WorkingSet->nWConstr[0] +
        WorkingSet->nWConstr[1];
      nmpcSine_B.idxStartIneq_f = nmpcSine_B.idxStartIneq_tmp + 1;
      nmpcSine_B.idx_Aineq = WorkingSet->nActiveConstr;
      for (nmpcSine_B.idxIneqOffset = nmpcSine_B.idxStartIneq_f;
           nmpcSine_B.idxIneqOffset <= nmpcSine_B.idx_Aineq;
           nmpcSine_B.idxIneqOffset++) {
        WorkingSet->isActiveConstr.data[(WorkingSet->isActiveIdx
          [WorkingSet->Wid.data[nmpcSine_B.idxIneqOffset - 1] - 1] +
          WorkingSet->Wlocalidx.data[nmpcSine_B.idxIneqOffset - 1]) - 2] = false;
      }

      WorkingSet->nWConstr[2] = 0;
      WorkingSet->nWConstr[3] = 0;
      WorkingSet->nWConstr[4] = 0;
      WorkingSet->nActiveConstr = nmpcSine_B.idxStartIneq_tmp;
      for (nmpcSine_B.idxIneqOffset = 0; nmpcSine_B.idxIneqOffset <
           nmpcSine_B.nWIneq_old; nmpcSine_B.idxIneqOffset++) {
        nmpcSine_addAineqConstr(WorkingSet, TrialState->
          workingset_old.data[nmpcSine_B.idxIneqOffset]);
      }

      for (nmpcSine_B.nWIneq_old = 0; nmpcSine_B.nWIneq_old <
           nmpcSine_B.nWLower_old; nmpcSine_B.nWIneq_old++) {
        nmpc_addBoundToActiveSetMatrix_(WorkingSet, 4,
          TrialState->workingset_old.data[nmpcSine_B.nWIneq_old +
          nmpcSine_B.nVar_j]);
      }

      for (nmpcSine_B.nWLower_old = 0; nmpcSine_B.nWLower_old <
           nmpcSine_B.nWUpper_old; nmpcSine_B.nWLower_old++) {
        nmpc_addBoundToActiveSetMatrix_(WorkingSet, 5,
          TrialState->workingset_old.data[(nmpcSine_B.nWLower_old +
          nmpcSine_B.nVar_j) + WorkingSet->sizes[3]]);
      }
    }
  }

  if (!success) {
    if (nmpcSine_B.mConstrMax_o - 1 >= 0) {
      memcpy(&TrialState->lambda.data[0], &TrialState->lambdaStopTest.data[0],
             static_cast<uint32_T>(nmpcSine_B.mConstrMax_o) * sizeof(real_T));
    }
  } else {
    nmpcSine_sortLambdaQP(TrialState->lambda.data, WorkingSet->nActiveConstr,
                          WorkingSet->sizes, WorkingSet->isActiveIdx,
                          WorkingSet->Wid.data, WorkingSet->Wlocalidx.data,
                          memspace->workspace_float.data);
  }

  return success;
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_maxConstraintViolation(const
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *obj, const real_T x_data[])
{
  real_T u1;
  real_T v;
  int32_T b;
  int32_T idx;
  if (obj->probType == 2) {
    nmpcSine_B.b_obj = *obj;
    v = maxConstraintViolation_AMats_re(&nmpcSine_B.b_obj, x_data);
  } else {
    nmpcSine_B.b_obj = *obj;
    v = maxConstraintViolation_AMats_no(&nmpcSine_B.b_obj, x_data);
  }

  if (obj->sizes[3] > 0) {
    b = static_cast<uint16_T>(obj->sizes[3]);
    for (idx = 0; idx < b; idx++) {
      u1 = -x_data[nmpcSine_B.b_obj.indexLB.data[idx] - 1] -
        nmpcSine_B.b_obj.lb.data[nmpcSine_B.b_obj.indexLB.data[idx] - 1];
      if ((!(v >= u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  if (obj->sizes[4] > 0) {
    b = static_cast<uint16_T>(obj->sizes[4]);
    for (idx = 0; idx < b; idx++) {
      u1 = x_data[nmpcSine_B.b_obj.indexUB.data[idx] - 1] -
        nmpcSine_B.b_obj.ub.data[nmpcSine_B.b_obj.indexUB.data[idx] - 1];
      if ((!(v >= u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  if (obj->sizes[0] > 0) {
    b = static_cast<uint16_T>(obj->sizes[0]);
    for (idx = 0; idx < b; idx++) {
      u1 = fabs(x_data[nmpcSine_B.b_obj.indexFixed.data[idx] - 1] -
                nmpcSine_B.b_obj.ub.data[nmpcSine_B.b_obj.indexFixed.data[idx] -
                1]);
      if ((!(v >= u1)) && (!rtIsNaN(u1))) {
        v = u1;
      }
    }
  }

  return v;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_normal(const real_T Hessian[15625], const real_T
  grad_data[], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
  sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction,
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *WorkingSet, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager,
  s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *
  QPObjective, const somzaGboVhDG7PNQS6E98jD_nmpcS_T *qpoptions,
  s7RdrPWkr8UPAUyTdDJkLaG_nmpcS_T *stepFlags)
{
  nmpcSine_driver_o(Hessian, grad_data, TrialState, memspace, WorkingSet,
                    QRManager, CholManager, QPObjective, qpoptions, qpoptions);
  nmpcSine_B.isEqAndIneqFeasible = (nmpcSine_maxConstraintViolation(WorkingSet,
    TrialState->xstar.data) <= 0.001);
  if ((TrialState->state > 0) || ((TrialState->state == 0) &&
       nmpcSine_B.isEqAndIneqFeasible)) {
    nmpcSine_B.penaltyParamTrial = MeritFunction->penaltyParam;
    nmpcSine_B.constrViolationEq = 0.0;
    for (nmpcSine_B.k_g = 0; nmpcSine_B.k_g < 120; nmpcSine_B.k_g++) {
      nmpcSine_B.constrViolationEq += fabs(TrialState->cEq[nmpcSine_B.k_g]);
    }

    nmpcSine_B.constrViolationIneq = 0.0;
    nmpcSine_B.b_m = static_cast<uint16_T>(WorkingSet->sizes[2]);
    for (nmpcSine_B.k_g = 0; nmpcSine_B.k_g < nmpcSine_B.b_m; nmpcSine_B.k_g++)
    {
      nmpcSine_B.TrialState_cIneq = TrialState->cIneq.data[nmpcSine_B.k_g];
      if (nmpcSine_B.TrialState_cIneq > 0.0) {
        nmpcSine_B.constrViolationIneq += nmpcSine_B.TrialState_cIneq;
      }
    }

    nmpcSine_B.constrViolationEq += nmpcSine_B.constrViolationIneq;
    nmpcSine_B.constrViolationIneq = MeritFunction->linearizedConstrViol;
    MeritFunction->linearizedConstrViol = 0.0;
    nmpcSine_B.constrViolationIneq += nmpcSine_B.constrViolationEq;
    if ((nmpcSine_B.constrViolationIneq > 2.2204460492503131E-16) &&
        (TrialState->fstar > 0.0)) {
      if (TrialState->sqpFval == 0.0) {
        nmpcSine_B.penaltyParamTrial = 1.0;
      } else {
        nmpcSine_B.penaltyParamTrial = 1.5;
      }

      nmpcSine_B.penaltyParamTrial = nmpcSine_B.penaltyParamTrial *
        TrialState->fstar / nmpcSine_B.constrViolationIneq;
    }

    if (nmpcSine_B.penaltyParamTrial < MeritFunction->penaltyParam) {
      MeritFunction->phi = nmpcSine_B.penaltyParamTrial *
        nmpcSine_B.constrViolationEq + TrialState->sqpFval;
      if (((MeritFunction->initConstrViolationEq +
            MeritFunction->initConstrViolationIneq) *
           nmpcSine_B.penaltyParamTrial + MeritFunction->initFval) -
          MeritFunction->phi > static_cast<real_T>
          (MeritFunction->nPenaltyDecreases) * MeritFunction->threshold) {
        MeritFunction->nPenaltyDecreases++;
        if ((MeritFunction->nPenaltyDecreases << 1) > TrialState->sqpIterations)
        {
          MeritFunction->threshold *= 10.0;
        }

        if (nmpcSine_B.penaltyParamTrial >= 1.0E-10) {
          MeritFunction->penaltyParam = nmpcSine_B.penaltyParamTrial;
        } else {
          MeritFunction->penaltyParam = 1.0E-10;
        }
      } else {
        MeritFunction->phi = MeritFunction->penaltyParam *
          nmpcSine_B.constrViolationEq + TrialState->sqpFval;
      }
    } else {
      if (nmpcSine_B.penaltyParamTrial >= 1.0E-10) {
        MeritFunction->penaltyParam = nmpcSine_B.penaltyParamTrial;
      } else {
        MeritFunction->penaltyParam = 1.0E-10;
      }

      MeritFunction->phi = MeritFunction->penaltyParam *
        nmpcSine_B.constrViolationEq + TrialState->sqpFval;
    }

    nmpcSine_B.penaltyParamTrial = TrialState->fstar -
      MeritFunction->penaltyParam * nmpcSine_B.constrViolationEq;
    if (nmpcSine_B.penaltyParamTrial <= 0.0) {
      MeritFunction->phiPrimePlus = nmpcSine_B.penaltyParamTrial;
    } else {
      MeritFunction->phiPrimePlus = 0.0;
    }
  } else if (TrialState->state != -6) {
    stepFlags->stepType = 2;
  }

  nmpcSine_sortLambdaQP(TrialState->lambda.data, WorkingSet->nActiveConstr,
                        WorkingSet->sizes, WorkingSet->isActiveIdx,
                        WorkingSet->Wid.data, WorkingSet->Wlocalidx.data,
                        memspace->workspace_float.data);
  nmpcSine_B.isEqAndIneqFeasible = (WorkingSet->mEqRemoved > 0);
  while ((WorkingSet->mEqRemoved > 0) && (WorkingSet->indexEqRemoved
          [WorkingSet->mEqRemoved - 1] >= 1)) {
    nmpcSine_addAeqConstr(WorkingSet, WorkingSet->indexEqRemoved
                          [WorkingSet->mEqRemoved - 1]);
    WorkingSet->mEqRemoved--;
  }

  if (nmpcSine_B.isEqAndIneqFeasible) {
    for (nmpcSine_B.k_g = 0; nmpcSine_B.k_g < 120; nmpcSine_B.k_g++) {
      WorkingSet->Wlocalidx.data[WorkingSet->sizes[0] + nmpcSine_B.k_g] =
        nmpcSine_B.k_g + 1;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_relaxed(const real_T Hessian[15625], const real_T
  grad_data[], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
  sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction,
  s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T
  *WorkingSet, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager,
  s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *
  QPObjective, somzaGboVhDG7PNQS6E98jD_nmpcS_T *qpoptions)
{
  nmpcSine_B.nVarOrig = WorkingSet->nVar - 1;
  nmpcSine_B.beta = 0.0;
  nmpcSine_B.idx_max = static_cast<uint16_T>(WorkingSet->nVar);
  for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig <
       nmpcSine_B.idx_max; nmpcSine_B.mFiniteLBOrig++) {
    nmpcSine_B.beta += Hessian[125 * nmpcSine_B.mFiniteLBOrig +
      nmpcSine_B.mFiniteLBOrig];
  }

  nmpcSine_B.beta /= static_cast<real_T>(WorkingSet->nVar);
  if (TrialState->sqpIterations <= 1) {
    nmpcSine_B.mIneq_n = QPObjective->nvar;
    if (QPObjective->nvar < 1) {
      nmpcSine_B.idx_max = 0;
    } else {
      nmpcSine_B.idx_max = 1;
      if (QPObjective->nvar > 1) {
        nmpcSine_B.smax_f = fabs(grad_data[0]);
        for (nmpcSine_B.mFiniteLBOrig = 2; nmpcSine_B.mFiniteLBOrig <=
             nmpcSine_B.mIneq_n; nmpcSine_B.mFiniteLBOrig++) {
          nmpcSine_B.s_a = fabs(grad_data[nmpcSine_B.mFiniteLBOrig - 1]);
          if (nmpcSine_B.s_a > nmpcSine_B.smax_f) {
            nmpcSine_B.idx_max = nmpcSine_B.mFiniteLBOrig;
            nmpcSine_B.smax_f = nmpcSine_B.s_a;
          }
        }
      }
    }

    nmpcSine_B.smax_f = fabs(grad_data[nmpcSine_B.idx_max - 1]);
    if ((nmpcSine_B.smax_f <= 1.0) || rtIsNaN(nmpcSine_B.smax_f)) {
      nmpcSine_B.smax_f = 1.0;
    }

    nmpcSine_B.smax_f *= 100.0;
  } else {
    nmpcSine_B.mIneq_n = WorkingSet->mConstr;
    if (WorkingSet->mConstr < 1) {
      nmpcSine_B.idx_max = 0;
    } else {
      nmpcSine_B.idx_max = 1;
      if (WorkingSet->mConstr > 1) {
        nmpcSine_B.smax_f = fabs(TrialState->lambdasqp.data[0]);
        for (nmpcSine_B.mFiniteLBOrig = 2; nmpcSine_B.mFiniteLBOrig <=
             nmpcSine_B.mIneq_n; nmpcSine_B.mFiniteLBOrig++) {
          nmpcSine_B.s_a = fabs(TrialState->
                                lambdasqp.data[nmpcSine_B.mFiniteLBOrig - 1]);
          if (nmpcSine_B.s_a > nmpcSine_B.smax_f) {
            nmpcSine_B.idx_max = nmpcSine_B.mFiniteLBOrig;
            nmpcSine_B.smax_f = nmpcSine_B.s_a;
          }
        }
      }
    }

    nmpcSine_B.smax_f = fabs(TrialState->lambdasqp.data[nmpcSine_B.idx_max - 1]);
  }

  QPObjective->nvar = WorkingSet->nVar;
  QPObjective->beta = nmpcSine_B.beta;
  QPObjective->rho = nmpcSine_B.smax_f;
  QPObjective->hasLinear = true;
  QPObjective->objtype = 4;
  nmpcSine_B.c_WorkingSet = *WorkingSet;
  nmpcSine_setProblemType(&nmpcSine_B.c_WorkingSet, 2);
  nmpcSine_B.mIneq_n = nmpcSine_B.c_WorkingSet.sizes[2] + 1;
  nmpcSine_B.mLBOrig = (nmpcSine_B.c_WorkingSet.sizes[3] -
                        nmpcSine_B.c_WorkingSet.sizes[2]) - 240;
  nmpcSine_B.idx_max = static_cast<uint16_T>(nmpcSine_B.c_WorkingSet.sizes[2]);
  if (nmpcSine_B.idx_max - 1 >= 0) {
    memcpy(&memspace->workspace_float.data[0],
           &nmpcSine_B.c_WorkingSet.bineq.data[0], static_cast<uint32_T>
           (nmpcSine_B.idx_max) * sizeof(real_T));
  }

  nmpcSine_xgemv_dffoyl(WorkingSet->nVar, nmpcSine_B.c_WorkingSet.sizes[2],
                        nmpcSine_B.c_WorkingSet.Aineq.data,
                        nmpcSine_B.c_WorkingSet.ldA, TrialState->xstar.data,
                        memspace->workspace_float.data);
  for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig <
       nmpcSine_B.idx_max; nmpcSine_B.mFiniteLBOrig++) {
    TrialState->xstar.data[(nmpcSine_B.nVarOrig + nmpcSine_B.mFiniteLBOrig) + 1]
      = static_cast<real_T>(memspace->
      workspace_float.data[nmpcSine_B.mFiniteLBOrig] > 0.0) *
      memspace->workspace_float.data[nmpcSine_B.mFiniteLBOrig];
  }

  memcpy(&memspace->workspace_float.data[0], &nmpcSine_B.c_WorkingSet.beq[0],
         120U * sizeof(real_T));
  nmpcSine_xgemv_dffoyl(WorkingSet->nVar, 120, nmpcSine_B.c_WorkingSet.Aeq.data,
                        nmpcSine_B.c_WorkingSet.ldA, TrialState->xstar.data,
                        memspace->workspace_float.data);
  for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig < 120;
       nmpcSine_B.mFiniteLBOrig++) {
    nmpcSine_B.idx_max = nmpcSine_B.mIneq_n + nmpcSine_B.mFiniteLBOrig;
    if (memspace->workspace_float.data[nmpcSine_B.mFiniteLBOrig] <= 0.0) {
      TrialState->xstar.data[nmpcSine_B.nVarOrig + nmpcSine_B.idx_max] = 0.0;
      TrialState->xstar.data[(nmpcSine_B.nVarOrig + nmpcSine_B.idx_max) + 120] =
        -memspace->workspace_float.data[nmpcSine_B.mFiniteLBOrig];
      nmpc_addBoundToActiveSetMatrix_(&nmpcSine_B.c_WorkingSet, 4,
        nmpcSine_B.mLBOrig + nmpcSine_B.idx_max);
      if (memspace->workspace_float.data[nmpcSine_B.mFiniteLBOrig] >= -0.001) {
        nmpc_addBoundToActiveSetMatrix_(&nmpcSine_B.c_WorkingSet, 4,
          (nmpcSine_B.mLBOrig + nmpcSine_B.idx_max) + 120);
      }
    } else {
      nmpcSine_B.mIneq_tmp = nmpcSine_B.nVarOrig + nmpcSine_B.idx_max;
      TrialState->xstar.data[nmpcSine_B.mIneq_tmp] =
        memspace->workspace_float.data[nmpcSine_B.mFiniteLBOrig];
      TrialState->xstar.data[nmpcSine_B.mIneq_tmp + 120] = 0.0;
      nmpc_addBoundToActiveSetMatrix_(&nmpcSine_B.c_WorkingSet, 4,
        (nmpcSine_B.mLBOrig + nmpcSine_B.idx_max) + 120);
      if (memspace->workspace_float.data[nmpcSine_B.mFiniteLBOrig] <= 0.001) {
        nmpc_addBoundToActiveSetMatrix_(&nmpcSine_B.c_WorkingSet, 4,
          nmpcSine_B.mLBOrig + nmpcSine_B.idx_max);
      }
    }
  }

  nmpcSine_B.nVarOrig = qpoptions->MaxIterations;
  qpoptions->MaxIterations = (qpoptions->MaxIterations +
    nmpcSine_B.c_WorkingSet.nVar) - WorkingSet->nVar;
  nmpcSine_driver_o(Hessian, grad_data, TrialState, memspace,
                    &nmpcSine_B.c_WorkingSet, QRManager, CholManager,
                    QPObjective, qpoptions, qpoptions);
  qpoptions->MaxIterations = nmpcSine_B.nVarOrig;
  nmpcSine_B.mIneq_n = nmpcSine_B.c_WorkingSet.sizes[3] - 241;
  nmpcSine_B.nVarOrig = 0;
  for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig < 120;
       nmpcSine_B.mFiniteLBOrig++) {
    nmpcSine_B.idx_max = (nmpcSine_B.c_WorkingSet.isActiveIdx[3] +
                          nmpcSine_B.mIneq_n) + nmpcSine_B.mFiniteLBOrig;
    nmpcSine_B.tf =
      nmpcSine_B.c_WorkingSet.isActiveConstr.data[nmpcSine_B.idx_max];
    nmpcSine_B.b_tf =
      nmpcSine_B.c_WorkingSet.isActiveConstr.data[nmpcSine_B.idx_max + 120];
    memspace->workspace_int.data[nmpcSine_B.mFiniteLBOrig] = nmpcSine_B.tf;
    memspace->workspace_int.data[nmpcSine_B.mFiniteLBOrig + 120] =
      nmpcSine_B.b_tf;
    nmpcSine_B.nVarOrig = (nmpcSine_B.nVarOrig + nmpcSine_B.tf) +
      nmpcSine_B.b_tf;
  }

  nmpcSine_B.mLBOrig = static_cast<uint16_T>(nmpcSine_B.c_WorkingSet.sizes[2]);
  for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig <
       nmpcSine_B.mLBOrig; nmpcSine_B.mFiniteLBOrig++) {
    nmpcSine_B.idx_max = nmpcSine_B.c_WorkingSet.isActiveConstr.data
      [((nmpcSine_B.c_WorkingSet.isActiveIdx[3] + nmpcSine_B.mIneq_n) -
        nmpcSine_B.c_WorkingSet.sizes[2]) + nmpcSine_B.mFiniteLBOrig];
    memspace->workspace_int.data[nmpcSine_B.mFiniteLBOrig + 240] =
      nmpcSine_B.idx_max;
    nmpcSine_B.nVarOrig += nmpcSine_B.idx_max;
  }

  if (TrialState->state != -6) {
    nmpcSine_B.idx_max = (nmpcSine_B.c_WorkingSet.nVarMax - WorkingSet->nVar) -
      1;
    nmpcSine_B.mIneq_tmp = WorkingSet->nVar + 1;
    nmpcSine_B.s_a = 0.0;
    nmpcSine_B.qpfvalQuadExcess = 0.0;
    if (nmpcSine_B.idx_max >= 1) {
      nmpcSine_B.mLBOrig = WorkingSet->nVar + nmpcSine_B.idx_max;
      for (nmpcSine_B.mFiniteLBOrig = nmpcSine_B.mIneq_tmp;
           nmpcSine_B.mFiniteLBOrig <= nmpcSine_B.mLBOrig;
           nmpcSine_B.mFiniteLBOrig++) {
        nmpcSine_B.s_a += fabs(TrialState->xstar.data[nmpcSine_B.mFiniteLBOrig -
          1]);
      }

      nmpcSine_B.idx_max = static_cast<uint16_T>(nmpcSine_B.idx_max);
      for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig <
           nmpcSine_B.idx_max; nmpcSine_B.mFiniteLBOrig++) {
        nmpcSine_B.qpfvalQuadExcess_tmp = TrialState->xstar.data
          [WorkingSet->nVar + nmpcSine_B.mFiniteLBOrig];
        nmpcSine_B.qpfvalQuadExcess += nmpcSine_B.qpfvalQuadExcess_tmp *
          nmpcSine_B.qpfvalQuadExcess_tmp;
      }
    }

    nmpcSine_B.beta = (TrialState->fstar - nmpcSine_B.smax_f * nmpcSine_B.s_a) -
      nmpcSine_B.beta / 2.0 * nmpcSine_B.qpfvalQuadExcess;
    nmpcSine_B.mIneq_n = (WorkingSet->nVarMax - WorkingSet->nVar) - 1;
    nmpcSine_B.smax_f = MeritFunction->penaltyParam;
    nmpcSine_B.s_a = 0.0;
    for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig < 120;
         nmpcSine_B.mFiniteLBOrig++) {
      nmpcSine_B.s_a += fabs(TrialState->cEq[nmpcSine_B.mFiniteLBOrig]);
    }

    nmpcSine_B.qpfvalQuadExcess = 0.0;
    nmpcSine_B.mLBOrig = static_cast<uint16_T>(WorkingSet->sizes[2]);
    for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig <
         nmpcSine_B.mLBOrig; nmpcSine_B.mFiniteLBOrig++) {
      nmpcSine_B.qpfvalQuadExcess_tmp = TrialState->
        cIneq.data[nmpcSine_B.mFiniteLBOrig];
      if (nmpcSine_B.qpfvalQuadExcess_tmp > 0.0) {
        nmpcSine_B.qpfvalQuadExcess += nmpcSine_B.qpfvalQuadExcess_tmp;
      }
    }

    nmpcSine_B.s_a += nmpcSine_B.qpfvalQuadExcess;
    nmpcSine_B.qpfvalQuadExcess = MeritFunction->linearizedConstrViol;
    nmpcSine_B.qpfvalQuadExcess_tmp = 0.0;
    if (nmpcSine_B.mIneq_n >= 1) {
      nmpcSine_B.mIneq_n += WorkingSet->nVar;
      for (nmpcSine_B.mFiniteLBOrig = nmpcSine_B.mIneq_tmp;
           nmpcSine_B.mFiniteLBOrig <= nmpcSine_B.mIneq_n;
           nmpcSine_B.mFiniteLBOrig++) {
        nmpcSine_B.qpfvalQuadExcess_tmp += fabs(TrialState->
          xstar.data[nmpcSine_B.mFiniteLBOrig - 1]);
      }
    }

    MeritFunction->linearizedConstrViol = nmpcSine_B.qpfvalQuadExcess_tmp;
    nmpcSine_B.qpfvalQuadExcess = (nmpcSine_B.s_a + nmpcSine_B.qpfvalQuadExcess)
      - nmpcSine_B.qpfvalQuadExcess_tmp;
    if ((nmpcSine_B.qpfvalQuadExcess > 2.2204460492503131E-16) &&
        (nmpcSine_B.beta > 0.0)) {
      if (TrialState->sqpFval == 0.0) {
        nmpcSine_B.smax_f = 1.0;
      } else {
        nmpcSine_B.smax_f = 1.5;
      }

      nmpcSine_B.smax_f = nmpcSine_B.smax_f * nmpcSine_B.beta /
        nmpcSine_B.qpfvalQuadExcess;
    }

    if (nmpcSine_B.smax_f < MeritFunction->penaltyParam) {
      MeritFunction->phi = nmpcSine_B.smax_f * nmpcSine_B.s_a +
        TrialState->sqpFval;
      if (((MeritFunction->initConstrViolationEq +
            MeritFunction->initConstrViolationIneq) * nmpcSine_B.smax_f +
           MeritFunction->initFval) - MeritFunction->phi > static_cast<real_T>
          (MeritFunction->nPenaltyDecreases) * MeritFunction->threshold) {
        MeritFunction->nPenaltyDecreases++;
        if ((MeritFunction->nPenaltyDecreases << 1) > TrialState->sqpIterations)
        {
          MeritFunction->threshold *= 10.0;
        }

        if (nmpcSine_B.smax_f >= 1.0E-10) {
          MeritFunction->penaltyParam = nmpcSine_B.smax_f;
        } else {
          MeritFunction->penaltyParam = 1.0E-10;
        }
      } else {
        MeritFunction->phi = MeritFunction->penaltyParam * nmpcSine_B.s_a +
          TrialState->sqpFval;
      }
    } else {
      if (nmpcSine_B.smax_f >= 1.0E-10) {
        MeritFunction->penaltyParam = nmpcSine_B.smax_f;
      } else {
        MeritFunction->penaltyParam = 1.0E-10;
      }

      MeritFunction->phi = MeritFunction->penaltyParam * nmpcSine_B.s_a +
        TrialState->sqpFval;
    }

    nmpcSine_B.beta -= MeritFunction->penaltyParam * nmpcSine_B.s_a;
    if (nmpcSine_B.beta <= 0.0) {
      MeritFunction->phiPrimePlus = nmpcSine_B.beta;
    } else {
      MeritFunction->phiPrimePlus = 0.0;
    }

    nmpcSine_B.mIneq_n = nmpcSine_B.c_WorkingSet.isActiveIdx[1] - 1;
    for (nmpcSine_B.mFiniteLBOrig = 0; nmpcSine_B.mFiniteLBOrig < 120;
         nmpcSine_B.mFiniteLBOrig++) {
      if (memspace->workspace_int.data[nmpcSine_B.mFiniteLBOrig] != 0) {
        nmpcSine_B.tf = (memspace->workspace_int.data[nmpcSine_B.mFiniteLBOrig +
                         120] != 0);
      } else {
        nmpcSine_B.tf = false;
      }

      nmpcSine_B.idx_max = nmpcSine_B.mIneq_n + nmpcSine_B.mFiniteLBOrig;
      TrialState->lambda.data[nmpcSine_B.idx_max] *= static_cast<real_T>
        (nmpcSine_B.tf);
    }

    nmpcSine_B.idx_max = nmpcSine_B.c_WorkingSet.isActiveIdx[2];
    nmpcSine_B.mIneq_n = nmpcSine_B.c_WorkingSet.nActiveConstr;
    for (nmpcSine_B.mFiniteLBOrig = nmpcSine_B.idx_max; nmpcSine_B.mFiniteLBOrig
         <= nmpcSine_B.mIneq_n; nmpcSine_B.mFiniteLBOrig++) {
      if (nmpcSine_B.c_WorkingSet.Wid.data[nmpcSine_B.mFiniteLBOrig - 1] == 3) {
        TrialState->lambda.data[nmpcSine_B.mFiniteLBOrig - 1] *=
          static_cast<real_T>(memspace->
                              workspace_int.data[nmpcSine_B.c_WorkingSet.Wlocalidx.data
                              [nmpcSine_B.mFiniteLBOrig - 1] + 239]);
      }
    }
  }

  nmpcSine_B.mFiniteLBOrig = (nmpcSine_B.c_WorkingSet.sizes[3] -
    nmpcSine_B.c_WorkingSet.sizes[2]) - 240;
  nmpcSine_B.idx_max = nmpcSine_B.c_WorkingSet.nActiveConstr;
  while ((nmpcSine_B.idx_max > nmpcSine_B.c_WorkingSet.sizes[0] + 120) &&
         (nmpcSine_B.nVarOrig > 0)) {
    if ((nmpcSine_B.c_WorkingSet.Wid.data[nmpcSine_B.idx_max - 1] == 4) &&
        (nmpcSine_B.c_WorkingSet.Wlocalidx.data[nmpcSine_B.idx_max - 1] >
         nmpcSine_B.mFiniteLBOrig)) {
      nmpcSine_B.beta = TrialState->
        lambda.data[nmpcSine_B.c_WorkingSet.nActiveConstr - 1];
      TrialState->lambda.data[nmpcSine_B.c_WorkingSet.nActiveConstr - 1] = 0.0;
      TrialState->lambda.data[nmpcSine_B.idx_max - 1] = nmpcSine_B.beta;
      nmpcSine_removeConstr(&nmpcSine_B.c_WorkingSet, nmpcSine_B.idx_max);
      nmpcSine_B.nVarOrig--;
    }

    nmpcSine_B.idx_max--;
  }

  QPObjective->nvar = WorkingSet->nVar;
  QPObjective->hasLinear = true;
  QPObjective->objtype = 3;
  *WorkingSet = nmpcSine_B.c_WorkingSet;
  nmpcSine_setProblemType(WorkingSet, 3);
  nmpcSine_sortLambdaQP(TrialState->lambda.data, WorkingSet->nActiveConstr,
                        WorkingSet->sizes, WorkingSet->isActiveIdx,
                        WorkingSet->Wid.data, WorkingSet->Wlocalidx.data,
                        memspace->workspace_float.data);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_step_k(s7RdrPWkr8UPAUyTdDJkLaG_nmpcS_T *stepFlags,
  real_T Hessian[15625], const real_T lb[125], const real_T ub[125],
  s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState, sG8JZ69axY52WWR6RKyApQC_nmpcS_T
  *MeritFunction, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T *memspace,
  s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet, s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T
  *QRManager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T *CholManager,
  s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *QPObjective, somzaGboVhDG7PNQS6E98jD_nmpcS_T *
  qpoptions)
{
  __m128d tmp;
  __m128d tmp_0;
  int32_T exitg1;
  boolean_T guard1;
  stepFlags->stepAccepted = true;
  nmpcSine_B.checkBoundViolation = true;
  nmpcSine_B.nVar = WorkingSet->nVar - 1;
  if (stepFlags->stepType != 3) {
    nmpcSine_B.idxStartIneq = static_cast<uint16_T>(WorkingSet->nVar);
    memcpy(&TrialState->xstar.data[0], &TrialState->xstarsqp[0],
           static_cast<uint32_T>(nmpcSine_B.idxStartIneq) * sizeof(real_T));
  } else if (nmpcSine_B.nVar >= 0) {
    memcpy(&TrialState->searchDir.data[0], &TrialState->xstar.data[0],
           static_cast<uint32_T>(nmpcSine_B.nVar + 1) * sizeof(real_T));
  }

  do {
    exitg1 = 0;
    guard1 = false;
    switch (stepFlags->stepType) {
     case 1:
      nmpcSine_normal(Hessian, TrialState->grad.data, TrialState, MeritFunction,
                      memspace, WorkingSet, QRManager, CholManager, QPObjective,
                      qpoptions, stepFlags);
      if (stepFlags->stepType == 2) {
      } else {
        if (nmpcSine_B.nVar >= 0) {
          memcpy(&TrialState->delta_x.data[0], &TrialState->xstar.data[0],
                 static_cast<uint32_T>(nmpcSine_B.nVar + 1) * sizeof(real_T));
        }

        guard1 = true;
      }
      break;

     case 2:
      nmpcSine_B.vectorUB_f = WorkingSet->nWConstr[0] + WorkingSet->nWConstr[1];
      nmpcSine_B.idxStartIneq = nmpcSine_B.vectorUB_f + 1;
      nmpcSine_B.idxEndIneq = WorkingSet->nActiveConstr;
      for (nmpcSine_B.k_i = nmpcSine_B.idxStartIneq; nmpcSine_B.k_i <=
           nmpcSine_B.idxEndIneq; nmpcSine_B.k_i++) {
        WorkingSet->isActiveConstr.data[(WorkingSet->isActiveIdx
          [WorkingSet->Wid.data[nmpcSine_B.k_i - 1] - 1] +
          WorkingSet->Wlocalidx.data[nmpcSine_B.k_i - 1]) - 2] = false;
      }

      WorkingSet->nWConstr[2] = 0;
      WorkingSet->nWConstr[3] = 0;
      WorkingSet->nWConstr[4] = 0;
      WorkingSet->nActiveConstr = nmpcSine_B.vectorUB_f;
      nmpcSine_B.idxStartIneq = TrialState->xstar.size[0];
      nmpcSine_B.idxEndIneq = TrialState->xstar.size[0];
      if (nmpcSine_B.idxEndIneq - 1 >= 0) {
        memcpy(&nmpcSine_B.b_data[0], &TrialState->xstar.data[0],
               static_cast<uint32_T>(nmpcSine_B.idxEndIneq) * sizeof(real_T));
      }

      nmpcSine_B.idxEndIneq = static_cast<uint16_T>(WorkingSet->sizes[3]);
      for (nmpcSine_B.k_i = 0; nmpcSine_B.k_i < nmpcSine_B.idxEndIneq;
           nmpcSine_B.k_i++) {
        nmpcSine_B.nrmGradInf = WorkingSet->lb.data[WorkingSet->
          indexLB.data[nmpcSine_B.k_i] - 1];
        if (-nmpcSine_B.b_data[WorkingSet->indexLB.data[nmpcSine_B.k_i] - 1] >
            nmpcSine_B.nrmGradInf) {
          if (rtIsInf(ub[WorkingSet->indexLB.data[nmpcSine_B.k_i] - 1])) {
            nmpcSine_B.b_data[WorkingSet->indexLB.data[nmpcSine_B.k_i] - 1] =
              -nmpcSine_B.nrmGradInf + fabs(nmpcSine_B.nrmGradInf);
          } else {
            nmpcSine_B.b_data[WorkingSet->indexLB.data[nmpcSine_B.k_i] - 1] =
              (WorkingSet->ub.data[WorkingSet->indexLB.data[nmpcSine_B.k_i] - 1]
               - nmpcSine_B.nrmGradInf) / 2.0;
          }
        }
      }

      nmpcSine_B.idxEndIneq = static_cast<uint16_T>(WorkingSet->sizes[4]);
      for (nmpcSine_B.k_i = 0; nmpcSine_B.k_i < nmpcSine_B.idxEndIneq;
           nmpcSine_B.k_i++) {
        nmpcSine_B.nrmGradInf = WorkingSet->ub.data[WorkingSet->
          indexUB.data[nmpcSine_B.k_i] - 1];
        if (nmpcSine_B.b_data[WorkingSet->indexUB.data[nmpcSine_B.k_i] - 1] >
            nmpcSine_B.nrmGradInf) {
          if (rtIsInf(lb[WorkingSet->indexUB.data[nmpcSine_B.k_i] - 1])) {
            nmpcSine_B.b_data[WorkingSet->indexUB.data[nmpcSine_B.k_i] - 1] =
              nmpcSine_B.nrmGradInf - fabs(nmpcSine_B.nrmGradInf);
          } else {
            nmpcSine_B.b_data[WorkingSet->indexUB.data[nmpcSine_B.k_i] - 1] =
              (nmpcSine_B.nrmGradInf - WorkingSet->lb.data
               [WorkingSet->indexUB.data[nmpcSine_B.k_i] - 1]) / 2.0;
          }
        }
      }

      if (nmpcSine_B.idxStartIneq - 1 >= 0) {
        memcpy(&TrialState->xstar.data[0], &nmpcSine_B.b_data[0],
               static_cast<uint32_T>(nmpcSine_B.idxStartIneq) * sizeof(real_T));
      }

      nmpcSine_relaxed(Hessian, TrialState->grad.data, TrialState, MeritFunction,
                       memspace, WorkingSet, QRManager, CholManager, QPObjective,
                       qpoptions);
      if (nmpcSine_B.nVar >= 0) {
        memcpy(&TrialState->delta_x.data[0], &TrialState->xstar.data[0],
               static_cast<uint32_T>(nmpcSine_B.nVar + 1) * sizeof(real_T));
      }

      guard1 = true;
      break;

     default:
      nmpcSine_B.checkBoundViolation = nmpcSine_soc(Hessian,
        TrialState->grad.data, TrialState, memspace, WorkingSet, QRManager,
        CholManager, QPObjective, qpoptions);
      stepFlags->stepAccepted = nmpcSine_B.checkBoundViolation;
      if (stepFlags->stepAccepted && (TrialState->state != -6)) {
        nmpcSine_B.idxStartIneq = static_cast<uint16_T>(nmpcSine_B.nVar + 1);
        nmpcSine_B.idxEndIneq = (static_cast<uint16_T>(nmpcSine_B.nVar + 1) / 2)
          << 1;
        nmpcSine_B.vectorUB_f = nmpcSine_B.idxEndIneq - 2;
        for (nmpcSine_B.k_i = 0; nmpcSine_B.k_i <= nmpcSine_B.vectorUB_f;
             nmpcSine_B.k_i += 2) {
          tmp = _mm_loadu_pd(&TrialState->xstar.data[nmpcSine_B.k_i]);
          tmp_0 = _mm_loadu_pd(&TrialState->socDirection.data[nmpcSine_B.k_i]);
          _mm_storeu_pd(&TrialState->delta_x.data[nmpcSine_B.k_i], _mm_add_pd
                        (tmp, tmp_0));
        }

        for (nmpcSine_B.k_i = nmpcSine_B.idxEndIneq; nmpcSine_B.k_i <
             nmpcSine_B.idxStartIneq; nmpcSine_B.k_i++) {
          TrialState->delta_x.data[nmpcSine_B.k_i] = TrialState->
            xstar.data[nmpcSine_B.k_i] + TrialState->
            socDirection.data[nmpcSine_B.k_i];
        }
      }

      guard1 = true;
      break;
    }

    if (guard1) {
      if (TrialState->state != -6) {
        exitg1 = 1;
      } else {
        nmpcSine_B.nrmGradInf = 0.0;
        nmpcSine_B.nrmDirInf = 1.0;
        for (nmpcSine_B.k_i = 0; nmpcSine_B.k_i < 125; nmpcSine_B.k_i++) {
          nmpcSine_B.u1 = fabs(TrialState->grad.data[nmpcSine_B.k_i]);
          if ((!(nmpcSine_B.nrmGradInf >= nmpcSine_B.u1)) && (!rtIsNaN
               (nmpcSine_B.u1))) {
            nmpcSine_B.nrmGradInf = nmpcSine_B.u1;
          }

          nmpcSine_B.u1 = fabs(TrialState->xstar.data[nmpcSine_B.k_i]);
          if ((!(nmpcSine_B.nrmDirInf >= nmpcSine_B.u1)) && (!rtIsNaN
               (nmpcSine_B.u1))) {
            nmpcSine_B.nrmDirInf = nmpcSine_B.u1;
          }
        }

        nmpcSine_B.nrmGradInf /= nmpcSine_B.nrmDirInf;
        if ((nmpcSine_B.nrmGradInf <= 2.2204460492503131E-16) || rtIsNaN
            (nmpcSine_B.nrmGradInf)) {
          nmpcSine_B.nrmGradInf = 2.2204460492503131E-16;
        }

        for (nmpcSine_B.k_i = 0; nmpcSine_B.k_i < 125; nmpcSine_B.k_i++) {
          nmpcSine_B.idxEndIneq = 125 * nmpcSine_B.k_i;
          for (nmpcSine_B.idxStartIneq = 0; nmpcSine_B.idxStartIneq <
               nmpcSine_B.k_i; nmpcSine_B.idxStartIneq++) {
            Hessian[nmpcSine_B.idxEndIneq + nmpcSine_B.idxStartIneq] = 0.0;
          }

          nmpcSine_B.idxEndIneq = 125 * nmpcSine_B.k_i + nmpcSine_B.k_i;
          Hessian[nmpcSine_B.idxEndIneq] = nmpcSine_B.nrmGradInf;
          nmpcSine_B.vectorUB_f = 123 - nmpcSine_B.k_i;
          if (nmpcSine_B.vectorUB_f >= 0) {
            memset(&Hessian[nmpcSine_B.idxEndIneq + 1], 0, static_cast<uint32_T>
                   (nmpcSine_B.vectorUB_f + 1) * sizeof(real_T));
          }
        }
      }
    }
  } while (exitg1 == 0);

  if (nmpcSine_B.checkBoundViolation) {
    nmpcSine_B.idxStartIneq = TrialState->delta_x.size[0];
    nmpcSine_B.idxEndIneq = TrialState->delta_x.size[0];
    if (nmpcSine_B.idxEndIneq - 1 >= 0) {
      memcpy(&nmpcSine_B.b_data[0], &TrialState->delta_x.data[0],
             static_cast<uint32_T>(nmpcSine_B.idxEndIneq) * sizeof(real_T));
    }

    nmpcSine_B.k_i = static_cast<uint16_T>(WorkingSet->sizes[3]);
    for (nmpcSine_B.nVar = 0; nmpcSine_B.nVar < nmpcSine_B.k_i; nmpcSine_B.nVar
         ++) {
      nmpcSine_B.nrmDirInf = nmpcSine_B.b_data[WorkingSet->
        indexLB.data[nmpcSine_B.nVar] - 1];
      nmpcSine_B.nrmGradInf = (TrialState->xstarsqp[WorkingSet->
        indexLB.data[nmpcSine_B.nVar] - 1] + nmpcSine_B.nrmDirInf) -
        lb[WorkingSet->indexLB.data[nmpcSine_B.nVar] - 1];
      if (nmpcSine_B.nrmGradInf < 0.0) {
        _mm_storeu_pd(&nmpcSine_B.dv8[0], _mm_sub_pd(_mm_set_pd
          (TrialState->xstar.data[WorkingSet->indexLB.data[nmpcSine_B.nVar] - 1],
           nmpcSine_B.nrmDirInf), _mm_set1_pd(nmpcSine_B.nrmGradInf)));
        nmpcSine_B.b_data[WorkingSet->indexLB.data[nmpcSine_B.nVar] - 1] =
          nmpcSine_B.dv8[0];
        TrialState->xstar.data[WorkingSet->indexLB.data[nmpcSine_B.nVar] - 1] =
          nmpcSine_B.dv8[1];
      }
    }

    nmpcSine_B.k_i = static_cast<uint16_T>(WorkingSet->sizes[4]);
    for (nmpcSine_B.nVar = 0; nmpcSine_B.nVar < nmpcSine_B.k_i; nmpcSine_B.nVar
         ++) {
      nmpcSine_B.nrmDirInf = nmpcSine_B.b_data[WorkingSet->
        indexUB.data[nmpcSine_B.nVar] - 1];
      nmpcSine_B.nrmGradInf = (ub[WorkingSet->indexUB.data[nmpcSine_B.nVar] - 1]
        - TrialState->xstarsqp[WorkingSet->indexUB.data[nmpcSine_B.nVar] - 1]) -
        nmpcSine_B.nrmDirInf;
      if (nmpcSine_B.nrmGradInf < 0.0) {
        _mm_storeu_pd(&nmpcSine_B.dv8[0], _mm_add_pd(_mm_set_pd
          (TrialState->xstar.data[WorkingSet->indexUB.data[nmpcSine_B.nVar] - 1],
           nmpcSine_B.nrmDirInf), _mm_set1_pd(nmpcSine_B.nrmGradInf)));
        nmpcSine_B.b_data[WorkingSet->indexUB.data[nmpcSine_B.nVar] - 1] =
          nmpcSine_B.dv8[0];
        TrialState->xstar.data[WorkingSet->indexUB.data[nmpcSine_B.nVar] - 1] =
          nmpcSine_B.dv8[1];
      }
    }

    if (nmpcSine_B.idxStartIneq - 1 >= 0) {
      memcpy(&TrialState->delta_x.data[0], &nmpcSine_B.b_data[0],
             static_cast<uint32_T>(nmpcSine_B.idxStartIneq) * sizeof(real_T));
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_evalObjAndConstr(int32_T obj_next_next_next_next_next_b_,
  const s_jex761Cl1dvQqVqRqjms8C_nmpc_T *obj_next_next_next_next_next_ne, const
  s_I4XPpWw7d7shktLagLlNtD_nmpc_T *obj_next_next_next_next_next__0, const real_T
  x[125], real_T Cineq_workspace_data[], int32_T ineq0, real_T Ceq_workspace[120],
  real_T *fval, int32_T *status)
{
  __m128i tmp;
  boolean_T exitg1;
  nmpcSine_getXUe(x, obj_next_next_next_next_next__0->runtimedata.x,
                  nmpcSine_B.X_g, nmpcSine_B.U_n, &nmpcSine_B.e_p);
  *fval = nmpcSine_costFcn_d(nmpcSine_B.X_g, nmpcSine_B.U_n,
    obj_next_next_next_next_next__0->userdata.References,
    obj_next_next_next_next_next__0->runtimedata.Parameters[0],
    obj_next_next_next_next_next__0->runtimedata.Parameters[1],
    obj_next_next_next_next_next__0->runtimedata.Parameters[2]);
  *status = 1;
  nmpcSine_B.y_fn = rtIsNaN(*fval);
  if (rtIsInf(*fval) || nmpcSine_B.y_fn) {
    if (nmpcSine_B.y_fn) {
      *status = -3;
    } else if (*fval < 0.0) {
      *status = -1;
    } else {
      *status = -2;
    }
  }

  if (*status == 1) {
    if (obj_next_next_next_next_next_b_ - 1 < 0) {
      nmpcSine_B.n = 0;
    } else {
      nmpcSine_B.n = static_cast<uint8_T>(obj_next_next_next_next_next_b_ - 1) +
        1;
    }

    nmpcSine_B.ineqRange_size_idx_1 = nmpcSine_B.n;
    if (nmpcSine_B.n > 0) {
      nmpcSine_B.ineqRange_data[0] = 0;
      nmpcSine_B.yk = 0;
      for (nmpcSine_B.k_a = 2; nmpcSine_B.k_a <= nmpcSine_B.n; nmpcSine_B.k_a++)
      {
        nmpcSine_B.yk++;
        nmpcSine_B.ineqRange_data[nmpcSine_B.k_a - 1] = nmpcSine_B.yk;
      }
    }

    nmpcSine_B.yk = nmpcSine_B.n - 1;
    nmpcSine_B.n = (nmpcSine_B.n / 4) << 2;
    nmpcSine_B.vectorUB_o = nmpcSine_B.n - 4;
    for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a <= nmpcSine_B.vectorUB_o;
         nmpcSine_B.k_a += 4) {
      tmp = _mm_loadu_si128((const __m128i *)
                            &nmpcSine_B.ineqRange_data[nmpcSine_B.k_a]);
      _mm_storeu_si128((__m128i *)&nmpcSine_B.ineqRange_data[nmpcSine_B.k_a],
                       _mm_add_epi32(tmp, _mm_set1_epi32(ineq0)));
    }

    for (nmpcSine_B.k_a = nmpcSine_B.n; nmpcSine_B.k_a <= nmpcSine_B.yk;
         nmpcSine_B.k_a++) {
      nmpcSine_B.ineqRange_data[nmpcSine_B.k_a] += ineq0;
    }

    nmpcSine_getXUe(x, obj_next_next_next_next_next_ne->x, nmpcSine_B.X_g,
                    nmpcSine_B.U_n, &nmpcSine_B.e_p);
    memset(&nmpcSine_B.reshapes_f1[0], 0, 120U * sizeof(real_T));
    for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < 6; nmpcSine_B.k_a++) {
      nmpcSine_B.ic_a[nmpcSine_B.k_a] = static_cast<real_T>(nmpcSine_B.k_a) +
        1.0;
    }

    for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < 21; nmpcSine_B.k_a++) {
      nmpcSine_B.yk = nmpcSine_B.k_a << 2;
      nmpcSine_B.b_U_b[nmpcSine_B.yk] = nmpcSine_B.U_n[nmpcSine_B.k_a];
      nmpcSine_B.b_U_b[nmpcSine_B.yk + 1] = nmpcSine_B.U_n[nmpcSine_B.k_a + 21];
      nmpcSine_B.b_U_b[nmpcSine_B.yk + 2] = nmpcSine_B.U_n[nmpcSine_B.k_a + 42];
      nmpcSine_B.b_U_b[nmpcSine_B.yk + 3] = nmpcSine_B.U_n[nmpcSine_B.k_a + 63];
    }

    for (nmpcSine_B.n = 0; nmpcSine_B.n < 6; nmpcSine_B.n++) {
      for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < 21; nmpcSine_B.k_a++) {
        nmpcSine_B.b_X_l[nmpcSine_B.n + 6 * nmpcSine_B.k_a] = nmpcSine_B.X_g[21 *
          nmpcSine_B.n + nmpcSine_B.k_a];
      }
    }

    for (nmpcSine_B.n = 0; nmpcSine_B.n < 20; nmpcSine_B.n++) {
      nmpcSine_B.b_ic_idx_0 = nmpcSine_B.b_X_l[6 * nmpcSine_B.n + 2];
      nmpcSine_B.t2_p = cos(nmpcSine_B.b_ic_idx_0);
      nmpcSine_B.t3_p = sin(nmpcSine_B.b_ic_idx_0);
      nmpcSine_B.vectorUB_o = (nmpcSine_B.n + 1) * 6;
      nmpcSine_B.b_ic_idx_0 = nmpcSine_B.b_X_l[nmpcSine_B.vectorUB_o + 2];
      nmpcSine_B.b_t2_a = cos(nmpcSine_B.b_ic_idx_0);
      nmpcSine_B.b_t3_j = sin(nmpcSine_B.b_ic_idx_0);
      nmpcSine_B.b_ic_idx_0 = nmpcSine_B.b_X_l[6 * nmpcSine_B.n + 4];
      nmpcSine_B.b_X_ek = nmpcSine_B.b_X_l[6 * nmpcSine_B.n + 3];
      nmpcSine_B.b_X_e[0] = -nmpcSine_B.b_ic_idx_0 * nmpcSine_B.t3_p +
        nmpcSine_B.t2_p * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_X_e[1] = nmpcSine_B.b_ic_idx_0 * nmpcSine_B.t2_p +
        nmpcSine_B.t3_p * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_X_o = nmpcSine_B.b_X_l[6 * nmpcSine_B.n + 5];
      nmpcSine_B.b_X_e[2] = nmpcSine_B.b_X_o;
      nmpcSine_B.yk = nmpcSine_B.n << 2;
      nmpcSine_B.b_U_bb = nmpcSine_B.b_U_b[nmpcSine_B.yk + 2];
      nmpcSine_B.b_U_a = nmpcSine_B.b_U_b[nmpcSine_B.yk + 3];
      nmpcSine_B.b_U_g = nmpcSine_B.b_U_b[nmpcSine_B.yk];
      nmpcSine_B.b_U_ex = nmpcSine_B.b_U_b[nmpcSine_B.yk + 1];
      nmpcSine_B.b_X_e[3] = (((((nmpcSine_B.t2_p * nmpcSine_B.b_U_bb +
        nmpcSine_B.t3_p * nmpcSine_B.b_U_a) - 0.80982825268897241 *
        nmpcSine_B.b_X_ek) + 0.079415032681564857 * nmpcSine_B.b_U_g) +
        0.079415032681564857 * nmpcSine_B.b_U_ex) + 63.0144902620752 *
        nmpcSine_B.b_ic_idx_0 * nmpcSine_B.b_X_o * 0.022897834317090127) +
        nmpcSine_B.b_X_o * nmpcSine_B.b_X_o * 63.0144902620752 *
        0.022897834317090127 * 0.699999982041001;
      nmpcSine_B.b_X_tmp_f = 125.558050313511 * nmpcSine_B.b_X_o *
        0.03333333342744111;
      nmpcSine_B.b_yk_idx_0 = 43.6722524126937 * nmpcSine_B.b_ic_idx_0 *
        0.03333333342744111 * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_yk_idx_1 = 63.0144902620752 * nmpcSine_B.b_ic_idx_0 *
        0.03333333342744111 * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_yk_idx_2 = 0.11560777874744173 * nmpcSine_B.b_U_g *
        0.699999982041001 * 0.415;
      nmpcSine_B.b_yk_idx_3 = 0.11560777874744173 * nmpcSine_B.b_U_ex *
        0.699999982041001 * 0.415;
      nmpcSine_B.b_X_e[4] = ((((((((((((nmpcSine_B.t2_p * nmpcSine_B.b_U_a -
        nmpcSine_B.t3_p * nmpcSine_B.b_U_bb) - 74.8832173268989 *
        nmpcSine_B.b_ic_idx_0 * 0.015869365852854363) - 74.8832173268989 *
        nmpcSine_B.b_ic_idx_0 * 0.48999997485740171 * 0.03333333342744111) +
        nmpcSine_B.b_X_tmp_f * 0.699999982041001) - 43.6722524126937 *
        nmpcSine_B.b_X_o * 0.015869365852854363 * nmpcSine_B.b_X_ek) -
        nmpcSine_B.b_yk_idx_2) + nmpcSine_B.b_yk_idx_3) - 43.6722524126937 *
        nmpcSine_B.b_X_o * 0.48999997485740171 * 0.03333333342744111 *
        nmpcSine_B.b_X_ek) + 63.0144902620752 * nmpcSine_B.b_X_o *
        0.48999997485740171 * 0.03333333342744111 * nmpcSine_B.b_X_ek) -
        nmpcSine_B.b_yk_idx_0 * 0.699999982041001) + nmpcSine_B.b_yk_idx_1 *
        0.699999982041001) + 63.0144902620752 * nmpcSine_B.t2_p *
        0.48999997485740171 * 0.03333333342744111 * nmpcSine_B.b_U_a) -
        63.0144902620752 * nmpcSine_B.t3_p * 0.48999997485740171 *
        0.03333333342744111 * nmpcSine_B.b_U_bb;
      nmpcSine_B.b_X_tmp_h22 = 0.11560777874744173 * nmpcSine_B.b_U_g * 0.415;
      nmpcSine_B.b_X_tmp_e = 0.11560777874744173 * nmpcSine_B.b_U_ex * 0.415;
      nmpcSine_B.b_X_e[5] = ((((((74.8832173268989 * nmpcSine_B.b_ic_idx_0 *
        0.03333333342744111 * 0.699999982041001 + ((-nmpcSine_B.b_X_tmp_f +
        nmpcSine_B.b_X_tmp_h22) - nmpcSine_B.b_X_tmp_e)) + nmpcSine_B.b_yk_idx_0)
        - nmpcSine_B.b_yk_idx_1) + 43.6722524126937 * nmpcSine_B.b_X_o *
        0.03333333342744111 * nmpcSine_B.b_X_ek * 0.699999982041001) -
        63.0144902620752 * nmpcSine_B.b_X_o * 0.03333333342744111 *
        nmpcSine_B.b_X_ek * 0.699999982041001) - 63.0144902620752 *
        nmpcSine_B.t2_p * 0.03333333342744111 * nmpcSine_B.b_U_a *
        0.699999982041001) + 63.0144902620752 * nmpcSine_B.t3_p *
        0.03333333342744111 * nmpcSine_B.b_U_bb * 0.699999982041001;
      nmpcSine_B.b_ic_idx_0 = nmpcSine_B.b_X_l[nmpcSine_B.vectorUB_o + 4];
      nmpcSine_B.b_X_ek = nmpcSine_B.b_X_l[nmpcSine_B.vectorUB_o + 3];
      nmpcSine_B.b_X_a[0] = -nmpcSine_B.b_ic_idx_0 * nmpcSine_B.b_t3_j +
        nmpcSine_B.b_t2_a * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_X_a[1] = nmpcSine_B.b_ic_idx_0 * nmpcSine_B.b_t2_a +
        nmpcSine_B.b_t3_j * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_X_o = nmpcSine_B.b_X_l[nmpcSine_B.vectorUB_o + 5];
      nmpcSine_B.b_X_a[2] = nmpcSine_B.b_X_o;
      nmpcSine_B.b_X_a[3] = (((((nmpcSine_B.b_t2_a * nmpcSine_B.b_U_bb +
        nmpcSine_B.b_t3_j * nmpcSine_B.b_U_a) - 0.80982825268897241 *
        nmpcSine_B.b_X_ek) + 0.079415032681564857 * nmpcSine_B.b_U_g) +
        0.079415032681564857 * nmpcSine_B.b_U_ex) + 63.0144902620752 *
        nmpcSine_B.b_ic_idx_0 * nmpcSine_B.b_X_o * 0.022897834317090127) +
        nmpcSine_B.b_X_o * nmpcSine_B.b_X_o * 63.0144902620752 *
        0.022897834317090127 * 0.699999982041001;
      nmpcSine_B.b_X_tmp_f = 125.558050313511 * nmpcSine_B.b_X_o *
        0.03333333342744111;
      nmpcSine_B.b_yk_idx_0 = 43.6722524126937 * nmpcSine_B.b_ic_idx_0 *
        0.03333333342744111 * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_yk_idx_1 = 63.0144902620752 * nmpcSine_B.b_ic_idx_0 *
        0.03333333342744111 * nmpcSine_B.b_X_ek;
      nmpcSine_B.b_X_a[4] = ((((((((((((nmpcSine_B.b_t2_a * nmpcSine_B.b_U_a -
        nmpcSine_B.b_t3_j * nmpcSine_B.b_U_bb) - 74.8832173268989 *
        nmpcSine_B.b_ic_idx_0 * 0.015869365852854363) - 74.8832173268989 *
        nmpcSine_B.b_ic_idx_0 * 0.48999997485740171 * 0.03333333342744111) +
        nmpcSine_B.b_X_tmp_f * 0.699999982041001) - 43.6722524126937 *
        nmpcSine_B.b_X_o * 0.015869365852854363 * nmpcSine_B.b_X_ek) -
        nmpcSine_B.b_yk_idx_2) + nmpcSine_B.b_yk_idx_3) - 43.6722524126937 *
        nmpcSine_B.b_X_o * 0.48999997485740171 * 0.03333333342744111 *
        nmpcSine_B.b_X_ek) + 63.0144902620752 * nmpcSine_B.b_X_o *
        0.48999997485740171 * 0.03333333342744111 * nmpcSine_B.b_X_ek) -
        nmpcSine_B.b_yk_idx_0 * 0.699999982041001) + nmpcSine_B.b_yk_idx_1 *
        0.699999982041001) + 63.0144902620752 * nmpcSine_B.b_t2_a *
        0.48999997485740171 * 0.03333333342744111 * nmpcSine_B.b_U_a) -
        63.0144902620752 * nmpcSine_B.b_t3_j * 0.48999997485740171 *
        0.03333333342744111 * nmpcSine_B.b_U_bb;
      nmpcSine_B.b_X_a[5] = ((((((74.8832173268989 * nmpcSine_B.b_ic_idx_0 *
        0.03333333342744111 * 0.699999982041001 + ((-nmpcSine_B.b_X_tmp_f +
        nmpcSine_B.b_X_tmp_h22) - nmpcSine_B.b_X_tmp_e)) + nmpcSine_B.b_yk_idx_0)
        - nmpcSine_B.b_yk_idx_1) + 43.6722524126937 * nmpcSine_B.b_X_o *
        0.03333333342744111 * nmpcSine_B.b_X_ek * 0.699999982041001) -
        63.0144902620752 * nmpcSine_B.b_X_o * 0.03333333342744111 *
        nmpcSine_B.b_X_ek * 0.699999982041001) - 63.0144902620752 *
        nmpcSine_B.b_t2_a * 0.03333333342744111 * nmpcSine_B.b_U_a *
        0.699999982041001) + 63.0144902620752 * nmpcSine_B.b_t3_j *
        0.03333333342744111 * nmpcSine_B.b_U_bb * 0.699999982041001;
      for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < 6; nmpcSine_B.k_a++) {
        nmpcSine_B.b_ic_idx_0 = nmpcSine_B.ic_a[nmpcSine_B.k_a];
        nmpcSine_B.reshapes_f1[static_cast<int32_T>(nmpcSine_B.b_ic_idx_0) - 1] =
          (nmpcSine_B.b_X_l[6 * nmpcSine_B.n + nmpcSine_B.k_a] +
           (nmpcSine_B.b_X_e[nmpcSine_B.k_a] + nmpcSine_B.b_X_a[nmpcSine_B.k_a])
           * 0.05) - nmpcSine_B.b_X_l[nmpcSine_B.vectorUB_o + nmpcSine_B.k_a];
        nmpcSine_B.ic_a[nmpcSine_B.k_a] = nmpcSine_B.b_ic_idx_0 + 6.0;
      }
    }

    for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < 80; nmpcSine_B.k_a++) {
      nmpcSine_B.bv[nmpcSine_B.k_a] = rtIsInf
        (obj_next_next_next_next_next_ne->OutputMin[nmpcSine_B.k_a]);
    }

    nmpcSine_all(nmpcSine_B.bv, nmpcSine_B.b_x_b);
    nmpcSine_B.y_fn = true;
    nmpcSine_B.k_a = 0;
    exitg1 = false;
    while ((!exitg1) && (nmpcSine_B.k_a < 4)) {
      if (!nmpcSine_B.b_x_b[nmpcSine_B.k_a]) {
        nmpcSine_B.y_fn = false;
        exitg1 = true;
      } else {
        nmpcSine_B.k_a++;
      }
    }

    if (nmpcSine_B.y_fn) {
      for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < 80; nmpcSine_B.k_a++) {
        nmpcSine_B.bv[nmpcSine_B.k_a] = rtIsInf
          (obj_next_next_next_next_next_ne->OutputMax[nmpcSine_B.k_a]);
      }

      nmpcSine_all(nmpcSine_B.bv, nmpcSine_B.b_x_b);
      nmpcSine_B.y_fn = true;
      nmpcSine_B.k_a = 0;
      exitg1 = false;
      while ((!exitg1) && (nmpcSine_B.k_a < 4)) {
        if (!nmpcSine_B.b_x_b[nmpcSine_B.k_a]) {
          nmpcSine_B.y_fn = false;
          exitg1 = true;
        } else {
          nmpcSine_B.k_a++;
        }
      }
    } else {
      nmpcSine_B.y_fn = false;
    }

    if (nmpcSine_B.y_fn) {
      nmpcSine_B.yk = 0;
      nmpcSine_B.n = 0;
    } else {
      for (nmpcSine_B.n = 0; nmpcSine_B.n < 160; nmpcSine_B.n++) {
        nmpcSine_B.c[nmpcSine_B.n] = 0.0;
        nmpcSine_B.icf[nmpcSine_B.n] = true;
      }

      nmpcSine_B.b_ic_idx_0 = 1.0;
      nmpcSine_B.t2_p = 2.0;
      nmpcSine_B.t3_p = 3.0;
      nmpcSine_B.b_t2_a = 4.0;
      for (nmpcSine_B.n = 0; nmpcSine_B.n < 20; nmpcSine_B.n++) {
        nmpcSine_B.b_t3_j = obj_next_next_next_next_next_ne->
          OutputMin[nmpcSine_B.n];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.b_ic_idx_0) - 1] =
          ((!rtIsInf(nmpcSine_B.b_t3_j)) && (!rtIsNaN(nmpcSine_B.b_t3_j)));
        nmpcSine_B.b_X_ek = obj_next_next_next_next_next_ne->
          OutputMin[nmpcSine_B.n + 20];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.t2_p) - 1] = ((!rtIsInf
          (nmpcSine_B.b_X_ek)) && (!rtIsNaN(nmpcSine_B.b_X_ek)));
        nmpcSine_B.b_X_o = obj_next_next_next_next_next_ne->
          OutputMin[nmpcSine_B.n + 40];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.t3_p) - 1] = ((!rtIsInf
          (nmpcSine_B.b_X_o)) && (!rtIsNaN(nmpcSine_B.b_X_o)));
        nmpcSine_B.b_U_bb = obj_next_next_next_next_next_ne->
          OutputMin[nmpcSine_B.n + 60];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.b_t2_a) - 1] = ((!rtIsInf
          (nmpcSine_B.b_U_bb)) && (!rtIsNaN(nmpcSine_B.b_U_bb)));
        nmpcSine_B.b_U_a = obj_next_next_next_next_next_ne->
          OutputMax[nmpcSine_B.n];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.b_ic_idx_0 + 4.0) - 1] =
          ((!rtIsInf(nmpcSine_B.b_U_a)) && (!rtIsNaN(nmpcSine_B.b_U_a)));
        nmpcSine_B.icf_tmp[0] = static_cast<int32_T>(nmpcSine_B.b_ic_idx_0) - 1;
        nmpcSine_B.icf_tmp[4] = static_cast<int32_T>(nmpcSine_B.b_ic_idx_0 + 4.0)
          - 1;
        nmpcSine_B.b_U_g = obj_next_next_next_next_next_ne->
          OutputMax[nmpcSine_B.n + 20];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.t2_p + 4.0) - 1] =
          ((!rtIsInf(nmpcSine_B.b_U_g)) && (!rtIsNaN(nmpcSine_B.b_U_g)));
        nmpcSine_B.icf_tmp[1] = static_cast<int32_T>(nmpcSine_B.t2_p) - 1;
        nmpcSine_B.icf_tmp[5] = static_cast<int32_T>(nmpcSine_B.t2_p + 4.0) - 1;
        nmpcSine_B.b_U_ex = obj_next_next_next_next_next_ne->
          OutputMax[nmpcSine_B.n + 40];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.t3_p + 4.0) - 1] =
          ((!rtIsInf(nmpcSine_B.b_U_ex)) && (!rtIsNaN(nmpcSine_B.b_U_ex)));
        nmpcSine_B.icf_tmp[2] = static_cast<int32_T>(nmpcSine_B.t3_p) - 1;
        nmpcSine_B.icf_tmp[6] = static_cast<int32_T>(nmpcSine_B.t3_p + 4.0) - 1;
        nmpcSine_B.b_X_tmp_f = obj_next_next_next_next_next_ne->
          OutputMax[nmpcSine_B.n + 60];
        nmpcSine_B.icf[static_cast<int32_T>(nmpcSine_B.b_t2_a + 4.0) - 1] =
          ((!rtIsInf(nmpcSine_B.b_X_tmp_f)) && (!rtIsNaN(nmpcSine_B.b_X_tmp_f)));
        nmpcSine_B.icf_tmp[3] = static_cast<int32_T>(nmpcSine_B.b_t2_a) - 1;
        nmpcSine_B.icf_tmp[7] = static_cast<int32_T>(nmpcSine_B.b_t2_a + 4.0) -
          1;
        for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < 8; nmpcSine_B.k_a++) {
          nmpcSine_B.icf_o[nmpcSine_B.k_a] =
            nmpcSine_B.icf[nmpcSine_B.icf_tmp[nmpcSine_B.k_a]];
        }

        if (nmpcSine_any(nmpcSine_B.icf_o)) {
          nmpcSine_B.b_yk_idx_0 = nmpcSine_B.X_g[nmpcSine_B.n + 1];
          nmpcSine_B.b_yk_idx_1 = nmpcSine_B.X_g[nmpcSine_B.n + 22];
          nmpcSine_B.b_yk_idx_2 = nmpcSine_B.X_g[nmpcSine_B.n + 43];
          nmpcSine_B.b_yk_idx_3 = nmpcSine_B.X_g[nmpcSine_B.n + 106];
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.b_ic_idx_0) - 1] =
            (nmpcSine_B.b_t3_j - nmpcSine_B.e_p) - nmpcSine_B.b_yk_idx_0;
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.t2_p) - 1] =
            (nmpcSine_B.b_X_ek - nmpcSine_B.e_p) - nmpcSine_B.b_yk_idx_1;
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.t3_p) - 1] =
            (nmpcSine_B.b_X_o - nmpcSine_B.e_p) - nmpcSine_B.b_yk_idx_2;
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.b_t2_a) - 1] =
            (nmpcSine_B.b_U_bb - nmpcSine_B.e_p) - nmpcSine_B.b_yk_idx_3;
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.b_ic_idx_0 + 4.0) - 1] =
            (nmpcSine_B.b_yk_idx_0 - nmpcSine_B.b_U_a) - nmpcSine_B.e_p;
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.t2_p + 4.0) - 1] =
            (nmpcSine_B.b_yk_idx_1 - nmpcSine_B.b_U_g) - nmpcSine_B.e_p;
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.t3_p + 4.0) - 1] =
            (nmpcSine_B.b_yk_idx_2 - nmpcSine_B.b_U_ex) - nmpcSine_B.e_p;
          nmpcSine_B.c[static_cast<int32_T>(nmpcSine_B.b_t2_a + 4.0) - 1] =
            (nmpcSine_B.b_yk_idx_3 - nmpcSine_B.b_X_tmp_f) - nmpcSine_B.e_p;
        }

        nmpcSine_B.b_ic_idx_0 += 8.0;
        nmpcSine_B.t2_p += 8.0;
        nmpcSine_B.t3_p += 8.0;
        nmpcSine_B.b_t2_a += 8.0;
      }

      nmpcSine_B.k_a = 0;
      for (nmpcSine_B.n = 0; nmpcSine_B.n < 160; nmpcSine_B.n++) {
        if (nmpcSine_B.icf[nmpcSine_B.n]) {
          nmpcSine_B.k_a++;
        }
      }

      nmpcSine_B.yk = nmpcSine_B.k_a;
      nmpcSine_B.k_a = 0;
      for (nmpcSine_B.n = 0; nmpcSine_B.n < 160; nmpcSine_B.n++) {
        if (nmpcSine_B.icf[nmpcSine_B.n]) {
          nmpcSine_B.tmp_data_n[nmpcSine_B.k_a] = static_cast<uint8_T>
            (nmpcSine_B.n);
          nmpcSine_B.k_a++;
        }
      }

      nmpcSine_B.n = 1;
      for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < nmpcSine_B.yk; nmpcSine_B.k_a++)
      {
        nmpcSine_B.varargin_1_data_j[nmpcSine_B.k_a] =
          nmpcSine_B.c[nmpcSine_B.tmp_data_n[nmpcSine_B.k_a]];
      }
    }

    nmpcSine_B.y_fn = ((nmpcSine_B.yk != 0) && (nmpcSine_B.n != 0));
    if (!nmpcSine_B.y_fn) {
      nmpcSine_B.sizes_g[0] = static_cast<uint8_T>(nmpcSine_B.yk);
    } else if (nmpcSine_B.y_fn) {
      nmpcSine_B.sizes_g[0] = static_cast<uint8_T>(nmpcSine_B.yk);
    } else {
      nmpcSine_B.sizes_g[0] = 0U;
    }

    nmpcSine_B.n = nmpcSine_B.sizes_g[0];
    nmpcSine_B.yk = nmpcSine_B.y_fn;
    for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < nmpcSine_B.yk; nmpcSine_B.k_a++) {
      if (nmpcSine_B.n - 1 >= 0) {
        memcpy(&nmpcSine_B.varargin_1_data[0], &nmpcSine_B.varargin_1_data_j[0],
               static_cast<uint32_T>(nmpcSine_B.n) * sizeof(real_T));
      }
    }

    for (nmpcSine_B.k_a = 0; nmpcSine_B.k_a < nmpcSine_B.ineqRange_size_idx_1;
         nmpcSine_B.k_a++) {
      Cineq_workspace_data[nmpcSine_B.ineqRange_data[nmpcSine_B.k_a] - 1] =
        nmpcSine_B.varargin_1_data[nmpcSine_B.k_a];
    }

    memcpy(&Ceq_workspace[0], &nmpcSine_B.reshapes_f1[0], 120U * sizeof(real_T));
    *status = nmpcSine_checkVectorNonFinite(obj_next_next_next_next_next_b_,
      Cineq_workspace_data, ineq0);
    if (*status == 1) {
      *status = nmpcSine_checkVectorNonFinite_n(Ceq_workspace);
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_computeLinearResiduals(const real_T x[125], int32_T nVar,
  real_T workspaceIneq_data[], const int32_T workspaceIneq_size[1], int32_T
  mLinIneq, const real_T AineqT_data[], const real_T bineq_data[], int32_T ldAi)
{
  if (mLinIneq > 0) {
    int32_T k;
    int32_T loop_ub;
    int32_T scalarLB;
    int32_T vectorUB;
    loop_ub = workspaceIneq_size[0];
    if (loop_ub - 1 >= 0) {
      memcpy(&nmpcSine_B.y_data_f[0], &workspaceIneq_data[0],
             static_cast<uint32_T>(loop_ub) * sizeof(real_T));
    }

    memcpy(&nmpcSine_B.y_data_f[0], &bineq_data[0], static_cast<uint32_T>
           (mLinIneq) * sizeof(real_T));
    if (loop_ub - 1 >= 0) {
      memcpy(&workspaceIneq_data[0], &nmpcSine_B.y_data_f[0], static_cast<
             uint32_T>(loop_ub) * sizeof(real_T));
    }

    k = static_cast<uint16_T>(mLinIneq);
    scalarLB = (static_cast<uint16_T>(mLinIneq) / 2) << 1;
    vectorUB = scalarLB - 2;
    for (loop_ub = 0; loop_ub <= vectorUB; loop_ub += 2) {
      __m128d tmp;
      tmp = _mm_loadu_pd(&workspaceIneq_data[loop_ub]);
      _mm_storeu_pd(&workspaceIneq_data[loop_ub], _mm_mul_pd(tmp, _mm_set1_pd
        (-1.0)));
    }

    for (loop_ub = scalarLB; loop_ub < k; loop_ub++) {
      workspaceIneq_data[loop_ub] = -workspaceIneq_data[loop_ub];
    }

    scalarLB = 0;
    vectorUB = (mLinIneq - 1) * ldAi + 1;
    for (loop_ub = 1; ldAi < 0 ? loop_ub >= vectorUB : loop_ub <= vectorUB;
         loop_ub += ldAi) {
      real_T c;
      int32_T e;
      c = 0.0;
      e = (loop_ub + nVar) - 1;
      for (k = loop_ub; k <= e; k++) {
        c += AineqT_data[k - 1] * x[k - loop_ub];
      }

      workspaceIneq_data[scalarLB] += c;
      scalarLB++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
real_T nmpcSine::nmpcSine_computeMeritFcn(real_T obj_penaltyParam, real_T fval,
  const real_T Cineq_workspace_data[], int32_T mIneq, const real_T
  Ceq_workspace[120], boolean_T evalWellDefined)
{
  real_T val;
  if (evalWellDefined) {
    real_T constrViolationEq;
    real_T constrViolationIneq;
    int32_T k;
    constrViolationEq = 0.0;
    for (k = 0; k < 120; k++) {
      constrViolationEq += fabs(Ceq_workspace[k]);
    }

    constrViolationIneq = 0.0;
    k = static_cast<uint16_T>(mIneq);
    for (int32_T idx = 0; idx < k; idx++) {
      real_T Cineq_workspace;
      Cineq_workspace = Cineq_workspace_data[idx];
      if (Cineq_workspace > 0.0) {
        constrViolationIneq += Cineq_workspace;
      }
    }

    val = (constrViolationEq + constrViolationIneq) * obj_penaltyParam + fval;
  } else {
    val = (rtInf);
  }

  return val;
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_linesearch(boolean_T *evalWellDefined, const real_T
  bineq_data[], int32_T WorkingSet_nVar, int32_T WorkingSet_ldA, const real_T
  WorkingSet_Aineq_data[], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState, real_T
  MeritFunction_penaltyParam, real_T MeritFunction_phi, real_T
  MeritFunction_phiPrimePlus, real_T MeritFunction_phiFullStep, int32_T
  FcnEvaluator_next_next_next_nex, const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
  *FcnEvaluator_next_next_next_n_0, const s_I4XPpWw7d7shktLagLlNtD_nmpc_T
  *FcnEvaluator_next_next_next_n_1, boolean_T socTaken, real_T *alpha, int32_T
  *exitflag)
{
  __m128d tmp;
  __m128d tmp_0;
  int32_T exitg1;
  boolean_T exitg2;
  nmpcSine_B.mLinIneq_c = TrialState->mIneq - TrialState->mNonlinIneq;
  *alpha = 1.0;
  *exitflag = 1;
  nmpcSine_B.phi_alpha = MeritFunction_phiFullStep;
  if (WorkingSet_nVar - 1 >= 0) {
    memcpy(&TrialState->searchDir.data[0], &TrialState->delta_x.data[0],
           static_cast<uint32_T>(WorkingSet_nVar) * sizeof(real_T));
  }

  do {
    exitg1 = 0;
    if (TrialState->FunctionEvaluations < 12500) {
      if ((*evalWellDefined) && (nmpcSine_B.phi_alpha <= *alpha * 0.0001 *
           MeritFunction_phiPrimePlus + MeritFunction_phi)) {
        exitg1 = 1;
      } else {
        *alpha *= 0.7;
        nmpcSine_B.k_p = static_cast<uint16_T>(WorkingSet_nVar);
        nmpcSine_B.scalarLB_tmp = (static_cast<uint16_T>(WorkingSet_nVar) / 2) <<
          1;
        nmpcSine_B.vectorUB_tmp = nmpcSine_B.scalarLB_tmp - 2;
        for (nmpcSine_B.idx = 0; nmpcSine_B.idx <= nmpcSine_B.vectorUB_tmp;
             nmpcSine_B.idx += 2) {
          tmp_0 = _mm_loadu_pd(&TrialState->xstar.data[nmpcSine_B.idx]);
          _mm_storeu_pd(&TrialState->delta_x.data[nmpcSine_B.idx], _mm_mul_pd
                        (_mm_set1_pd(*alpha), tmp_0));
        }

        for (nmpcSine_B.idx = nmpcSine_B.scalarLB_tmp; nmpcSine_B.idx <
             nmpcSine_B.k_p; nmpcSine_B.idx++) {
          TrialState->delta_x.data[nmpcSine_B.idx] = *alpha *
            TrialState->xstar.data[nmpcSine_B.idx];
        }

        if (socTaken) {
          nmpcSine_B.phi_alpha = *alpha * *alpha;
          if ((WorkingSet_nVar >= 1) && (!(nmpcSine_B.phi_alpha == 0.0))) {
            nmpcSine_B.scalarLB = (WorkingSet_nVar / 2) << 1;
            nmpcSine_B.vectorUB_p = nmpcSine_B.scalarLB - 2;
            for (nmpcSine_B.idx = 0; nmpcSine_B.idx <= nmpcSine_B.vectorUB_p;
                 nmpcSine_B.idx += 2) {
              tmp_0 = _mm_loadu_pd(&TrialState->socDirection.data[nmpcSine_B.idx]);
              tmp = _mm_loadu_pd(&TrialState->delta_x.data[nmpcSine_B.idx]);
              _mm_storeu_pd(&TrialState->delta_x.data[nmpcSine_B.idx],
                            _mm_add_pd(_mm_mul_pd(_mm_set1_pd
                (nmpcSine_B.phi_alpha), tmp_0), tmp));
            }

            for (nmpcSine_B.idx = nmpcSine_B.scalarLB; nmpcSine_B.idx <
                 WorkingSet_nVar; nmpcSine_B.idx++) {
              TrialState->delta_x.data[nmpcSine_B.idx] += nmpcSine_B.phi_alpha *
                TrialState->socDirection.data[nmpcSine_B.idx];
            }
          }
        }

        nmpcSine_B.tooSmallX = true;
        nmpcSine_B.idx = 0;
        exitg2 = false;
        while ((!exitg2) && (nmpcSine_B.idx <= static_cast<uint16_T>
                             (WorkingSet_nVar) - 1)) {
          nmpcSine_B.phi_alpha = fabs(TrialState->xstarsqp[nmpcSine_B.idx]);
          if ((nmpcSine_B.phi_alpha <= 1.0) || rtIsNaN(nmpcSine_B.phi_alpha)) {
            nmpcSine_B.phi_alpha = 1.0;
          }

          if (0.001 * nmpcSine_B.phi_alpha <= fabs(TrialState->
               delta_x.data[nmpcSine_B.idx])) {
            nmpcSine_B.tooSmallX = false;
            exitg2 = true;
          } else {
            nmpcSine_B.idx++;
          }
        }

        if (nmpcSine_B.tooSmallX) {
          *exitflag = -2;
          exitg1 = 1;
        } else {
          for (nmpcSine_B.idx = 0; nmpcSine_B.idx <= nmpcSine_B.vectorUB_tmp;
               nmpcSine_B.idx += 2) {
            tmp_0 = _mm_loadu_pd(&TrialState->xstarsqp_old[nmpcSine_B.idx]);
            tmp = _mm_loadu_pd(&TrialState->delta_x.data[nmpcSine_B.idx]);
            _mm_storeu_pd(&TrialState->xstarsqp[nmpcSine_B.idx], _mm_add_pd
                          (tmp_0, tmp));
          }

          for (nmpcSine_B.idx = nmpcSine_B.scalarLB_tmp; nmpcSine_B.idx <
               nmpcSine_B.k_p; nmpcSine_B.idx++) {
            TrialState->xstarsqp[nmpcSine_B.idx] = TrialState->
              xstarsqp_old[nmpcSine_B.idx] + TrialState->
              delta_x.data[nmpcSine_B.idx];
          }

          nmpcSine_evalObjAndConstr(FcnEvaluator_next_next_next_nex,
            FcnEvaluator_next_next_next_n_0, FcnEvaluator_next_next_next_n_1,
            TrialState->xstarsqp, TrialState->cIneq.data, TrialState->iNonIneq0,
            TrialState->cEq, &TrialState->sqpFval, &nmpcSine_B.k_p);
          nmpcSine_computeLinearResiduals(TrialState->xstarsqp, WorkingSet_nVar,
            TrialState->cIneq.data, TrialState->cIneq.size,
            nmpcSine_B.mLinIneq_c, WorkingSet_Aineq_data, bineq_data,
            WorkingSet_ldA);
          TrialState->FunctionEvaluations++;
          *evalWellDefined = (nmpcSine_B.k_p == 1);
          nmpcSine_B.phi_alpha = nmpcSine_computeMeritFcn
            (MeritFunction_penaltyParam, TrialState->sqpFval,
             TrialState->cIneq.data, TrialState->mIneq, TrialState->cEq,
             *evalWellDefined);
        }
      }
    } else {
      *exitflag = 0;
      exitg1 = 1;
    }
  } while (exitg1 == 0);
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_driver(const real_T bineq_data[], const real_T lb[125],
  const real_T ub[125], s_2COE1uYisQtyPYvPjrXP9G_nmpc_T *TrialState,
  sG8JZ69axY52WWR6RKyApQC_nmpcS_T *MeritFunction, const
  coder_internal_stickyStruct_2_T *FcnEvaluator, s_kmYqIq13KlaOrGCTq3ShMG_nmpc_T
  *memspace, s_OouGZKqwdkE6H2b5xBmmyD_nmpc_T *WorkingSet,
  s_0RmwrXfzGd5lqbHvgKQe2_nmpcS_T *QRManager, s_mDApTYzDBpxuvxemclsuEF_nmpc_T
  *CholManager, s_xtSBzQGTZuMYOTjcuMqLQH_nmpc_T *QPObjective, const int32_T
  fscales_lineq_constraint_size[1], const int32_T fscales_cineq_constraint_size
  [1], real_T Hessian[15625])
{
  __m128d tmp;
  __m128d tmp_0;
  static const int8_T s[15625] = { 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 1 };

  static const char_T r[7] = { 'f', 'm', 'i', 'n', 'c', 'o', 'n' };

  for (nmpcSine_B.k = 0; nmpcSine_B.k < 15625; nmpcSine_B.k++) {
    Hessian[nmpcSine_B.k] = s[nmpcSine_B.k];
  }

  nmpcSine_B.nVar_tmp_tmp = WorkingSet->nVar;
  nmpcSine_B.mFixed_d = WorkingSet->sizes[0];
  nmpcSine_B.mIneq_i = WorkingSet->sizes[2];
  nmpcSine_B.mLB_g = WorkingSet->sizes[3];
  nmpcSine_B.mUB_n = WorkingSet->sizes[4];
  nmpcSine_B.mConstr = (((WorkingSet->sizes[0] + WorkingSet->sizes[2]) +
    WorkingSet->sizes[3]) + WorkingSet->sizes[4]) + 120;
  nmpcSine_B.mLinIneq = WorkingSet->sizes[2] - TrialState->mNonlinIneq;
  nmpcSine_B.u1_n = ((WorkingSet->sizes[2] + WorkingSet->sizes[3]) +
                     WorkingSet->sizes[4]) + (WorkingSet->sizes[0] << 1);
  if (WorkingSet->nVar >= nmpcSine_B.u1_n) {
    nmpcSine_B.u1_n = WorkingSet->nVar;
  }

  nmpcSine_B.u1_n *= 10;
  TrialState->steplength = 1.0;
  nmpcSine_test_exit(MeritFunction, WorkingSet, TrialState, lb, ub,
                     &nmpcSine_B.Flags.gradOK, &nmpcSine_B.Flags.fevalOK,
                     &nmpcSine_B.Flags.done, &nmpcSine_B.Flags.stepAccepted,
                     &nmpcSine_B.Flags.failedLineSearch,
                     &nmpcSine_B.Flags.stepType);
  nmpcSine_saveJacobian(TrialState, WorkingSet->nVar, WorkingSet->sizes[2],
                        WorkingSet->Aineq.data, TrialState->iNonIneq0,
                        WorkingSet->Aeq.data, WorkingSet->ldA);
  TrialState->sqpFval_old = TrialState->sqpFval;
  for (nmpcSine_B.k = 0; nmpcSine_B.k < 125; nmpcSine_B.k++) {
    TrialState->xstarsqp_old[nmpcSine_B.k] = TrialState->xstarsqp[nmpcSine_B.k];
    TrialState->grad_old.data[nmpcSine_B.k] = TrialState->grad.data[nmpcSine_B.k];
  }

  nmpcSine_B.d_ix_tmp = TrialState->mIneq;
  nmpcSine_B.loop_ub_p = TrialState->cIneq_old.size[0];
  nmpcSine_B.loop_ub_c = TrialState->cIneq_old.size[0];
  if (nmpcSine_B.loop_ub_c - 1 >= 0) {
    memcpy(&nmpcSine_B.y_data_c[0], &TrialState->cIneq_old.data[0], static_cast<
           uint32_T>(nmpcSine_B.loop_ub_c) * sizeof(real_T));
  }

  if (nmpcSine_B.d_ix_tmp - 1 >= 0) {
    memcpy(&nmpcSine_B.y_data_c[0], &TrialState->cIneq.data[0],
           static_cast<uint32_T>(nmpcSine_B.d_ix_tmp) * sizeof(real_T));
  }

  if (nmpcSine_B.loop_ub_p - 1 >= 0) {
    memcpy(&TrialState->cIneq_old.data[0], &nmpcSine_B.y_data_c[0],
           static_cast<uint32_T>(nmpcSine_B.loop_ub_p) * sizeof(real_T));
  }

  memcpy(&TrialState->cEq_old[0], &TrialState->cEq[0], 120U * sizeof(real_T));
  if (!nmpcSine_B.Flags.done) {
    TrialState->sqpIterations = 1;
  }

  while (!nmpcSine_B.Flags.done) {
    if ((!nmpcSine_B.Flags.stepAccepted) && (!nmpcSine_B.Flags.failedLineSearch))
    {
      nmpcSine_B.expl_temp_j.IterDisplayQP = false;
      nmpcSine_B.expl_temp_j.RemainFeasible = false;
      nmpcSine_B.expl_temp_j.ProbRelTolFactor = 1.0;
      nmpcSine_B.expl_temp_j.ConstrRelTolFactor = 1.0;
      nmpcSine_B.expl_temp_j.PricingTolerance = 0.0;
      nmpcSine_B.expl_temp_j.ObjectiveLimit = (rtMinusInf);
      nmpcSine_B.expl_temp_j.ConstraintTolerance = 0.001;
      nmpcSine_B.expl_temp_j.OptimalityTolerance = 2.2204460492503131E-14;
      nmpcSine_B.expl_temp_j.StepTolerance = 1.0E-6;
      nmpcSine_B.expl_temp_j.MaxIterations = nmpcSine_B.u1_n;
      for (nmpcSine_B.k = 0; nmpcSine_B.k < 7; nmpcSine_B.k++) {
        nmpcSine_B.expl_temp_j.SolverName[nmpcSine_B.k] = r[nmpcSine_B.k];
      }
    }

    while ((!nmpcSine_B.Flags.stepAccepted) &&
           (!nmpcSine_B.Flags.failedLineSearch)) {
      if (nmpcSine_B.Flags.stepType != 3) {
        nmpcSi_updateWorkingSetForNewQP(TrialState->xstarsqp, WorkingSet,
          nmpcSine_B.mIneq_i, TrialState->mNonlinIneq, TrialState->cIneq.data,
          TrialState->cEq, nmpcSine_B.mLB_g, lb, nmpcSine_B.mUB_n, ub,
          nmpcSine_B.mFixed_d);
      }

      nmpcSine_B.expl_temp_h = nmpcSine_B.expl_temp_j;
      nmpcSine_step_k(&nmpcSine_B.Flags, Hessian, lb, ub, TrialState,
                      MeritFunction, memspace, WorkingSet, QRManager,
                      CholManager, QPObjective, &nmpcSine_B.expl_temp_h);
      if (nmpcSine_B.Flags.stepAccepted) {
        nmpcSine_B.loop_ub_c = static_cast<uint16_T>(nmpcSine_B.nVar_tmp_tmp);
        nmpcSine_B.loop_ub_p = (static_cast<uint16_T>(nmpcSine_B.nVar_tmp_tmp) /
          2) << 1;
        nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p - 2;
        for (nmpcSine_B.k = 0; nmpcSine_B.k <= nmpcSine_B.vectorUB; nmpcSine_B.k
             += 2) {
          tmp = _mm_loadu_pd(&TrialState->xstarsqp[nmpcSine_B.k]);
          tmp_0 = _mm_loadu_pd(&TrialState->delta_x.data[nmpcSine_B.k]);
          _mm_storeu_pd(&TrialState->xstarsqp[nmpcSine_B.k], _mm_add_pd(tmp,
            tmp_0));
        }

        for (nmpcSine_B.k = nmpcSine_B.loop_ub_p; nmpcSine_B.k <
             nmpcSine_B.loop_ub_c; nmpcSine_B.k++) {
          TrialState->xstarsqp[nmpcSine_B.k] += TrialState->
            delta_x.data[nmpcSine_B.k];
        }

        nmpcSine_evalObjAndConstr(FcnEvaluator->next.next.next.next.next.b_value,
          &FcnEvaluator->next.next.next.next.next.next.next.b_value.workspace.runtimedata,
          &FcnEvaluator->next.next.next.next.next.next.next.next.b_value.workspace,
          TrialState->xstarsqp, TrialState->cIneq.data, TrialState->iNonIneq0,
          TrialState->cEq, &TrialState->sqpFval, &nmpcSine_B.k);
        nmpcSine_B.Flags.fevalOK = (nmpcSine_B.k == 1);
        TrialState->FunctionEvaluations++;
        nmpcSine_computeLinearResiduals(TrialState->xstarsqp,
          nmpcSine_B.nVar_tmp_tmp, TrialState->cIneq.data,
          TrialState->cIneq.size, nmpcSine_B.mLinIneq, WorkingSet->Aineq.data,
          bineq_data, WorkingSet->ldA);
        MeritFunction->phiFullStep = nmpcSine_computeMeritFcn
          (MeritFunction->penaltyParam, TrialState->sqpFval,
           TrialState->cIneq.data, nmpcSine_B.mIneq_i, TrialState->cEq,
           nmpcSine_B.Flags.fevalOK);
      }

      if ((nmpcSine_B.Flags.stepType == 1) && nmpcSine_B.Flags.stepAccepted &&
          nmpcSine_B.Flags.fevalOK && (MeritFunction->phi <
           MeritFunction->phiFullStep) && (TrialState->sqpFval <
           TrialState->sqpFval_old)) {
        nmpcSine_B.Flags.stepType = 3;
        nmpcSine_B.Flags.stepAccepted = false;
      } else {
        nmpcSine_linesearch(&nmpcSine_B.Flags.fevalOK, bineq_data,
                            WorkingSet->nVar, WorkingSet->ldA,
                            WorkingSet->Aineq.data, TrialState,
                            MeritFunction->penaltyParam, MeritFunction->phi,
                            MeritFunction->phiPrimePlus,
                            MeritFunction->phiFullStep,
                            FcnEvaluator->next.next.next.next.next.b_value,
                            &FcnEvaluator->next.next.next.next.next.next.next.b_value.workspace.runtimedata,
                            &FcnEvaluator->next.next.next.next.next.next.next.next.b_value.workspace,
                            ((nmpcSine_B.Flags.stepType == 3) &&
                             nmpcSine_B.Flags.stepAccepted),
                            &nmpcSine_B.TrialState_lambdasqp, &nmpcSine_B.k);
        TrialState->steplength = nmpcSine_B.TrialState_lambdasqp;
        if (nmpcSine_B.k > 0) {
          nmpcSine_B.Flags.stepAccepted = true;
        } else {
          nmpcSine_B.Flags.failedLineSearch = true;
        }
      }
    }

    if (nmpcSine_B.Flags.stepAccepted && (!nmpcSine_B.Flags.failedLineSearch)) {
      nmpcSine_B.loop_ub_c = static_cast<uint16_T>(nmpcSine_B.nVar_tmp_tmp);
      nmpcSine_B.loop_ub_p = (static_cast<uint16_T>(nmpcSine_B.nVar_tmp_tmp) / 2)
        << 1;
      nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p - 2;
      for (nmpcSine_B.k = 0; nmpcSine_B.k <= nmpcSine_B.vectorUB; nmpcSine_B.k +=
           2) {
        tmp = _mm_loadu_pd(&TrialState->xstarsqp_old[nmpcSine_B.k]);
        tmp_0 = _mm_loadu_pd(&TrialState->delta_x.data[nmpcSine_B.k]);
        _mm_storeu_pd(&TrialState->xstarsqp[nmpcSine_B.k], _mm_add_pd(tmp, tmp_0));
      }

      for (nmpcSine_B.k = nmpcSine_B.loop_ub_p; nmpcSine_B.k <
           nmpcSine_B.loop_ub_c; nmpcSine_B.k++) {
        TrialState->xstarsqp[nmpcSine_B.k] = TrialState->
          xstarsqp_old[nmpcSine_B.k] + TrialState->delta_x.data[nmpcSine_B.k];
      }

      nmpcSine_B.loop_ub_c = static_cast<uint16_T>(nmpcSine_B.mConstr);
      nmpcSine_B.loop_ub_p = (static_cast<uint16_T>(nmpcSine_B.mConstr) / 2) <<
        1;
      nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p - 2;
      for (nmpcSine_B.k = 0; nmpcSine_B.k <= nmpcSine_B.vectorUB; nmpcSine_B.k +=
           2) {
        tmp = _mm_loadu_pd(&TrialState->lambda.data[nmpcSine_B.k]);
        tmp_0 = _mm_loadu_pd(&TrialState->lambdasqp.data[nmpcSine_B.k]);
        _mm_storeu_pd(&TrialState->lambdasqp.data[nmpcSine_B.k], _mm_add_pd
                      (_mm_mul_pd(_mm_sub_pd(tmp, tmp_0), _mm_set1_pd
          (TrialState->steplength)), tmp_0));
      }

      for (nmpcSine_B.k = nmpcSine_B.loop_ub_p; nmpcSine_B.k <
           nmpcSine_B.loop_ub_c; nmpcSine_B.k++) {
        nmpcSine_B.TrialState_lambdasqp = TrialState->
          lambdasqp.data[nmpcSine_B.k];
        TrialState->lambdasqp.data[nmpcSine_B.k] = (TrialState->
          lambda.data[nmpcSine_B.k] - nmpcSine_B.TrialState_lambdasqp) *
          TrialState->steplength + nmpcSine_B.TrialState_lambdasqp;
      }

      TrialState->sqpFval_old = TrialState->sqpFval;
      for (nmpcSine_B.k = 0; nmpcSine_B.k < 125; nmpcSine_B.k++) {
        TrialState->xstarsqp_old[nmpcSine_B.k] = TrialState->
          xstarsqp[nmpcSine_B.k];
        TrialState->grad_old.data[nmpcSine_B.k] = TrialState->
          grad.data[nmpcSine_B.k];
      }

      nmpcSine_B.loop_ub_p = TrialState->cIneq_old.size[0];
      nmpcSine_B.loop_ub_c = TrialState->cIneq_old.size[0];
      if (nmpcSine_B.loop_ub_c - 1 >= 0) {
        memcpy(&nmpcSine_B.y_data_c[0], &TrialState->cIneq_old.data[0],
               static_cast<uint32_T>(nmpcSine_B.loop_ub_c) * sizeof(real_T));
      }

      if (nmpcSine_B.d_ix_tmp - 1 >= 0) {
        memcpy(&nmpcSine_B.y_data_c[0], &TrialState->cIneq.data[0], static_cast<
               uint32_T>(nmpcSine_B.d_ix_tmp) * sizeof(real_T));
      }

      if (nmpcSine_B.loop_ub_p - 1 >= 0) {
        memcpy(&TrialState->cIneq_old.data[0], &nmpcSine_B.y_data_c[0],
               static_cast<uint32_T>(nmpcSine_B.loop_ub_p) * sizeof(real_T));
      }

      memcpy(&TrialState->cEq_old[0], &TrialState->cEq[0], 120U * sizeof(real_T));
      nmpcSine_B.Flags.gradOK = true;
      evalObjAndConstrAndDerivatives
        (FcnEvaluator->next.next.next.next.next.b_value,
         &FcnEvaluator->next.next.next.next.next.next.next.b_value.workspace.runtimedata,
         &FcnEvaluator->next.next.next.next.next.next.next.next.b_value.workspace,
         TrialState->xstarsqp, TrialState->grad.data, TrialState->cIneq.data,
         TrialState->iNonIneq0, TrialState->cEq, WorkingSet->Aineq.data,
         TrialState->iNonIneq0, WorkingSet->ldA, WorkingSet->Aeq.data,
         WorkingSet->ldA, &TrialState->sqpFval, &nmpcSine_B.k);
      TrialState->FunctionEvaluations++;
      nmpcSine_B.Flags.fevalOK = (nmpcSine_B.k == 1);
    } else {
      TrialState->sqpFval = TrialState->sqpFval_old;
      memcpy(&TrialState->xstarsqp[0], &TrialState->xstarsqp_old[0], 125U *
             sizeof(real_T));
      nmpcSine_B.loop_ub_p = TrialState->cIneq.size[0];
      nmpcSine_B.loop_ub_c = TrialState->cIneq.size[0];
      if (nmpcSine_B.loop_ub_c - 1 >= 0) {
        memcpy(&nmpcSine_B.y_data_c[0], &TrialState->cIneq.data[0], static_cast<
               uint32_T>(nmpcSine_B.loop_ub_c) * sizeof(real_T));
      }

      if (nmpcSine_B.d_ix_tmp - 1 >= 0) {
        memcpy(&nmpcSine_B.y_data_c[0], &TrialState->cIneq_old.data[0],
               static_cast<uint32_T>(nmpcSine_B.d_ix_tmp) * sizeof(real_T));
      }

      if (nmpcSine_B.loop_ub_p - 1 >= 0) {
        memcpy(&TrialState->cIneq.data[0], &nmpcSine_B.y_data_c[0],
               static_cast<uint32_T>(nmpcSine_B.loop_ub_p) * sizeof(real_T));
      }

      memcpy(&TrialState->cEq[0], &TrialState->cEq_old[0], 120U * sizeof(real_T));
    }

    nmpcSine_test_exit_n(&nmpcSine_B.Flags, memspace, MeritFunction,
                         fscales_lineq_constraint_size,
                         fscales_cineq_constraint_size, WorkingSet, TrialState,
                         QRManager, lb, ub);
    if ((!nmpcSine_B.Flags.done) && nmpcSine_B.Flags.stepAccepted) {
      nmpcSine_B.Flags.stepAccepted = false;
      nmpcSine_B.Flags.stepType = 1;
      nmpcSine_B.Flags.failedLineSearch = false;
      nmpcSine_B.d_ix = (nmpcSine_B.mFixed_d + TrialState->iNonIneq0) + 119;
      nmpcSine_B.loop_ub_c = WorkingSet->ldA;
      nmpcSine_B.loop_ub_p = static_cast<uint16_T>(nmpcSine_B.nVar_tmp_tmp);
      memcpy(&TrialState->delta_gradLag.data[0], &TrialState->grad.data[0],
             static_cast<uint32_T>(nmpcSine_B.loop_ub_p) * sizeof(real_T));
      nmpcSine_B.y_size_idx_0 = TrialState->delta_gradLag.size[0];
      nmpcSine_B.loop_ub_p = TrialState->delta_gradLag.size[0];
      if (nmpcSine_B.loop_ub_p - 1 >= 0) {
        memcpy(&nmpcSine_B.y_data_c[0], &TrialState->delta_gradLag.data[0],
               static_cast<uint32_T>(nmpcSine_B.loop_ub_p) * sizeof(real_T));
      }

      if (nmpcSine_B.nVar_tmp_tmp >= 1) {
        nmpcSine_B.loop_ub_p = (nmpcSine_B.nVar_tmp_tmp / 2) << 1;
        nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p - 2;
        for (nmpcSine_B.k = 0; nmpcSine_B.k <= nmpcSine_B.vectorUB; nmpcSine_B.k
             += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.y_data_c[nmpcSine_B.k]);
          tmp_0 = _mm_loadu_pd(&TrialState->grad_old.data[nmpcSine_B.k]);
          _mm_storeu_pd(&nmpcSine_B.y_data_c[nmpcSine_B.k], _mm_sub_pd(tmp,
            tmp_0));
        }

        for (nmpcSine_B.k = nmpcSine_B.loop_ub_p; nmpcSine_B.k <
             nmpcSine_B.nVar_tmp_tmp; nmpcSine_B.k++) {
          nmpcSine_B.y_data_c[nmpcSine_B.k] -= TrialState->
            grad_old.data[nmpcSine_B.k];
        }
      }

      nmpcSine_B.ix = nmpcSine_B.mFixed_d;
      if (nmpcSine_B.y_size_idx_0 - 1 >= 0) {
        memcpy(&TrialState->delta_gradLag.data[0], &nmpcSine_B.y_data_c[0],
               static_cast<uint32_T>(nmpcSine_B.y_size_idx_0) * sizeof(real_T));
      }

      nmpcSine_B.y_size_idx_0 = WorkingSet->ldA * 119 + 1;
      for (nmpcSine_B.loop_ub_p = 1; nmpcSine_B.loop_ub_c < 0 ?
           nmpcSine_B.loop_ub_p >= nmpcSine_B.y_size_idx_0 :
           nmpcSine_B.loop_ub_p <= nmpcSine_B.y_size_idx_0; nmpcSine_B.loop_ub_p
           += nmpcSine_B.loop_ub_c) {
        nmpcSine_B.h = (nmpcSine_B.loop_ub_p + nmpcSine_B.nVar_tmp_tmp) - 1;
        for (nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p; nmpcSine_B.vectorUB <=
             nmpcSine_B.h; nmpcSine_B.vectorUB++) {
          nmpcSine_B.k = nmpcSine_B.vectorUB - nmpcSine_B.loop_ub_p;
          TrialState->delta_gradLag.data[nmpcSine_B.k] += WorkingSet->
            Aeq.data[nmpcSine_B.vectorUB - 1] * TrialState->
            lambdasqp.data[nmpcSine_B.ix];
        }

        nmpcSine_B.ix++;
      }

      nmpcSine_B.ix = nmpcSine_B.mFixed_d;
      for (nmpcSine_B.loop_ub_p = 1; nmpcSine_B.loop_ub_c < 0 ?
           nmpcSine_B.loop_ub_p >= nmpcSine_B.y_size_idx_0 :
           nmpcSine_B.loop_ub_p <= nmpcSine_B.y_size_idx_0; nmpcSine_B.loop_ub_p
           += nmpcSine_B.loop_ub_c) {
        nmpcSine_B.h = (nmpcSine_B.loop_ub_p + nmpcSine_B.nVar_tmp_tmp) - 1;
        for (nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p; nmpcSine_B.vectorUB <=
             nmpcSine_B.h; nmpcSine_B.vectorUB++) {
          nmpcSine_B.k = nmpcSine_B.vectorUB - nmpcSine_B.loop_ub_p;
          TrialState->delta_gradLag.data[nmpcSine_B.k] +=
            TrialState->JacCeqTrans_old.data[nmpcSine_B.vectorUB - 1] *
            -TrialState->lambdasqp.data[nmpcSine_B.ix];
        }

        nmpcSine_B.ix++;
      }

      if (TrialState->mNonlinIneq > 0) {
        nmpcSine_B.y_size_idx_0 = (TrialState->iNonIneq0 - 1) * WorkingSet->ldA
          + 1;
        nmpcSine_B.ix = nmpcSine_B.d_ix;
        nmpcSine_B.h_tmp = (TrialState->mNonlinIneq - 1) * WorkingSet->ldA;
        nmpcSine_B.h = nmpcSine_B.h_tmp + nmpcSine_B.y_size_idx_0;
        for (nmpcSine_B.loop_ub_p = nmpcSine_B.y_size_idx_0;
             nmpcSine_B.loop_ub_c < 0 ? nmpcSine_B.loop_ub_p >= nmpcSine_B.h :
             nmpcSine_B.loop_ub_p <= nmpcSine_B.h; nmpcSine_B.loop_ub_p +=
             nmpcSine_B.loop_ub_c) {
          nmpcSine_B.o_l = (nmpcSine_B.loop_ub_p + nmpcSine_B.nVar_tmp_tmp) - 1;
          for (nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p; nmpcSine_B.vectorUB <=
               nmpcSine_B.o_l; nmpcSine_B.vectorUB++) {
            nmpcSine_B.k = nmpcSine_B.vectorUB - nmpcSine_B.loop_ub_p;
            TrialState->delta_gradLag.data[nmpcSine_B.k] +=
              WorkingSet->Aineq.data[nmpcSine_B.vectorUB - 1] *
              TrialState->lambdasqp.data[nmpcSine_B.ix];
          }

          nmpcSine_B.ix++;
        }

        nmpcSine_B.y_size_idx_0 = nmpcSine_B.h_tmp + 1;
        for (nmpcSine_B.loop_ub_p = 1; nmpcSine_B.loop_ub_c < 0 ?
             nmpcSine_B.loop_ub_p >= nmpcSine_B.y_size_idx_0 :
             nmpcSine_B.loop_ub_p <= nmpcSine_B.y_size_idx_0;
             nmpcSine_B.loop_ub_p += nmpcSine_B.loop_ub_c) {
          nmpcSine_B.ix = (nmpcSine_B.loop_ub_p + nmpcSine_B.nVar_tmp_tmp) - 1;
          for (nmpcSine_B.vectorUB = nmpcSine_B.loop_ub_p; nmpcSine_B.vectorUB <=
               nmpcSine_B.ix; nmpcSine_B.vectorUB++) {
            nmpcSine_B.k = nmpcSine_B.vectorUB - nmpcSine_B.loop_ub_p;
            TrialState->delta_gradLag.data[nmpcSine_B.k] +=
              TrialState->JacCineqTrans_old.data[nmpcSine_B.vectorUB - 1] *
              -TrialState->lambdasqp.data[nmpcSine_B.d_ix];
          }

          nmpcSine_B.d_ix++;
        }
      }

      nmpcSine_saveJacobian(TrialState, nmpcSine_B.nVar_tmp_tmp,
                            nmpcSine_B.mIneq_i, WorkingSet->Aineq.data,
                            TrialState->iNonIneq0, WorkingSet->Aeq.data,
                            WorkingSet->ldA);
      nmpcSine_BFGSUpdate(nmpcSine_B.nVar_tmp_tmp, Hessian,
                          TrialState->delta_x.data,
                          TrialState->delta_gradLag.data,
                          memspace->workspace_float.data);
      TrialState->sqpIterations++;
    }
  }
}

// Function for MATLAB Function: '<S24>/NLMPC'
void nmpcSine::nmpcSine_fmincon(const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
  *fun_workspace_runtimedata, const sAc4bxvmmjmxjQV9i3feLrE_nmpcS_T
  *fun_workspace_userdata, const real_T x0[125], const real_T Aineq_data[],
  const real_T bineq_data[], const int32_T bineq_size[1], const real_T lb[125],
  const real_T ub[125], const s_jex761Cl1dvQqVqRqjms8C_nmpc_T
  *nonlcon_workspace_runtimedata, real_T x[125], real_T *fval, real_T *exitflag,
  sttYSJM5GCi2c1Eu0R50efC_nmpcS_T *output)
{
  __m128d tmp;
  boolean_T guard1;
  nmpcSine_c4_mpclib_anonFcn2(nonlcon_workspace_runtimedata->x,
    nonlcon_workspace_runtimedata->OutputMin,
    nonlcon_workspace_runtimedata->OutputMax, x0, nmpcSine_B.Cineq_data,
    nmpcSine_B.Cineq_size, nmpcSine_B.Ceq, nmpcSine_B.JacCineqTrans_data,
    nmpcSine_B.JacCineqTrans_size, nmpcSine_B.JacCeqTrans);
  nmpcSine_B.i_k = nmpcSine_B.Cineq_size[0] * nmpcSine_B.Cineq_size[1];
  nmpcSine_B.mLinIneq_tmp = bineq_size[0];
  nmpcSine_B.mIneq = bineq_size[0] + nmpcSine_B.i_k;
  nmpcSine_B.mConstrMax = (nmpcSine_B.mIneq + nmpcSine_B.mIneq) + 611;
  if (nmpcSine_B.mIneq + 366 >= nmpcSine_B.mConstrMax) {
    nmpcSine_B.maxDims = nmpcSine_B.mIneq + 366;
  } else {
    nmpcSine_B.maxDims = nmpcSine_B.mConstrMax;
  }

  nmpcSine_factoryConstruct(nmpcSine_B.mIneq + 366, nmpcSine_B.mConstrMax,
    nmpcSine_B.mIneq, nmpcSine_B.i_k, &nmpcSine_B.TrialState);
  nmpcSine_B.FcnEvaluator.next.next.next.next.next.b_value = nmpcSine_B.i_k;
  nmpcSine_B.FcnEvaluator.next.next.next.next.next.next.next.b_value.workspace.runtimedata
    = *nonlcon_workspace_runtimedata;
  nmpcSine_B.FcnEvaluator.next.next.next.next.next.next.next.next.b_value.workspace.runtimedata
    = *fun_workspace_runtimedata;
  nmpcSine_B.FcnEvaluator.next.next.next.next.next.next.next.next.b_value.workspace.userdata
    = *fun_workspace_userdata;
  nmpcSine_factoryConstruct_h53b(nmpcSine_B.mIneq + 366,
    nmpcSine_B.QPObjective.grad.size, nmpcSine_B.QPObjective.Hx.size,
    &nmpcSine_B.QPObjective.hasLinear, &nmpcSine_B.QPObjective.nvar,
    &nmpcSine_B.QPObjective.maxVar, &nmpcSine_B.QPObjective.beta,
    &nmpcSine_B.QPObjective.rho, &nmpcSine_B.QPObjective.objtype,
    &nmpcSine_B.QPObjective.prev_objtype, &nmpcSine_B.QPObjective.prev_nvar,
    &nmpcSine_B.QPObjective.prev_hasLinear, &nmpcSine_B.QPObjective.gammaScalar);
  nmpcSine_B.QPObjective.nvar = 125;
  nmpcSine_B.QPObjective.hasLinear = true;
  nmpcSine_B.QPObjective.objtype = 3;
  nmpcSine_B.memspace.workspace_float.size[0] = nmpcSine_B.maxDims;
  nmpcSine_B.memspace.workspace_float.size[1] = nmpcSine_B.mIneq + 366;
  nmpcSine_B.memspace.workspace_int.size[0] = nmpcSine_B.maxDims;
  nmpcSine_B.memspace.workspace_sort.size[0] = nmpcSine_B.maxDims;
  nmpcSine_factoryConstruct_h53bm(nmpcSine_B.mIneq, nmpcSine_B.mIneq + 366,
    nmpcSine_B.mConstrMax, &nmpcSine_B.WorkingSet);
  nmpcSine_B.mLB = 0;
  nmpcSine_B.mUB = 0;
  nmpcSine_B.mFixed = 0;
  for (nmpcSine_B.iw0 = 0; nmpcSine_B.iw0 < 125; nmpcSine_B.iw0++) {
    nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.iw0] = x0[nmpcSine_B.iw0];
    nmpcSine_B.b_c_o = lb[nmpcSine_B.iw0];
    guard1 = false;
    if ((!rtIsInf(nmpcSine_B.b_c_o)) && (!rtIsNaN(nmpcSine_B.b_c_o))) {
      if (fabs(nmpcSine_B.b_c_o - ub[nmpcSine_B.iw0]) < 0.001) {
        nmpcSine_B.mFixed++;
        nmpcSine_B.WorkingSet.indexFixed.data[nmpcSine_B.mFixed - 1] =
          nmpcSine_B.iw0 + 1;
      } else {
        nmpcSine_B.mLB++;
        nmpcSine_B.WorkingSet.indexLB.data[nmpcSine_B.mLB - 1] = nmpcSine_B.iw0
          + 1;
        guard1 = true;
      }
    } else {
      guard1 = true;
    }

    if (guard1) {
      nmpcSine_B.b_c_o = ub[nmpcSine_B.iw0];
      if ((!rtIsInf(nmpcSine_B.b_c_o)) && (!rtIsNaN(nmpcSine_B.b_c_o))) {
        nmpcSine_B.mUB++;
        nmpcSine_B.WorkingSet.indexUB.data[nmpcSine_B.mUB - 1] = nmpcSine_B.iw0
          + 1;
      }
    }
  }

  nmpcSine_B.WorkingSet.mConstrMax = nmpcSine_B.mConstrMax;
  nmpcSine_B.mConstrMax = nmpcSine_B.mIneq + nmpcSine_B.mLB;
  nmpcSine_B.iw0 = ((nmpcSine_B.mConstrMax + nmpcSine_B.mUB) + nmpcSine_B.mFixed)
    + 120;
  nmpcSine_B.WorkingSet.mConstr = nmpcSine_B.iw0;
  nmpcSine_B.WorkingSet.mConstrOrig = nmpcSine_B.iw0;
  nmpcSine_B.WorkingSet.sizes[0] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.sizes[1] = 120;
  nmpcSine_B.WorkingSet.sizes[2] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.sizes[3] = nmpcSine_B.mLB;
  nmpcSine_B.WorkingSet.sizes[4] = nmpcSine_B.mUB;
  for (nmpcSine_B.iEq0 = 0; nmpcSine_B.iEq0 < 5; nmpcSine_B.iEq0++) {
    nmpcSine_B.WorkingSet.sizesNormal[nmpcSine_B.iEq0] =
      nmpcSine_B.WorkingSet.sizes[nmpcSine_B.iEq0];
  }

  nmpcSine_B.WorkingSet.sizesPhaseOne[0] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.sizesPhaseOne[1] = 120;
  nmpcSine_B.WorkingSet.sizesPhaseOne[2] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.sizesPhaseOne[3] = nmpcSine_B.mLB + 1;
  nmpcSine_B.WorkingSet.sizesPhaseOne[4] = nmpcSine_B.mUB;
  nmpcSine_B.WorkingSet.sizesRegularized[0] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.sizesRegularized[1] = 120;
  nmpcSine_B.WorkingSet.sizesRegularized[2] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.sizesRegularized[3] = nmpcSine_B.mConstrMax + 240;
  nmpcSine_B.WorkingSet.sizesRegularized[4] = nmpcSine_B.mUB;
  nmpcSine_B.WorkingSet.sizesRegPhaseOne[0] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.sizesRegPhaseOne[1] = 120;
  nmpcSine_B.WorkingSet.sizesRegPhaseOne[2] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.sizesRegPhaseOne[3] = nmpcSine_B.mConstrMax + 241;
  nmpcSine_B.WorkingSet.sizesRegPhaseOne[4] = nmpcSine_B.mUB;
  nmpcSine_B.WorkingSet.isActiveIdxNormal[0] = 1;
  nmpcSine_B.WorkingSet.isActiveIdxNormal[1] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.isActiveIdxNormal[2] = 120;
  nmpcSine_B.WorkingSet.isActiveIdxNormal[3] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.isActiveIdxNormal[4] = nmpcSine_B.mLB;
  nmpcSine_B.WorkingSet.isActiveIdxNormal[5] = nmpcSine_B.mUB;
  for (nmpcSine_B.iEq0 = 0; nmpcSine_B.iEq0 < 6; nmpcSine_B.iEq0++) {
    nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iEq0] =
      nmpcSine_B.WorkingSet.isActiveIdxNormal[nmpcSine_B.iEq0];
  }

  for (nmpcSine_B.iw0 = 0; nmpcSine_B.iw0 < 5; nmpcSine_B.iw0++) {
    nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iw0 + 1] +=
      nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iw0];
  }

  for (nmpcSine_B.iEq0 = 0; nmpcSine_B.iEq0 < 6; nmpcSine_B.iEq0++) {
    nmpcSine_B.WorkingSet.isActiveIdx[nmpcSine_B.iEq0] =
      nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iEq0];
  }

  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[0] = 1;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[1] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[2] = 120;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[3] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[4] = nmpcSine_B.mLB + 1;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[5] = nmpcSine_B.mUB;
  for (nmpcSine_B.iw0 = 0; nmpcSine_B.iw0 < 5; nmpcSine_B.iw0++) {
    nmpcSine_B.WorkingSet.isActiveIdxNormal[nmpcSine_B.iw0 + 1] +=
      nmpcSine_B.WorkingSet.isActiveIdxNormal[nmpcSine_B.iw0];
    nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iw0 + 1] +=
      nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iw0];
  }

  for (nmpcSine_B.iEq0 = 0; nmpcSine_B.iEq0 < 6; nmpcSine_B.iEq0++) {
    nmpcSine_B.WorkingSet.isActiveIdxPhaseOne[nmpcSine_B.iEq0] =
      nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iEq0];
  }

  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[0] = 1;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[1] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[2] = 120;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[3] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[4] = nmpcSine_B.mConstrMax + 240;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[5] = nmpcSine_B.mUB;
  for (nmpcSine_B.iw0 = 0; nmpcSine_B.iw0 < 5; nmpcSine_B.iw0++) {
    nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iw0 + 1] +=
      nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iw0];
  }

  for (nmpcSine_B.iEq0 = 0; nmpcSine_B.iEq0 < 6; nmpcSine_B.iEq0++) {
    nmpcSine_B.WorkingSet.isActiveIdxRegularized[nmpcSine_B.iEq0] =
      nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.iEq0];
  }

  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[0] = 1;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[1] = nmpcSine_B.mFixed;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[2] = 120;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[3] = nmpcSine_B.mIneq;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[4] = nmpcSine_B.mConstrMax + 241;
  nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[5] = nmpcSine_B.mUB;
  for (nmpcSine_B.mConstrMax = 0; nmpcSine_B.mConstrMax < 5;
       nmpcSine_B.mConstrMax++) {
    nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.mConstrMax + 1] +=
      nmpcSine_B.WorkingSet.isActiveIdxRegPhaseOne[nmpcSine_B.mConstrMax];
  }

  if (nmpcSine_B.mIneq > 0) {
    for (nmpcSine_B.mConstrMax = 0; nmpcSine_B.mConstrMax <
         nmpcSine_B.mLinIneq_tmp; nmpcSine_B.mConstrMax++) {
      for (nmpcSine_B.iw0 = 0; nmpcSine_B.iw0 < 125; nmpcSine_B.iw0++) {
        nmpcSine_B.WorkingSet.Aineq.data[nmpcSine_B.iw0 +
          nmpcSine_B.WorkingSet.ldA * nmpcSine_B.mConstrMax] =
          Aineq_data[nmpcSine_B.mLinIneq_tmp * nmpcSine_B.iw0 +
          nmpcSine_B.mConstrMax];
      }
    }
  }

  nmpcSine_B.mLB = static_cast<uint16_T>(nmpcSine_B.mLB);
  for (nmpcSine_B.mConstrMax = 0; nmpcSine_B.mConstrMax < nmpcSine_B.mLB;
       nmpcSine_B.mConstrMax++) {
    nmpcSine_B.iw0 = nmpcSine_B.WorkingSet.indexLB.data[nmpcSine_B.mConstrMax];
    nmpcSine_B.b_c_o = lb[nmpcSine_B.iw0 - 1];
    if ((nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.iw0 - 1] >= nmpcSine_B.b_c_o)
        || rtIsNaN(nmpcSine_B.b_c_o)) {
    } else {
      nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.iw0 - 1] = nmpcSine_B.b_c_o;
    }
  }

  nmpcSine_B.mUB = static_cast<uint16_T>(nmpcSine_B.mUB);
  for (nmpcSine_B.mConstrMax = 0; nmpcSine_B.mConstrMax < nmpcSine_B.mUB;
       nmpcSine_B.mConstrMax++) {
    nmpcSine_B.iw0 = nmpcSine_B.WorkingSet.indexUB.data[nmpcSine_B.mConstrMax];
    nmpcSine_B.b_c_o = ub[nmpcSine_B.iw0 - 1];
    if ((nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.iw0 - 1] <= nmpcSine_B.b_c_o)
        || rtIsNaN(nmpcSine_B.b_c_o)) {
    } else {
      nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.iw0 - 1] = nmpcSine_B.b_c_o;
    }
  }

  nmpcSine_B.mConstrMax = static_cast<uint16_T>(nmpcSine_B.mFixed);
  for (nmpcSine_B.iw0 = 0; nmpcSine_B.iw0 < nmpcSine_B.mConstrMax;
       nmpcSine_B.iw0++) {
    nmpcSine_B.iEq0 = nmpcSine_B.WorkingSet.indexFixed.data[nmpcSine_B.iw0];
    nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.iEq0 - 1] = ub[nmpcSine_B.iEq0 - 1];
  }

  evalObjAndConstrAndDerivatives(nmpcSine_B.i_k, nonlcon_workspace_runtimedata,
    &nmpcSine_B.FcnEvaluator.next.next.next.next.next.next.next.next.b_value.workspace,
    nmpcSine_B.TrialState.xstarsqp, nmpcSine_B.TrialState.grad.data,
    nmpcSine_B.TrialState.cIneq.data, nmpcSine_B.TrialState.iNonIneq0,
    nmpcSine_B.TrialState.cEq, nmpcSine_B.WorkingSet.Aineq.data,
    nmpcSine_B.TrialState.iNonIneq0, nmpcSine_B.WorkingSet.ldA,
    nmpcSine_B.WorkingSet.Aeq.data, nmpcSine_B.WorkingSet.ldA,
    &nmpcSine_B.TrialState.sqpFval, &nmpcSine_B.iw0);
  nmpcSine_B.TrialState.FunctionEvaluations = 1;
  nmpcSine_B.iw0 = nmpcSine_B.WorkingSet.ldA;
  if (bineq_size[0] > 0) {
    nmpcSine_B.loop_ub = nmpcSine_B.TrialState.cIneq.size[0];
    if (nmpcSine_B.loop_ub - 1 >= 0) {
      memcpy(&nmpcSine_B.y_data_k[0], &nmpcSine_B.TrialState.cIneq.data[0],
             static_cast<uint32_T>(nmpcSine_B.loop_ub) * sizeof(real_T));
    }

    if (nmpcSine_B.mLinIneq_tmp - 1 >= 0) {
      memcpy(&nmpcSine_B.y_data_k[0], &bineq_data[0], static_cast<uint32_T>
             (nmpcSine_B.mLinIneq_tmp) * sizeof(real_T));
    }

    if (nmpcSine_B.loop_ub - 1 >= 0) {
      memcpy(&nmpcSine_B.TrialState.cIneq.data[0], &nmpcSine_B.y_data_k[0],
             static_cast<uint32_T>(nmpcSine_B.loop_ub) * sizeof(real_T));
    }

    nmpcSine_B.iEq0 = (bineq_size[0] / 2) << 1;
    nmpcSine_B.loop_ub = nmpcSine_B.iEq0 - 2;
    for (nmpcSine_B.b_iy = 0; nmpcSine_B.b_iy <= nmpcSine_B.loop_ub;
         nmpcSine_B.b_iy += 2) {
      tmp = _mm_loadu_pd(&nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.b_iy]);
      _mm_storeu_pd(&nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.b_iy],
                    _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
    }

    for (nmpcSine_B.b_iy = nmpcSine_B.iEq0; nmpcSine_B.b_iy <
         nmpcSine_B.mLinIneq_tmp; nmpcSine_B.b_iy++) {
      nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.b_iy] =
        -nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.b_iy];
    }

    nmpcSine_B.iEq0 = 0;
    nmpcSine_B.loop_ub = (bineq_size[0] - 1) * nmpcSine_B.WorkingSet.ldA + 1;
    for (nmpcSine_B.b_iy = 1; nmpcSine_B.iw0 < 0 ? nmpcSine_B.b_iy >=
         nmpcSine_B.loop_ub : nmpcSine_B.b_iy <= nmpcSine_B.loop_ub;
         nmpcSine_B.b_iy += nmpcSine_B.iw0) {
      nmpcSine_B.b_c_o = 0.0;
      for (nmpcSine_B.mLinIneq_tmp = nmpcSine_B.b_iy; nmpcSine_B.mLinIneq_tmp <=
           nmpcSine_B.b_iy + 124; nmpcSine_B.mLinIneq_tmp++) {
        nmpcSine_B.b_c_o +=
          nmpcSine_B.WorkingSet.Aineq.data[nmpcSine_B.mLinIneq_tmp - 1] *
          nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.mLinIneq_tmp -
          nmpcSine_B.b_iy];
      }

      nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.iEq0] += nmpcSine_B.b_c_o;
      nmpcSine_B.iEq0++;
    }
  }

  nmpcSine_B.iw0 = nmpcSine_B.WorkingSet.ldA * nmpcSine_B.mFixed;
  nmpcSine_B.iEq0 = 0;
  for (nmpcSine_B.loop_ub = 0; nmpcSine_B.loop_ub < 120; nmpcSine_B.loop_ub++) {
    nmpcSine_B.b_c_o = -nmpcSine_B.TrialState.cEq[nmpcSine_B.loop_ub];
    nmpcSine_B.WorkingSet.beq[nmpcSine_B.loop_ub] = nmpcSine_B.b_c_o;
    nmpcSine_B.WorkingSet.bwset.data[nmpcSine_B.mFixed + nmpcSine_B.loop_ub] =
      nmpcSine_B.b_c_o;
    memcpy(&nmpcSine_B.WorkingSet.ATwset.data[nmpcSine_B.iw0],
           &nmpcSine_B.WorkingSet.Aeq.data[nmpcSine_B.iEq0], 125U * sizeof
           (real_T));
    nmpcSine_B.iw0 += nmpcSine_B.WorkingSet.ldA;
    nmpcSine_B.iEq0 += nmpcSine_B.WorkingSet.ldA;
  }

  nmpcSine_B.iEq0 = (nmpcSine_B.mIneq / 2) << 1;
  nmpcSine_B.loop_ub = nmpcSine_B.iEq0 - 2;
  for (nmpcSine_B.mFixed = 0; nmpcSine_B.mFixed <= nmpcSine_B.loop_ub;
       nmpcSine_B.mFixed += 2) {
    tmp = _mm_loadu_pd(&nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.mFixed]);
    _mm_storeu_pd(&nmpcSine_B.WorkingSet.bineq.data[nmpcSine_B.mFixed],
                  _mm_mul_pd(tmp, _mm_set1_pd(-1.0)));
  }

  for (nmpcSine_B.mFixed = nmpcSine_B.iEq0; nmpcSine_B.mFixed < nmpcSine_B.mIneq;
       nmpcSine_B.mFixed++) {
    nmpcSine_B.WorkingSet.bineq.data[nmpcSine_B.mFixed] =
      -nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.mFixed];
  }

  for (nmpcSine_B.mFixed = 0; nmpcSine_B.mFixed < nmpcSine_B.mLB;
       nmpcSine_B.mFixed++) {
    nmpcSine_B.WorkingSet.lb.data[nmpcSine_B.WorkingSet.indexLB.data[nmpcSine_B.mFixed]
      - 1] = -lb[nmpcSine_B.WorkingSet.indexLB.data[nmpcSine_B.mFixed] - 1] +
      x0[nmpcSine_B.WorkingSet.indexLB.data[nmpcSine_B.mFixed] - 1];
  }

  for (nmpcSine_B.mFixed = 0; nmpcSine_B.mFixed < nmpcSine_B.mUB;
       nmpcSine_B.mFixed++) {
    nmpcSine_B.WorkingSet.ub.data[nmpcSine_B.WorkingSet.indexUB.data[nmpcSine_B.mFixed]
      - 1] = ub[nmpcSine_B.WorkingSet.indexUB.data[nmpcSine_B.mFixed] - 1] -
      x0[nmpcSine_B.WorkingSet.indexUB.data[nmpcSine_B.mFixed] - 1];
  }

  for (nmpcSine_B.mFixed = 0; nmpcSine_B.mFixed < nmpcSine_B.mConstrMax;
       nmpcSine_B.mFixed++) {
    _mm_storeu_pd(&nmpcSine_B.dv5[0], _mm_sub_pd(_mm_set1_pd
      (ub[nmpcSine_B.WorkingSet.indexFixed.data[nmpcSine_B.mFixed] - 1]),
      _mm_set1_pd(x0[nmpcSine_B.WorkingSet.indexFixed.data[nmpcSine_B.mFixed] -
                  1])));
    nmpcSine_B.WorkingSet.ub.data[nmpcSine_B.WorkingSet.indexFixed.data[nmpcSine_B.mFixed]
      - 1] = nmpcSine_B.dv5[0];
    nmpcSine_B.WorkingSet.bwset.data[nmpcSine_B.mFixed] = nmpcSine_B.dv5[1];
  }

  nmpcSine_initActiveSet(&nmpcSine_B.WorkingSet);
  nmpcSine_B.MeritFunction.initFval = nmpcSine_B.TrialState.sqpFval;
  nmpcSine_B.MeritFunction.penaltyParam = 1.0;
  nmpcSine_B.MeritFunction.threshold = 0.0001;
  nmpcSine_B.MeritFunction.nPenaltyDecreases = 0;
  nmpcSine_B.MeritFunction.linearizedConstrViol = 0.0;
  nmpcSine_B.b_c_o = 0.0;
  for (nmpcSine_B.mFixed = 0; nmpcSine_B.mFixed < 120; nmpcSine_B.mFixed++) {
    nmpcSine_B.b_c_o += fabs(nmpcSine_B.TrialState.cEq[nmpcSine_B.mFixed]);
  }

  nmpcSine_B.MeritFunction.initConstrViolationEq = nmpcSine_B.b_c_o;
  nmpcSine_B.b_c_o = 0.0;
  for (nmpcSine_B.mFixed = 0; nmpcSine_B.mFixed < nmpcSine_B.mIneq;
       nmpcSine_B.mFixed++) {
    nmpcSine_B.scale = nmpcSine_B.TrialState.cIneq.data[nmpcSine_B.mFixed];
    if (nmpcSine_B.scale > 0.0) {
      nmpcSine_B.b_c_o += nmpcSine_B.scale;
    }
  }

  nmpcSine_B.MeritFunction.initConstrViolationIneq = nmpcSine_B.b_c_o;
  nmpcSine_B.MeritFunction.phi = 0.0;
  nmpcSine_B.MeritFunction.phiPrimePlus = 0.0;
  nmpcSine_B.MeritFunction.phiFullStep = 0.0;
  nmpcSine_B.MeritFunction.feasRelativeFactor = 0.0;
  nmpcSine_B.MeritFunction.nlpPrimalFeasError = 0.0;
  nmpcSine_B.MeritFunction.nlpDualFeasError = 0.0;
  nmpcSine_B.MeritFunction.nlpComplError = 0.0;
  nmpcSine_B.MeritFunction.firstOrderOpt = 0.0;
  nmpcSine_factoryConstruct_h5(nmpcSine_B.maxDims, nmpcSine_B.maxDims,
    &nmpcSine_B.QRManager.ldq, nmpcSine_B.QRManager.QR.size,
    nmpcSine_B.QRManager.Q.data, nmpcSine_B.QRManager.Q.size,
    nmpcSine_B.QRManager.jpvt.data, nmpcSine_B.QRManager.jpvt.size,
    &nmpcSine_B.QRManager.mrows, &nmpcSine_B.QRManager.ncols,
    nmpcSine_B.QRManager.tau.size, &nmpcSine_B.QRManager.minRowCol,
    &nmpcSine_B.QRManager.usedPivoting);
  nmpcSine_factoryConstruct_h53(nmpcSine_B.maxDims,
    nmpcSine_B.CholManager.FMat.size, &nmpcSine_B.CholManager.ldm,
    &nmpcSine_B.CholManager.ndims, &nmpcSine_B.CholManager.info,
    &nmpcSine_B.CholManager.scaleFactor, &nmpcSine_B.CholManager.ConvexCheck,
    &nmpcSine_B.CholManager.regTol_, &nmpcSine_B.CholManager.workspace_,
    &nmpcSine_B.CholManager.workspace2_);
  nmpcSine_B.tmp_size_o[0] = bineq_size[0];
  nmpcSine_B.tmp_size_c[0] = nmpcSine_B.i_k;
  nmpcSine_driver(bineq_data, lb, ub, &nmpcSine_B.TrialState,
                  &nmpcSine_B.MeritFunction, &nmpcSine_B.FcnEvaluator,
                  &nmpcSine_B.memspace, &nmpcSine_B.WorkingSet,
                  &nmpcSine_B.QRManager, &nmpcSine_B.CholManager,
                  &nmpcSine_B.QPObjective, nmpcSine_B.tmp_size_o,
                  nmpcSine_B.tmp_size_c, nmpcSine_B.unusedExpr);
  *fval = nmpcSine_B.TrialState.sqpFval;
  *exitflag = nmpcSine_B.TrialState.sqpExitFlag;
  output->iterations = nmpcSine_B.TrialState.sqpIterations;
  output->funcCount = nmpcSine_B.TrialState.FunctionEvaluations;
  output->algorithm[0] = 's';
  output->algorithm[1] = 'q';
  output->algorithm[2] = 'p';
  output->constrviolation = nmpcSine_B.MeritFunction.nlpPrimalFeasError;
  nmpcSine_B.b_c_o = 0.0;
  nmpcSine_B.scale = 3.3121686421112381E-170;
  for (nmpcSine_B.i_k = 0; nmpcSine_B.i_k < 125; nmpcSine_B.i_k++) {
    x[nmpcSine_B.i_k] = nmpcSine_B.TrialState.xstarsqp[nmpcSine_B.i_k];
    nmpcSine_B.absxk = fabs(nmpcSine_B.TrialState.delta_x.data[nmpcSine_B.i_k]);
    if (nmpcSine_B.absxk > nmpcSine_B.scale) {
      nmpcSine_B.t = nmpcSine_B.scale / nmpcSine_B.absxk;
      nmpcSine_B.b_c_o = nmpcSine_B.b_c_o * nmpcSine_B.t * nmpcSine_B.t + 1.0;
      nmpcSine_B.scale = nmpcSine_B.absxk;
    } else {
      nmpcSine_B.t = nmpcSine_B.absxk / nmpcSine_B.scale;
      nmpcSine_B.b_c_o += nmpcSine_B.t * nmpcSine_B.t;
    }
  }

  output->stepsize = nmpcSine_B.scale * sqrt(nmpcSine_B.b_c_o);
  output->lssteplength = nmpcSine_B.TrialState.steplength;
  output->firstorderopt = nmpcSine_B.MeritFunction.firstOrderOpt;
}

void nmpcSine::nmpcSine_stateTransitionFcnDT(const real_T xk[6], const real_T u
  [3], real_T xk1[6])
{
  __m128d tmp;
  real_T k2_tmp;

  // Start for MATLABSystem: '<S22>/MATLAB System'
  //  Discrete State Transition Function
  //  Control input
  //  Time step
  //  Compute RK4 intermediate steps
  //  auto-generated state function of nonlinear grey box
  //  Load the parameters
  //    SystemDynamics
  //     F = SystemDynamics(D11,D22,D33,K,M11,M22,M33,NU,PSI,R,U1,U2,UPSILON,W1,W2,X_G) 
  //     This function was generated by the Symbolic Math Toolbox version 24.2.
  //     15-Oct-2024 00:22:47
  nmpcSine_B.t2_o = cos(xk[2]);
  nmpcSine_B.t3_g = sin(xk[2]);

  //  Compute derivative of states
  _mm_storeu_pd(&nmpcSine_B.k1[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd(xk[4], -xk[4]),
    _mm_set_pd(nmpcSine_B.t2_o, nmpcSine_B.t3_g)), _mm_mul_pd(_mm_set_pd
    (nmpcSine_B.t3_g, nmpcSine_B.t2_o), _mm_set1_pd(xk[3]))));
  nmpcSine_B.k1[2] = xk[5];
  nmpcSine_B.k1[3] = (((((nmpcSine_B.t2_o * 0.0 + nmpcSine_B.t3_g * 0.0) -
    0.80982825268897241 * xk[3]) + 0.079415032681564857 * u[0]) +
                       0.079415032681564857 * u[1]) + 63.0144902620752 * xk[4] *
                      xk[5] * 0.022897834317090127) + xk[5] * xk[5] *
    63.0144902620752 * 0.022897834317090127 * 0.699999982041001;
  nmpcSine_B.k1_tmp = 125.558050313511 * xk[5] * 0.03333333342744111;
  nmpcSine_B.k1_tmp_c = 43.6722524126937 * xk[4] * 0.03333333342744111 * xk[3];
  nmpcSine_B.k1_tmp_cj = 63.0144902620752 * xk[4] * 0.03333333342744111 * xk[3];
  nmpcSine_B.k1_tmp_m = 0.11560777874744173 * u[0] * 0.699999982041001 * 0.415;
  nmpcSine_B.k1_tmp_j = 0.11560777874744173 * u[1] * 0.699999982041001 * 0.415;
  nmpcSine_B.k1[4] = ((((((((((((nmpcSine_B.t2_o * 0.0 - nmpcSine_B.t3_g * 0.0)
    - 74.8832173268989 * xk[4] * 0.015869365852854363) - 74.8832173268989 * xk[4]
    * 0.48999997485740171 * 0.03333333342744111) + nmpcSine_B.k1_tmp *
    0.699999982041001) - 43.6722524126937 * xk[5] * 0.015869365852854363 * xk[3])
    - nmpcSine_B.k1_tmp_m) + nmpcSine_B.k1_tmp_j) - 43.6722524126937 * xk[5] *
    0.48999997485740171 * 0.03333333342744111 * xk[3]) + 63.0144902620752 * xk[5]
    * 0.48999997485740171 * 0.03333333342744111 * xk[3]) - nmpcSine_B.k1_tmp_c *
                        0.699999982041001) + nmpcSine_B.k1_tmp_cj *
                       0.699999982041001) + 63.0144902620752 * nmpcSine_B.t2_o *
                      0.48999997485740171 * 0.03333333342744111 * 0.0) -
    63.0144902620752 * nmpcSine_B.t3_g * 0.48999997485740171 *
    0.03333333342744111 * 0.0;
  nmpcSine_B.k1_tmp_k = 0.11560777874744173 * u[0] * 0.415;
  nmpcSine_B.k1_tmp_mx = 0.11560777874744173 * u[1] * 0.415;
  nmpcSine_B.k1[5] = ((((((74.8832173268989 * xk[4] * 0.03333333342744111 *
    0.699999982041001 + ((-nmpcSine_B.k1_tmp + nmpcSine_B.k1_tmp_k) -
    nmpcSine_B.k1_tmp_mx)) + nmpcSine_B.k1_tmp_c) - nmpcSine_B.k1_tmp_cj) +
                        43.6722524126937 * xk[5] * 0.03333333342744111 * xk[3] *
                        0.699999982041001) - 63.0144902620752 * xk[5] *
                       0.03333333342744111 * xk[3] * 0.699999982041001) -
                      63.0144902620752 * nmpcSine_B.t2_o * 0.03333333342744111 *
                      0.0 * 0.699999982041001) + 63.0144902620752 *
    nmpcSine_B.t3_g * 0.03333333342744111 * 0.0 * 0.699999982041001;

  //  Output Function
  nmpcSine_B.k1_tmp = 0.5 * u[2];
  for (int32_T i = 0; i <= 4; i += 2) {
    // Start for MATLABSystem: '<S22>/MATLAB System'
    tmp = _mm_loadu_pd(&nmpcSine_B.k1[i]);
    _mm_storeu_pd(&nmpcSine_B.x_l[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
      (nmpcSine_B.k1_tmp), tmp), _mm_loadu_pd(&xk[i])));
  }

  // Start for MATLABSystem: '<S22>/MATLAB System'
  //  auto-generated state function of nonlinear grey box
  //  Load the parameters
  //    SystemDynamics
  //     F = SystemDynamics(D11,D22,D33,K,M11,M22,M33,NU,PSI,R,U1,U2,UPSILON,W1,W2,X_G) 
  //     This function was generated by the Symbolic Math Toolbox version 24.2.
  //     15-Oct-2024 00:22:47
  nmpcSine_B.t2_o = cos(nmpcSine_B.x_l[2]);
  nmpcSine_B.t3_g = sin(nmpcSine_B.x_l[2]);

  //  Compute derivative of states
  _mm_storeu_pd(&nmpcSine_B.k2[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
    (nmpcSine_B.x_l[4], -nmpcSine_B.x_l[4]), _mm_set_pd(nmpcSine_B.t2_o,
    nmpcSine_B.t3_g)), _mm_mul_pd(_mm_set_pd(nmpcSine_B.t3_g, nmpcSine_B.t2_o),
    _mm_set1_pd(nmpcSine_B.x_l[3]))));
  nmpcSine_B.k2[2] = nmpcSine_B.x_l[5];
  nmpcSine_B.k2[3] = (((((nmpcSine_B.t2_o * 0.0 + nmpcSine_B.t3_g * 0.0) -
    0.80982825268897241 * nmpcSine_B.x_l[3]) + 0.079415032681564857 * u[0]) +
                       0.079415032681564857 * u[1]) + 63.0144902620752 *
                      nmpcSine_B.x_l[4] * nmpcSine_B.x_l[5] *
                      0.022897834317090127) + nmpcSine_B.x_l[5] *
    nmpcSine_B.x_l[5] * 63.0144902620752 * 0.022897834317090127 *
    0.699999982041001;
  nmpcSine_B.k1_tmp_c = 125.558050313511 * nmpcSine_B.x_l[5] *
    0.03333333342744111;
  nmpcSine_B.k1_tmp_cj = 43.6722524126937 * nmpcSine_B.x_l[4] *
    0.03333333342744111 * nmpcSine_B.x_l[3];
  k2_tmp = 63.0144902620752 * nmpcSine_B.x_l[4] * 0.03333333342744111 *
    nmpcSine_B.x_l[3];
  nmpcSine_B.k2[4] = ((((((((((((nmpcSine_B.t2_o * 0.0 - nmpcSine_B.t3_g * 0.0)
    - 74.8832173268989 * nmpcSine_B.x_l[4] * 0.015869365852854363) -
    74.8832173268989 * nmpcSine_B.x_l[4] * 0.48999997485740171 *
    0.03333333342744111) + nmpcSine_B.k1_tmp_c * 0.699999982041001) -
    43.6722524126937 * nmpcSine_B.x_l[5] * 0.015869365852854363 *
    nmpcSine_B.x_l[3]) - nmpcSine_B.k1_tmp_m) + nmpcSine_B.k1_tmp_j) -
    43.6722524126937 * nmpcSine_B.x_l[5] * 0.48999997485740171 *
    0.03333333342744111 * nmpcSine_B.x_l[3]) + 63.0144902620752 *
    nmpcSine_B.x_l[5] * 0.48999997485740171 * 0.03333333342744111 *
    nmpcSine_B.x_l[3]) - nmpcSine_B.k1_tmp_cj * 0.699999982041001) + k2_tmp *
                       0.699999982041001) + 63.0144902620752 * nmpcSine_B.t2_o *
                      0.48999997485740171 * 0.03333333342744111 * 0.0) -
    63.0144902620752 * nmpcSine_B.t3_g * 0.48999997485740171 *
    0.03333333342744111 * 0.0;
  nmpcSine_B.k2[5] = ((((((74.8832173268989 * nmpcSine_B.x_l[4] *
    0.03333333342744111 * 0.699999982041001 + ((-nmpcSine_B.k1_tmp_c +
    nmpcSine_B.k1_tmp_k) - nmpcSine_B.k1_tmp_mx)) + nmpcSine_B.k1_tmp_cj) -
    k2_tmp) + 43.6722524126937 * nmpcSine_B.x_l[5] * 0.03333333342744111 *
                        nmpcSine_B.x_l[3] * 0.699999982041001) -
                       63.0144902620752 * nmpcSine_B.x_l[5] *
                       0.03333333342744111 * nmpcSine_B.x_l[3] *
                       0.699999982041001) - 63.0144902620752 * nmpcSine_B.t2_o *
                      0.03333333342744111 * 0.0 * 0.699999982041001) +
    63.0144902620752 * nmpcSine_B.t3_g * 0.03333333342744111 * 0.0 *
    0.699999982041001;

  //  Output Function
  for (int32_T i = 0; i <= 4; i += 2) {
    // Start for MATLABSystem: '<S22>/MATLAB System'
    tmp = _mm_loadu_pd(&nmpcSine_B.k2[i]);
    _mm_storeu_pd(&nmpcSine_B.x_l[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd
      (nmpcSine_B.k1_tmp), tmp), _mm_loadu_pd(&xk[i])));
  }

  // Start for MATLABSystem: '<S22>/MATLAB System'
  //  auto-generated state function of nonlinear grey box
  //  Load the parameters
  //    SystemDynamics
  //     F = SystemDynamics(D11,D22,D33,K,M11,M22,M33,NU,PSI,R,U1,U2,UPSILON,W1,W2,X_G) 
  //     This function was generated by the Symbolic Math Toolbox version 24.2.
  //     15-Oct-2024 00:22:47
  nmpcSine_B.t2_o = cos(nmpcSine_B.x_l[2]);
  nmpcSine_B.t3_g = sin(nmpcSine_B.x_l[2]);

  //  Compute derivative of states
  _mm_storeu_pd(&nmpcSine_B.k3[0], _mm_add_pd(_mm_mul_pd(_mm_set_pd
    (nmpcSine_B.x_l[4], -nmpcSine_B.x_l[4]), _mm_set_pd(nmpcSine_B.t2_o,
    nmpcSine_B.t3_g)), _mm_mul_pd(_mm_set_pd(nmpcSine_B.t3_g, nmpcSine_B.t2_o),
    _mm_set1_pd(nmpcSine_B.x_l[3]))));
  nmpcSine_B.k3[2] = nmpcSine_B.x_l[5];
  nmpcSine_B.k3[3] = (((((nmpcSine_B.t2_o * 0.0 + nmpcSine_B.t3_g * 0.0) -
    0.80982825268897241 * nmpcSine_B.x_l[3]) + 0.079415032681564857 * u[0]) +
                       0.079415032681564857 * u[1]) + 63.0144902620752 *
                      nmpcSine_B.x_l[4] * nmpcSine_B.x_l[5] *
                      0.022897834317090127) + nmpcSine_B.x_l[5] *
    nmpcSine_B.x_l[5] * 63.0144902620752 * 0.022897834317090127 *
    0.699999982041001;
  nmpcSine_B.k1_tmp = 125.558050313511 * nmpcSine_B.x_l[5] * 0.03333333342744111;
  nmpcSine_B.k1_tmp_c = 43.6722524126937 * nmpcSine_B.x_l[4] *
    0.03333333342744111 * nmpcSine_B.x_l[3];
  nmpcSine_B.k1_tmp_cj = 63.0144902620752 * nmpcSine_B.x_l[4] *
    0.03333333342744111 * nmpcSine_B.x_l[3];
  nmpcSine_B.k3[4] = ((((((((((((nmpcSine_B.t2_o * 0.0 - nmpcSine_B.t3_g * 0.0)
    - 74.8832173268989 * nmpcSine_B.x_l[4] * 0.015869365852854363) -
    74.8832173268989 * nmpcSine_B.x_l[4] * 0.48999997485740171 *
    0.03333333342744111) + nmpcSine_B.k1_tmp * 0.699999982041001) -
    43.6722524126937 * nmpcSine_B.x_l[5] * 0.015869365852854363 *
    nmpcSine_B.x_l[3]) - nmpcSine_B.k1_tmp_m) + nmpcSine_B.k1_tmp_j) -
    43.6722524126937 * nmpcSine_B.x_l[5] * 0.48999997485740171 *
    0.03333333342744111 * nmpcSine_B.x_l[3]) + 63.0144902620752 *
    nmpcSine_B.x_l[5] * 0.48999997485740171 * 0.03333333342744111 *
    nmpcSine_B.x_l[3]) - nmpcSine_B.k1_tmp_c * 0.699999982041001) +
                       nmpcSine_B.k1_tmp_cj * 0.699999982041001) +
                      63.0144902620752 * nmpcSine_B.t2_o * 0.48999997485740171 *
                      0.03333333342744111 * 0.0) - 63.0144902620752 *
    nmpcSine_B.t3_g * 0.48999997485740171 * 0.03333333342744111 * 0.0;
  nmpcSine_B.k3[5] = ((((((74.8832173268989 * nmpcSine_B.x_l[4] *
    0.03333333342744111 * 0.699999982041001 + ((-nmpcSine_B.k1_tmp +
    nmpcSine_B.k1_tmp_k) - nmpcSine_B.k1_tmp_mx)) + nmpcSine_B.k1_tmp_c) -
    nmpcSine_B.k1_tmp_cj) + 43.6722524126937 * nmpcSine_B.x_l[5] *
                        0.03333333342744111 * nmpcSine_B.x_l[3] *
                        0.699999982041001) - 63.0144902620752 * nmpcSine_B.x_l[5]
                       * 0.03333333342744111 * nmpcSine_B.x_l[3] *
                       0.699999982041001) - 63.0144902620752 * nmpcSine_B.t2_o *
                      0.03333333342744111 * 0.0 * 0.699999982041001) +
    63.0144902620752 * nmpcSine_B.t3_g * 0.03333333342744111 * 0.0 *
    0.699999982041001;

  //  Output Function
  for (int32_T i = 0; i <= 4; i += 2) {
    // Start for MATLABSystem: '<S22>/MATLAB System'
    tmp = _mm_loadu_pd(&nmpcSine_B.k3[i]);
    _mm_storeu_pd(&nmpcSine_B.x_l[i], _mm_add_pd(_mm_mul_pd(_mm_set1_pd(u[2]),
      tmp), _mm_loadu_pd(&xk[i])));
  }

  // Start for MATLABSystem: '<S22>/MATLAB System'
  //  auto-generated state function of nonlinear grey box
  //  Load the parameters
  //    SystemDynamics
  //     F = SystemDynamics(D11,D22,D33,K,M11,M22,M33,NU,PSI,R,U1,U2,UPSILON,W1,W2,X_G) 
  //     This function was generated by the Symbolic Math Toolbox version 24.2.
  //     15-Oct-2024 00:22:47
  nmpcSine_B.t2_o = cos(nmpcSine_B.x_l[2]);
  nmpcSine_B.t3_g = sin(nmpcSine_B.x_l[2]);

  //  Compute derivative of states
  //  Output Function
  //  Compute next state using RK4 formula
  nmpcSine_B.k1_tmp = u[2] / 6.0;
  tmp = _mm_set1_pd(2.0);
  tmp = _mm_add_pd(_mm_mul_pd(_mm_add_pd(_mm_add_pd(_mm_add_pd(_mm_mul_pd(tmp,
    _mm_loadu_pd(&nmpcSine_B.k2[0])), _mm_loadu_pd(&nmpcSine_B.k1[0])),
    _mm_mul_pd(tmp, _mm_loadu_pd(&nmpcSine_B.k3[0]))), _mm_add_pd(_mm_mul_pd
    (_mm_set_pd(nmpcSine_B.x_l[4], -nmpcSine_B.x_l[4]), _mm_set_pd
     (nmpcSine_B.t2_o, nmpcSine_B.t3_g)), _mm_mul_pd(_mm_set_pd(nmpcSine_B.t3_g,
    nmpcSine_B.t2_o), _mm_set1_pd(nmpcSine_B.x_l[3])))), _mm_set1_pd
    (nmpcSine_B.k1_tmp)), _mm_loadu_pd(&xk[0]));
  _mm_storeu_pd(&xk1[0], tmp);

  // Start for MATLABSystem: '<S22>/MATLAB System'
  xk1[2] = (((2.0 * nmpcSine_B.k2[2] + xk[5]) + 2.0 * nmpcSine_B.k3[2]) +
            nmpcSine_B.x_l[5]) * nmpcSine_B.k1_tmp + xk[2];
  xk1[3] = (((((((nmpcSine_B.t2_o * 0.0 + nmpcSine_B.t3_g * 0.0) -
                 0.80982825268897241 * nmpcSine_B.x_l[3]) + 0.079415032681564857
                * u[0]) + 0.079415032681564857 * u[1]) + 63.0144902620752 *
              nmpcSine_B.x_l[4] * nmpcSine_B.x_l[5] * 0.022897834317090127) +
             nmpcSine_B.x_l[5] * nmpcSine_B.x_l[5] * 63.0144902620752 *
             0.022897834317090127 * 0.699999982041001) + ((2.0 * nmpcSine_B.k2[3]
              + nmpcSine_B.k1[3]) + 2.0 * nmpcSine_B.k3[3])) * nmpcSine_B.k1_tmp
    + xk[3];
  nmpcSine_B.k1_tmp_c = 125.558050313511 * nmpcSine_B.x_l[5] *
    0.03333333342744111;
  nmpcSine_B.k1_tmp_cj = 43.6722524126937 * nmpcSine_B.x_l[4] *
    0.03333333342744111 * nmpcSine_B.x_l[3];
  k2_tmp = 63.0144902620752 * nmpcSine_B.x_l[4] * 0.03333333342744111 *
    nmpcSine_B.x_l[3];
  xk1[4] = ((((((((((((((nmpcSine_B.t2_o * 0.0 - nmpcSine_B.t3_g * 0.0) -
                        74.8832173268989 * nmpcSine_B.x_l[4] *
                        0.015869365852854363) - 74.8832173268989 *
                       nmpcSine_B.x_l[4] * 0.48999997485740171 *
                       0.03333333342744111) + nmpcSine_B.k1_tmp_c *
                      0.699999982041001) - 43.6722524126937 * nmpcSine_B.x_l[5] *
                     0.015869365852854363 * nmpcSine_B.x_l[3]) -
                    nmpcSine_B.k1_tmp_m) + nmpcSine_B.k1_tmp_j) -
                  43.6722524126937 * nmpcSine_B.x_l[5] * 0.48999997485740171 *
                  0.03333333342744111 * nmpcSine_B.x_l[3]) + 63.0144902620752 *
                 nmpcSine_B.x_l[5] * 0.48999997485740171 * 0.03333333342744111 *
                 nmpcSine_B.x_l[3]) - nmpcSine_B.k1_tmp_cj * 0.699999982041001)
               + k2_tmp * 0.699999982041001) + 63.0144902620752 *
              nmpcSine_B.t2_o * 0.48999997485740171 * 0.03333333342744111 * 0.0)
             - 63.0144902620752 * nmpcSine_B.t3_g * 0.48999997485740171 *
             0.03333333342744111 * 0.0) + ((2.0 * nmpcSine_B.k2[4] +
              nmpcSine_B.k1[4]) + 2.0 * nmpcSine_B.k3[4])) * nmpcSine_B.k1_tmp +
    xk[4];
  xk1[5] = ((((((((74.8832173268989 * nmpcSine_B.x_l[4] * 0.03333333342744111 *
                   0.699999982041001 + ((-nmpcSine_B.k1_tmp_c +
    nmpcSine_B.k1_tmp_k) - nmpcSine_B.k1_tmp_mx)) + nmpcSine_B.k1_tmp_cj) -
                 k2_tmp) + 43.6722524126937 * nmpcSine_B.x_l[5] *
                0.03333333342744111 * nmpcSine_B.x_l[3] * 0.699999982041001) -
               63.0144902620752 * nmpcSine_B.x_l[5] * 0.03333333342744111 *
               nmpcSine_B.x_l[3] * 0.699999982041001) - 63.0144902620752 *
              nmpcSine_B.t2_o * 0.03333333342744111 * 0.0 * 0.699999982041001) +
             63.0144902620752 * nmpcSine_B.t3_g * 0.03333333342744111 * 0.0 *
             0.699999982041001) + ((2.0 * nmpcSine_B.k2[5] + nmpcSine_B.k1[5]) +
             2.0 * nmpcSine_B.k3[5])) * nmpcSine_B.k1_tmp + xk[5];
}

real_T nmpcSine::nmpcSine_xnrm2_m32uc(int32_T n, const real_T x[72], int32_T ix0)
{
  real_T y;

  // Start for MATLABSystem: '<S22>/MATLAB System'
  y = 0.0;
  if (n >= 1) {
    if (n == 1) {
      y = fabs(x[ix0 - 1]);
    } else {
      real_T scale;
      int32_T kend;
      scale = 3.3121686421112381E-170;
      kend = ix0 + n;
      for (int32_T k = ix0; k < kend; k++) {
        real_T absxk;
        absxk = fabs(x[k - 1]);
        if (absxk > scale) {
          real_T t;
          t = scale / absxk;
          y = y * t * t + 1.0;
          scale = absxk;
        } else {
          real_T t;
          t = absxk / scale;
          y += t * t;
        }
      }

      y = scale * sqrt(y);
    }
  }

  // End of Start for MATLABSystem: '<S22>/MATLAB System'
  return y;
}

void nmpcSine::nmpcSine_qr_m32(const real_T A[72], real_T Q[72], real_T R[36])
{
  __m128d tmp;
  real_T b_tau[6];
  real_T work[6];
  real_T atmp;
  real_T beta1;
  real_T c;
  int32_T d;
  int32_T exitg1;
  int32_T i;
  int32_T iac;
  int32_T iaii;
  int32_T ii;
  int32_T jA;
  int32_T knt;
  int32_T lastc;
  int32_T lastv;
  int32_T scalarLB;
  boolean_T exitg2;
  for (i = 0; i < 6; i++) {
    // Start for MATLABSystem: '<S22>/MATLAB System'
    b_tau[i] = 0.0;
  }

  // Start for MATLABSystem: '<S22>/MATLAB System'
  memcpy(&Q[0], &A[0], 72U * sizeof(real_T));
  for (i = 0; i < 6; i++) {
    // Start for MATLABSystem: '<S22>/MATLAB System'
    work[i] = 0.0;
  }

  // Start for MATLABSystem: '<S22>/MATLAB System'
  for (iaii = 0; iaii < 6; iaii++) {
    ii = iaii * 12 + iaii;
    i = ii + 2;
    atmp = Q[ii];
    b_tau[iaii] = 0.0;
    beta1 = nmpcSine_xnrm2_m32uc(11 - iaii, Q, ii + 2);
    if (beta1 != 0.0) {
      c = Q[ii];
      beta1 = nmpcSine_rt_hypotd_snf(c, beta1);
      if (c >= 0.0) {
        beta1 = -beta1;
      }

      if (fabs(beta1) < 1.0020841800044864E-292) {
        knt = -1;
        do {
          knt++;
          scalarLB = ii - iaii;
          lastc = (((((scalarLB - ii) + 11) / 2) << 1) + ii) + 2;
          jA = lastc - 2;
          for (lastv = i; lastv <= jA; lastv += 2) {
            tmp = _mm_loadu_pd(&Q[lastv - 1]);
            _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd
              (9.9792015476736E+291)));
          }

          for (lastv = lastc; lastv <= scalarLB + 12; lastv++) {
            Q[lastv - 1] *= 9.9792015476736E+291;
          }

          beta1 *= 9.9792015476736E+291;
          atmp *= 9.9792015476736E+291;
        } while ((fabs(beta1) < 1.0020841800044864E-292) && (knt + 1 < 20));

        beta1 = nmpcSine_rt_hypotd_snf(atmp, nmpcSine_xnrm2_m32uc(11 - iaii, Q,
          ii + 2));
        if (atmp >= 0.0) {
          beta1 = -beta1;
        }

        b_tau[iaii] = (beta1 - atmp) / beta1;
        atmp = 1.0 / (atmp - beta1);
        for (lastv = i; lastv <= jA; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = lastc; lastv <= scalarLB + 12; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        for (lastv = 0; lastv <= knt; lastv++) {
          beta1 *= 1.0020841800044864E-292;
        }

        atmp = beta1;
      } else {
        b_tau[iaii] = (beta1 - c) / beta1;
        atmp = 1.0 / (c - beta1);
        knt = ii - iaii;
        scalarLB = (((((knt - ii) + 11) / 2) << 1) + ii) + 2;
        lastc = scalarLB - 2;
        for (lastv = i; lastv <= lastc; lastv += 2) {
          tmp = _mm_loadu_pd(&Q[lastv - 1]);
          _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(atmp)));
        }

        for (lastv = scalarLB; lastv <= knt + 12; lastv++) {
          Q[lastv - 1] *= atmp;
        }

        atmp = beta1;
      }
    }

    Q[ii] = atmp;
    if (iaii + 1 < 6) {
      Q[ii] = 1.0;
      scalarLB = ii + 13;
      if (b_tau[iaii] != 0.0) {
        lastv = 12 - iaii;
        i = ii - iaii;
        while ((lastv > 0) && (Q[i + 11] == 0.0)) {
          lastv--;
          i--;
        }

        knt = 5 - iaii;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          lastc = (knt - 1) * 12 + ii;
          jA = lastc + 13;
          do {
            exitg1 = 0;
            if (jA <= (lastc + lastv) + 12) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            memset(&work[0], 0, static_cast<uint32_T>(lastc + 1) * sizeof(real_T));
          }

          knt = (12 * lastc + ii) + 13;
          for (iac = scalarLB; iac <= knt; iac += 12) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[(ii + jA) - iac] * Q[jA - 1];
            }

            jA = div_nde_s32_floor((iac - ii) - 13, 12);
            work[jA] += c;
          }
        }

        if (!(-b_tau[iaii] == 0.0)) {
          jA = ii;
          for (iac = 0; iac <= lastc; iac++) {
            c = work[iac];
            if (c != 0.0) {
              c *= -b_tau[iaii];
              knt = jA + 13;
              scalarLB = (lastv + jA) + 12;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((ii + d) - jA) - 13] * c;
              }
            }

            jA += 12;
          }
        }
      }

      Q[ii] = atmp;
    }
  }

  for (ii = 0; ii < 6; ii++) {
    for (iaii = 0; iaii <= ii; iaii++) {
      // Start for MATLABSystem: '<S22>/MATLAB System'
      R[iaii + 6 * ii] = Q[12 * ii + iaii];
    }

    for (iaii = ii + 2; iaii < 7; iaii++) {
      R[(iaii + 6 * ii) - 1] = 0.0;
    }

    // Start for MATLABSystem: '<S22>/MATLAB System'
    work[ii] = 0.0;
  }

  // Start for MATLABSystem: '<S22>/MATLAB System'
  for (i = 5; i >= 0; i--) {
    iaii = (i * 12 + i) + 12;
    if (i + 1 < 6) {
      Q[iaii - 12] = 1.0;
      scalarLB = iaii + 1;
      if (b_tau[i] != 0.0) {
        lastv = 12 - i;
        ii = (iaii - i) - 1;
        while ((lastv > 0) && (Q[ii] == 0.0)) {
          lastv--;
          ii--;
        }

        knt = 5 - i;
        exitg2 = false;
        while ((!exitg2) && (knt > 0)) {
          lastc = (knt - 1) * 12 + iaii;
          jA = lastc + 1;
          do {
            exitg1 = 0;
            if (jA <= lastc + lastv) {
              if (Q[jA - 1] != 0.0) {
                exitg1 = 1;
              } else {
                jA++;
              }
            } else {
              knt--;
              exitg1 = 2;
            }
          } while (exitg1 == 0);

          if (exitg1 == 1) {
            exitg2 = true;
          }
        }

        lastc = knt - 1;
      } else {
        lastv = 0;
        lastc = -1;
      }

      if (lastv > 0) {
        if (lastc + 1 != 0) {
          if (lastc >= 0) {
            memset(&work[0], 0, static_cast<uint32_T>(lastc + 1) * sizeof(real_T));
          }

          knt = (12 * lastc + iaii) + 1;
          for (iac = scalarLB; iac <= knt; iac += 12) {
            c = 0.0;
            d = iac + lastv;
            for (jA = iac; jA < d; jA++) {
              c += Q[((iaii + jA) - iac) - 12] * Q[jA - 1];
            }

            jA = div_nde_s32_floor((iac - iaii) - 1, 12);
            work[jA] += c;
          }
        }

        if (!(-b_tau[i] == 0.0)) {
          jA = iaii;
          for (iac = 0; iac <= lastc; iac++) {
            c = work[iac];
            if (c != 0.0) {
              c *= -b_tau[i];
              knt = jA + 1;
              scalarLB = lastv + jA;
              for (d = knt; d <= scalarLB; d++) {
                Q[d - 1] += Q[((iaii + d) - jA) - 13] * c;
              }
            }

            jA += 12;
          }
        }
      }
    }

    knt = iaii - i;
    scalarLB = (((((knt - iaii) + 11) / 2) << 1) + iaii) - 10;
    lastc = scalarLB - 2;
    for (lastv = iaii - 10; lastv <= lastc; lastv += 2) {
      tmp = _mm_loadu_pd(&Q[lastv - 1]);
      _mm_storeu_pd(&Q[lastv - 1], _mm_mul_pd(tmp, _mm_set1_pd(-b_tau[i])));
    }

    for (lastv = scalarLB; lastv <= knt; lastv++) {
      Q[lastv - 1] *= -b_tau[i];
    }

    Q[iaii - 12] = 1.0 - b_tau[i];
    for (ii = 0; ii < i; ii++) {
      Q[(iaii - ii) - 13] = 0.0;
    }
  }
}

void nmpcSine::nmpcSine_Publisher_setupImpl_m(const
  ros_slros2_internal_block_Pub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[23];
  static const char_T b_zeroDelimTopic_0[23] = "/controller/opt_status";
  qos_profile = rmw_qos_profile_default;

  // Start for MATLABSystem: '<S48>/SinkBlock'
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 23; i++) {
    // Start for MATLABSystem: '<S48>/SinkBlock'
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Pub_nmpcSine_264.createPublisher(&b_zeroDelimTopic[0], qos_profile);
}

void nmpcSine::nmpcSine_Publisher_setupImpl_m3(const
  ros_slros2_internal_block_Pub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[21];
  static const char_T b_zeroDelimTopic_0[21] = "/controller/ref_path";
  qos_profile = rmw_qos_profile_default;

  // Start for MATLABSystem: '<S50>/SinkBlock'
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 21; i++) {
    // Start for MATLABSystem: '<S50>/SinkBlock'
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Pub_nmpcSine_305.createPublisher(&b_zeroDelimTopic[0], qos_profile);
}

void nmpcSine::nmpcSin_Publisher_setupImpl_m32(const
  ros_slros2_internal_block_Pub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[19];
  static const char_T b_zeroDelimTopic_0[19] = "/controller/states";
  qos_profile = rmw_qos_profile_default;

  // Start for MATLABSystem: '<S52>/SinkBlock'
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 19; i++) {
    // Start for MATLABSystem: '<S52>/SinkBlock'
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Pub_nmpcSine_253.createPublisher(&b_zeroDelimTopic[0], qos_profile);
}

void nmpcSine::nmpcSine_Subscriber_setupImpl(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[23];
  static const char_T b_zeroDelimTopic_0[23] = "/sensors/emlid_gps_fix";
  qos_profile = rmw_qos_profile_default;

  // Start for MATLABSystem: '<S4>/SourceBlock'
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 23; i++) {
    // Start for MATLABSystem: '<S4>/SourceBlock'
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Sub_nmpcSine_272.createSubscriber(&b_zeroDelimTopic[0], qos_profile);
}

void nmpcSine::nmpcSine_Subscriber_setupImpl_m(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[21];
  static const char_T b_zeroDelimTopic_0[21] = "/sensors/mag_heading";
  qos_profile = rmw_qos_profile_default;

  // Start for MATLABSystem: '<S5>/SourceBlock'
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 21; i++) {
    // Start for MATLABSystem: '<S5>/SourceBlock'
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Sub_nmpcSine_273.createSubscriber(&b_zeroDelimTopic[0], qos_profile);
}

void nmpcSine::nmpcSin_Subscriber_setupImpl_m3(const
  ros_slros2_internal_block_Sub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[20];
  static const char_T b_zeroDelimTopic_0[20] = "/sensors/ouster/imu";
  qos_profile = rmw_qos_profile_default;

  // Start for MATLABSystem: '<S6>/SourceBlock'
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 20; i++) {
    // Start for MATLABSystem: '<S6>/SourceBlock'
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Sub_nmpcSine_274.createSubscriber(&b_zeroDelimTopic[0], qos_profile);
}

void nmpcSine::nmpcSine_Publisher_setupImpl(const
  ros_slros2_internal_block_Pub_T *obj)
{
  rmw_qos_profile_t qos_profile;
  sJ4ih70VmKcvCeguWN0mNVF deadline;
  sJ4ih70VmKcvCeguWN0mNVF lifespan;
  sJ4ih70VmKcvCeguWN0mNVF liveliness_lease_duration;
  char_T b_zeroDelimTopic[18];
  static const char_T b_zeroDelimTopic_0[18] = "/platform/cmd_vel";
  qos_profile = rmw_qos_profile_default;

  // Start for MATLABSystem: '<S9>/SinkBlock'
  deadline.sec = 0.0;
  deadline.nsec = 0.0;
  lifespan.sec = 0.0;
  lifespan.nsec = 0.0;
  liveliness_lease_duration.sec = 0.0;
  liveliness_lease_duration.nsec = 0.0;
  SET_QOS_VALUES(qos_profile, RMW_QOS_POLICY_HISTORY_KEEP_LAST, (size_t)1.0,
                 RMW_QOS_POLICY_DURABILITY_VOLATILE,
                 RMW_QOS_POLICY_RELIABILITY_BEST_EFFORT, deadline, lifespan,
                 RMW_QOS_POLICY_LIVELINESS_AUTOMATIC, liveliness_lease_duration,
                 (bool)obj->QOSAvoidROSNamespaceConventions);
  for (int32_T i = 0; i < 18; i++) {
    // Start for MATLABSystem: '<S9>/SinkBlock'
    b_zeroDelimTopic[i] = b_zeroDelimTopic_0[i];
  }

  Pub_nmpcSine_168.createPublisher(&b_zeroDelimTopic[0], qos_profile);
}

// Model step function
void nmpcSine::step()
{
  __m128d tmp;
  __m128d tmp_0;
  static const int8_T a[160] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  static const int8_T d[80] = { 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 };

  if ((&nmpcSine_M)->Timing.TaskCounters.TID[2] == 0) {
    // MATLABSystem: '<S4>/SourceBlock'
    nmpcSine_B.SourceBlock_o1_o = Sub_nmpcSine_272.getLatestMessage
      (&nmpcSine_B.rtb_SourceBlock_o2_e_o);

    // Outputs for Enabled SubSystem: '<S4>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S53>/Enable'

    if (nmpcSine_B.SourceBlock_o1_o) {
      // SignalConversion generated from: '<S53>/In1' incorporates:
      //   MATLABSystem: '<S4>/SourceBlock'

      nmpcSine_B.In1_j = nmpcSine_B.rtb_SourceBlock_o2_e_o;
    }

    // End of Outputs for SubSystem: '<S4>/Enabled Subsystem'
  }

  nmpcSine_B.unnamed_idx_0 = ((&nmpcSine_M)->Timing.TaskCounters.TID[1] == 0);
  if (nmpcSine_B.unnamed_idx_0) {
    // MATLABSystem: '<S5>/SourceBlock'
    nmpcSine_B.SourceBlock_o1 = Sub_nmpcSine_273.getLatestMessage
      (&nmpcSine_B.SourceBlock_o2_b);
  }

  // MATLABSystem: '<S6>/SourceBlock'
  nmpcSine_B.b_varargout_1 = Sub_nmpcSine_274.getLatestMessage
    (&nmpcSine_B.rtb_SourceBlock_o2_h);

  // MATLABSystem: '<Root>/Get Parameter'
  ParamGet_nmpcSine_289.getParameter(3U,
    &nmpcSine_B.TmpSignalConversionAtMATLAB[0], &nmpcSine_B.len);

  // Chart: '<Root>/Control Logic' incorporates:
  //   MATLABSystem: '<Root>/Get Parameter'
  //   MATLABSystem: '<S6>/SourceBlock'
  //   SignalConversion generated from: '<S2>/ SFunction '
  //
  if (nmpcSine_DW.temporalCounter_i1 < 2047) {
    nmpcSine_DW.temporalCounter_i1 = static_cast<uint16_T>
      (nmpcSine_DW.temporalCounter_i1 + 1);
  }

  if (nmpcSine_DW.is_active_c4_nmpcSine == 0) {
    nmpcSine_DW.is_active_c4_nmpcSine = 1U;
    nmpcSine_DW.is_c4_nmpcSine = nmpcSine_IN_Wait;
    nmpcSine_B.controlOn = 0;
    nmpcSine_DW.en1 = 0.0;
    nmpcSine_DW.en2 = 0.0;
    nmpcSine_DW.en3 = 0.0;
    nmpcSine_B.gps0[0] = 0.0;
    nmpcSine_B.gps0[1] = 0.0;
    nmpcSine_B.gps0[2] = 0.0;
  } else {
    switch (nmpcSine_DW.is_c4_nmpcSine) {
     case nmpcSine_IN_Control:
      nmpcSine_B.controlOn = 1;
      if (nmpcSine_DW.temporalCounter_i1 >= 1500) {
        nmpcSine_DW.is_c4_nmpcSine = nmpcSine_IN_Stop;
        nmpcSine_B.controlOn = 0;
      }
      break;

     case nmpcSine_IN_Stop:
      nmpcSine_B.controlOn = 0;
      break;

     default:
      // case IN_Wait:
      nmpcSine_B.controlOn = 0;
      if ((nmpcSine_DW.en1 != 0.0) && (nmpcSine_DW.en2 != 0.0) &&
          (nmpcSine_DW.en3 != 0.0) && (nmpcSine_B.len == 3U)) {
        nmpcSine_DW.temporalCounter_i1 = 0U;
        nmpcSine_DW.is_c4_nmpcSine = nmpcSine_IN_Control;
        nmpcSine_B.gps0[0] = nmpcSine_B.In1_j.latitude;
        nmpcSine_B.gps0[1] = nmpcSine_B.In1_j.longitude;
        nmpcSine_B.gps0[2] = nmpcSine_B.In1_j.altitude;
        nmpcSine_B.params[0] = nmpcSine_B.TmpSignalConversionAtMATLAB[0];
        nmpcSine_B.params[1] = nmpcSine_B.TmpSignalConversionAtMATLAB[1];
        nmpcSine_B.params[2] = nmpcSine_B.TmpSignalConversionAtMATLAB[2];
        nmpcSine_B.controlOn = 1;
      } else {
        if (nmpcSine_B.SourceBlock_o1_o) {
          nmpcSine_DW.en1 = 1.0;
        }

        if (nmpcSine_B.b_varargout_1) {
          nmpcSine_DW.en2 = 1.0;
        }

        if (nmpcSine_B.SourceBlock_o1) {
          nmpcSine_DW.en3 = 1.0;
        }
      }
      break;
    }
  }

  // End of Chart: '<Root>/Control Logic'
  if (nmpcSine_B.unnamed_idx_0) {
    // Outputs for Enabled SubSystem: '<S5>/Enabled Subsystem' incorporates:
    //   EnablePort: '<S54>/Enable'

    if (nmpcSine_B.SourceBlock_o1) {
      // SignalConversion generated from: '<S54>/In1'
      nmpcSine_B.In1_c = nmpcSine_B.SourceBlock_o2_b;
    }

    // End of Outputs for SubSystem: '<S5>/Enabled Subsystem'
  }

  // Outputs for Enabled SubSystem: '<S6>/Enabled Subsystem' incorporates:
  //   EnablePort: '<S55>/Enable'

  // Start for MATLABSystem: '<S6>/SourceBlock'
  if (nmpcSine_B.b_varargout_1) {
    // SignalConversion generated from: '<S55>/In1'
    nmpcSine_B.In1 = nmpcSine_B.rtb_SourceBlock_o2_h;
  }

  // End of Outputs for SubSystem: '<S6>/Enabled Subsystem'

  // Outputs for Enabled SubSystem: '<Root>/NMPC Controller' incorporates:
  //   EnablePort: '<S3>/Enable'

  if (nmpcSine_B.controlOn > 0) {
    nmpcSine_B.unnamed_idx_0 = ((&nmpcSine_M)->Timing.TaskCounters.TID[1] == 0);
    if (nmpcSine_B.unnamed_idx_0) {
      // MATLAB Function: '<S3>/HeadingToPsi'
      if (rtIsNaN(nmpcSine_B.In1_c.data - 26.566667)) {
        nmpcSine_B.psi = (rtNaN);
      } else if (rtIsInf(nmpcSine_B.In1_c.data - 26.566667)) {
        nmpcSine_B.psi = (rtNaN);
      } else if (nmpcSine_B.In1_c.data - 26.566667 == 0.0) {
        nmpcSine_B.psi = 0.0;
      } else {
        nmpcSine_B.psi = fmod(nmpcSine_B.In1_c.data - 26.566667, 360.0);
        if (nmpcSine_B.psi == 0.0) {
          nmpcSine_B.psi = 0.0;
        } else if (nmpcSine_B.psi < 0.0) {
          nmpcSine_B.psi += 360.0;
        }
      }

      if ((nmpcSine_B.psi == 0.0) && (nmpcSine_B.In1_c.data - 26.566667 > 0.0))
      {
        nmpcSine_B.psi = 360.0;
      }

      nmpcSine_B.psi *= 0.017453292519943295;

      // End of MATLAB Function: '<S3>/HeadingToPsi'
    }

    // MATLAB Function: '<S3>/LatLonToXY'
    nmpcSine_B.s1 = nmpcSine_B.gps0[0];
    nmpcSine_sind(&nmpcSine_B.s1);
    nmpcSine_B.c1 = nmpcSine_B.gps0[0];
    nmpcSine_cosd(&nmpcSine_B.c1);
    nmpcSine_B.s2 = nmpcSine_B.In1_j.latitude;
    nmpcSine_sind(&nmpcSine_B.s2);
    nmpcSine_B.c2 = nmpcSine_B.In1_j.latitude;
    nmpcSine_cosd(&nmpcSine_B.c2);
    nmpcSine_B.b_g = nmpcSine_B.gps0[1];
    nmpcSine_cosd(&nmpcSine_B.b_g);
    nmpcSine_B.c_c = nmpcSine_B.In1_j.longitude;
    nmpcSine_cosd(&nmpcSine_B.c_c);
    nmpcSine_B.epsilon = nmpcSine_B.gps0[1];
    nmpcSine_sind(&nmpcSine_B.epsilon);
    nmpcSine_B.e = nmpcSine_B.In1_j.longitude;
    nmpcSine_sind(&nmpcSine_B.e);
    nmpcSine_B.w1 = 1.0 / sqrt(1.0 - nmpcSine_B.s1 * nmpcSine_B.s1 *
      0.0066943799901413165);
    nmpcSine_B.w2 = 1.0 / sqrt(1.0 - nmpcSine_B.s2 * nmpcSine_B.s2 *
      0.0066943799901413165);
    tmp = _mm_set_pd(nmpcSine_B.c2 * nmpcSine_B.e, nmpcSine_B.c2 *
                     nmpcSine_B.c_c);
    tmp_0 = _mm_set_pd(nmpcSine_B.c1 * nmpcSine_B.epsilon, nmpcSine_B.c1 *
                       nmpcSine_B.b_g);
    _mm_storeu_pd(&nmpcSine_B.dv4[0], _mm_add_pd(_mm_mul_pd(_mm_sub_pd
      (_mm_mul_pd(tmp, _mm_set1_pd(nmpcSine_B.w2)), _mm_mul_pd(tmp_0,
      _mm_set1_pd(nmpcSine_B.w1))), _mm_set1_pd(6.378137E+6)), _mm_sub_pd
      (_mm_mul_pd(_mm_set1_pd(nmpcSine_B.In1_j.altitude), tmp), _mm_mul_pd
       (_mm_set1_pd(nmpcSine_B.gps0[2]), tmp_0))));
    nmpcSine_B.c1 = nmpcSine_B.gps0[1];
    nmpcSine_cosd(&nmpcSine_B.c1);
    nmpcSine_B.c2 = nmpcSine_B.gps0[1];
    nmpcSine_sind(&nmpcSine_B.c2);
    nmpcSine_B.b_g = nmpcSine_B.gps0[0];
    nmpcSine_sind(&nmpcSine_B.b_g);
    nmpcSine_B.c_c = nmpcSine_B.gps0[0];
    nmpcSine_cosd(&nmpcSine_B.c_c);

    // ZeroOrderHold: '<S3>/Zero-Order Hold'
    if ((&nmpcSine_M)->Timing.TaskCounters.TID[2] == 0) {
      // Outputs for Enabled SubSystem: '<S10>/Correct1' incorporates:
      //   EnablePort: '<S18>/Enable'

      // MATLABSystem: '<S18>/MATLAB System' incorporates:
      //   Constant: '<S10>/BlockOrdering'

      nmpcSine_B.MATLABSystem_o3_i = true;

      // MATLABSystem: '<S18>/MATLAB System' incorporates:
      //   DataStoreRead: '<S18>/Data Store ReadP'
      //   DataStoreRead: '<S18>/Data Store ReadX'

      //  GPS Measurment Function
      //  Position is measured
      for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
        for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
          nmpcSine_B.imvec[nmpcSine_B.i] = nmpcSine_DW.x[nmpcSine_B.i];
        }

        // Start for MATLABSystem: '<S18>/MATLAB System' incorporates:
        //   DataStoreRead: '<S18>/Data Store ReadX'

        nmpcSine_B.epsilon = 1.4901161193847656E-8 * fabs
          (nmpcSine_DW.x[nmpcSine_B.b_j]);
        if ((nmpcSine_B.epsilon <= 1.4901161193847656E-8) || rtIsNaN
            (nmpcSine_B.epsilon)) {
          nmpcSine_B.epsilon = 1.4901161193847656E-8;
        }

        nmpcSine_B.imvec[nmpcSine_B.b_j] = nmpcSine_DW.x[nmpcSine_B.b_j] +
          nmpcSine_B.epsilon;

        //  GPS Measurment Function
        //  Position is measured
        tmp = _mm_div_pd(_mm_sub_pd(_mm_loadu_pd(&nmpcSine_B.imvec[0]),
          _mm_loadu_pd(&nmpcSine_DW.x[0])), _mm_set1_pd(nmpcSine_B.epsilon));
        _mm_storeu_pd(&nmpcSine_B.b_dHdx[nmpcSine_B.b_j << 1], tmp);
      }

      //  GPS Measurment Function
      //  Position is measured
      for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 2; nmpcSine_B.b_j++) {
        nmpcSine_B.coffset = nmpcSine_B.b_j * 6 - 1;
        for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
          nmpcSine_B.aoffset = nmpcSine_B.i * 6 - 1;
          nmpcSine_B.epsilon = 0.0;
          for (nmpcSine_B.b_k = 0; nmpcSine_B.b_k < 6; nmpcSine_B.b_k++) {
            nmpcSine_B.epsilon += nmpcSine_DW.P[(nmpcSine_B.aoffset +
              nmpcSine_B.b_k) + 1] * nmpcSine_B.b_dHdx[(nmpcSine_B.b_k << 1) +
              nmpcSine_B.b_j];
          }

          nmpcSine_B.K_c[(nmpcSine_B.coffset + nmpcSine_B.i) + 1] =
            nmpcSine_B.epsilon;
        }
      }

      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        // Start for MATLABSystem: '<S18>/MATLAB System'
        nmpcSine_B.K[nmpcSine_B.i] = nmpcSine_B.K_c[nmpcSine_B.i];
        nmpcSine_B.K[nmpcSine_B.i + 8] = nmpcSine_B.K_c[nmpcSine_B.i + 6];
      }

      // Start for MATLABSystem: '<S18>/MATLAB System' incorporates:
      //   Constant: '<S10>/R1'

      nmpcSine_B.K[6] = 0.007;
      nmpcSine_B.K[7] = 0.0;
      nmpcSine_B.K[14] = 0.0;
      nmpcSine_B.K[15] = 0.007;

      // MATLABSystem: '<S18>/MATLAB System'
      nmpcSine_qr(nmpcSine_B.K, nmpcSine_B.a__1_m, nmpcSine_B.R);
      nmpcSine_B.Sy[0] = nmpcSine_B.R[0];
      nmpcSine_B.Sy[1] = nmpcSine_B.R[2];
      nmpcSine_B.Sy[2] = nmpcSine_B.R[1];
      nmpcSine_B.Sy[3] = nmpcSine_B.R[3];
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset
             ++) {
          // Start for MATLABSystem: '<S18>/MATLAB System' incorporates:
          //   DataStoreRead: '<S18>/Data Store ReadP'

          nmpcSine_B.epsilon = 0.0;
          for (nmpcSine_B.aoffset = 0; nmpcSine_B.aoffset < 6;
               nmpcSine_B.aoffset++) {
            nmpcSine_B.epsilon += nmpcSine_DW.P[6 * nmpcSine_B.aoffset +
              nmpcSine_B.i] * nmpcSine_DW.P[6 * nmpcSine_B.aoffset +
              nmpcSine_B.coffset];
          }

          nmpcSine_B.A[nmpcSine_B.i + 6 * nmpcSine_B.coffset] =
            nmpcSine_B.epsilon;
        }
      }

      for (nmpcSine_B.i = 0; nmpcSine_B.i < 2; nmpcSine_B.i++) {
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset
             ++) {
          // Start for MATLABSystem: '<S18>/MATLAB System'
          nmpcSine_B.epsilon = 0.0;
          for (nmpcSine_B.aoffset = 0; nmpcSine_B.aoffset < 6;
               nmpcSine_B.aoffset++) {
            nmpcSine_B.epsilon += nmpcSine_B.A[6 * nmpcSine_B.aoffset +
              nmpcSine_B.coffset] * nmpcSine_B.b_dHdx[(nmpcSine_B.aoffset << 1)
              + nmpcSine_B.i];
          }

          nmpcSine_B.K_c[nmpcSine_B.i + (nmpcSine_B.coffset << 1)] =
            nmpcSine_B.epsilon;
        }
      }

      // MATLABSystem: '<S18>/MATLAB System'
      for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
        // Start for MATLABSystem: '<S18>/MATLAB System'
        nmpcSine_B.i = nmpcSine_B.b_j << 1;
        nmpcSine_B.C[nmpcSine_B.i] = nmpcSine_B.K_c[nmpcSine_B.i];
        nmpcSine_B.C[nmpcSine_B.i + 1] = nmpcSine_B.K_c[nmpcSine_B.i + 1];
      }

      // Start for MATLABSystem: '<S18>/MATLAB System'
      nmpcSine_trisolve(nmpcSine_B.Sy, nmpcSine_B.C);

      // MATLABSystem: '<S18>/MATLAB System'
      for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
        // Start for MATLABSystem: '<S18>/MATLAB System'
        nmpcSine_B.i = nmpcSine_B.b_j << 1;
        nmpcSine_B.C_m[nmpcSine_B.i] = nmpcSine_B.C[nmpcSine_B.i];
        nmpcSine_B.C_m[nmpcSine_B.i + 1] = nmpcSine_B.C[nmpcSine_B.i + 1];
      }

      // Start for MATLABSystem: '<S18>/MATLAB System'
      nmpcSine_B.Sy[0] = nmpcSine_B.R[0];
      nmpcSine_B.Sy[1] = nmpcSine_B.R[1];
      nmpcSine_B.Sy[2] = nmpcSine_B.R[2];
      nmpcSine_B.Sy[3] = nmpcSine_B.R[3];
      nmpcSine_trisolve_m(nmpcSine_B.Sy, nmpcSine_B.C_m);

      // MATLABSystem: '<S18>/MATLAB System'
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        nmpcSine_B.b_j = nmpcSine_B.i << 1;
        nmpcSine_B.K_c[nmpcSine_B.i] = nmpcSine_B.C_m[nmpcSine_B.b_j];
        nmpcSine_B.K_c[nmpcSine_B.i + 6] = nmpcSine_B.C_m[nmpcSine_B.b_j + 1];
      }

      for (nmpcSine_B.i = 0; nmpcSine_B.i <= 10; nmpcSine_B.i += 2) {
        // Start for MATLABSystem: '<S18>/MATLAB System'
        tmp = _mm_loadu_pd(&nmpcSine_B.K_c[nmpcSine_B.i]);
        _mm_storeu_pd(&nmpcSine_B.C[nmpcSine_B.i], _mm_mul_pd(tmp, _mm_set1_pd
          (-1.0)));
      }

      // MATLABSystem: '<S18>/MATLAB System' incorporates:
      //   DataStoreRead: '<S18>/Data Store ReadP'

      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        // Start for MATLABSystem: '<S18>/MATLAB System'
        nmpcSine_B.b_j = nmpcSine_B.i << 1;
        nmpcSine_B.epsilon = nmpcSine_B.b_dHdx[nmpcSine_B.b_j + 1];
        nmpcSine_B.e = nmpcSine_B.b_dHdx[nmpcSine_B.b_j];
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset <= 4; nmpcSine_B.coffset
             += 2) {
          tmp = _mm_loadu_pd(&nmpcSine_B.C[nmpcSine_B.coffset + 6]);
          tmp_0 = _mm_loadu_pd(&nmpcSine_B.C[nmpcSine_B.coffset]);
          _mm_storeu_pd(&nmpcSine_B.A[nmpcSine_B.coffset + 6 * nmpcSine_B.i],
                        _mm_add_pd(_mm_mul_pd(_mm_set1_pd(nmpcSine_B.epsilon),
            tmp), _mm_mul_pd(_mm_set1_pd(nmpcSine_B.e), tmp_0)));
        }
      }

      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        nmpcSine_B.b_j = 6 * nmpcSine_B.i + nmpcSine_B.i;
        nmpcSine_B.A[nmpcSine_B.b_j]++;
      }

      for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
        nmpcSine_B.coffset = nmpcSine_B.b_j * 6 - 1;
        for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
          nmpcSine_B.aoffset = nmpcSine_B.i * 6 - 1;
          nmpcSine_B.epsilon = 0.0;
          for (nmpcSine_B.b_k = 0; nmpcSine_B.b_k < 6; nmpcSine_B.b_k++) {
            nmpcSine_B.epsilon += nmpcSine_DW.P[(nmpcSine_B.aoffset +
              nmpcSine_B.b_k) + 1] * nmpcSine_B.A[nmpcSine_B.b_k * 6 +
              nmpcSine_B.b_j];
          }

          nmpcSine_B.y_ju[(nmpcSine_B.coffset + nmpcSine_B.i) + 1] =
            nmpcSine_B.epsilon;
        }

        // Start for MATLABSystem: '<S18>/MATLAB System' incorporates:
        //   Constant: '<S10>/R1'
        //   DataStoreRead: '<S18>/Data Store ReadP'

        nmpcSine_B.epsilon = nmpcSine_B.K_c[nmpcSine_B.b_j + 6];
        nmpcSine_B.e = nmpcSine_B.K_c[nmpcSine_B.b_j];
        for (nmpcSine_B.i = 0; nmpcSine_B.i <= 0; nmpcSine_B.i += 2) {
          nmpcSine_B.coffset = (nmpcSine_B.i + 1) << 1;
          nmpcSine_B.aoffset = nmpcSine_B.i << 1;
          _mm_storeu_pd(&nmpcSine_B.C[nmpcSine_B.i + (nmpcSine_B.b_j << 1)],
                        _mm_add_pd(_mm_mul_pd(_mm_set_pd
            (nmpcSine_ConstP.R1_Value[nmpcSine_B.coffset + 1],
             nmpcSine_ConstP.R1_Value[nmpcSine_B.aoffset + 1]), _mm_set1_pd
            (nmpcSine_B.epsilon)), _mm_mul_pd(_mm_set_pd
            (nmpcSine_ConstP.R1_Value[nmpcSine_B.coffset],
             nmpcSine_ConstP.R1_Value[nmpcSine_B.aoffset]), _mm_set1_pd
            (nmpcSine_B.e))));
        }
      }

      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset
             ++) {
          // Start for MATLABSystem: '<S18>/MATLAB System'
          nmpcSine_B.y_d[nmpcSine_B.coffset + (nmpcSine_B.i << 3)] =
            nmpcSine_B.y_ju[6 * nmpcSine_B.i + nmpcSine_B.coffset];
        }

        // Start for MATLABSystem: '<S18>/MATLAB System'
        nmpcSine_B.coffset = nmpcSine_B.i << 1;
        nmpcSine_B.b_j = nmpcSine_B.i << 3;
        nmpcSine_B.y_d[nmpcSine_B.b_j + 6] = nmpcSine_B.C[nmpcSine_B.coffset];
        nmpcSine_B.y_d[nmpcSine_B.b_j + 7] = nmpcSine_B.C[nmpcSine_B.coffset + 1];
      }

      // MATLABSystem: '<S18>/MATLAB System'
      nmpcSine_qr_m(nmpcSine_B.y_d, nmpcSine_B.a__1_b, nmpcSine_B.A);

      // MATLAB Function: '<S3>/LatLonToXY' incorporates:
      //   DataStoreRead: '<S18>/Data Store ReadX'
      //   MATLABSystem: '<S18>/MATLAB System'
      //
      nmpcSine_B.s1 = (((nmpcSine_B.s2 * nmpcSine_B.w2 - nmpcSine_B.s1 *
                         nmpcSine_B.w1) * 6.3354393272928195E+6 +
                        (nmpcSine_B.In1_j.altitude * nmpcSine_B.s2 -
                         nmpcSine_B.gps0[2] * nmpcSine_B.s1)) * nmpcSine_B.c_c +
                       (nmpcSine_B.c1 * nmpcSine_B.dv4[0] + nmpcSine_B.c2 *
                        nmpcSine_B.dv4[1]) * -nmpcSine_B.b_g) - nmpcSine_DW.x[0];
      nmpcSine_B.s2 = (-nmpcSine_B.c2 * nmpcSine_B.dv4[0] + nmpcSine_B.c1 *
                       nmpcSine_B.dv4[1]) - nmpcSine_DW.x[1];
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        // DataStoreWrite: '<S18>/Data Store WriteP' incorporates:
        //   MATLABSystem: '<S18>/MATLAB System'
        //
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset
             ++) {
          nmpcSine_DW.P[nmpcSine_B.coffset + 6 * nmpcSine_B.i] = nmpcSine_B.A[6 *
            nmpcSine_B.coffset + nmpcSine_B.i];
        }

        // End of DataStoreWrite: '<S18>/Data Store WriteP'

        // DataStoreWrite: '<S18>/Data Store WriteX' incorporates:
        //   DataStoreRead: '<S18>/Data Store ReadX'
        //   MATLABSystem: '<S18>/MATLAB System'
        //
        nmpcSine_DW.x[nmpcSine_B.i] += nmpcSine_B.K_c[nmpcSine_B.i + 6] *
          nmpcSine_B.s2 + nmpcSine_B.K_c[nmpcSine_B.i] * nmpcSine_B.s1;
      }

      // End of Outputs for SubSystem: '<S10>/Correct1'
    }

    // End of ZeroOrderHold: '<S3>/Zero-Order Hold'

    // BusAssignment: '<S16>/Bus Assignment'
    memset(&nmpcSine_B.BusAssignment_h, 0, sizeof(SL_Bus_geometry_msgs_Point));
    if (nmpcSine_B.unnamed_idx_0) {
      // Outputs for Enabled SubSystem: '<S10>/Correct2' incorporates:
      //   EnablePort: '<S19>/Enable'

      // MATLABSystem: '<S19>/MATLAB System' incorporates:
      //   Constant: '<S10>/R2'
      //   DataStoreRead: '<S19>/Data Store ReadP'
      //   DataStoreRead: '<S19>/Data Store ReadX'

      //  Magnetometer Measurement Function
      //  Heading is measured
      //  Heading is bounded between 0 and 360 degrees
      for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
        for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
          nmpcSine_B.imvec[nmpcSine_B.i] = nmpcSine_DW.x[nmpcSine_B.i];
        }

        nmpcSine_B.epsilon = 1.4901161193847656E-8 * fabs
          (nmpcSine_DW.x[nmpcSine_B.b_j]);
        if ((nmpcSine_B.epsilon <= 1.4901161193847656E-8) || rtIsNaN
            (nmpcSine_B.epsilon)) {
          nmpcSine_B.epsilon = 1.4901161193847656E-8;
        }

        nmpcSine_B.imvec[nmpcSine_B.b_j] = nmpcSine_DW.x[nmpcSine_B.b_j] +
          nmpcSine_B.epsilon;

        //  Magnetometer Measurement Function
        //  Heading is measured
        //  Heading is bounded between 0 and 360 degrees
        nmpcSine_B.z_c[nmpcSine_B.b_j] = (nmpcSine_B.imvec[2] - nmpcSine_DW.x[2])
          / nmpcSine_B.epsilon;
      }

      //  Magnetometer Measurement Function
      //  Heading is measured
      //  Heading is bounded between 0 and 360 degrees
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        nmpcSine_B.aoffset = nmpcSine_B.i * 6 - 1;
        nmpcSine_B.epsilon = 0.0;
        for (nmpcSine_B.b_k = 0; nmpcSine_B.b_k < 6; nmpcSine_B.b_k++) {
          nmpcSine_B.epsilon += nmpcSine_DW.P[(nmpcSine_B.aoffset +
            nmpcSine_B.b_k) + 1] * nmpcSine_B.z_c[nmpcSine_B.b_k];
        }

        nmpcSine_B.c_A_c[nmpcSine_B.i] = nmpcSine_B.epsilon;
      }

      nmpcSine_B.c_A_c[6] = 0.017453292519943295;
      nmpcSine_B.psi -= nmpcSine_DW.x[2];
      nmpcSine_B.ref[0] = 0.0;
      nmpcSine_B.ref[1] = 6.2831853071795862;

      // End of Outputs for SubSystem: '<S10>/Correct2'
      nmpcSine_B.b_j = 0;

      // Outputs for Enabled SubSystem: '<S10>/Correct2' incorporates:
      //   EnablePort: '<S19>/Enable'

      // MATLABSystem: '<S19>/MATLAB System'
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 1; nmpcSine_B.i++) {
        nmpcSine_B.s1 = nmpcSine_B.c_A_c[0];
        nmpcSine_B.s2 = 0.0;
        nmpcSine_B.w1 = 3.3121686421112381E-170;
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset
             ++) {
          nmpcSine_B.w2 = fabs(nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1]);
          if (nmpcSine_B.w2 > nmpcSine_B.w1) {
            nmpcSine_B.c1 = nmpcSine_B.w1 / nmpcSine_B.w2;
            nmpcSine_B.s2 = nmpcSine_B.s2 * nmpcSine_B.c1 * nmpcSine_B.c1 + 1.0;
            nmpcSine_B.w1 = nmpcSine_B.w2;
          } else {
            nmpcSine_B.c1 = nmpcSine_B.w2 / nmpcSine_B.w1;
            nmpcSine_B.s2 += nmpcSine_B.c1 * nmpcSine_B.c1;
          }
        }

        nmpcSine_B.s2 = nmpcSine_B.w1 * sqrt(nmpcSine_B.s2);
        if (nmpcSine_B.s2 != 0.0) {
          nmpcSine_B.s2 = nmpcSine_rt_hypotd_snf(nmpcSine_B.c_A_c[0],
            nmpcSine_B.s2);
          if (nmpcSine_B.c_A_c[0] >= 0.0) {
            nmpcSine_B.s2 = -nmpcSine_B.s2;
          }

          if (fabs(nmpcSine_B.s2) < 1.0020841800044864E-292) {
            nmpcSine_B.aoffset = -1;
            do {
              nmpcSine_B.aoffset++;
              for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset <= 4;
                   nmpcSine_B.coffset += 2) {
                tmp = _mm_loadu_pd(&nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1]);
                _mm_storeu_pd(&nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1],
                              _mm_mul_pd(tmp, _mm_set1_pd(9.9792015476736E+291)));
              }

              nmpcSine_B.s2 *= 9.9792015476736E+291;
              nmpcSine_B.s1 *= 9.9792015476736E+291;
            } while ((fabs(nmpcSine_B.s2) < 1.0020841800044864E-292) &&
                     (nmpcSine_B.aoffset + 1 < 20));

            nmpcSine_B.s2 = 0.0;
            nmpcSine_B.w1 = 3.3121686421112381E-170;
            for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6;
                 nmpcSine_B.coffset++) {
              nmpcSine_B.w2 = fabs(nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1]);
              if (nmpcSine_B.w2 > nmpcSine_B.w1) {
                nmpcSine_B.c1 = nmpcSine_B.w1 / nmpcSine_B.w2;
                nmpcSine_B.s2 = nmpcSine_B.s2 * nmpcSine_B.c1 * nmpcSine_B.c1 +
                  1.0;
                nmpcSine_B.w1 = nmpcSine_B.w2;
              } else {
                nmpcSine_B.c1 = nmpcSine_B.w2 / nmpcSine_B.w1;
                nmpcSine_B.s2 += nmpcSine_B.c1 * nmpcSine_B.c1;
              }
            }

            nmpcSine_B.s2 = nmpcSine_rt_hypotd_snf(nmpcSine_B.s1, nmpcSine_B.w1 *
              sqrt(nmpcSine_B.s2));
            if (nmpcSine_B.s1 >= 0.0) {
              nmpcSine_B.s2 = -nmpcSine_B.s2;
            }

            nmpcSine_B.s1 = 1.0 / (nmpcSine_B.s1 - nmpcSine_B.s2);
            for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset <= 4;
                 nmpcSine_B.coffset += 2) {
              tmp = _mm_loadu_pd(&nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1]);
              _mm_storeu_pd(&nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1],
                            _mm_mul_pd(tmp, _mm_set1_pd(nmpcSine_B.s1)));
            }

            for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset <=
                 nmpcSine_B.aoffset; nmpcSine_B.coffset++) {
              nmpcSine_B.s2 *= 1.0020841800044864E-292;
            }

            nmpcSine_B.s1 = nmpcSine_B.s2;
          } else {
            nmpcSine_B.s1 = 1.0 / (nmpcSine_B.c_A_c[0] - nmpcSine_B.s2);
            for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset <= 4;
                 nmpcSine_B.coffset += 2) {
              tmp = _mm_loadu_pd(&nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1]);
              _mm_storeu_pd(&nmpcSine_B.c_A_c[nmpcSine_B.coffset + 1],
                            _mm_mul_pd(tmp, _mm_set1_pd(nmpcSine_B.s1)));
            }

            nmpcSine_B.s1 = nmpcSine_B.s2;
          }
        }

        nmpcSine_B.c_A_c[0] = nmpcSine_B.s1;
        nmpcSine_B.b_j++;
      }

      if (nmpcSine_B.b_j - 1 >= 0) {
        for (nmpcSine_B.i = 0; nmpcSine_B.i < 2; nmpcSine_B.i++) {
          nmpcSine_B.ref[nmpcSine_B.i] = 6.2831853071795862 * static_cast<real_T>
            (nmpcSine_B.i) - 3.1415926535897931;
        }
      }

      // End of Outputs for SubSystem: '<S10>/Correct2'
      nmpcSine_B.b_j = 0;

      // Outputs for Enabled SubSystem: '<S10>/Correct2' incorporates:
      //   EnablePort: '<S19>/Enable'

      // Start for MATLABSystem: '<S19>/MATLAB System'
      if (fabs(fabs(nmpcSine_B.psi) / nmpcSine_B.ref[1]) > 0.001) {
        for (nmpcSine_B.i = 0; nmpcSine_B.i < 1; nmpcSine_B.i++) {
          nmpcSine_B.b_j++;
        }
      }

      // MATLABSystem: '<S19>/MATLAB System'
      if (nmpcSine_B.b_j - 1 >= 0) {
        _mm_storeu_pd(&nmpcSine_B.dv4[0], _mm_sub_pd(_mm_set_pd(nmpcSine_B.ref[1],
          nmpcSine_B.psi), _mm_set1_pd(nmpcSine_B.ref[0])));
        nmpcSine_B.resToWrap_data_idx_0 = nmpcSine_B.dv4[0];
        nmpcSine_B.y_data_idx_0 = nmpcSine_B.dv4[1];
      }

      if (nmpcSine_B.b_j - 1 >= 0) {
        nmpcSine_B.l_data = nmpcSine_mod(nmpcSine_B.resToWrap_data_idx_0,
          nmpcSine_B.y_data_idx_0);
      }

      nmpcSine_B.o = nmpcSine_B.psi;
      if (nmpcSine_B.b_j - 1 >= 0) {
        nmpcSine_B.o = nmpcSine_B.l_data + nmpcSine_B.ref[0];
      }

      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        // Start for MATLABSystem: '<S19>/MATLAB System' incorporates:
        //   DataStoreRead: '<S19>/Data Store ReadP'

        nmpcSine_B.Pxy[nmpcSine_B.i] = 0.0;
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset
             ++) {
          nmpcSine_B.epsilon = 0.0;
          for (nmpcSine_B.aoffset = 0; nmpcSine_B.aoffset < 6;
               nmpcSine_B.aoffset++) {
            nmpcSine_B.epsilon += nmpcSine_DW.P[6 * nmpcSine_B.aoffset +
              nmpcSine_B.i] * nmpcSine_DW.P[6 * nmpcSine_B.aoffset +
              nmpcSine_B.coffset];
          }

          nmpcSine_B.Pxy[nmpcSine_B.i] += nmpcSine_B.epsilon *
            nmpcSine_B.z_c[nmpcSine_B.coffset];
        }
      }

      for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
        // Start for MATLABSystem: '<S19>/MATLAB System' incorporates:
        //   DataStoreRead: '<S19>/Data Store ReadX'

        nmpcSine_B.imvec[nmpcSine_B.b_j] = nmpcSine_DW.x[nmpcSine_B.b_j];
      }

      // DataStoreRead: '<S19>/Data Store ReadP'
      memcpy(&nmpcSine_B.A[0], &nmpcSine_DW.P[0], 36U * sizeof(real_T));

      // DataStoreWrite: '<S19>/Data Store WriteX' incorporates:
      //   Constant: '<S10>/R2'
      //   MATLABSystem: '<S19>/MATLAB System'
      //
      EKFCorrector_correctStateAndSqr(nmpcSine_B.imvec, nmpcSine_B.A,
        nmpcSine_B.o, nmpcSine_B.Pxy, nmpcSine_B.c_A_c[0], nmpcSine_B.z_c,
        0.017453292519943295, nmpcSine_DW.x, nmpcSine_DW.P);

      // End of Outputs for SubSystem: '<S10>/Correct2'
    }

    // Outputs for Enabled SubSystem: '<S10>/Correct3' incorporates:
    //   EnablePort: '<S20>/Enable'

    // MATLABSystem: '<S20>/MATLAB System' incorporates:
    //   Constant: '<S10>/R3'
    //   DataStoreRead: '<S20>/Data Store ReadP'
    //   DataStoreRead: '<S20>/Data Store ReadX'
    //   Gain: '<S3>/Gain'

    EKFCorrectorAdditive_getMeasure(0.024494897427831779, nmpcSine_DW.x,
      nmpcSine_DW.P, &nmpcSine_B.s2, nmpcSine_B.Pxy, &nmpcSine_B.s1,
      nmpcSine_B.imvec, &nmpcSine_B.psi);
    nmpcSine_B.s2 = -nmpcSine_B.In1.angular_velocity.z - nmpcSine_B.s2;
    for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
      nmpcSine_B.z_c[nmpcSine_B.b_j] = nmpcSine_B.Pxy[nmpcSine_B.b_j];
    }

    // Start for MATLABSystem: '<S20>/MATLAB System'
    nmpcSine_trisolve_m3(nmpcSine_B.s1, nmpcSine_B.z_c);

    // MATLABSystem: '<S20>/MATLAB System'
    for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
      nmpcSine_B.Pxy[nmpcSine_B.b_j] = nmpcSine_B.z_c[nmpcSine_B.b_j];
    }

    // Start for MATLABSystem: '<S20>/MATLAB System'
    nmpcSine_trisolve_m3(nmpcSine_B.s1, nmpcSine_B.Pxy);
    for (nmpcSine_B.i = 0; nmpcSine_B.i <= 4; nmpcSine_B.i += 2) {
      tmp = _mm_loadu_pd(&nmpcSine_B.Pxy[nmpcSine_B.i]);
      _mm_storeu_pd(&nmpcSine_B.z_c[nmpcSine_B.i], _mm_mul_pd(tmp, _mm_set1_pd
        (-1.0)));
    }

    // MATLABSystem: '<S20>/MATLAB System'
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset <= 4; nmpcSine_B.coffset +=
           2) {
        tmp = _mm_loadu_pd(&nmpcSine_B.z_c[nmpcSine_B.coffset]);
        _mm_storeu_pd(&nmpcSine_B.A[nmpcSine_B.coffset + 6 * nmpcSine_B.i],
                      _mm_mul_pd(tmp, _mm_set1_pd(nmpcSine_B.imvec[nmpcSine_B.i])));
      }
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      nmpcSine_B.b_j = 6 * nmpcSine_B.i + nmpcSine_B.i;
      nmpcSine_B.A[nmpcSine_B.b_j]++;
    }

    for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
      // Start for MATLABSystem: '<S20>/MATLAB System' incorporates:
      //   DataStoreRead: '<S20>/Data Store ReadP'

      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        nmpcSine_B.aoffset = nmpcSine_B.i * 6 - 1;
        nmpcSine_B.epsilon = 0.0;
        for (nmpcSine_B.b_k = 0; nmpcSine_B.b_k < 6; nmpcSine_B.b_k++) {
          nmpcSine_B.epsilon += nmpcSine_DW.P[(nmpcSine_B.aoffset +
            nmpcSine_B.b_k) + 1] * nmpcSine_B.A[nmpcSine_B.b_k * 6 +
            nmpcSine_B.b_j];
        }

        nmpcSine_B.y_b[nmpcSine_B.i + 7 * nmpcSine_B.b_j] = nmpcSine_B.epsilon;
      }

      nmpcSine_B.y_b[7 * nmpcSine_B.b_j + 6] = nmpcSine_B.Pxy[nmpcSine_B.b_j] *
        nmpcSine_B.psi;
    }

    nmpcSine_qr_m3(nmpcSine_B.y_b, nmpcSine_B.a__1_e, nmpcSine_B.A);
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      // DataStoreWrite: '<S20>/Data Store WriteP' incorporates:
      //   MATLABSystem: '<S20>/MATLAB System'
      //
      for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset++)
      {
        nmpcSine_DW.P[nmpcSine_B.coffset + 6 * nmpcSine_B.i] = nmpcSine_B.A[6 *
          nmpcSine_B.coffset + nmpcSine_B.i];
      }

      // End of DataStoreWrite: '<S20>/Data Store WriteP'

      // DataStoreWrite: '<S20>/Data Store WriteX' incorporates:
      //   DataStoreRead: '<S20>/Data Store ReadX'
      //   MATLABSystem: '<S20>/MATLAB System'
      //
      nmpcSine_DW.x[nmpcSine_B.i] += nmpcSine_B.Pxy[nmpcSine_B.i] *
        nmpcSine_B.s2;
    }

    // End of Outputs for SubSystem: '<S10>/Correct3'

    // MATLAB Function: '<S3>/MATLAB Function1' incorporates:
    //   Concatenate: '<S3>/Vector Concatenate1'
    //   Constant: '<S3>/Constant5'

    nmpcSine_B.ref[0] = 0.0;
    nmpcSine_B.ref[1] = 0.0;
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 40; nmpcSine_B.i++) {
      nmpcSine_B.VectorConcatenate1[nmpcSine_B.i] = 0.0;
      nmpcSine_B.VectorConcatenate1[nmpcSine_B.i + 40] = 0.0;
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 20; nmpcSine_B.i++) {
      nmpcSine_B.psi = (nmpcSine_DW.k + static_cast<real_T>(nmpcSine_B.i)) - 1.0;
      if (rtIsNaN(nmpcSine_B.psi)) {
        nmpcSine_B.psi = (rtNaN);
      } else if (rtIsInf(nmpcSine_B.psi)) {
        nmpcSine_B.psi = (rtNaN);
      } else if (nmpcSine_B.psi == 0.0) {
        nmpcSine_B.psi = 0.0;
      } else {
        nmpcSine_B.psi = fmod(nmpcSine_B.psi, 1501.0);
        if (nmpcSine_B.psi == 0.0) {
          nmpcSine_B.psi = 0.0;
        } else if (nmpcSine_B.psi < 0.0) {
          nmpcSine_B.psi += 1501.0;
        }
      }

      if (nmpcSine_B.i == 0) {
        nmpcSine_B.ref[0] = nmpcSine_ConstP.Constant5_Value[static_cast<int32_T>
          (nmpcSine_B.psi + 1.0) - 1];
        nmpcSine_B.ref[1] = nmpcSine_ConstP.Constant5_Value[static_cast<int32_T>
          (nmpcSine_B.psi + 1.0) + 1500];
      }

      nmpcSine_B.VectorConcatenate1[nmpcSine_B.i] =
        nmpcSine_ConstP.Constant5_Value[static_cast<int32_T>(nmpcSine_B.psi +
        1.0) - 1];
      nmpcSine_B.VectorConcatenate1[nmpcSine_B.i + 20] =
        nmpcSine_ConstP.Constant5_Value[static_cast<int32_T>(nmpcSine_B.psi +
        1.0) + 1500];
    }

    nmpcSine_DW.k++;

    // End of MATLAB Function: '<S3>/MATLAB Function1'

    // Delay: '<S25>/mv_Delay' incorporates:
    //   Constant: '<S25>/ones'
    //   Product: '<S25>/Product'
    //   UnitDelay: '<S3>/Unit Delay'

    if (nmpcSine_DW.icLoad) {
      // Product: '<S25>/Product' incorporates:
      //   UnitDelay: '<S3>/Unit Delay'

      nmpcSine_B.psi = nmpcSine_DW.UnitDelay_DSTATE[0];
      nmpcSine_B.s1 = nmpcSine_DW.UnitDelay_DSTATE[1];
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 21; nmpcSine_B.i++) {
        nmpcSine_DW.mv_Delay_DSTATE[nmpcSine_B.i] = nmpcSine_B.psi;
        nmpcSine_DW.mv_Delay_DSTATE[nmpcSine_B.i + 21] = nmpcSine_B.s1;
      }
    }

    memcpy(&nmpcSine_B.x_Delay[0], &nmpcSine_DW.mv_Delay_DSTATE[0], 42U * sizeof
           (real_T));

    // End of Delay: '<S25>/mv_Delay'
    for (nmpcSine_B.i = 0; nmpcSine_B.i <= 16; nmpcSine_B.i += 2) {
      // Selector: '<S25>/Selector1' incorporates:
      //   Constant: '<S25>/Constant1'

      nmpcSine_B.b_j = static_cast<int32_T>
        (nmpcSine_ConstP.Constant1_Value[nmpcSine_B.i]);

      // Selector: '<S25>/Selector1' incorporates:
      //   Constant: '<S25>/Constant1'

      nmpcSine_B.Selector1_a[nmpcSine_B.i] = nmpcSine_B.x_Delay[nmpcSine_B.b_j -
        1];

      // Selector: '<S25>/Selector1' incorporates:
      //   Constant: '<S25>/Constant1'

      nmpcSine_B.coffset = static_cast<int32_T>
        (nmpcSine_ConstP.Constant1_Value[nmpcSine_B.i + 1]);

      // Selector: '<S25>/Selector1' incorporates:
      //   Constant: '<S25>/Constant1'

      nmpcSine_B.Selector1_a[nmpcSine_B.i + 1] =
        nmpcSine_B.x_Delay[nmpcSine_B.coffset - 1];
      nmpcSine_B.Selector1_a[nmpcSine_B.i + 19] =
        nmpcSine_B.x_Delay[nmpcSine_B.b_j + 20];
      nmpcSine_B.Selector1_a[nmpcSine_B.i + 20] =
        nmpcSine_B.x_Delay[nmpcSine_B.coffset + 20];
    }

    // Selector: '<S25>/Selector1' incorporates:
    //   Constant: '<S25>/Constant1'

    for (nmpcSine_B.i = 18; nmpcSine_B.i < 19; nmpcSine_B.i++) {
      nmpcSine_B.b_j = static_cast<int32_T>
        (nmpcSine_ConstP.Constant1_Value[nmpcSine_B.i]);
      nmpcSine_B.Selector1_a[nmpcSine_B.i] = nmpcSine_B.x_Delay[nmpcSine_B.b_j -
        1];
      nmpcSine_B.Selector1_a[nmpcSine_B.i + 19] =
        nmpcSine_B.x_Delay[nmpcSine_B.b_j + 20];
    }

    // Delay: '<S25>/x_Delay' incorporates:
    //   Constant: '<S25>/ones'
    //   DataStoreRead: '<S21>/Data Store Read'
    //   Product: '<S25>/Product1'

    if (nmpcSine_DW.icLoad_j) {
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 21; nmpcSine_B.coffset
             ++) {
          nmpcSine_DW.x_Delay_DSTATE[nmpcSine_B.coffset + 21 * nmpcSine_B.i] =
            nmpcSine_DW.x[nmpcSine_B.i];
        }
      }
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset <= 16; nmpcSine_B.coffset +=
           2) {
        // Selector: '<S25>/Selector' incorporates:
        //   Constant: '<S25>/Constant'

        nmpcSine_B.b_j = 19 * nmpcSine_B.i + nmpcSine_B.coffset;

        // Selector: '<S25>/Selector' incorporates:
        //   Constant: '<S25>/Constant'
        //   Delay: '<S25>/x_Delay'

        nmpcSine_B.Selector_b[nmpcSine_B.b_j] = nmpcSine_DW.x_Delay_DSTATE[(21 *
          nmpcSine_B.i + static_cast<int32_T>
          (nmpcSine_ConstP.Constant_Value_a[nmpcSine_B.coffset])) - 1];
        nmpcSine_B.Selector_b[nmpcSine_B.b_j + 1] = nmpcSine_DW.x_Delay_DSTATE
          [(21 * nmpcSine_B.i + static_cast<int32_T>
            (nmpcSine_ConstP.Constant_Value_a[nmpcSine_B.coffset + 1])) - 1];
      }

      // Selector: '<S25>/Selector' incorporates:
      //   Constant: '<S25>/Constant'
      //   Delay: '<S25>/x_Delay'

      for (nmpcSine_B.coffset = 18; nmpcSine_B.coffset < 19; nmpcSine_B.coffset
           ++) {
        nmpcSine_B.Selector_b[nmpcSine_B.coffset + 19 * nmpcSine_B.i] =
          nmpcSine_DW.x_Delay_DSTATE[(21 * nmpcSine_B.i + static_cast<int32_T>
          (nmpcSine_ConstP.Constant_Value_a[nmpcSine_B.coffset])) - 1];
      }
    }

    // Delay: '<S25>/slack_delay' incorporates:
    //   Constant: '<S14>/e.init_zero'

    if (nmpcSine_DW.icLoad_m) {
      nmpcSine_DW.slack_delay_DSTATE = 0.0;
    }

    // MATLAB Function: '<S24>/NLMPC' incorporates:
    //   BusCreator: '<S3>/ParamBus'
    //   Concatenate: '<S3>/Vector Concatenate1'
    //   DataStoreRead: '<S21>/Data Store Read'
    //   Delay: '<S25>/slack_delay'
    //   Selector: '<S25>/Selector'
    //   Selector: '<S25>/Selector1'
    //   UnitDelay: '<S3>/Unit Delay'

    nmpcSine_B.expl_temp.Parameters[0] = nmpcSine_B.params[0];
    nmpcSine_B.expl_temp.Parameters[1] = nmpcSine_B.params[1];
    nmpcSine_B.expl_temp.Parameters[2] = nmpcSine_B.params[2];
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 19; nmpcSine_B.coffset++)
      {
        nmpcSine_B.Selector[nmpcSine_B.i + 6 * nmpcSine_B.coffset] =
          nmpcSine_B.Selector_b[19 * nmpcSine_B.i + nmpcSine_B.coffset];
      }
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      nmpcSine_B.Selector[nmpcSine_B.i + 114] = nmpcSine_B.Selector_b[19 *
        nmpcSine_B.i + 18];
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 19; nmpcSine_B.i++) {
      nmpcSine_B.b_j = nmpcSine_B.i << 1;
      nmpcSine_B.Selector1[nmpcSine_B.b_j] = nmpcSine_B.Selector1_a[nmpcSine_B.i];
      nmpcSine_B.Selector1[nmpcSine_B.b_j + 1] =
        nmpcSine_B.Selector1_a[nmpcSine_B.i + 19];
    }

    nmpcSine_B.Selector1[38] = nmpcSine_B.Selector1_a[18];
    nmpcSine_B.Selector1[39] = nmpcSine_B.Selector1_a[37];
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 4; nmpcSine_B.i++) {
      nmpcSine_B.psi = 0.0;
      for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 40; nmpcSine_B.coffset++)
      {
        nmpcSine_B.psi += static_cast<real_T>(a[(nmpcSine_B.coffset << 2) +
          nmpcSine_B.i]) * nmpcSine_B.Selector1[nmpcSine_B.coffset];
      }

      nmpcSine_B.R[nmpcSine_B.i] = nmpcSine_B.psi;
    }

    memcpy(&nmpcSine_B.z0[0], &nmpcSine_B.Selector[0], 120U * sizeof(real_T));
    nmpcSine_B.z0[120] = nmpcSine_B.R[0];
    nmpcSine_B.z0[121] = nmpcSine_B.R[1];
    nmpcSine_B.z0[122] = nmpcSine_B.R[2];
    nmpcSine_B.z0[123] = nmpcSine_B.R[3];
    nmpcSine_B.z0[124] = nmpcSine_DW.slack_delay_DSTATE;
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 80; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp.OutputWeights[nmpcSine_B.i] = d[nmpcSine_B.i];
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 120; nmpcSine_B.i++) {
      nmpcSine_B.zUB[nmpcSine_B.i] = (rtInf);
    }

    nmpcSine_B.zUB[120] = (rtInf);
    nmpcSine_B.zUB[121] = (rtInf);
    nmpcSine_B.zUB[122] = (rtInf);
    nmpcSine_B.zUB[123] = (rtInf);
    nmpcSine_B.zUB[124] = (rtInf);
    nmpcSine_getXUe(nmpcSine_B.z0, nmpcSine_DW.x, nmpcSine_DW.x_Delay_DSTATE,
                    nmpcSine_B.U0, &nmpcSine_B.psi);
    nmpcSine_B.psi = nmpcSine_costFcn(nmpcSine_DW.x_Delay_DSTATE, nmpcSine_B.U0,
      nmpcSine_B.VectorConcatenate1, nmpcSine_B.params[0], nmpcSine_B.params[1],
      nmpcSine_B.params[2]);
    if (nmpcSine_B.psi <= nmpcSine_B.psi) {
      nmpcSine_B.zUB[124] = 0.0;
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 40; nmpcSine_B.i++) {
      nmpcSine_B.Selector1[nmpcSine_B.i] = -22.5;
      nmpcSine_B.dv1[nmpcSine_B.i] = 20.3;
      nmpcSine_B.dv2[nmpcSine_B.i] = -10.0;
      nmpcSine_B.dv3[nmpcSine_B.i] = 10.0;
      nmpcSine_B.expl_temp.MVScaledTarget[nmpcSine_B.i] = 0.0;
      nmpcSine_B.expl_temp.MVRateMax[nmpcSine_B.i] = 10.0;
      nmpcSine_B.expl_temp.MVRateMin[nmpcSine_B.i] = -10.0;
      nmpcSine_B.expl_temp.MVMax[nmpcSine_B.i] = 20.3;
      nmpcSine_B.expl_temp.MVMin[nmpcSine_B.i] = -22.5;
    }

    nmpcSine_getUBounds(nmpcSine_DW.UnitDelay_DSTATE, nmpcSine_B.Selector1,
                        nmpcSine_B.dv1, nmpcSine_B.dv2, nmpcSine_B.dv3,
                        nmpcSine_B.A_data, nmpcSine_B.A_size,
                        nmpcSine_B.B_data_p, nmpcSine_B.B_size);
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 120; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp.StateMax[nmpcSine_B.i] = (rtInf);
      nmpcSine_B.expl_temp.StateMin[nmpcSine_B.i] = (rtMinusInf);
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 80; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp.OutputMax[nmpcSine_B.i] = (rtInf);
      nmpcSine_B.expl_temp.OutputMin[nmpcSine_B.i] = (rtMinusInf);
    }

    nmpcSine_B.expl_temp.ECRWeight = 100000.0;
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 40; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp.MVRateWeights[nmpcSine_B.i] = 0.1;
      nmpcSine_B.expl_temp.MVWeights[nmpcSine_B.i] = 0.0;
    }

    memcpy(&nmpcSine_B.expl_temp.ref[0], &nmpcSine_B.VectorConcatenate1[0], 80U *
           sizeof(real_T));
    nmpcSine_B.expl_temp.lastMV[0] = nmpcSine_DW.UnitDelay_DSTATE[0];
    nmpcSine_B.expl_temp.lastMV[1] = nmpcSine_DW.UnitDelay_DSTATE[1];
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp.x[nmpcSine_B.i] = nmpcSine_DW.x[nmpcSine_B.i];
    }

    nmpcSine_B.expl_temp_d.PassivityUsePredictedX = true;
    nmpcSine_B.expl_temp_d.OutputPassivityIndex = 0.1;
    nmpcSine_B.expl_temp_d.InputPassivityIndex = 0.0;
    nmpcSine_B.expl_temp_d.UDIndex[0] = 3.0;
    nmpcSine_B.expl_temp_d.MVIndex[0] = 1.0;
    nmpcSine_B.expl_temp_d.UDIndex[1] = 4.0;
    nmpcSine_B.expl_temp_d.MVIndex[1] = 2.0;
    nmpcSine_B.expl_temp_d.NumOfInputs = 4.0;
    nmpcSine_B.expl_temp_d.NumOfOutputs = 4.0;
    nmpcSine_B.expl_temp_d.NumOfStates = 6.0;
    nmpcSine_B.expl_temp_d.PredictionHorizon = 20.0;
    memset(&nmpcSine_B.expl_temp_d.MVTarget[0], 0, 40U * sizeof(real_T));
    memcpy(&nmpcSine_B.expl_temp_d.References[0],
           &nmpcSine_B.VectorConcatenate1[0], 80U * sizeof(real_T));
    nmpcSine_B.expl_temp_d.LastMV[0] = nmpcSine_DW.UnitDelay_DSTATE[0];
    nmpcSine_B.expl_temp_d.LastMV[1] = nmpcSine_DW.UnitDelay_DSTATE[1];
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp_d.CurrentStates[nmpcSine_B.i] =
        nmpcSine_DW.x[nmpcSine_B.i];
    }

    nmpcSine_B.expl_temp_d.Ts = 0.1;
    nmpcSine_B.expl_temp_g.Parameters[0] = nmpcSine_B.params[0];
    nmpcSine_B.expl_temp_g.Parameters[1] = nmpcSine_B.params[1];
    nmpcSine_B.expl_temp_g.Parameters[2] = nmpcSine_B.params[2];
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 40; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp_g.MVScaledTarget[nmpcSine_B.i] = 0.0;
      nmpcSine_B.expl_temp_g.MVRateMax[nmpcSine_B.i] = 10.0;
      nmpcSine_B.expl_temp_g.MVRateMin[nmpcSine_B.i] = -10.0;
      nmpcSine_B.expl_temp_g.MVMax[nmpcSine_B.i] = 20.3;
      nmpcSine_B.expl_temp_g.MVMin[nmpcSine_B.i] = -22.5;
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 120; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp_g.StateMax[nmpcSine_B.i] = (rtInf);
      nmpcSine_B.expl_temp_g.StateMin[nmpcSine_B.i] = (rtMinusInf);
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 80; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp_g.OutputMax[nmpcSine_B.i] = (rtInf);
      nmpcSine_B.expl_temp_g.OutputMin[nmpcSine_B.i] = (rtMinusInf);
    }

    nmpcSine_B.expl_temp_g.ECRWeight = 100000.0;
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 40; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp_g.MVRateWeights[nmpcSine_B.i] = 0.1;
      nmpcSine_B.expl_temp_g.MVWeights[nmpcSine_B.i] = 0.0;
    }

    memcpy(&nmpcSine_B.expl_temp_g.OutputWeights[0],
           &nmpcSine_B.expl_temp.OutputWeights[0], 80U * sizeof(real_T));
    memcpy(&nmpcSine_B.expl_temp_g.ref[0], &nmpcSine_B.VectorConcatenate1[0],
           80U * sizeof(real_T));
    nmpcSine_B.expl_temp_g.lastMV[0] = nmpcSine_DW.UnitDelay_DSTATE[0];
    nmpcSine_B.expl_temp_g.lastMV[1] = nmpcSine_DW.UnitDelay_DSTATE[1];
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      nmpcSine_B.expl_temp_g.x[nmpcSine_B.i] = nmpcSine_DW.x[nmpcSine_B.i];
    }

    for (nmpcSine_B.i = 0; nmpcSine_B.i < 120; nmpcSine_B.i++) {
      nmpcSine_B.dv[nmpcSine_B.i] = (rtMinusInf);
    }

    nmpcSine_B.dv[120] = (rtMinusInf);
    nmpcSine_B.dv[121] = (rtMinusInf);
    nmpcSine_B.dv[122] = (rtMinusInf);
    nmpcSine_B.dv[123] = (rtMinusInf);
    nmpcSine_B.dv[124] = 0.0;
    nmpcSine_fmincon(&nmpcSine_B.expl_temp, &nmpcSine_B.expl_temp_d,
                     nmpcSine_B.z0, nmpcSine_B.A_data, nmpcSine_B.B_data_p,
                     nmpcSine_B.B_size, nmpcSine_B.dv, nmpcSine_B.zUB,
                     &nmpcSine_B.expl_temp_g, nmpcSine_B.z, &nmpcSine_B.s1,
                     &nmpcSine_B.psi, &nmpcSine_B.Out);
    if ((nmpcSine_B.psi == 0.0) && (nmpcSine_B.Out.constrviolation > 0.001)) {
      nmpcSine_B.psi = -2.0;
    }

    // Update for Delay: '<S25>/slack_delay' incorporates:
    //   DataStoreRead: '<S21>/Data Store Read'
    //   MATLAB Function: '<S24>/NLMPC'

    nmpcSine_getXUe(nmpcSine_B.z, nmpcSine_DW.x, nmpcSine_DW.x_Delay_DSTATE,
                    nmpcSine_B.U0, &nmpcSine_DW.slack_delay_DSTATE);

    // MATLAB Function: '<S24>/NLMPC' incorporates:
    //   UnitDelay: '<S3>/Unit Delay'

    if ((nmpcSine_B.psi > 0.0) || (nmpcSine_B.psi == 0.0)) {
      nmpcSine_B.mv[0] = nmpcSine_B.U0[0];
      nmpcSine_B.mv[1] = nmpcSine_B.U0[21];
    } else {
      nmpcSine_B.mv[0] = nmpcSine_DW.UnitDelay_DSTATE[0];
      nmpcSine_B.mv[1] = nmpcSine_DW.UnitDelay_DSTATE[1];
    }

    // DataTypeConversion: '<S15>/Data Type Conversion' incorporates:
    //   MATLAB Function: '<S24>/NLMPC'

    nmpcSine_B.epsilon = floor(nmpcSine_B.psi);
    if (nmpcSine_B.epsilon < 2.147483648E+9) {
      if (nmpcSine_B.epsilon >= -2.147483648E+9) {
        // BusAssignment: '<S15>/Bus Assignment'
        nmpcSine_B.BusAssignment_cw.data = static_cast<int32_T>
          (nmpcSine_B.epsilon);
      } else {
        // BusAssignment: '<S15>/Bus Assignment'
        nmpcSine_B.BusAssignment_cw.data = MIN_int32_T;
      }
    } else {
      // BusAssignment: '<S15>/Bus Assignment'
      nmpcSine_B.BusAssignment_cw.data = MAX_int32_T;
    }

    // End of DataTypeConversion: '<S15>/Data Type Conversion'

    // MATLABSystem: '<S48>/SinkBlock'
    Pub_nmpcSine_264.publish(&nmpcSine_B.BusAssignment_cw);

    // BusAssignment: '<S16>/Bus Assignment'
    nmpcSine_B.BusAssignment_h.x = nmpcSine_B.ref[0];
    nmpcSine_B.BusAssignment_h.y = nmpcSine_B.ref[1];

    // MATLABSystem: '<S50>/SinkBlock'
    Pub_nmpcSine_305.publish(&nmpcSine_B.BusAssignment_h);

    // BusAssignment: '<S17>/Bus Assignment' incorporates:
    //   DataStoreRead: '<S21>/Data Store Read'
    //   DataStoreRead: '<S21>/Data Store Read1'
    //   MATLAB Function: '<S21>/MATLAB Function'

    nmpcSine_B.BusAssignment.x = nmpcSine_DW.x[0];
    nmpcSine_B.BusAssignment.y = nmpcSine_DW.x[1];
    nmpcSine_B.BusAssignment.psi = nmpcSine_DW.x[2];
    nmpcSine_B.BusAssignment.u = nmpcSine_DW.x[3];
    nmpcSine_B.BusAssignment.v = nmpcSine_DW.x[4];
    nmpcSine_B.BusAssignment.r = nmpcSine_DW.x[5];
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset++)
      {
        // Outputs for Atomic SubSystem: '<S10>/Output'
        // MATLAB Function: '<S21>/MATLAB Function'
        nmpcSine_B.psi = 0.0;

        // End of Outputs for SubSystem: '<S10>/Output'
        for (nmpcSine_B.aoffset = 0; nmpcSine_B.aoffset < 6; nmpcSine_B.aoffset
             ++) {
          // Outputs for Atomic SubSystem: '<S10>/Output'
          nmpcSine_B.psi += nmpcSine_DW.P[6 * nmpcSine_B.aoffset + nmpcSine_B.i]
            * nmpcSine_DW.P[6 * nmpcSine_B.aoffset + nmpcSine_B.coffset];

          // End of Outputs for SubSystem: '<S10>/Output'
        }

        // Outputs for Atomic SubSystem: '<S10>/Output'
        // MATLAB Function: '<S21>/MATLAB Function' incorporates:
        //   DataStoreRead: '<S21>/Data Store Read1'

        nmpcSine_B.BusAssignment.covariance[nmpcSine_B.i + 6 *
          nmpcSine_B.coffset] = nmpcSine_B.psi;

        // End of Outputs for SubSystem: '<S10>/Output'
      }
    }

    // End of BusAssignment: '<S17>/Bus Assignment'

    // MATLABSystem: '<S52>/SinkBlock'
    Pub_nmpcSine_253.publish(&nmpcSine_B.BusAssignment);

    // Outputs for Atomic SubSystem: '<S10>/Predict'
    // SignalConversion generated from: '<S22>/MATLAB System' incorporates:
    //   Constant: '<S3>/Constant4'

    nmpcSine_B.TmpSignalConversionAtMATLAB[0] = nmpcSine_B.mv[0];
    nmpcSine_B.TmpSignalConversionAtMATLAB[1] = nmpcSine_B.mv[1];
    nmpcSine_B.TmpSignalConversionAtMATLAB[2] = 0.1;

    // MATLABSystem: '<S22>/MATLAB System' incorporates:
    //   Constant: '<S10>/Q'
    //   DataStoreRead: '<S22>/Data Store ReadP'
    //   DataStoreRead: '<S22>/Data Store ReadX'

    nmpcSine_stateTransitionFcnDT(nmpcSine_DW.x,
      nmpcSine_B.TmpSignalConversionAtMATLAB, nmpcSine_B.z_c);
    for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        nmpcSine_B.imvec[nmpcSine_B.i] = nmpcSine_DW.x[nmpcSine_B.i];
      }

      nmpcSine_B.epsilon = 1.4901161193847656E-8 * fabs
        (nmpcSine_DW.x[nmpcSine_B.b_j]);
      if ((nmpcSine_B.epsilon <= 1.4901161193847656E-8) || rtIsNaN
          (nmpcSine_B.epsilon)) {
        nmpcSine_B.epsilon = 1.4901161193847656E-8;
      }

      nmpcSine_B.imvec[nmpcSine_B.b_j] = nmpcSine_DW.x[nmpcSine_B.b_j] +
        nmpcSine_B.epsilon;
      nmpcSine_stateTransitionFcnDT(nmpcSine_B.imvec,
        nmpcSine_B.TmpSignalConversionAtMATLAB, nmpcSine_B.Pxy);
      for (nmpcSine_B.i = 0; nmpcSine_B.i <= 4; nmpcSine_B.i += 2) {
        tmp = _mm_loadu_pd(&nmpcSine_B.Pxy[nmpcSine_B.i]);
        tmp_0 = _mm_loadu_pd(&nmpcSine_B.z_c[nmpcSine_B.i]);
        _mm_storeu_pd(&nmpcSine_B.A[nmpcSine_B.i + 6 * nmpcSine_B.b_j],
                      _mm_div_pd(_mm_sub_pd(tmp, tmp_0), _mm_set1_pd
          (nmpcSine_B.epsilon)));
      }
    }

    for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
      for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
        nmpcSine_B.aoffset = nmpcSine_B.i * 6 - 1;
        nmpcSine_B.epsilon = 0.0;
        for (nmpcSine_B.b_k = 0; nmpcSine_B.b_k < 6; nmpcSine_B.b_k++) {
          nmpcSine_B.epsilon += nmpcSine_DW.P[(nmpcSine_B.aoffset +
            nmpcSine_B.b_k) + 1] * nmpcSine_B.A[nmpcSine_B.b_k * 6 +
            nmpcSine_B.b_j];
        }

        nmpcSine_B.coffset = 12 * nmpcSine_B.b_j + nmpcSine_B.i;
        nmpcSine_B.y[nmpcSine_B.coffset] = nmpcSine_B.epsilon;
        nmpcSine_B.y[nmpcSine_B.coffset + 6] = nmpcSine_ConstP.Q_Value[6 *
          nmpcSine_B.i + nmpcSine_B.b_j];
      }
    }

    nmpcSine_qr_m32(nmpcSine_B.y, nmpcSine_B.a__1, nmpcSine_B.A);

    // DataStoreWrite: '<S22>/Data Store WriteP' incorporates:
    //   MATLABSystem: '<S22>/MATLAB System'
    //
    for (nmpcSine_B.i = 0; nmpcSine_B.i < 6; nmpcSine_B.i++) {
      for (nmpcSine_B.coffset = 0; nmpcSine_B.coffset < 6; nmpcSine_B.coffset++)
      {
        nmpcSine_DW.P[nmpcSine_B.coffset + 6 * nmpcSine_B.i] = nmpcSine_B.A[6 *
          nmpcSine_B.coffset + nmpcSine_B.i];
      }
    }

    // End of DataStoreWrite: '<S22>/Data Store WriteP'
    for (nmpcSine_B.b_j = 0; nmpcSine_B.b_j < 6; nmpcSine_B.b_j++) {
      // Start for MATLABSystem: '<S22>/MATLAB System' incorporates:
      //   DataStoreRead: '<S22>/Data Store ReadX'

      nmpcSine_B.imvec[nmpcSine_B.b_j] = nmpcSine_DW.x[nmpcSine_B.b_j];
    }

    // Start for MATLABSystem: '<S22>/MATLAB System' incorporates:
    //   DataStoreWrite: '<S22>/Data Store WriteX'

    nmpcSine_stateTransitionFcnDT(nmpcSine_B.imvec,
      nmpcSine_B.TmpSignalConversionAtMATLAB, nmpcSine_DW.x);

    // End of Outputs for SubSystem: '<S10>/Predict'

    // Update for UnitDelay: '<S3>/Unit Delay'
    nmpcSine_DW.UnitDelay_DSTATE[0] = nmpcSine_B.mv[0];
    nmpcSine_DW.UnitDelay_DSTATE[1] = nmpcSine_B.mv[1];

    // Update for Delay: '<S25>/mv_Delay' incorporates:
    //   MATLAB Function: '<S24>/NLMPC'

    nmpcSine_DW.icLoad = false;
    memcpy(&nmpcSine_DW.mv_Delay_DSTATE[0], &nmpcSine_B.U0[0], 42U * sizeof
           (real_T));

    // Update for Delay: '<S25>/x_Delay'
    nmpcSine_DW.icLoad_j = false;

    // Update for Delay: '<S25>/slack_delay'
    nmpcSine_DW.icLoad_m = false;
  }

  // End of Outputs for SubSystem: '<Root>/NMPC Controller'

  // BusAssignment: '<S1>/Bus Assignment'
  memset(&nmpcSine_B.BusAssignment_c, 0, sizeof(SL_Bus_geometry_msgs_Twist));

  // MATLAB Function: '<S1>/ConvertToCmdVel' incorporates:
  //   Constant: '<Root>/Coefficients'

  nmpcSine_B.psi = 5.9653387259350676 * nmpcSine_B.mv[0] + 134.05331619330624;
  if ((nmpcSine_B.psi <= 0.0) || rtIsNaN(nmpcSine_B.psi)) {
    nmpcSine_B.psi = 0.0;
  }

  if (nmpcSine_B.psi >= 255.0) {
    nmpcSine_B.ref[0] = 255.0;
  } else {
    nmpcSine_B.ref[0] = nmpcSine_B.psi;
  }

  nmpcSine_B.psi = 5.9653387259350676 * nmpcSine_B.mv[1] + 134.05331619330624;
  if ((nmpcSine_B.psi <= 0.0) || rtIsNaN(nmpcSine_B.psi)) {
    nmpcSine_B.psi = 0.0;
  }

  if (nmpcSine_B.psi >= 255.0) {
    nmpcSine_B.psi = 255.0;
  }

  nmpcSine_B.psi = 2.0 - (255.0 - nmpcSine_B.psi) * 4.0 / 255.0;
  nmpcSine_B.s1 = 2.0 - (255.0 - nmpcSine_B.ref[0]) * 4.0 / 255.0;

  // Switch: '<S1>/Switch' incorporates:
  //   Switch: '<S1>/Switch1'

  if (nmpcSine_B.controlOn > 0) {
    // BusAssignment: '<S1>/Bus Assignment' incorporates:
    //   Gain: '<S1>/Gain'
    //   MATLAB Function: '<S1>/ConvertToCmdVel'

    nmpcSine_B.BusAssignment_c.linear.x = (nmpcSine_B.psi + nmpcSine_B.s1) / 2.0;
    nmpcSine_B.BusAssignment_c.angular.z = -((nmpcSine_B.s1 - nmpcSine_B.psi) /
      0.84 * 0.42);
  } else {
    // BusAssignment: '<S1>/Bus Assignment' incorporates:
    //   Constant: '<S1>/Constant'
    //   Constant: '<S1>/Constant1'
    //   Gain: '<S1>/Gain'

    nmpcSine_B.BusAssignment_c.linear.x = 0.0;
    nmpcSine_B.BusAssignment_c.angular.z = -0.0;
  }

  // End of Switch: '<S1>/Switch'

  // MATLABSystem: '<S9>/SinkBlock'
  Pub_nmpcSine_168.publish(&nmpcSine_B.BusAssignment_c);
  rate_scheduler((&nmpcSine_M));
}

// Model initialize function
void nmpcSine::initialize()
{
  // Registration code

  // initialize non-finites
  rt_InitInfAndNaN(sizeof(real_T));

  {
    static const char_T prmName[8] = "weights";
    real_T varargin_1;

    // SystemInitialize for Enabled SubSystem: '<Root>/NMPC Controller'
    // Start for DataStoreMemory: '<S10>/DataStoreMemory - P'
    memcpy(&nmpcSine_DW.P[0], &nmpcSine_ConstP.DataStoreMemoryP_InitialValue[0],
           36U * sizeof(real_T));

    // InitializeConditions for Delay: '<S25>/mv_Delay'
    nmpcSine_DW.icLoad = true;

    // InitializeConditions for Delay: '<S25>/x_Delay'
    nmpcSine_DW.icLoad_j = true;

    // InitializeConditions for Delay: '<S25>/slack_delay'
    nmpcSine_DW.icLoad_m = true;

    // SystemInitialize for MATLAB Function: '<S3>/MATLAB Function1'
    nmpcSine_DW.k = 1.0;

    // Start for MATLABSystem: '<S48>/SinkBlock'
    nmpcSine_DW.obj_k.QOSAvoidROSNamespaceConventions = false;
    nmpcSine_DW.obj_k.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj_k.isSetupComplete = false;
    nmpcSine_DW.obj_k.isInitialized = 1;
    nmpcSine_Publisher_setupImpl_m(&nmpcSine_DW.obj_k);
    nmpcSine_DW.obj_k.isSetupComplete = true;

    // Start for MATLABSystem: '<S50>/SinkBlock'
    nmpcSine_DW.obj_i.QOSAvoidROSNamespaceConventions = false;
    nmpcSine_DW.obj_i.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj_i.isSetupComplete = false;
    nmpcSine_DW.obj_i.isInitialized = 1;
    nmpcSine_Publisher_setupImpl_m3(&nmpcSine_DW.obj_i);
    nmpcSine_DW.obj_i.isSetupComplete = true;

    // Start for MATLABSystem: '<S52>/SinkBlock'
    nmpcSine_DW.obj_o.QOSAvoidROSNamespaceConventions = false;
    nmpcSine_DW.obj_o.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj_o.isSetupComplete = false;
    nmpcSine_DW.obj_o.isInitialized = 1;
    nmpcSin_Publisher_setupImpl_m32(&nmpcSine_DW.obj_o);
    nmpcSine_DW.obj_o.isSetupComplete = true;

    // End of SystemInitialize for SubSystem: '<Root>/NMPC Controller'

    // Start for MATLABSystem: '<S4>/SourceBlock'
    nmpcSine_DW.obj_km.QOSAvoidROSNamespaceConventions = false;
    nmpcSine_DW.obj_km.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj_km.isSetupComplete = false;
    nmpcSine_DW.obj_km.isInitialized = 1;
    nmpcSine_Subscriber_setupImpl(&nmpcSine_DW.obj_km);
    nmpcSine_DW.obj_km.isSetupComplete = true;

    // Start for MATLABSystem: '<S5>/SourceBlock'
    nmpcSine_DW.obj_l.QOSAvoidROSNamespaceConventions = false;
    nmpcSine_DW.obj_l.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj_l.isSetupComplete = false;
    nmpcSine_DW.obj_l.isInitialized = 1;
    nmpcSine_Subscriber_setupImpl_m(&nmpcSine_DW.obj_l);
    nmpcSine_DW.obj_l.isSetupComplete = true;

    // Start for MATLABSystem: '<S6>/SourceBlock'
    nmpcSine_DW.obj_p.QOSAvoidROSNamespaceConventions = false;
    nmpcSine_DW.obj_p.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj_p.isSetupComplete = false;
    nmpcSine_DW.obj_p.isInitialized = 1;
    nmpcSin_Subscriber_setupImpl_m3(&nmpcSine_DW.obj_p);
    nmpcSine_DW.obj_p.isSetupComplete = true;

    // Start for MATLABSystem: '<Root>/Get Parameter'
    nmpcSine_DW.obj.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj.isInitialized = 1;
    varargin_1 = -1.0;
    ParamGet_nmpcSine_289.initParam(&prmName[0]);
    ParamGet_nmpcSine_289.setInitialValue(&varargin_1, 1U);
    nmpcSine_DW.obj.isSetupComplete = true;

    // Start for MATLABSystem: '<S9>/SinkBlock'
    nmpcSine_DW.obj_c.QOSAvoidROSNamespaceConventions = false;
    nmpcSine_DW.obj_c.matlabCodegenIsDeleted = false;
    nmpcSine_DW.obj_c.isSetupComplete = false;
    nmpcSine_DW.obj_c.isInitialized = 1;
    nmpcSine_Publisher_setupImpl(&nmpcSine_DW.obj_c);
    nmpcSine_DW.obj_c.isSetupComplete = true;
  }
}

// Model terminate function
void nmpcSine::terminate()
{
  // Terminate for MATLABSystem: '<S4>/SourceBlock'
  if (!nmpcSine_DW.obj_km.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj_km.matlabCodegenIsDeleted = true;
    if ((nmpcSine_DW.obj_km.isInitialized == 1) &&
        nmpcSine_DW.obj_km.isSetupComplete) {
      Sub_nmpcSine_272.resetSubscriberPtr();//();
    }
  }

  // End of Terminate for MATLABSystem: '<S4>/SourceBlock'

  // Terminate for MATLABSystem: '<S5>/SourceBlock'
  if (!nmpcSine_DW.obj_l.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj_l.matlabCodegenIsDeleted = true;
    if ((nmpcSine_DW.obj_l.isInitialized == 1) &&
        nmpcSine_DW.obj_l.isSetupComplete) {
      Sub_nmpcSine_273.resetSubscriberPtr();//();
    }
  }

  // End of Terminate for MATLABSystem: '<S5>/SourceBlock'

  // Terminate for MATLABSystem: '<S6>/SourceBlock'
  if (!nmpcSine_DW.obj_p.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj_p.matlabCodegenIsDeleted = true;
    if ((nmpcSine_DW.obj_p.isInitialized == 1) &&
        nmpcSine_DW.obj_p.isSetupComplete) {
      Sub_nmpcSine_274.resetSubscriberPtr();//();
    }
  }

  // End of Terminate for MATLABSystem: '<S6>/SourceBlock'

  // Terminate for MATLABSystem: '<Root>/Get Parameter'
  if (!nmpcSine_DW.obj.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj.matlabCodegenIsDeleted = true;
  }

  // End of Terminate for MATLABSystem: '<Root>/Get Parameter'

  // Terminate for Enabled SubSystem: '<Root>/NMPC Controller'
  // Terminate for MATLABSystem: '<S48>/SinkBlock'
  if (!nmpcSine_DW.obj_k.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj_k.matlabCodegenIsDeleted = true;
    if ((nmpcSine_DW.obj_k.isInitialized == 1) &&
        nmpcSine_DW.obj_k.isSetupComplete) {
      Pub_nmpcSine_264.resetPublisherPtr();//();
    }
  }

  // End of Terminate for MATLABSystem: '<S48>/SinkBlock'

  // Terminate for MATLABSystem: '<S50>/SinkBlock'
  if (!nmpcSine_DW.obj_i.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj_i.matlabCodegenIsDeleted = true;
    if ((nmpcSine_DW.obj_i.isInitialized == 1) &&
        nmpcSine_DW.obj_i.isSetupComplete) {
      Pub_nmpcSine_305.resetPublisherPtr();//();
    }
  }

  // End of Terminate for MATLABSystem: '<S50>/SinkBlock'

  // Terminate for MATLABSystem: '<S52>/SinkBlock'
  if (!nmpcSine_DW.obj_o.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj_o.matlabCodegenIsDeleted = true;
    if ((nmpcSine_DW.obj_o.isInitialized == 1) &&
        nmpcSine_DW.obj_o.isSetupComplete) {
      Pub_nmpcSine_253.resetPublisherPtr();//();
    }
  }

  // End of Terminate for MATLABSystem: '<S52>/SinkBlock'
  // End of Terminate for SubSystem: '<Root>/NMPC Controller'

  // Terminate for MATLABSystem: '<S9>/SinkBlock'
  if (!nmpcSine_DW.obj_c.matlabCodegenIsDeleted) {
    nmpcSine_DW.obj_c.matlabCodegenIsDeleted = true;
    if ((nmpcSine_DW.obj_c.isInitialized == 1) &&
        nmpcSine_DW.obj_c.isSetupComplete) {
      Pub_nmpcSine_168.resetPublisherPtr();//();
    }
  }

  // End of Terminate for MATLABSystem: '<S9>/SinkBlock'
}

// Constructor
nmpcSine::nmpcSine() :
  nmpcSine_B(),
  nmpcSine_DW(),
  nmpcSine_M()
{
  // Currently there is no constructor body generated.
}

// Destructor
nmpcSine::~nmpcSine()
{
  // Currently there is no destructor body generated.
}

// Real-Time Model get method
RT_MODEL_nmpcSine_T * nmpcSine::getRTM()
{
  return (&nmpcSine_M);
}

const char_T* RT_MODEL_nmpcSine_T::getErrorStatus() const
{
  return (errorStatus);
}

void RT_MODEL_nmpcSine_T::setErrorStatus(const char_T* const volatile
  aErrorStatus)
{
  (errorStatus = aErrorStatus);
}

//
// File trailer for generated code.
//
// [EOF]
//
