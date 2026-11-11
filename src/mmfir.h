//     mmfir.h: MiniMax FIR filter design per McClellan, Parks & Rabiner
//      Copyright (c) 2014/15/21 Rob Sykes <robs@users.sourceforge.net>
//
// mmfir is free software; you can redistribute it and/or modify it under
// the terms of the GNU General Public License as published by the Free Software
// Foundation; either version 3 of the License, or (at your option) any later
// version.
//
// mmfir is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR
// A PARTICULAR PURPOSE. See the GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License along with
// this program; if not, see <http://www.gnu.org/licenses/>.



#if !defined mmfir_included
#define mmfir_included

#include <stdarg.h>

#if defined __cplusplus
extern "C" {
#endif



typedef enum {
  MmfirSymmetric,          // Type I/II (symmetric impulse-response).
  MmfirAntiSymmetric,      // Type III/IV (anti-symmetric impulse-response).
  MmfirDifferentiator      // Ditto, but with passband error-weighting *= 2/f
} MmfirFilterClass;



typedef struct {
  double freqL, freqR;     // Frequency at left & right band-edges. In [0,1]
  double ampL, ampR;       // Desired amplitude of response at given edges.
  double weightL, weightR; // Error weighting (multiplier) at given edges.

  // Internal use only:
  int ends, weightF, endP, hadRemoval;
  double eAv, ePk;
  double portion0, portion;
} MmfirBandSpec;



typedef enum {
  MmfirSuccess,

  MmfirWarning,            // Section marker
  MmfirGaveUp = MmfirWarning,
  MmfirGaveUp1,            // During first pass.

// No filter will be returned:
  MmfirMissedTarget,
  MmfirError,              // Section marker
  MmfirTooFewPeaks = MmfirError,
  MmfirTooManyPeaks,
  MmfirNumericalError,

  MmfirInvocationError,    // Section marker
  MmfirOutOfMemory = MmfirInvocationError,
  MmfirInvalidNumBands,
  MmfirInvalidOrder,
  MmfirInvalidFrequency,
  MmfirInvalidWeight,
  MmfirInvalidAmplitude,
  MmfirInvalidDcAmplitude,
  MmfirInvalidNyquistAmplitude,

// Internal use only:
  MmfirRuntimeError,
  MmfirOngoing
} MmfirResult;



typedef struct {
  // Useful stuff:
  MmfirResult result;      // As defined above.
  char const * text;       // Textual version of result.
  double minimax;          // Final approximated minimax weighted error.

  // Academic interest:
  int FEs;                 // # of Lagrange interp. fn. evaluations made.
  int iterations;          // # of Remez iterations made.

  double controls[3];      // Allows client to know if these were clamped.
} MmfirReport;



typedef struct {
  double f;    // Frequency in [0,1] of this point
  double a;    // Desired amplitude response @ this point.
  double w;    // Error weight @ this point.
  double e;    // Internal use only.
} MmfirPoint;



typedef struct {
  MmfirPoint * space;
  MmfirPoint * peaks;
  int          spaceLength;
  int          peaksLength;
} MmfirState;



typedef MmfirResult (* MmfirRespFn)(
  MmfirBandSpec     const * bandSpecs,       // As bandSpecs below.
  int               numPoints,
  MmfirPoint        points[/* numPoints */], // Provide .a & .w for each .f
  va_list           additionalArgs);



MmfirReport mmfir(
  double            h[], // In which to store the N filter coeffiecients.
  MmfirFilterClass  filterClass,
  int               N,          // = filter order + 1
  int               numBands,
  MmfirBandSpec     bandSpecs[/* numBands */], // N.B. destroyed.
  double            accuracy,
  double            persistence,
  double            robustness,
  double            target,
  int               flags,
  MmfirState        * reserved, // Set to 0.
  MmfirRespFn       respFn,     // Set to 0 if resp. is given in bandSpecs.
  ...);             // Additional args for (non-zero) respFn.



#if defined __cplusplus
}
#endif

#endif
