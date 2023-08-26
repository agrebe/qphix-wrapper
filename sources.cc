#include "sources.h"
#include "qphix/abs_solver.h"

using namespace std;
using namespace QPhiX;

void point_source(double ** qphix_ferm,
                  int x, int y, int z, int t,
                  int spin, int col,
                  int nx, int nt) {
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) qphix_ferm;
  int vol = (nx * nx * nx * nt);      // number of sites
  int max_index = vol / (2 * OUTER_SOALEN); // number of vecs per checkerboard
  // zero out the source
  for (int cb = 0; cb < 2; cb ++) // checkerboard
    for (int index = 0; index < max_index; index ++) // volume
      for (int c = 0; c < 3; c ++) // color
        for (int s = 0; s < 4; s ++) // spin
          for (int r = 0; r < 2; r ++) // real vs imag
            for (int v = 0; v < OUTER_SOALEN; v ++) // index in vector
              psi_s[cb][index][c][s][r][v] = 0;

  // add a 1 to the desired location

  // find index of position
  // first find checkerboard
  int cb = (x + y + z + t) % 2;
  x /= 2;
  // then find vector index
  int v = x % OUTER_SOALEN;
  x /= OUTER_SOALEN;
  // finally get index in lattice
  // x runs fastest, t is slowest
  int index = ((t * nx + z) * nx + y) * nx / (OUTER_SOALEN * 2) + x;
  psi_s[cb][index][col][spin][0][v] = 1;
}

void wall_source(double ** qphix_ferm,
                 int t,
                 int spin, int col,
                 int nx, int nt) {
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) qphix_ferm;
  int vol = (nx * nx * nx * nt);      // number of sites
  int max_index = vol / (2 * OUTER_SOALEN); // number of vecs per checkerboard
  // zero out the source
  for (int cb = 0; cb < 2; cb ++) // checkerboard
    for (int index = 0; index < max_index; index ++) // volume
      for (int c = 0; c < 3; c ++) // color
        for (int s = 0; s < 4; s ++) // spin
          for (int r = 0; r < 2; r ++) // real vs imag
            for (int v = 0; v < OUTER_SOALEN; v ++) // index in vector
              psi_s[cb][index][c][s][r][v] = 0;

  // add a 1 to the desired timeslice

  for (int z = 0; z < nx; z ++) {
    for (int y = 0; y < nx; y ++) {
      for (int x_full = 0; x_full < nx; x_full ++) {
        int x = x_full;
        // find index of position
        // first find checkerboard
        int cb = (x + y + z + t) % 2;
        x /= 2;
        // then find vector index
        int v = x % OUTER_SOALEN;
        x /= OUTER_SOALEN;
        // finally get index in lattice
        // x runs fastest, t is slowest
        int index = ((t * nx + z) * nx + y) * nx / (OUTER_SOALEN * 2) + x;
        psi_s[cb][index][col][spin][0][v] = 1;
      }
    }
  }
}

// parity projectors for sources
void project_positive(double ** qphix_ferm,
                      int nx, int nt) {
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) qphix_ferm;
  int vol = (nx * nx * nx * nt);      // number of sites
  int max_index = vol / (2 * OUTER_SOALEN); // number of vecs per checkerboard
  // zero out the source
  for (int cb = 0; cb < 2; cb ++) // checkerboard
    for (int index = 0; index < max_index; index ++) // volume
      for (int c = 0; c < 3; c ++) // color
          for (int r = 0; r < 2; r ++) // real vs imag
            for (int v = 0; v < OUTER_SOALEN; v ++) { // index in vector
              // upper spin
              double temp = psi_s[cb][index][c][0][r][v] + psi_s[cb][index][c][2][r][v];
              temp /= 2;
              psi_s[cb][index][c][0][r][v] = psi_s[cb][index][c][2][r][v] = temp;
              temp = psi_s[cb][index][c][1][r][v] + psi_s[cb][index][c][3][r][v];
              temp /= 2;
              psi_s[cb][index][c][1][r][v] = psi_s[cb][index][c][3][r][v] = temp;
            }
}

void project_negative(double ** qphix_ferm,
                      int nx, int nt) {
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) qphix_ferm;
  int vol = (nx * nx * nx * nt);      // number of sites
  int max_index = vol / (2 * OUTER_SOALEN); // number of vecs per checkerboard
  // zero out the source
  for (int cb = 0; cb < 2; cb ++) // checkerboard
    for (int index = 0; index < max_index; index ++) // volume
      for (int c = 0; c < 3; c ++) // color
          for (int r = 0; r < 2; r ++) // real vs imag
            for (int v = 0; v < OUTER_SOALEN; v ++) { // index in vector
              // upper spin
              double temp = psi_s[cb][index][c][0][r][v] - psi_s[cb][index][c][2][r][v];
              temp /= 2;
              psi_s[cb][index][c][0][r][v] = temp;
              psi_s[cb][index][c][2][r][v] = -temp;
              temp = psi_s[cb][index][c][1][r][v] - psi_s[cb][index][c][3][r][v];
              temp /= 2;
              psi_s[cb][index][c][1][r][v] = temp;
              psi_s[cb][index][c][3][r][v] = -temp;
            }
}

// convert to spin_mat format
// this will be a tensor with indices (t, z, y, x, c2, c1, s2, s1)
// s1, c1 are source and s2, c2 are sink
// this matches what is written out by QIO after calling chroma
void to_spin_mat(double * output,
    double ** qphix_ferm,
    int s1, int c1,
    int nx, int nt) {
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) qphix_ferm;
  int nvecs = nx / (2 * OUTER_SOALEN);
  int pxy = nvecs * nx;
  int pxyz = pxy * nx;
  // adapted from qdp_packer_parscalar.h
  for (int cb = 0; cb < 2; cb ++) {
    for (int64_t t = 0; t < nt; t++) {
      for (int64_t z = 0; z < nx; z++) {
        for (int64_t y = 0; y < nx; y++) {
          for (int64_t s_index = 0; s_index < nvecs; s_index++) {
            for (int c2 = 0; c2 < 3; c2++) {
              for (int s2 = 0; s2 < 4; s2++) {
                for (int x = 0; x < OUTER_SOALEN; x++) {

                  int ind = t * pxyz + z * pxy + y * nvecs + s_index; //((t*Nz+z)*Ny+y)*nvecs+s;
                  int x_coord = s_index * OUTER_SOALEN + x;

                  // parity bit: cb = (x ^ y ^ z ^ t) & 1, where x is actual x-coordinate
                  // inverting gives (x % 2) = (cb ^ y ^ z ^ t) & 1
                  int pos_index = ((t * nx + z) * nx + y) * nx + x_coord * 2 + ((cb ^ y ^ z ^ t) & 1);
                  for (int i = 0; i < 2; i ++) { // real versus imag
                    int final_index = ((((pos_index * 3 + c2) * 3 + c1) * 4 + s2) * 4 + s1) * 2 + i;
                    output[final_index] = psi_s[cb][ind][c2][s2][i][x];
                  }
                }
              }
            }
          } // s
        } // y
      } // z
    } // t
  } // cb
}

// convert to qc format
// this will be a fermion with indices (t, z, y, x, c, s)
// spin (and real/imag) run fastest
void to_qc_single(float * output,
    double ** qphix_ferm,
    int nx, int nt) {
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) qphix_ferm;
  int nvecs = nx / (2 * OUTER_SOALEN);
  int pxy = nvecs * nx;
  int pxyz = pxy * nx;
  // adapted from qdp_packer_parscalar.h
  for (int cb = 0; cb < 2; cb ++) {
    for (int64_t t = 0; t < nt; t++) {
      for (int64_t z = 0; z < nx; z++) {
        for (int64_t y = 0; y < nx; y++) {
          for (int64_t s_index = 0; s_index < nvecs; s_index++) {
            for (int c = 0; c < 3; c++) {
              for (int s = 0; s < 4; s++) {
                for (int x = 0; x < OUTER_SOALEN; x++) {

                  int ind = t * pxyz + z * pxy + y * nvecs + s_index; //((t*Nz+z)*Ny+y)*nvecs+s;
                  int x_coord = s_index * OUTER_SOALEN + x;

                  // parity bit: cb = (x ^ y ^ z ^ t) & 1, where x is actual x-coordinate
                  // inverting gives (x % 2) = (cb ^ y ^ z ^ t) & 1
                  int pos_index = ((t * nx + z) * nx + y) * nx + x_coord * 2 + ((cb ^ y ^ z ^ t) & 1);
                  for (int i = 0; i < 2; i ++) { // real versus imag
                    int qc_index = ((pos_index * 3 + c) * 4 + s) * 2 + i;
                    output[qc_index] = psi_s[cb][ind][c][s][i][x];
                  }
                }
              }
            }
          } // s
        } // y
      } // z
    } // t
  } // cb
}

// convert to qc format
// this will be a fermion with indices (t, z, y, x, c, s)
// spin (and real/imag) run fastest
void from_qc_single(float * input,
    double ** qphix_ferm,
    int nx, int nt) {
  QPHIX_FERM(psi_s) = (double (**)[3][4][2][OUTER_SOALEN]) qphix_ferm;
  int nvecs = nx / (2 * OUTER_SOALEN);
  int pxy = nvecs * nx;
  int pxyz = pxy * nx;
  // adapted from qdp_packer_parscalar.h
  for (int cb = 0; cb < 2; cb ++) {
    for (int64_t t = 0; t < nt; t++) {
      for (int64_t z = 0; z < nx; z++) {
        for (int64_t y = 0; y < nx; y++) {
          for (int64_t s_index = 0; s_index < nvecs; s_index++) {
            for (int c = 0; c < 3; c++) {
              for (int s = 0; s < 4; s++) {
                for (int x = 0; x < OUTER_SOALEN; x++) {

                  int ind = t * pxyz + z * pxy + y * nvecs + s_index; //((t*Nz+z)*Ny+y)*nvecs+s;
                  int x_coord = s_index * OUTER_SOALEN + x;

                  // parity bit: cb = (x ^ y ^ z ^ t) & 1, where x is actual x-coordinate
                  // inverting gives (x % 2) = (cb ^ y ^ z ^ t) & 1
                  int pos_index = ((t * nx + z) * nx + y) * nx + x_coord * 2 + ((cb ^ y ^ z ^ t) & 1);
                  for (int i = 0; i < 2; i ++) { // real versus imag
                    int qc_index = ((pos_index * 3 + c) * 4 + s) * 2 + i;
                    psi_s[cb][ind][c][s][i][x] = input[qc_index];
                  }
                }
              }
            }
          } // s
        } // y
      } // z
    } // t
  } // cb
}
