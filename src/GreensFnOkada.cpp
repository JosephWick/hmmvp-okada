/* GreensFnOkada
 * Joseph Wick, 7/30/2021
 *
 * Uses the fortran script provided in dc3dm/external to fill in an h-matrix with
 * Okada solutions
 *
*/


#include <iostream>

using namespace std;
extern "C"{
void dc3d_(char* SPACE, double* ALPHA, double* X, double* Y, double* Z,
              double* DEPTH, double* DIP, double* AL1, double* AL2, double* AW1,
              double* AW2, double* DISL1, double* DISL2, double* DISL3, double* UX,
              double* UY, double* UZ, double* UXX, double* UYX, double* UZX,
              double* UXY, double* UYY, double* UZY, double* UXZ, double* UYZ,
              double* UZZ);
}


class GreensFnOkada : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) {return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  //geometry
  Matd _x;
  UInt _order;
  double _delta;

  // halfspace
  char _h;

  // elastic properties
  double _mu;
  double _nu;
  double _alpha;

  // size of blocks
  double _dz;

  // angle of fault
  double _dip;

  double _depth;

  // length/width of fault
  double _L;
  double _W;

  // dislocations
  double _d1;
  double _d2;
  double _d3;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnOkada::Eval (UInt i, UInt j) const {
  //find individual values to the hmatrix here
  // i is the reeiver, j is the source

  // first lets find dimensions of mesh
  int meshL = _L/_dz;
  int meshW = _W/_dz;

  double zL = 0;
  double zW = 0;

  // getting parameters
  // for observer position; obs/rec is measured from center
  double obsx = _x(1, i)+(0.5*_dz);
  double obsy = _x(2, i)-(0.5*_dz);
  double obsz = _x(3, i);

  printf("_x(1,i): %d, _x(2,i): %d", _x(1,i), _x(2,i));

  // for source depth; source measured from top left
  double srcdepth = _x(2,j);

  // pointer business
  double tmp1 = 0.5*_L;
  double tmp2 = 0.5*_W;
  char * ph = const_cast<char*>(&_h);
  double * pAlpha = const_cast<double*>(&_alpha);
  double * pDip = const_cast<double*>(&_dip);
  double * pL = const_cast<double*>(&_L);
  double * pW = const_cast<double*>(&_W);
  double * pD1 = const_cast<double*>(&_d1);
  double * pD2 = const_cast<double*>(&_d2);
  double * pD3 = const_cast<double*>(&_d3);

  // outputs
  double ux;
  double uy;
  double uz;
  double uxx;
  double uyx;
  double uzx;
  double uxy;
  double uyy;
  double uzy;
  double uxz;
  double uyz;
  double uzz;

  dc3d_(ph, pAlpha, &obsx, &obsy, &obsz, &srcdepth, pDip, &zL, pL, &zW, pW, pD1, pD2, pD3, &ux, &uy, &uz, &uxx, &uyx, &uzx, &uxy, &uyy, &uzy, &uxz, &uyz, &uzz);
  printf("i: %d, j: %d, obsx: %d, obsy: %d uxy: %d, uyx: %d\n", i, j, obsx, obsy, uxy, uyx);
  return _mu*(uxy + uyx);
}

void GreensFnOkada::Init(const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  double d;
  double tmp;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (kvf->GetDouble("order", d)) _order = (UInt) d;

  kvf->GetDouble("delta", _delta);
  if (_delta < 0) throw Exception("delta must be >= 0.");
  printf("delta: %f\n", _delta);

  _h = 'f';
  kvf->GetDouble("halfspace", tmp);
  if (tmp == 1.0) _h = 'h';
  printf("h: %c\n", _h);

  kvf->GetDouble("mu", _mu);
  if (_mu < 1) throw Exception("mu must be greater than 0.");
  printf("mu: %f\n", _mu);

  kvf->GetDouble("nu", _nu);
  if (_nu <=0) throw Exception("nu must be greater than 0.");
  printf("nu: %f\n", _nu);

  tmp = (2*_mu*_nu)/(1-2*_nu);
  _alpha = (tmp+_mu)/(tmp+2*_mu);
  printf("alpha: %f\n", _alpha);

  kvf->GetDouble("dz", _dz);
  if (_dz <=0) throw Exception("dz must be greater than 0.");
  printf("dz: %f\n", _dz);

  kvf->GetDouble("dip", _dip);
  if (_dip <0 || _dip>360) throw Exception("dip must be 0<dip<360.");
  printf("dip: %f\n", _dip);

  kvf->GetDouble("depth", _depth);
  printf("depth: %f\n", _depth);

  kvf->GetDouble("L", _L);
  if (_L<0) throw Exception("L must be greater than 0.");
  printf("L: %f\n", _L);

  kvf->GetDouble("W", _W);
  if (_W < 0 ) throw Exception("W must be greater than 0.");
  printf("W: %f\n", _W);

  kvf->GetDouble("d1", _d1);
  kvf->GetDouble("d2", _d2);
  kvf->GetDouble("d3", _d3);
  printf("d1: %f, d2: %f, d3: %f", _d1, _d2, _d3);

}

bool GreensFnOkada::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
