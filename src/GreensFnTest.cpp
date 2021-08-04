
#include <iostream>

using namespace std;
extern "C"{
void * dc3d_(char *SPACE, double *ALPHA, double *X, double *Y, double *Z,
              double *DEPTH, double *DIP, double *AL1, double *AL2, double *AW1,
              double *AW2, double *DISL1, double *DISL2, double *DISL3, double *ux,
              double *uy, double *uz, double *uxx, double *uyx, double* uzx,
              double *uxy, double *uyy, double *uzy, double *uxz, double *uyz,
              double *uzz);
}

class GreensFnTest : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
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

  //length/width of fault
  double _len;
  double _wid;

  // dislocations
  double _d1;
  double _d2;
  double _d3;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnTest::Eval (UInt i, UInt j) const {
  //find individual values to the hmatrix here
  // i is the reeiver, j is the source

  // first lets find dimensions of mesh
  //printf("\n__%d, %d, %d__\n", _len, _wid, _dz);
  int meshL = _len/_dz;
  int meshW = _wid/_dz;

  //printf("%d, %d", meshL, meshW);

  // getting parameters
  // for observer position
  //double obsx = ((i % meshL) + 0.5) * _dz;
  //double obsy = ((i % meshW) + 0.5) * _dz;
  //double obsz = 0.0;

  return i+j;
}

void GreensFnTest::Init (const KeyValueFile* kvf) throw (Exception) {
  double* d;
  const Matd* m;
  double tmp;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (kvf->GetDouble("order", d)) _order = (UInt) d;
  printf("order: %d\n", d);

  kvf->GetDouble("delta", _delta);
  if (_delta < 0) throw Exception("delta must be >= 0.");
  printf("delta: %d\n", _delta);

  kvf->GetDouble("dz", d);
  _dz = d;
  printf("\nd: , %d; dz: %d\n", d, _dz);
  if (_dz <=0) throw Exception("dz must be greater than 0.");



}

bool GreensFnTest::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
