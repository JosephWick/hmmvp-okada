
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
  double _L;
  double _W;

  // dislocations
  double _d1;
  double _d2;
  double _d3;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnTest::Eval (UInt i, UInt j) const {

  int meshL = _L/_dz;
  int meshW = _W/_dz;

  return i+j;
}

void GreensFnTest::Init (const KeyValueFile* kvf) throw (Exception) {
  double d;
  const Matd* m;
  double tmp;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (kvf->GetDouble("order", d)) _order = (UInt) d;

  kvf->GetDouble("delta", _delta);
  if (_delta < 0) throw Exception("delta must be >= 0.");

  _h = 'f';
  kvf->GetDouble("halfspace", tmp);
  if (tmp == 1.0) _h = 'h';

  kvf->GetDouble("mu", _mu);
  if (_mu < 1) throw Exception("mu must be greater than 0.");

  kvf->GetDouble("nu", _nu);
  if (_nu <=0) throw Exception("nu must be greater than 0.");

  tmp = (2*_mu*_nu)/(1-2*_nu);
  _alpha = (tmp+_mu)/(tmp+2*_mu);

  kvf->GetDouble("dz", _dz);
  if (_dz <=0) throw Exception("dz must be greater than 0.");

  kvf->GetDouble("dip", _dip);
  if (_dip <0 || _dip>360) throw Exception("dip must be 0<dip<360.");

  kvf->GetDouble("depth", _depth);

  kvf->GetDouble("L", _L);
  if (_L<0) throw Exception("L must be greater than 0.");

  kvf->GetDouble("W", _W);
  if (_W < 0 ) throw Exception("W must be greater than 0.");

  kvf->GetDouble("d1", _d1);
  kvf->GetDouble("d2", _d2);
  kvf->GetDouble("d3", _d3);


}

bool GreensFnTest::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
