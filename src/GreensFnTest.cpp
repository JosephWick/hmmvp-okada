

class GreensFnTest : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_z, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // geometry
  Matd _x;
  Matd _y;
  Matd _z; // hm sizing only

  double _Ny;
  double _Nz;

  // mesh sizing
  Matd _L;
  Matd _W;

  // rigidity
  double _G;

  // depth offset
  double _trans;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnTest::Eval (UInt i, UInt j) const {
  // i is the reiver, j is the source; both start at 1
  // keep in mind that i/j are the cell number not location
  //printf("ij: %d, %d\n", i, j);

  // declaration

  // inputs for kernel equation
  double x2; // receiver
  double x3;
  double y2; // src
  double y3;

  double L; // block len x2
  double W; // block width x3

  double D; // src depth

  int srcy;
  int srcz;

  srcy = (j%(int)_Ny);
  srcz = (int)(j/(int)_Ny) +1;
  if (srcy == 0) {
    srcy = _Ny;
    srcz -= 1;
  }

  y2 = (double)_y(2, srcy);
  y3 = (double)_y(3, srcz);
  x2 = (double)_x(2, i);
  x3 = (double)_x(3, i);

  return pow( pow(x2-y2,2) + pow(x3-y3,2), 0.5);

}

void GreensFnTest::Init (const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  const Matd* n;
  const Matd* o;

  const Matd* l;
  const Matd* w;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (!kvf->GetMatd("Y", n)) throw Exception("Missing Y.");
  _y = *n;
  if (_y.Size(1) != 3) throw Exception("Y must be 3xN.");

  if (!kvf->GetMatd("Z", o)) throw Exception("Missing Z.");
  _z = *o;
  if (_z.Size(1) != 3) throw Exception("Z must be 3xN.");

  if (!kvf->GetMatd("L", l)) throw Exception("Missing L.");
  _L = *l;

  if (!kvf->GetMatd("W", w)) throw Exception("Missing W.");
  _W = *w;

  kvf->GetDouble("G", _G);
  if (_G <=0) throw Exception("G must be greater than 0.");

  kvf->GetDouble("Ny", _Ny);
  if (_Ny <= 0) throw Exception("Ny must be greater than 0.");

  kvf->GetDouble("Nz", _Nz);
  if (_Nz <= 0) throw Exception("Nz must be greater than 0.");
}

bool GreensFnTest::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
