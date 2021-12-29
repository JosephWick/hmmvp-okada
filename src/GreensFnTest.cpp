

class GreensFnTest : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) { return NewHd(_x, _y, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  Matd _x;
  Matd _y;
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
  // find individual values to the hmatrix here
  // i is the receiver, j is the source

  double z = (double)_y(3,i);

  return z;
}

void GreensFnTest::Init (const KeyValueFile* kvf) throw (Exception) {
  double d;
  const Matd* m;
  const Matd* n;
  UInt tmp;
  UInt tmp2;
  double d2;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (!kvf->GetMatd("Y", n)) throw Exception("Missing Y.");
  _y = *n;
  if (_y.Size(1) != 3) throw Exception("Y must be 3xN.");

  if (kvf->GetDouble("order", d)) _order = (UInt) d;
  printf("order: %d\n", _order);

  kvf->GetDouble("delta", _delta);
  printf("delta: %f\n", _delta);

  printf(" %d\n", d2);

}

bool GreensFnTest::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
