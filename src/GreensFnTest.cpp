
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
  return i+j;
}

void GreensFnTest::Init (const KeyValueFile* kvf) throw (Exception) {
  double d;
  const Matd* m;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (kvf->GetDouble("order", d)) _order = (UInt) d;

  kvf->GetDouble("delta", _delta);
  if (_delta < 0) throw Exception("delta must be >= 0.");

}

bool GreensFnTest::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k = 0, ic = 0; ic < cs.size(); ic++)
    for (UInt ir = 0; ir < rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
