/* GreensFnShear1312
 * Joseph Wick 8/16/2021
 *
 * Uses the kernel written by Valere Lambert in ViscoelasticCycles
 * and applies it to hmmvp
*/

class GreensFnShear1312 : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) {return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // geometry
  Matd _x;

  // size of blocks
  double _dz;

  // size of mesh
  double _L;
  double _W;

  // rigidity
  double _G;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnShear1312::Eval (UInt i, UInt j) const {
  // i is the reeiver, j is the source

  double x2 = (double)_x(2,i) - (double)_x(2,j);
  double x3 = (double)_x(3,i) - (double)_x(3,j);

  double s1312 = log(10);

}

void GreensFnShear1312::Init(const KeyValueFile* kvf) throw (Exception) {
  const matd* m;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  kvf->GetDouble("dz", _dz);
  if (_dz <=0) throw Exception("dz must be greater than 0.");
  printf("dz: %f\n", _dz);

  kvf->GetDouble("G", _G);
  if (_G <=0) throw Exception("G must be greater than 0.");

  kvf->GetDouble("W", _W);
  if (_W <= 0) throw Exception("W must be greater than 0.");

  kvf->GetDouble("L", _L);
  if (_L <= 0) throw Exception("L must be greater than 0.");

}

bool GreensFnShear1312::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k=0, ic=0; ic<cs.size(); ic++)
    for (UInt ir = 0; ir<rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
