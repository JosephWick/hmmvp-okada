/* GreensFnOkadaS12
 * Joseph Wick, 8/5/2021
 *
 * Uses the Okada solution for s12 written by Valere Lambert in
 * QDBIM2D and applies it to hmmvp
*/

class GreensFnOkadaS12 : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) {return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // geometry of fault
  Matd _x;

  // size of blocks
  double _dz;

 //rigidity
 double _G;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnOkadaS12::Eval (UInt i, UInt j) const {

  // args
  double x2 = 0.0; //(double)_x(1,i) + 0.5*_dz;
  double x3 = (double)_x(2,i) - 0.5*_dz;

  double y2 = 0.0; //_x(1,j);
  double y3 = _x(2,j);

  double W = _dz;

  double s12 = (_G/(2*M_PI))*( -(x3-y3)/(pow((x2-y2),2) + pow((x3-y3),2))
                              +(x3+y3)/(pow((x2-y2),2) + pow((x3+y3),2))
                              +(x3-y3-W)/(pow((x2-y2),2) + pow((x3-y3-W),2))
                              -(x3+y3+W)/(pow((x2-y2),2) + pow((x3+y3+W),2)) );

  return s12;
}

void GreensFnOkadaS12::Init(const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  kvf->GetDouble("dz", _dz);
  if (_dz <=0) throw Exception("dz must be greater than 0.");
  printf("dz: %f\n", _dz);

  kvf->GetDouble("G", _G);
  if (_G <=0) throw Exception("G must be greater than 0.");

  double t;
  kvf->GetDouble("tol", t);
  printf("tol: %f", t);

}

bool GreensFnOkadaS12::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k=0, ic=0; ic<cs.size(); ic++)
    for (UInt ir = 0; ir<rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
