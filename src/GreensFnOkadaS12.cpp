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

  // mesh size
  double _L; //x2 direction
  double _W; //x3 direction

 // rigidity
 double _G;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnOkadaS12::Eval (UInt i, UInt j) const {
  // i is the reeiver, j is the source
  // keep in mind that i/j are the cell number not location

  // args
  int idz = i*_dz;
  int jdz = j*_dz;

  double x2loc = (int)_L%(idz);
  double x3loc = (int)_W%(idz);
  double x2 = (double)_x(2,x2loc) - 0.5*_dz;
  double x3 = (double)_x(3,x3loc) - 0.5*_dz;

  double y2loc = _L%(jdz);
  double y3loc = _W%(jdz);
  double y2 = _x(2,y2loc);
  double y3 = _x(3,y3loc);

  double W = _dz;

  double s12 = (_G/(2*M_PI))*( -(x3-y3)/(pow((x2-y2),2) + pow((x3-y3),2))
                              +(x3+y3)/(pow((x2-y2),2) + pow((x3+y3),2))
                              +(x3-y3-W)/(pow((x2-y2),2) + pow((x3-y3-W),2))
                              -(x3+y3+W)/(pow((x2-y2),2) + pow((x3+y3+W),2)) );

  printf("x2: %f, x3: %f, y2: %f, y3: %f, W: %f, s: %f\n", x2, x3, y2, y3, W, s12);
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

  kvf->GetDouble("W", _W);
  if (_W <= 0) throw Exception("W must be greater than 0.");

  kvf->GetDouble("L", _L);
  if (_L <= 0) throw Exception("L must be greater than 0.");

}

bool GreensFnOkadaS12::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k=0, ic=0; ic<cs.size(); ic++)
    for (UInt ir = 0; ir<rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
