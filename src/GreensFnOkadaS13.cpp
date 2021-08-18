/* GreensFnS13
 * Joseph Wick, 8/16/2021
 *
 * Uses the Okada solution for s12 written by Valere Lambert in
 * ViscoelasticCycles and applies it to hmmvp
*/

class GreensFnOkadaS13 : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) {return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // geomtry
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

inline double GreensFnOkadaS13::Eval (UInt i, UInt j) const {
  // i is the reiver, j is the source
  // keep in mind that i/j are the cell number not location
  //printf("ij: %d, %d\n", i, j);

  // declaration

  // inputs for kernel equation
  double x2; // receiver
  double x3;
  double y2; // src
  double y3;

  // for kernel; receiver relative to src
  y2 = (double)_x(2,i);
  y3 = (double)_x(3,i);
  x2 = (double)_x(2,j) - 0.5*_dz;
  x3 = (double)_x(3,j) + 0.5*_dz;

  double W = _dz;
  double s13 = (_G/(2*M_PI))*( (x2-y2)/(pow((x2-y2),2) + pow((x3-y3),2))
                              -(x2-y2)/(pow((x2-y2),2) - pow((x3+y3),2))
                              -(x2-y2)/(pow((x2-y2),2) + pow((x3-y3-W),2))
                              +(x2-y2)/(pow((x2-y2),2) + pow((x3+y3+W),2))
                              );

  printf("x2: %f, x3: %f, y2: %f, y3: %f, W: %f, s: %f\n", x2, x3, y2, y3, _dz, s13);

  return s13;
}

void GreensFnOkadaS13::Init(const KeyValueFile* kvf) throw (Exception) {
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

bool GreensFnOkadaS13::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k=0, ic=0; ic<cs.size(); ic++)
    for (UInt ir=0; ir<rs.size(); ir++,k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
