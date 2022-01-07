/* GreensFnOkadaS12
 * Joseph Wick, 8/5/2021
 *
 * Uses the Okada solution for s12 written by Valere Lambert in
 * QDBIM2D and applies it to hmmvp
*/

class GreensFnOkadaS12 : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) {return NewHd(_y, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // geometry of fault
  Matd _x;
  Matd _y;

  // mesh sizing
  Matd _L;
  Matd _W;

  // rigidity
  double _G;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnOkadaS12::Eval (UInt i, UInt j) const {
  // i is the receiver, j is the source
  // keep in mind that i/j are the cell number not location
  //printf("ij: %d, %d\n", i, j);

  // declaration

  // inputs for kernel equation
  double x2; // receiver
  double x3;
  double y2; // src
  double y3;

  double L;
  double W;

  // for kernel; receiver relative to src
  y2 = (double)_y(2,j);
  y3 = (double)_y(3,j);
  x2 = (double)_x(2,i);
  x3 = (double)_x(3,i);

  L = _L(1, j);
  W = _W(1, j);;

  double s12 = (_G/(2*M_PI))*(-(x3-y3)/(pow((x2-y2),2)     + pow((x3-y3),2))
                              +(x3+y3)/(pow((x2-y2),2)     + pow((x3+y3),2))
                              +(x3-y3-W)/(pow((x2-y2),2) + pow((x3-y3-W),2))
                              -(x3+y3+W)/(pow((x2-y2),2) + pow((x3+y3+W),2)) );

  //printf("x2: %f, x3: %f, y2: %f, y3: %f, W: %f, s: %f\n", x2, x3, y2, y3, _dz, s12);

  //if (isinf(s12)){
    //s12 = -999;
    /*if (i<400 && j>400)*/
    //printf("inf at i=%d, j=%d\n", i, j);
  //}

  return s12;
}

void GreensFnOkadaS12::Init(const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (!kvf->GetMatd("Y", m)) throw Exception("Missing Y.");
  _y = *m;
  if (_y.Size(1) != 3) throw Exception("Y must be 3xN.");

  if (!kvf->GetMatd("L", l)) throw Exception("Missing L.");
  _L = *l;

  if (!kvf->GetMatd("W", w)) throw Exception("Missing W.");
  _W = *w;

  kvf->GetDouble("G", _G);
  if (_G <=0) throw Exception("G must be greater than 0.");

}

bool GreensFnOkadaS12::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k=0, ic=0; ic<cs.size(); ic++)
    for (UInt ir = 0; ir<rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
