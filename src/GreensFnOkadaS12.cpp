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
  //printf("ij: %d, %d\n", i, j);
  // args
  int cellsL = _L/_dz;
  int cellsW = _W/_dz;

  printf("cellsL: %d, cellsW: %d\n", cellsL, cellsW);

  int x2loc = (i%cellsL);
  if (x2loc == 0){ // for when remainer is zero, location is at far end
    x2loc = cellsL;
  }
  if (x2loc < 1){
    x2loc = 1;
  }
  int x3loc;
  if (cellsW > 0)
    x3loc = ceil(i/cellsW);
  else
    x3loc = 1;

  printf("x2loc: %d, x3loc: %d\n", x2loc, x3loc);

  double x2 = (double)_x(2,x2loc) - 0.5*_dz;
  double x3 = (double)_x(3,x3loc) + 0.5*_dz;

  double y2loc = (j%cellsL);
  if (y2loc == 0){
    y2loc = cellsL;
  }
  if (y2loc < 1){
    y2loc = 1;
  }
  double y3loc;
  if (cellsW > 0)
    y3loc = ceil(j/cellsW);
  else
    y3loc = 1;

  printf("y2loc: %d, y3loc: %d\n", y2loc, y3loc);
  double y2 = _x(2,y2loc);
  double y3 = _x(3,y3loc);

  // W replaced here with _dz
  double s12 = (_G/(2*M_PI))*( -(x3-y3)/(pow((x2-y2),2) + pow((x3-y3),2))
                              +(x3+y3)/(pow((x2-y2),2) + pow((x3+y3),2))
                              +(x3-y3-_dz)/(pow((x2-y2),2) + pow((x3-y3-_dz),2))
                              -(x3+y3+_dz)/(pow((x2-y2),2) + pow((x3+y3+_dz),2)) );

  printf("x2: %f, x3: %f, y2: %f, y3: %f, W: %f, s: %f\n", x2, x3, y2, y3, _dz, s12);
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
