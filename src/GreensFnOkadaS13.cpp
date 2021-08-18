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

  // number of cells in either dimension
  int cellsL;
  int cellsW;

  // positions within mesh
  int x2loc;
  int x3loc;
  int y2loc;
  int y3loc;

  // inputs for kernel equation
  double x2;
  double x3;
  double y2;
  double y3;

  // get num cells
  cellsL = _L/_dz;
  cellsW = _W/_dz;

  //printf("i: %d, j: %d, cellsL: %d, cellsW: %d\n", i, j, cellsL, cellsW);

  // receiver loc
  x2loc;
  if (cellsL > 0) {
    x2loc = (i%cellsL);
    if (x2loc == 0) x2loc = cellsL;
  } else {
    x2loc = 1;
  }
   x3loc;
  if (cellsL > 0)
    x3loc = ceil((double)i/cellsL);
  else
    x3loc = i;

  //printf("x2loc: %d, x3loc: %d\n", x2loc, x3loc);

  // src loc
  y2loc;
  if (cellsL > 0) {
    y2loc = j%cellsL;
    if (y2loc == 0) y2loc = cellsL;
  } else {
    y2loc = 1;
  }
  y3loc;
  if (cellsL > 0)
    y3loc = ceil((double)j/cellsL);
  else
    y3loc = j;

  //printf("y2loc: %d, y3loc: %d\n", y2loc, y3loc);

  // for kernel
  x2 = (double)_x(2,x2loc) - 0.5*_dz;
  x3 = (double)_x(3,x3loc) + 0.5*_dz;
  y2 = (double)_x(2,y2loc);
  y3 = (double)_x(3,y3loc);

  double W = _dz;
  double s13 = (_G/(2*M_PI))*( (x2-y2)/(pow((x2-y2),2) + pow((x3-y3),2))
                              -(x2-y2)/(pow((x2-y2),2) - pow((x3+y3),2))
                              -(x2-y2)/(pow((x2-y2),2) + pow((x3-y3-W),2))
                              +(x2-y2)/(pow((x2-y2),2) + pow((x3+y3+W),2))
                              );

  printf("s13: %f\n", s13);

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
