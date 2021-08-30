/* GreensFnShear1212
 * Joseph Wick, 8/16/2021
 *
 * Uses the shear kernel written by Valere Lambert in
 * ViscoelasticCycles and applies it to hmmvp
*/

class GreensFnShear1212 : public ImplGreensFn {
public:
  virtual void Init(const KeyValueFile* kvf) throw (Exception);
  virtual Hd* ComputeHd (double eta) {return NewHd(_x, _x, NULL, eta); }
  virtual bool Call(const CompressBlockInfo& cbi, const vector<UInt>& rs,
                    const vector<UInt>& cs, double* B) const;

private:
  // geometry
  Matd _x;
  Matd _y;

  // size of blocks
  double _dz;

  // size of mesh
  double _L;
  double _W;

  // rigidity
  double _G;

  // depth offset
  double _trans;

  double Eval(UInt i, UInt j) const;
};

inline double GreensFnShear1212::Eval(UInt i, UInt j) const {
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
  y2 = (double)_y(2,i);
  y3 = (double)_y(3,i);
  x2 = (double)_x(2,j) - y2;
  x3 = (double)_x(3,j) - y3;

  double D = (double)_x(3,i) + _trans;

  double s1212 = (_G/M_PI)*( atan((x3-D)/(x2+_L/2))
                           -atan((x3-D)/(x2-_L/2))
                           +atan((x3-D-_W)/(x2-_L/2))
                           -atan((x3-D-_W)/(x2+_L/2))
                           -atan((x3+D+_W)/(x2-_L/2))
                           -atan((x3+D)/(x2+_L/2))
                           +atan((x3+D)/(x2+_L/2))
                           +atan((x3+D+_W)/(x2+_L/2)));

   double p = (x3-(2*D+_W)/2)/_W;
   double bc = -2*_G*((x2/_L +0.5 >= 0)-(x2/_L -0.5 <= 0))*((p+0.5>=0)-(p-0.5>=0));

   //printf("i: %d, j: %d, D: %f, x2: %f, x3: %f, s: %f\n", i,j, D, x2, x3, s1212);

   return s1212+bc;

}

void GreensFnShear1212::Init(const KeyValueFile* kvf) throw (Exception) {
  const Matd* m;
  const Matd* mm;

  if (!kvf->GetMatd("X", m)) throw Exception("Missing X.");
  _x = *m;
  if (_x.Size(1) != 3) throw Exception("X must be 3xN.");

  if (!kvf->GetMatd("Y", mm)) throw Exception("Missing Y.");
  _y = *mm;
  if (_y.Size(1) != 3) throw Exception("Y must be 3xN.");

  kvf->GetDouble("dz", _dz);
  if (_dz <=0) throw Exception("dz must be greater than 0.");
  printf("dz: %f\n", _dz);

  kvf->GetDouble("G", _G);
  if (_G <=0) throw Exception("G must be greater than 0.");

  kvf->GetDouble("W", _W);
  if (_W <= 0) throw Exception("W must be greater than 0.");

  kvf->GetDouble("L", _L);
  if (_L <= 0) throw Exception("L must be greater than 0.");

  kvf->GetDouble("transition", _trans);
  if (_trans < 0) throw Exception("transition depth should be positive");

}

bool GreensFnShear1212::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k=0, ic=0; ic<cs.size(); ic++)
    for (UInt ir = 0; ir<rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
