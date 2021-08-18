/* GreensFnShear1313
 * Joseph Wick, 8/17/2021
 *
 * Uses the shear kernel written by Valere Lambert in
 * ViscoelasticCycles and applies it to hmmvp
*/

class GreensFnShear1313 : public ImplGreensFn {
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

inline double GreensFnShear1313::Eval(UInt i, UInt j) const {
  // i is the reiver, j is the source
  // keep in mind that i/j are the cell number not location
  //printf("ij: %d, %d\n", i, j);

  // declaration

  // inputs for kernel equation
  double x2; // receiver
  double x3;
  double y2; // src
  double y3;

  double D; // receiver depth

  // for kernel; receiver relative to src
  y2 = (double)_x(2,i);
  y3 = (double)_x(3,i);
  x2 = (double)_x(2,j) - 0.5*_dz - y2;
  x3 = (double)_x(3,j) + 0.5*_dz - y3;

  D = (double)_x(3,i); - 0.5*_dz;

  double s1313 = (_G/M_PI)*( atan((x2+_L/2)/(x3-D))
                           -atan((x2-_L/2)/(x3-D))
                           -atan((x2+_L/2)/(x3-D-_W))
                           +atan((x2-_L/2)/(x3-D-_W))
                           +atan((x2+_L/2)/(x3+D))
                           -atan((x2-_L/2)/(x3+D))
                           -atan((x2+_L/2)/(x3+D+_W))
                           +atan((x2-_L/2)/(x3+D+_W)) );

   double p = (x3-(2*D+_W)/2)/_W;
   double bc = -2*_G*((x2/_L +0.5 >= 0)-(x2/_L -0.5 <= 0))*((p+0.5>=0)-(p-0.5>=0));

   printf("D: %f, L: %f, W: %f, x2: %f, x3: %f, s: %f\n", D, _L, _W, x2, x3, s1313+bc);

   return s1313+bc;

}

void GreensFnShear1313::Init(const KeyValueFile* kvf) throw (Exception) {
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

bool GreensFnShear1313::
Call (const CompressBlockInfo& cbi, const vector<UInt>& rs,
      const vector<UInt>& cs, double* B) const {
  for (UInt k=0, ic=0; ic<cs.size(); ic++)
    for (UInt ir = 0; ir<rs.size(); ir++, k++)
      B[k] = Eval(rs[ir], cs[ic]);
  return true;
}
