namespace vcg {

// Very simple volume class.
// just an example of the interface that the trivial walker expects

template <class VOX_TYPE>
class MySimpleVolume : public BasicGrid<float>
{
public:
  typedef VOX_TYPE VoxelType;

  const Point3i &ISize() {return siz;}   /// Dimensioni griglia come numero di celle per lato

  float Val(const int &x,const int &y,const int &z) const {
    return cV(x,y,z).V();
    //else return numeric_limits<float>::quiet_NaN( );
  }

  float &Val(const int &x,const int &y,const int &z) {
    return V(x,y,z).V();
    //else return numeric_limits<float>::quiet_NaN( );
  }

  VOX_TYPE &V(const int &x,const int &y,const int &z) {
    return Vol[x+y*siz[0]+z*siz[0]*siz[1]];
  }

  VOX_TYPE &V(const Point3i &pi) {
    return Vol[ pi[0] + pi[1]*siz[0] + pi[2]*siz[0]*siz[1] ];
  }

  const VOX_TYPE &cV(const int &x,const int &y,const int &z) const {
    return Vol[x+y*siz[0]+z*siz[0]*siz[1]];
  }
  bool ValidCell(const Point3i & /*p0*/, const Point3i & /*p1*/) const { return true;}

  template < class VertexPointerType >
  void GetXIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,XAxis>(p1,p2,v,thr); }

  template < class VertexPointerType >
  void GetYIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,YAxis>(p1,p2,v,thr); }

  template < class VertexPointerType >
  void GetZIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  { GetIntercept<VertexPointerType,ZAxis>(p1,p2,v,thr); }

  /// The following members/methods are just for this particular case.
  /// The above one are the one required by the marching cube interface.

  std::vector<VoxelType> Vol;

  typedef enum { XAxis=0,YAxis=1,ZAxis=2} VolumeAxis;

  template < class VertexPointerType,  VolumeAxis AxisVal >
  void GetIntercept(const vcg::Point3i &p1, const vcg::Point3i &p2, VertexPointerType &v, const float thr)
  {
    float f1 = V(p1).V()-thr;
    float f2 = V(p2).V()-thr;
    float u = (float) f1/(f1-f2);
    if(AxisVal==XAxis) v->P().X() = (float) p1.X()*(1-u) + u*p2.X();
    else v->P().X() = (float) p1.X();
    if(AxisVal==YAxis) v->P().Y() = (float) p1.Y()*(1-u) + u*p2.Y();
    else v->P().Y() = (float) p1.Y();
    if(AxisVal==ZAxis) v->P().Z() = (float) p1.Z()*(1-u) + u*p2.Z();
    else v->P().Z() = (float) p1.Z();
    //this->IPfToPf(v->P(),v->P());
    if(VoxelType::HasNormal()) v->N() = V(p1).N()*(1-u) + V(p2).N()*u;
  }



  void Init(Point3i _sz)
  {
    siz=_sz;
    Vol.resize(siz[0]*siz[1]*siz[2]);
  }



};

class MySimpleVoxel
{
private:
  float _v;
public:
  float &V() {return _v;}
  float V() const {return _v;}
  static bool HasNormal() {return false;}
  vcg::Point3f N() const {return Point3f(0,0,0);}
  vcg::Point3f &N()  { static Point3f _p(0,0,0); return _p;}
};
}
