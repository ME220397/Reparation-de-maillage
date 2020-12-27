#ifndef HOLE_FILLING_H
#define HOLE_FILLING_H
#include <iostream>
#include <cmath>
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>
#include "MeshReconstruction.h"
#include "IO.h"
#include "mainwindow.h"

using namespace std ;
using namespace MeshReconstruction;

// ******************************
// Function computing polynomial P_i (degree 2 polynomials in 3 variables)
inline float myp(int i, MyMesh::Point X) {
    float x = X[0] ;
    float y = X[1] ;
    float z = X[2] ;
    switch (i) {
    case 0:
        return x*x ;
        break;
    case 1:
        return y*y ;
        break;
    case 2:
        return z*z ;
        break;
    case 3:
        return x*y ;
        break;
    case 4:
        return y*z ;
        break;
    case 5:
        return x*z ;
        break;
    case 6:
        return x ;
        break;
    case 7:
        return y ;
        break;
    case 8:
        return z;
        break;
    case 9:
        return 1;
        break;
    default:
        throw std::runtime_error("Error on indice i : unknown polynomial P.");
    }
}

// ******************************
// Function computing thin spline RBF
inline float myphi(float r) { if (r == 0) return 0 ; else return r*r*log(r) ; }

// ******************************
// Class encoding an implicit function built from the basis of RBF
// f(X) = Sum_i=0^n-1 alpha_i * phi(|| X-c_i ||) + Sum_j=0^9 beta_j * P_j(X)
class Implicit_RBF {
private:
    int _n ;
    int _d ;
    vector<float> _alpha, _beta ;
    vector<MyMesh::Point> _center ;

public:
    Implicit_RBF (vector<float> alpha, vector<float> beta, vector<MyMesh::Point> center) ;
    float val(MyMesh::Point) const ;
};

// ******************************
// Class encoding all the functions required for implicit hole filling
class Hole_Filling
{
private:
    MyMesh &_mesh ;
    OpenMesh::VPropHandleT<bool> _vprop ;

    // Constantes servant au rendu
    // TODO : à mettre dans des boutons ...
    const float _scale = 4 ; // Boîte englobante de rendu = _scale * BB des centres
    const float _discr = 1./40. ; // BB discrétisée en 1/_discr voxels
public:
    Hole_Filling(MyMesh &mesh);
    inline ~Hole_Filling() { std::cout << "Ending hole filling" << std::endl ; }

    // Computation of boundary and its neighborhood
    MyMesh::HalfedgeHandle find_boundary_edge() ;
    vector<MyMesh::VertexHandle> find_boundary(MyMesh::HalfedgeHandle heh) ;
    void init_mark_boundary(const vector<MyMesh::VertexHandle> & bnd) ;
    vector<MyMesh::VertexHandle> next_neighbors(const vector<MyMesh::VertexHandle> & bnd) ;

    // Computation of RBF
    pair<pair<Eigen::MatrixXd &,Eigen::VectorXd &>,vector<MyMesh::Point> &> compute_approx_mat(vector<MyMesh::VertexHandle> vlist) ;
    pair<vector<float>&, vector<float>&> solve_approx(const pair<Eigen::MatrixXd &, Eigen::VectorXd &> &p, int n, int d) ;

    // Hole filling
    void fill_hole(string out);

    //Colors
    void set_color_boundary();

    // IO
    void colorize_prop() ;
    void colorize_verts(const vector<MyMesh::VertexHandle> &vlist) ;
    Rect3 estimate_BB(const vector<MyMesh::VertexHandle> &vlist) ;
    void poly_n_out (const Implicit_RBF & implicit, Rect3 domain, std::string filename) ;
};

// Other ....

inline MyMesh::Point min (MyMesh::Point X, MyMesh::Point Y) { MyMesh::Point P ; for (int i=0; i<3; i++) P[i] = min(X[i], Y[i]) ; return P ;}
inline MyMesh::Point max (MyMesh::Point X, MyMesh::Point Y) { MyMesh::Point P ; for (int i=0; i<3; i++) P[i] = max(X[i], Y[i]) ; return P ;}

inline ostream & operator<< (ostream & out, Vec3 v) { out << v.x << ", " << v.y << ", " << v.z ; return out ;}

#endif // HOLE_FILLING_H
