#include <iostream>
#include "courbures.h"
#include <Eigen/Core>
#include <Eigen/Geometry>
#include <Eigen/Dense>


void Courbures::set_fixed_colors()
{
    MyMesh::VertexIter v_it ;

    int i = 0 ;
    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        i++ ;
        if (i % 2)
            _mesh.set_color(*v_it, MyMesh::Color(0, 255, 0)) ;
        else
            _mesh.set_color(*v_it, MyMesh::Color(0, 0, 255)) ;
    }
}

void Courbures::normales_locales()
{
    _mesh.request_vertex_normals() ;

    int i ;
    MyMesh::Normal n ;

    for (MyMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        n = MyMesh::Normal(0,0,0) ;
        i = 0 ;
        for(MyMesh::VertexFaceIter vf_it = _mesh.vf_iter(*v_it);vf_it.is_valid();++vf_it)
        {
            i++ ;
            n += _mesh.calc_face_normal(*vf_it) ;
        }
        if (i!=0)
        {
            n /= i ;
            n = n / norm(n) ;
        }
        _mesh.set_normal(*v_it, n) ;
    }
}

std::vector<MyMesh::VertexHandle> Courbures::get_two_neighborhood(const MyMesh::VertexHandle vh)
{
    OpenMesh::VPropHandleT<bool>    vprop_flag;
    _mesh.add_property(vprop_flag,  "vprop_flag");

    // Initialisation
    for (MyMesh::VertexIter v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        _mesh.property(vprop_flag, *v_it) = false ;
    }
    // Circulateur sur le premier cercle
    std::vector<MyMesh::VertexHandle> neigh, neigh2 ;

	_mesh.property(vprop_flag, vh) = true ;
    for(MyMesh::VertexVertexIter vv_it = _mesh.vv_iter(vh);vv_it.is_valid();++vv_it)
    {
        neigh.push_back(*vv_it) ; // ajout du point à la liste
        _mesh.property(vprop_flag, *vv_it) = true ;
    }

    // Parcours du premier cercle et ajout du second cercle par circulateurs
    for (int i=0; i<neigh.size(); i++)
    {
        MyMesh::VertexHandle vh = neigh.at(i) ;
        for(MyMesh::VertexVertexIter vv_it = _mesh.vv_iter(vh);vv_it.is_valid();++vv_it)
        {
            if (!_mesh.property(vprop_flag, *vv_it)) // sommet non encore rencontré
            neigh2.push_back(*vv_it) ; // ajout du point à la liste
            _mesh.property(vprop_flag, *vv_it) = true ;
        }
    }

    // Concaténation des deux cercles
    neigh.insert(neigh.end(), neigh2.begin(), neigh2.end());

    _mesh.remove_property(vprop_flag);

    return neigh ;
}

MyQuad Courbures::fit_quad(MyMesh::VertexHandle vh)
{
    MyQuad q ;

    std::vector<MyMesh::VertexHandle> neigh = get_two_neighborhood(vh) ;

//    std::cout << "-- Neighborhood" << std::endl ;
    for (int i = 0; i< neigh.size(); i++)
    {
        MyMesh::Point p = _mesh.point(neigh.at(i)) ;
//        std::cout << p[0] << ", " << p[1] << ", " << p[2] << " ; " << std::endl ;
    }

    // Calcul de la matrice de changement de base
    MyMesh::Normal n = _mesh.normal(vh) ;
    MyMesh::Point p = _mesh.point(vh);
//    std::cout << "-- normale" << std::endl ;
//    std::cout << n[0] << ", " << n[1] << ", " << n[2] << std::endl ;
//    std::cout << "-- point" << std::endl ;
//    std::cout << p[0] << ", " << p[1] << ", " << p[2] << std::endl ;

    Eigen::Vector3d ne(n[0], n[1], n[2]), Oz(0,0,1), axis;
    Eigen::Vector3d p_e(p[0], p[1], p[2]), pi_e ;

    axis = ne.cross(Oz) ;
    double sina = axis.norm(), cosa = ne.dot(Oz), angle ;
    if (sina >= 0)
        angle = acos(cosa) ;
    else
        angle = -acos(cosa) ;
    axis = axis.normalized() ;

    Eigen::AngleAxisd r(angle, axis) ;
//    std::cout << "-- rotation" << std:: endl ;
//    std::cout << r.matrix()(0,0) << ", " << r.matrix()(0, 1) << ", " << r.matrix()(0,2) << "; " ;
//    std::cout << r.matrix()(1,0) << ", " << r.matrix()(1, 1) << ", " << r.matrix()(1,2) << "; " ;
//    std::cout << r.matrix()(2,0) << ", " << r.matrix()(2, 1) << ", " << r.matrix()(2,2) << std::endl ;
    Eigen::Translation3d t(-p_e[0], -p_e[1], -p_e[2]) ;
    Eigen::Transform<double, 3, Eigen::Affine> ch_base = r * t ;


    // Calcul de la matrice / vecteur de moindres carrés linéaires (Eigen)

    if (neigh.size() >= 5)
    {
        int n(neigh.size()) ;
        Eigen::MatrixXd A(n,5);
        Eigen::VectorXd B(n);

//        std::cout << "-- Après changement de base" << std::endl ;
        for(int i=0; i<neigh.size(); i++)
        {
            MyMesh::Point p = _mesh.point(neigh.at(i)) ;
            pi_e << p[0], p[1], p[2] ;
            pi_e = ch_base * pi_e ; // Application du changement de base
//            std::cout << pi_e[0] << ", " << pi_e[1] << ", " << pi_e[2] << std::endl ;

		// TODO : initialisation de A / B
            A.row(i) << pi_e.x()*pi_e.x(), pi_e.x()*pi_e.y(), pi_e.y()*pi_e.y(), pi_e.x(), pi_e.y();
            B.row(i) << pi_e.z();
        }

        // Résolution aux moindres carrés par SVD
        Eigen::VectorXd coef(5) ;
        coef = A.bdcSvd(Eigen::ComputeThinU | Eigen::ComputeThinV).solve(B) ;
        MyQuad q(coef[0], coef[1], coef[2], coef[3], coef[4], r, t) ;
//        std::cout << "-- quad : " << std::endl ;
//        std::cout << q[0] << ", " << q[1] << ", " << q[2] << ", " << q[3] << ", " << q[4] << ", " << q[5] << ", " << std::endl ;
        return q ;
    }
    else
    {
        std::cout << "Quad fitting : not enough neighbors" ;
        throw "Quad fitting : not enough neighbors" ;
    }
}

void Courbures::compute_KH()
{
    MyMesh::VertexIter v_it ;
    MyQuad  q ;
    MyMesh::VertexHandle vh ;


    Eigen::Matrix2d Kp;
    Eigen::EigenSolver<Eigen::Matrix2d> solver;
    _mesh.add_property(vprop_K,  "vprop_K");
    _mesh.add_property(vprop_H,  "vprop_H");


    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        vh = *v_it ;
        q = fit_quad(vh) ;

        _mesh.property(vprop_K, vh) = 4 * q[0] * q[2] - q[1]*q[1] ; // K = det (K_P) = 4 a_0 a_2 - a_1 ^ 2
        _mesh.property(vprop_H, vh) = q[0] + q[2] ; // K = Tr (K_P) = a_0 + a_2
    }
}

OpenMesh::Vec3uc Courbures::color_scale_hot(double v){
    double value = v*255;
    value = abs(value);
    return OpenMesh::Vec3uc(255 ,0 + value, 0);
}

OpenMesh::Vec3uc Courbures::color_scale_cold(double v){
    double value = v*255;
    value = abs(value);
    return OpenMesh::Vec3uc(0 ,255 + value, 255);
}

OpenMesh::Vec3uc color_post_cold(double v){
    double value = v*255;
    value = abs(value);
    return OpenMesh::Vec3uc(0 ,255,0+value);
}

OpenMesh::Vec3uc color_pre_hot(double v){
    double value = v*255;
    value = abs(value);
    return OpenMesh::Vec3uc(0+value ,255, 0);
}

void Courbures::set_K_colors(int choice, bool approx)
{
    MyMesh::VertexIter v_it ;
    MyStats<double> dist ;
    OpenMesh::Vec3uc tmp_color ;
    double r = 255, g = 255, b = 255 ;
    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        if(choice == K){
            if(approx){
                dist.push_back(_mesh.property(vprop_K, *v_it));
            }
            else dist.push_back(_mesh.property(vprop_K_p, *v_it));
        }
        else if(choice == H){
            if(approx)
                dist.push_back(_mesh.property(vprop_H, *v_it));
            else dist.push_back(_mesh.property(vprop_H_p, *v_it));
        }
    }
    double m = dist.mean(), sigma = dist.stdev(m) ;
    double tmp, lambda ;
    std::cout << "K min : " << dist.min() << " - max : " << dist.max() << std::endl ;
    std::cout << "K mean : " << m << " - sigma : " << sigma << std::endl;
    for (v_it = _mesh.vertices_begin(); v_it != _mesh.vertices_end(); ++v_it)
    {
        if( choice == K)
            tmp  = _mesh.property(vprop_K, *v_it) - m ;
        else if(choice == H)
            tmp  = _mesh.property(vprop_H, *v_it) - m ;
		// Carte de couleur de tmp (pour les valeurs comprises entre
        // m-sigma / m+sigma
        tmp_color = OpenMesh::Vec3uc(150,150,150);
        if(tmp >= dist.min() && tmp <= m - sigma)
            tmp_color = color_scale_cold(tmp/dist.min());

            //tmp_color = OpenMesh::Vec3uc(0,0, 255);
        if(tmp >=m + sigma&& tmp <= dist.max())
            tmp_color = color_scale_hot(tmp/dist.max());
            //tmp_color = OpenMesh::Vec3uc(255,0,0);
        if(tmp > m-sigma && tmp < m ){
            tmp_color = color_post_cold(tmp/m-sigma);
        }
        if(tmp < m+sigma && tmp >= m ){
            tmp_color = color_pre_hot(tmp/m+sigma);
        }

        if(tmp == m){
            tmp_color = Vec3uc(255 ,255, 255);
        }

        _mesh.set_color(*v_it, MyMesh::Color(tmp_color)) ;
    }
}
