//
// Created by david on 2017-10-04.
//

#ifndef DMRG_CLASS_SVD_H
#define DMRG_CLASS_SVD_H

#include "nmspc_tensor_extra.h"
#include <Eigen/SVD>
#include <iomanip>
#include <Eigen/QR>


template<typename Scalar>
class class_SVD{
private:
    double SVDThreshold         = 1e-14;
    double truncation_error     = 0;
    int chi                     = 0;
    using MatrixType  = Eigen::Matrix<Scalar,Eigen::Dynamic, Eigen::Dynamic, Eigen::ColMajor>;
    Eigen::BDCSVD<MatrixType> SVD;
public:

    class_SVD(){
        setThreshold(SVDThreshold);
    }

    double get_truncation_error();
    void setThreshold(double newThreshold);
    Eigen::Tensor<Scalar, 2> pseudo_inverse(const Eigen::Tensor<Scalar,2> &tensor);

    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,2> &tensor);

    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,2> &tensor, const long chi_max);

    std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
    decompose(const Eigen::Tensor<Scalar,3> &tensor, const long rows,const long cols);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,2> &tensor);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,3> &tensor, long rows, long cols);


    std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,4> &tensor, long chi_max);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt  (const Eigen::Tensor<Scalar,4> &tensor);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>, double >
    schmidt_with_norm  (const Eigen::Tensor<Scalar,4> &tensor);

    std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
    schmidt_unnormalized (const Eigen::Tensor<Scalar,4> &tensor);




};


//
// Definitions
//


template<typename Scalar>
double class_SVD<Scalar>::get_truncation_error(){
    return truncation_error;
}

template<typename Scalar>
void class_SVD<Scalar>::setThreshold(double newThreshold) {
    SVD.setThreshold(newThreshold);
}

template<typename Scalar>
Eigen::Tensor<Scalar, 2>
class_SVD<Scalar>::pseudo_inverse(const Eigen::Tensor<Scalar, 2> &tensor){
    if (tensor.dimension(0) <= 0)  {throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(0)");}
    if (tensor.dimension(1) <= 0)  {throw std::runtime_error("pseudo_inverse error: Dimension is zero: tensor.dimension(1)");}
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    return Textra::Matrix_to_Tensor2(mat.completeOrthogonalDecomposition().pseudoInverse() );
}



template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD<Scalar>::decompose(const Eigen::Tensor<Scalar,2> &tensor) {
    if (tensor.dimension(0) <= 0)  {throw std::runtime_error("decompose error: Dimension is zero: tensor.dimension(0)");}
    if (tensor.dimension(1) <= 0)  {throw std::runtime_error("decompose error: Dimension is zero: tensor.dimension(1)");}
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).squaredNorm();
    long chi = SVD.rank();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(), SVD.matrixV().rows(), chi).shuffle(Textra::array2{1,0})
    );
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD<Scalar>::decompose(const Eigen::Tensor<Scalar,3> &tensor,const long rows,const long cols) {
    if (rows <= 0)  {throw std::runtime_error("decompose error: Dimension is zero: rows");}
    if (cols <= 0)  {throw std::runtime_error("decompose error: Dimension is zero: cols");}
    Eigen::Map<const MatrixType> mat (tensor.data(), rows, cols);
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi = SVD.rank();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(), SVD.matrixV().rows(), chi).shuffle(Textra::array2{1,0})
    );
}



template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 2> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 2> >
class_SVD<Scalar>::decompose(const Eigen::Tensor<Scalar,2> &tensor, const long chi_max) {
    if (tensor.dimension(0) <= 0)  {throw std::runtime_error("decompose error: Dimension is zero: tensor.dimension(0)");}
    if (tensor.dimension(1) <= 0)  {throw std::runtime_error("decompose error: Dimension is zero: tensor.dimension(1)");}

    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).squaredNorm();
    long chi = std::min(SVD.rank(), chi_max);
    return std::make_tuple
            (Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), SVD.matrixU().rows(), chi),
             Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>() , chi),
             Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate() ,  chi, SVD.matrixV().rows() ).shuffle(Textra::array2{1,0})
            );
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,2> &tensor, long d, long chiL, long chi_max, long chiR) {
    if (d <= 0)   {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 0, a.k.a. \"d\"");}
    if (chiL <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 1, a.k.a. \"chiL\"");}
    if (chiR <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 3, a.k.a. \"chiR\"");}

    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chi            = std::min(SVD.rank(),chi_max);

    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chi).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chi), d, chiL, chi),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chi).normalized().template cast<Scalar>(),  chi ),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chi).conjugate(),  d, chiR, chi ).shuffle(Textra::array3{ 0, 2, 1 })
    );
}


template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,2> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    if (dL <= 0)  {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 0, a.k.a. \"dL\"");}
    if (dR <= 0)  {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 2, a.k.a. \"dR\"");}
    if (chiL <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 1, a.k.a. \"chiL\"");}
    if (chiR <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 3, a.k.a. \"chiR\"");}

    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0), tensor.dimension(1));
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(),  chiC ),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC ).shuffle(Textra::array3{ 0, 2, 1 })
    );
}


template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,3> &tensor, long rows,long cols) {
    long d    = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long chiR = tensor.dimension(2);
    if (d <= 0)   {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 0, a.k.a. \"d\"");}
    if (chiL <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 1, a.k.a. \"chiL\"");}
    if (chiR <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 3, a.k.a. \"chiR\"");}
    Eigen::Map<const MatrixType> mat (tensor.data(), rows,cols);
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), d, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(),  chiC ),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  d, chiR, chiC ).shuffle(Textra::array3{ 0, 2, 1 })
    );
}



template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    if (dL <= 0)  {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 0, a.k.a. \"dL\"");}
    if (dR <= 0)  {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 2, a.k.a. \"dR\"");}
    if (chiL <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 1, a.k.a. \"chiL\"");}
    if (chiR <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 3, a.k.a. \"chiR\"");}
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC = SVD.rank();
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(), chiC),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1})
    );
}

template<typename Scalar>
std::tuple<Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3> >
class_SVD<Scalar>::schmidt(const Eigen::Tensor<Scalar,4> &tensor, long chi_max) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    if (dL <= 0)  {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 0, a.k.a. \"dL\"");}
    if (dR <= 0)  {throw std::runtime_error("schmidt error: Dimension is zero: theta dim 2, a.k.a. \"dR\"");}
    if (chiL <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 1, a.k.a. \"chiL\"");}
    if (chiR <= 0){throw std::runtime_error("schmidt error: Dimension is zero: theta dim 3, a.k.a. \"chiR\"");}

    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC            = std::min(SVD.rank(),chi_max);
    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(
            Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
            Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(),  chiC ),
            Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1})
    );
}


template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>,double  >
class_SVD<Scalar>::schmidt_with_norm(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    if (dL <= 0)  {throw std::runtime_error("schmidt_with_norm error: Dimension is zero: theta dim 0, a.k.a. \"dL\"");}
    if (dR <= 0)  {throw std::runtime_error("schmidt_with_norm error: Dimension is zero: theta dim 2, a.k.a. \"dR\"");}
    if (chiL <= 0){throw std::runtime_error("schmidt_with_norm error: Dimension is zero: theta dim 1, a.k.a. \"chiL\"");}
    if (chiR <= 0){throw std::runtime_error("schmidt_with_norm error: Dimension is zero: theta dim 3, a.k.a. \"chiR\"");}

    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC = SVD.rank();
    if (chiC <= 0){
        std::cerr << "M" << mat << std::endl;
        std::cerr << "U: \n" << SVD.matrixU() << std::endl;
        std::cerr << "S: \n" << SVD.singularValues() << std::endl;
        std::cerr << "V: \n" << SVD.matrixV() << std::endl;
        std::cerr << "rank: " << chiC << std::endl;
        throw std::runtime_error("schmidt_with_norm error: Dimension is zero: SVD rank a.k.a. \"chiC\"");
    }

    long num_tail = SVD.nonzeroSingularValues()-chiC;
    truncation_error = num_tail <= 0? 0 : SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
//    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).normalized().template cast<Scalar>(), chiC),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1}),
                           SVD.singularValues().head(chiC).norm()
    );
}

template<typename Scalar>
std::tuple <Eigen::Tensor<Scalar, 3> ,Eigen::Tensor<Scalar, 1>, Eigen::Tensor<Scalar, 3>  >
class_SVD<Scalar>::schmidt_unnormalized(const Eigen::Tensor<Scalar,4> &tensor) {
    long dL   = tensor.dimension(0);
    long chiL = tensor.dimension(1);
    long dR   = tensor.dimension(2);
    long chiR = tensor.dimension(3);
    Eigen::Map<const MatrixType> mat (tensor.data(), tensor.dimension(0)*tensor.dimension(1), tensor.dimension(2)*tensor.dimension(3));
    if(mat.isZero()){std::cout << "WARNING: SVD on zeroed matrix" << std::endl;}
    SVD.compute(mat, Eigen::ComputeThinU | Eigen::ComputeThinV);
    long chiC = SVD.rank();


    truncation_error = SVD.singularValues().tail(SVD.nonzeroSingularValues()-chiC).squaredNorm();
    return std::make_tuple(Textra::Matrix_to_Tensor(SVD.matrixU().leftCols(chiC), dL, chiL, chiC),
                           Textra::Matrix_to_Tensor(SVD.singularValues().head(chiC).template cast<Scalar>(), chiC),
                           Textra::Matrix_to_Tensor(SVD.matrixV().leftCols(chiC).conjugate(),  dR, chiR, chiC).shuffle(Textra::array3{0,2,1})
    );
}
// ============================ //
//   Explicit instantiations
// ============================ //
//
//template class class_SVD<double>;
//template class class_SVD<std::complex<double>>;
//






#endif //DMRG_CLASS_SVD_H
