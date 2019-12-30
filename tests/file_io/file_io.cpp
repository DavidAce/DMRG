//
// Created by david on 2018-12-16.
//


#include <general/nmspc_tensor_extra.h>
#include <general/nmspc_type_check.h>
#include <iomanip>
#include <iostream>
#include <h5pp/h5pp.h>


/*! \brief Prints the content of a vector nicely */
template<typename T>
std::ostream &operator<<(std::ostream &out, const std::vector<T> &v) {
    if (!v.empty()) {
        out << "[ ";
        std::copy(v.begin(), v.end(), std::ostream_iterator<T>(out, " "));
        out << "]";
    }
    return out;
}

using namespace std::complex_literals;



// Store some dummy data to an output file


int main()
{
    using Scalar = std::complex<double>;

    static_assert(TypeCheck::has_member_data<std::vector<double>>() and "Compile time type-checker failed. Could not properly detect class member data. Check that you are using a supported compiler!");
    static_assert(TypeCheck::has_member_data<std::vector<double>>() and "Compile time type-checker failed. Could not properly detect class member data. Check that you are using a supported compiler!");
    static_assert(TypeCheck::has_member_data<std::vector<double>>() and "Compile time type-checker failed. Could not properly detect class member data. Check that you are using a supported compiler!");

    // Generate dummy data
    std::vector<double> vectorD = { 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0,
                               0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                               -1.0, 0.0, 1.0, 0.0,-1.0, 0.0, 0.0,-1.0, 0.0, 1.0, 0.0, 1.0};
    std::vector<Scalar> vectorC = {-0.191154, 0.964728, -0.0351791,0.177544};

    Eigen::Tensor<Scalar,4> tensorC(3,3,3,3);
    tensorC.setRandom();

    std::string stringX = "Dummy string with spaces";


    std::string output_filename         = "Testing/h5pp_test.h5";
    h5pp::File h5ppFile (output_filename, h5pp::AccessMode::READWRITE, h5pp::CreateMode::TRUNCATE,0);
    std::cout << "Writing vectorD: " <<  vectorD  << std::endl;
    std::cout << "Writing vectorC: " <<  vectorC  << std::endl;
    std::cout << "Writing tensorC: " <<  tensorC  << std::endl;
    std::cout << "Writing stringX: " <<  stringX  << std::endl;

    h5ppFile.writeDataset(vectorD,"vectorD");
    h5ppFile.writeDataset(vectorC,"vectorC");
    h5ppFile.writeDataset(tensorC,"tensorC");
    h5ppFile.writeDataset(stringX,"stringX");
    h5ppFile.writeAttributeToLink(std::string("This is an attribute"), "TestAttr","vectorD");
    h5ppFile.writeAttributeToLink(std::string("This is an attribute"), "TestAttr","vectorC");
    h5ppFile.writeAttributeToLink(std::string("This is an attribute"), "TestAttr","tensorC");
    h5ppFile.writeAttributeToLink(std::string("This is an attribute"), "TestAttr","stringX");


    // Read the data back
    std::vector<double>     vectorD_read;
    std::vector<Scalar>     vectorC_read;
    Eigen::Tensor<Scalar,4> tensorC_read;
    std::string             stringX_read;


    h5ppFile.readDataset(vectorD_read,"vectorD");
    h5ppFile.readDataset(vectorC_read,"vectorC");
    h5ppFile.readDataset(tensorC_read,"tensorC");
    h5ppFile.readDataset(stringX_read,"stringX");
    std::cout << "Reading vectorD: " <<  vectorD_read  << std::endl;
    std::cout << "Reading vectorC: " <<  vectorC_read  << std::endl;
    std::cout << "Reading tensorC: " <<  tensorC_read  << std::endl;
    std::cout << "Reading stringX: " <<  stringX_read  << std::endl;

    if (vectorD != vectorD_read){throw std::runtime_error("vectorD != vectorD_read");}
    if (vectorC != vectorC_read){throw std::runtime_error("vectorC != vectorC_read");}
    if (stringX != stringX_read){throw std::runtime_error("stringX != stringX_read");}

    //Tensor comparison isn't as straightforward
    Eigen::Map<Eigen::VectorXcd> tensorMap(tensorC.data(), tensorC.size());
    Eigen::Map<Eigen::VectorXcd> tensorMapRead(tensorC_read.data(), tensorC_read.size());




//    for (int i = 0; i < tensorC.size(); i++ ){
//        std::cout << tensorC.coeff(i) << "    " << tensorC_read.coeff(i) << std::endl;
//    }

    for (int i = 0; i < tensorC.dimension(0); i++ ) {
        for (int j = 0; j < tensorC.dimension(1); j++ ) {
            for (int k = 0; k < tensorC.dimension(2); k++ ) {
                for (int l = 0; l < tensorC.dimension(3); l++ ) {
                    std::cout << "[ " << i << "  " << j << " " << k << " " << l << " ]: " << tensorC(i,j,k,l) << "   " << tensorC_read(i,j,k,l) << std::endl;
                }
                std::cout << std::endl;
            }
        }
    }



    if (tensorMap != tensorMapRead){throw std::runtime_error("tensorMap != tensorMapRead");}

    return 0;
}