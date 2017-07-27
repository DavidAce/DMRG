#include <class_MPS.h>
#include <iomanip>


using namespace std;


void class_MPS::initialize(const long local_dimension_){
    local_dimension = local_dimension_;
    GA.resize(array3{local_dimension,1,1});
    GA.setZero();
    GA(0,0,0) = 1;
    LA.resize(array1{1});
    LA.setConstant(1.0);

    GB.resize(array3{local_dimension,1,1});
    GB.setZero();
    GB(0,0,0) = 1;
    LB.resize(array1{1});
    LB.setConstant(1.0);

    L_tail = LA;

}


Tensor4 class_MPS::get_theta() const {
    return asDiagonal(L_tail) //whatever L_A was in the previous step
            .contract(GA, idxlist1{idx2(1,1)})
            .contract(asDiagonal(LA), idxlist1{idx2(2,0)})
            .contract(GB, idxlist1{idx2(2,1)})
            .contract(asDiagonal(LB), idxlist1{idx2(3,0)})
            .shuffle(array4{1,0,2,3});
}

void class_MPS::swap_AB(){
    swap = !swap;
    L_tail = LA;

    //Swap Gamma
    tmp3    = GA;
    GA      = GB;
    GB      = tmp3;

    //Swap Lambda
    tmp1    = LA;
    LA      = LB;
    LB      = tmp1;

//    L_tail = L.A;
//    G.swap_AB();
//    L.swap_AB();
}

void class_MPS::print_error(const Tensor4 &Hamiltonian) {
    cout << setprecision(16);
    Tensor4 LGLGL = asDiagonal(L_tail)
            .contract(GA, idxlist1{idx2(1,1)})
            .contract(asDiagonal(LA), idxlist1{idx2(2,0)})
            .contract(GB, idxlist1{idx2(2,1)})
            .contract(asDiagonal(LB), idxlist1{idx2(3,0)})
            .shuffle(array4{1,2,0,3});

    Tensor0 result = Hamiltonian.contract(LGLGL.conjugate(), idxlist2{idx2(0,0), idx2(1,1)})
            .contract(LGLGL, idxlist4{idx2(0,0),idx2(1,1),idx2(2,2),idx2(3,3)});

    cout << "E_iDMRG = " << result(0) << endl;
    cout << "E_exact = " << -1.063544409973372 << endl;
    cout << std::scientific;
    cout << "Error   = " << -1.063544409973372 - result(0) << endl;

}
