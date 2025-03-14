//
// Created by Hunter Whitlock on 3/8/25.
//

#include "main.hpp"

int main() {
    const std::__fs::filesystem::path DATA_DIRECTORY {"../Data"};
    struct XYdata *dataList = new struct XYdata;
    struct XYdata *testData = new struct XYdata;

    testData->XData = {0, 1, 2, 3, 4, 5, 6};
    testData->YData = {1, 2, 4, 8, 16, 32, 64};

    std::vector<double> coefficients;
    coefficients = {1,4,0};
    double optimizedValue = -1;

    importData(dataList, DATA_DIRECTORY);


    struct XYdata *tmpDataListPtr = dataList;
    while(tmpDataListPtr != nullptr) {
        moveXYVectorToOrigin(tmpDataListPtr->XData,tmpDataListPtr->YData);
        tmpDataListPtr = tmpDataListPtr->pNext;
    }

    std::cout << "typical least squares residual: " <<  standardQuadratic(coefficients, (void*)(testData)) << std::endl;

    nlopt::opt opt(nlopt::algorithm::LN_COBYLA, 3);
    opt.set_min_objective((nlopt::func) standardQuadratic, (void*)(testData));
    //opt.set_ftol_abs(.01);
    //opt.set_maxtime(15);
    std::vector<double> lb = {-10,-10,-10};
    std::vector<double> ub = {10,10,10};
//    opt.set_lower_bounds(lb);
//    opt.set_upper_bounds(ub);

    opt.set_maxeval(1000);
    //opt.add_equality_constraint(nonlinearConstraint, (void*)(dataList), .0001);

    try {
        nlopt::result result = opt.optimize(coefficients, optimizedValue);
    } catch(...) {
        //std::cout << "typical least squares residual: " <<  standardQuadratic(coefficients, (void*)(dataList)) << std::endl;
    }




    std::cout << "typical least squares residual: " <<  standardQuadratic(coefficients, (void*)(testData)) << std::endl;
        for(int i = 0; i < 3; i++) {
        std::cout << coefficients.at(i) << ' ';
    }
        std::cout << " num evals " << opt.get_numevals() << std::endl;

    return 0;
}