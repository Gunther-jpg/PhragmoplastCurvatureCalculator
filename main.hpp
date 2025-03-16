//
// Created by Hunter Whitlock on 3/8/25.
//
#pragma once
#define PHRAGMOPLASTCURVATURECALCULATOR_MAIN_H
#ifndef main_h

#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <filesystem>
#include <nlopt.hpp>
#include <math.h>

struct XYdata {
    std::string filename;
    std::vector<double> XData = {};
    std::vector<double> YData = {};
    std::vector<double> Coefficients = {};
    struct XYdata* pNext = nullptr;
};

void txtToXYData(std::ifstream &input, struct XYdata &XYdata) {
    std::string line;
    while(input.peek() != -1 and !input.eof()) {

        input >> line;
        XYdata.XData.push_back(std::stod(line));

        input >> line;
        XYdata.YData.push_back(std::stod(line));
    }
}

void importData(struct XYdata* &pHead, const std::filesystem::path &DATA_DIRECTORY) {
    struct XYdata *pCur = pHead;

    for (auto const& dir_entry : std::filesystem::directory_iterator{DATA_DIRECTORY}) {

        std::ifstream txtFile(dir_entry.path()); //opens the current file

        if(pCur != pHead) {
            struct XYdata *newEntry = new struct XYdata;
            pCur->pNext = newEntry; //sets the previous XYdata struct to point to the new one
            txtToXYData(txtFile, *pCur);
            pCur = newEntry; //updates pCure for next iteration

        } else {
            txtToXYData(txtFile, *pHead);
        }

    }
}

void moveXYVectorToOrigin(std::vector<double> &X, std::vector<double> &Y) {
    double X0 = X.front(), Y0 = Y.front();

    for(int i = 0; i < X.size(); i++) {
        X.at(i) = X.at(i) - X0;
        Y.at(i) = Y.at(i) - Y0;
    }

}

double standardQuadraticWithRotatedDataPoints(const std::vector<double> &x,std::vector<double> &grad, void* f_data) {
    struct XYdata *data = (struct XYdata *)f_data;
    double chi_squared = 0, observedValue = 0, expectedValue = 0;

    return 0.0;
}

double standardQuadratic(const std::vector<double> &x, void* f_data) { //need to add gradients back
    struct XYdata *data = (struct XYdata *)f_data;

    double chi_squared = 0, observedValue = 0, expectedValue = 0;

    for(int i = 0; i < data->XData.size(); i++) {
        observedValue = data->YData.at(i);
        expectedValue = x.at(0) * pow(data->XData.at(i),2) + x.at(1) * data->XData.at(i) + x.at(2); //ax^2 +bx + c

        chi_squared += abs(pow(observedValue - expectedValue, 2) / expectedValue); //abs is to encourage minimization of the magnitude of the residuals
    }

    return chi_squared;
}

double generalParabola(const std::vector<double> &x, void* f_data) {
    struct XYdata *data = (struct XYdata *)f_data;

    double leastSquareResidual = 0;
    leastSquareResidual = x.at(0);
    leastSquareResidual = 0;

    for(int i = 0; i < data->XData.size(); i++) {
        leastSquareResidual += pow(x.at(0) * pow(data->XData.at(i),2) + x.at(1) * data->YData.at(i) * data->XData.at(i) + x.at(2) * pow(data->YData.at(i),2)
                           + x.at(3) * data->XData.at(i) + x.at(4) * data->YData.at(i) + data->XData.at(i),2); //ax^2 +bxy+cy^2+dx+ey+f=0
    }
    leastSquareResidual /= data->XData.size();

    return abs(leastSquareResidual); //abs is to encourage minimization of the magnitude of the residuals
}
//need to add a restriction function or two

double nonlinearConstraint(const std::vector<double> &x, std::vector<double> &grad, void *data) {
    return pow(x.at(1),2) - 4*x.at(0)*x.at(2);
}



#endif //PHRAGMOPLASTCURVATURECALCULATOR_MAIN_H
