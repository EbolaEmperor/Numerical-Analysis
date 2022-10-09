#include <iostream>
#include <cmath>
#include "interpolation.h"

int main(){
    HermiteInterpolation poly;
    poly.addCondition(0, 0, 0);
    poly.addCondition(0, 1, 75);
    poly.addCondition(3, 0, 225);
    poly.addCondition(3, 1, 77);
    poly.addCondition(5, 0, 383);
    poly.addCondition(5, 1, 80);
    poly.addCondition(8, 0, 623);
    poly.addCondition(8, 1, 74);
    poly.addCondition(13, 0, 993);
    poly.addCondition(13, 1, 72);

    HermiteInterpolation::setOutput(HermiteInterpolation::OUTPUT_LATEX);
    std::cout << poly << std::endl;

    // 若需要输出 Tikz可识别 的格式， 请将下面的注释取消
    HermiteInterpolation::setOutput(HermiteInterpolation::OUTPUT_TIKZ);
    std::cout << poly << std::endl;

    return 0;
}