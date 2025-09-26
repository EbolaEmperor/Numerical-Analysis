#include "plot.h"
#include <vector>
#include <cmath>

int main() {
    using namespace plotcpp;
    
    // 生成测试数据
    std::vector<double> x, y1, y2, y3, y4;
    for (double i = 0; i <= 10; i += 0.1) {
        x.push_back(i);
        y1.push_back(sin(i));
        y2.push_back(sin(i + 1));
        y3.push_back(sin(i + 2));
        y4.push_back(sin(i + 3));
    }
    
    try {
        Plot plot(true);
        plot.SetTerminalRaw("pngcairo size 800,600");
        plot.SetOutput("linewidth_demo.png");
        plot.SetTitle("线宽演示");
        plot.SetXLabel("x");
        plot.SetYLabel("y");
        
        // 使用不同线宽和线型绘制多条曲线
        plot.Draw2D(
            Lines(x.begin(), x.end(), y1.begin(), "线宽 1.0", LineType::Solid, 1.0),
            Lines(x.begin(), x.end(), y2.begin(), "线宽 2.0", LineType::Dashed, 2.0),
            Lines(x.begin(), x.end(), y3.begin(), "线宽 3.0", LineType::Dotted, 3.0),
            Lines(x.begin(), x.end(), y4.begin(), "线宽 4.0", LineType::DashDot, 4.0)
        );
        
        plot.Flush();
        
        std::cout << "线宽演示图已保存为 linewidth_demo.png" << std::endl;
        
    } catch (const std::exception& e) {
        std::cerr << "绘图错误: " << e.what() << std::endl;
        return 1;
    }
    
    return 0;
}