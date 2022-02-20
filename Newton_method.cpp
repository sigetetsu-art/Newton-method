#include <iostream>
#include <cmath>
#include <fstream>
#include <vector>
using namespace std;

constexpr double initial_alpha = 1.0;
constexpr double initial_x = 1.0;
constexpr double initial_y = 1.0;
constexpr double threshold = 0.0001;
constexpr double raito = 0.40;//raito:係数αの減衰比率
constexpr double coeff = 0.001;//c1 傾きの減少量

struct grad{
    double fx;
    double fy;
    double fxy;//(1の定数)
    double fxx;//(2の定数)
    double fyy;//(2の定数)
};

inline auto func(double x, double y){
    return x * x +  y * y + x * y - 2 * x;
}

inline auto calc_grad(double x, double y)-> grad{
    return {2 * x + y - 2, 2 * y + x, 1.0, 2.0, 2.0};
}

inline auto calc_grad_norm(grad grad){
    return sqrt(grad.fx * grad.fx + grad.fy * grad.fy);
}

inline auto Armijo_line_serch(double x, double y, grad grad, double alpha, double x_element, double y_element, 
                                                                                            double coeff, double raito)
{
    for(int i = 0; i< 10; i++){
        if(func(x - grad.fx * alpha, y - grad.fy * alpha) <= func(x, y) - coeff * alpha * (grad.fx * x_element + grad.fy * y_element)) break;
        else alpha = raito * alpha;
    }
    return alpha;
}

int main(int ac, char *av[]){
    double x = initial_x;
    double y = initial_y;
    bool loop_flag = true;
    double alpha;
    double det = 0;
    double x_element = 0.0;
    double y_element = 0.0;
    grad grad;

    ofstream fout (av[1]);
    if(!fout){
        cout << "Can't open the file !" << endl;
        exit(8);
    }

    fout << "initial x = " << x << " initila y = " << y << " f(x, y) = " << func(x, y) <<endl;

    while(loop_flag){
        grad = calc_grad(x, y);
        if(calc_grad_norm(grad) < threshold){
            loop_flag = false;
        }
        det = grad.fxx * grad.fyy - grad.fxy * grad.fxy;
        if(det < 0.000001) det = 0.000001;
        x_element = grad.fyy * grad.fx - grad.fxy * grad.fy;
        y_element = -grad.fxy * grad.fx + grad.fxx * grad.fy;
        alpha = Armijo_line_serch(x, y, grad, initial_alpha, x_element, y_element, coeff, raito);
        x = x - alpha * x_element / det;
        y = y - alpha * y_element / det;
        fout << "x = " << x << " y = " << y << " f(x, y) = " << func(x, y) <<endl;
    }
    cout << "x = " << x << " y = " << y << " Extream value f(x, y) = " << func(x, y) <<endl;

    fout.close();
    return 0;
}


