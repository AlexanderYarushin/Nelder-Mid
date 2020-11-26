#include <iostream>
#include <cmath>
#include <vector>

using namespace std;


class Vector{
public:
    double x;
    double y;
    Vector(){};
    Vector(double _x, double _y) : x(_x), y(_y){};

    Vector operator +(Vector v){
        Vector n_vec(this->x + v.x, this->y + v.y);
        return n_vec;
    }

    Vector operator - (Vector v){
        Vector n_vec(this->x - v.x, this->y - v.y);
        return n_vec;
    }

    Vector operator * (double scale){
        Vector n_vec(this->x * scale, this->y * scale);
        return n_vec;
    }

    Vector operator /(double scale){
        Vector n_vec(this->x / scale, this->y / scale);
        return n_vec;
    }
};


double f(Vector point){
    double x = point.x;
    double y = point.y;

    return pow(x,2) + (x*y) + pow(y,2) - (6*x) - (9*y);
}

void sort(vector<Vector> &v){
    for(size_t i = 0; i < v.size(); ++i){
        for(size_t j = 0; j < v.size(); ++j){
            if(i != j){
                if(f(v[i]) < f(v[j])){
                    Vector tmp = v[i];
                    v[i] = v[j];
                    v[j] = tmp;
                }
            }
        }
    }
}

void log(Vector v){

    cout << "f(" << v.x << "," << v.y << ") = " << f(v) << endl;

}

Vector NelderMid(){
    double alpha = 1;
    double beta = 0.5;
    double gamma = 2;

    Vector v1(0,0);
    Vector v2(1,0);
    Vector v3(0,1);


    Vector b;
    Vector g;
    Vector w;

    int count = 0;
    while(true){
        vector<Vector> adict = {v1,v2,v3};

        sort(adict);

        b = adict[0];
        g = adict[1];
        w = adict[2];

        Vector mid = (g + b) / 2;


        //Отражение

        Vector xr = mid + ((mid - w) * alpha);

        if(f(xr) < f(g)){
            w = xr;
        }else{
            if(f(xr) < f(w)){
                w = xr;
            }
            Vector c = (w + mid) / 2;
            if(f(c) < f(w)){
                w = c;
            }
        }

        //Растяжение

        if(f(xr) < f(b)){
            Vector xe = mid + (xr - mid) * gamma;

            if(f(xe) < f(w)){
                w = xe;
            }else{
                w = xr;
            }
        }

        //Сжатие

        if(f(xr) > f(g)){
            Vector xc = mid + (w - mid) * beta;
            if(f(xc) < f(w)){
                w = xc;
            }
        }

        v1 = w;
        v2 = g;
        v3 = b;

        count ++;
        if(sqrt(pow(v1.x - v2.x,2) + pow(v1.y - v2.y,2)) < 0.001) break;

    }
    cout << "Iteration count: " << count << endl;
    return b;

}

int main()
{

    Vector result = NelderMid();
    log(result);
    return 0;
}
