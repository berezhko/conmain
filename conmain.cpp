#define _USE_MATH_DEFINES
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <strings.h>

#define M 11
 
class model {
    double mag;
    double ang;
    double step;
public:
    model(double modelledMagnitude, double modelledAngle, double modelledFrequency, double samplingFrequency):
        mag(modelledMagnitude),
        ang(modelledAngle),
        step(2.0 * M_PI * modelledFrequency / samplingFrequency)
    {}
    double next() {
        return mag*sin(ang += step);
    }
};


class data {
    double ti[M] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    double vi[M] = {0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0};
    int size = 0;
    double freq;

    double Mul(double t, double *tj, int *ex)
    {
        int j;
        int cont;
        double mul = 1;
    
        for (j = 0; j < size; j++){
            if (ex[j])
                continue;
            mul *= (t-tj[j]);
        }
        return mul;
    }
    
    double P(double t)
    {
        int i;
        double res;
        int *ex;
        
        ex = (int *) malloc(size*sizeof(int));
        bzero(ex, size*sizeof(int));
    
        res = 0;
    
        for (i = 0; i < size; i++){
            bzero(ex, size*sizeof(int));
            ex[i] = 1;
            res += vi[i] * (Mul(t, ti, ex) / Mul(ti[i], ti, ex));
        }
    
        free(ex);
        return res;
    }
    
    double P1(double t)
    {
        int i, k;
        double res, res2;
        int *ex;
        
        ex = (int *) malloc(size*sizeof(int));
        bzero(ex, size*sizeof(int));
    
        res = 0;
        res2 = 0;
    
        for (i = 0; i < size; i++){
            ex[i] = 1;
            res2 = 0;
            for (k = 0; k < size; k++){
                if (k == i)
                    continue;
                ex[k] = 1;
                res2 += Mul(t, ti, ex);
                ex[k] = 0;
            }
            res += vi[i] * (res2 / Mul(ti[i], ti, ex));
            ex[i] = 0;
        }
    
        free(ex);
        return res;
    }
    
    double P2(double t)
    {
        int i, k, p;
        double res, res2, res3;
        int *ex;
        
        ex = (int *) malloc(size*sizeof(int));
        bzero(ex, size*sizeof(int));
    
        res = 0;
        res2 = 0;
        res3 = 0;
    
        for (i = 0; i < size; i++){
            ex[i] = 1;
            res2 = 0;
            for (k = 0; k < size; k++){
                if (k == i)
                    continue;
                ex[k] = 1;
                res3 = 0;
                for (p = 0; p < size; p++){
                    if (p == i || p == k)
                        continue;
                    ex[p] = 1;
                    res3 += Mul(t, ti, ex);
                    ex[p] = 0;
                }
                res2 += res3;
                ex[k] = 0;
            }
            res += vi[i] * (res2 / Mul(ti[i], ti, ex));
            ex[i] = 0;
        }
    
        free(ex);
        return res;
    }
    
    double P3(double t)
    {
        int i, k, p, l;
        double res, res2, res3, res4;
        int *ex;
        
        ex = (int *) malloc(size*sizeof(int));
        bzero(ex, size*sizeof(int));
    
        res = 0;
        res2 = 0;
        res3 = 0;
        res4 = 0;
    
        for (i = 0; i < size; i++){
            ex[i] = 1;
            res2 = 0;
            for (k = 0; k < size; k++){
                if (k == i)
                    continue;
                ex[k] = 1;
                res3 = 0;
                for (p = 0; p < size; p++){
                    if (p == i || p == k)
                        continue;
                    ex[p] = 1;
                    res4 = 0;
                    for (l = 0; l < size; l++){
                        if (l == i || l == k || l == p)
                            continue;
                        ex[l] = 1;
                        res4 += Mul(t, ti, ex);
                        ex[l] = 0;
                    }
                    res3 += res4;
                    ex[p] = 0;
                }
                res2 += res3;
                ex[k] = 0;
            }
            res += vi[i] * (res2 / Mul(ti[i], ti, ex));
            ex[i] = 0;
        }
    
        free(ex);
        return res;
    }
    
    double P4(double t)
    {
        int i, k, p, l, m;
        double res, res2, res3, res4, res5;
        int *ex;
        
        ex = (int *) malloc(size*sizeof(int));
        bzero(ex, size*sizeof(int));
    
        res = 0;
        res2 = 0;
        res3 = 0;
        res4 = 0;
        res5 = 0;
    
        for (i = 0; i < size; i++){
            ex[i] = 1;
            res2 = 0;
            for (k = 0; k < size; k++){
                if (k == i)
                    continue;
                ex[k] = 1;
                res3 = 0;
                for (p = 0; p < size; p++){
                    if (p == i || p == k)
                        continue;
                    ex[p] = 1;
                    res4 = 0;
                    for (l = 0; l < size; l++){
                        if (l == i || l == k || l == p)
                            continue;
                        ex[l] = 1;
                        res5 = 0;
                        for (m = 0; m < size; m++){
                            if (m == i || m == k || m == p || m == l)
                                continue;
                            ex[m] = 1;
                            res5 += Mul(t, ti, ex);
                            ex[m] = 0;
                        }
                        res4 += res5;
                        ex[l] = 0;
                    }
                    res3 += res4;
                    ex[p] = 0;
                }
                res2 += res3;
                ex[k] = 0;
            }
            res += vi[i] * (res2 / Mul(ti[i], ti, ex));
            ex[i] = 0;
        }
    
        free(ex);
        return res;
    }

public:
    data(double samplingFrequency)
    {
        freq = 1/samplingFrequency;
    }

    void add(double v)
    {
        if (size < M){
            vi[size] = v;
            ti[size] = (size+1)*freq;
            size++;
        }else{
            for (int i = 0; i < M-1; i++){
                vi[i] = vi[i+1];
                ti[i] = ti[i+1];
            }
            vi[M-1] = v;
            ti[M-1] += freq;
        }
    }

    double computedMag()
    {
        if (size < M)
            return 0;

        double t = ti[size/2];
        double f = 1/2.0/M_PI*sqrt(-P2(t)/P(t));
        return sqrt((pow(2*M_PI*f*P3(t), 2) + pow(P4(t), 2))/pow(2*M_PI*f, 8));
    }

    double computedFreq()
    {
        if (size < M)
            return 0;

        double t = ti[size/2];
        return 1/2.0/M_PI*sqrt(-P2(t)/P(t));
    }

    double computedAng()
    {
        if (size < M)
            return 0;

        static double pi;
        static int tmp;
        double t = ti[size/2];
        double f = 1/2.0/M_PI*sqrt(-P2(t)/P(t));
        double v = sqrt((pow(2*M_PI*f*P3(t), 2) + pow(P4(t), 2))/pow(2*M_PI*f, 8));

        if (P2(t) < 0){
            tmp = 0;
            return pi + acos(P1(t)/v/2/M_PI/f) - 2*M_PI*f*t;
        }else{
            if (!tmp)
                pi += 2*M_PI;
            tmp = 1;
            return pi - acos(P1(t)/v/2/M_PI/f) - 2*M_PI*f*t;
        }
    }
};


int main(int argc, char **argv) {
    int c = 0;
    const double samplingFrequency = 4000.0;
    const double modelMag = 380;
    const double modelAng = 0.25 * M_PI;
    const double modelFreq = 350.123456789;
 
    model mdl(modelMag, modelAng, modelFreq, samplingFrequency);
    data dat(samplingFrequency);
    do {
        for (int i=0; i < 20; ++i) {
            double samp = mdl.next();
            dat.add(samp);
            double computedMag = dat.computedMag();
            double computedAng = dat.computedAng();
            double computedFreq = dat.computedFreq();
            printf("samp %f, mag %f, ang %f, freq %f\n", samp, computedMag, computedAng, computedFreq);
        }
        c++;
    } while (c < 30);
}
