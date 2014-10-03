#include "conmain.h"

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
    long double *alfa;
    int size = 0;
    double freq;

public:
    data(double samplingFrequency)
    {
        freq = 1/samplingFrequency;
        alfa = (long double *) malloc(M*sizeof(long double));
    }

    //Добавляем новый элемент в массивы данных
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
            free(alfa);
            //Расчет коэффициентов представления методом Гаусса
            if ((alfa = gauss(M, ti, vi)) == NULL){
                fprintf(stderr, "[ERROR GAUSS] Ошибка в определении коэфициентов\n")
                exit(10);
            }
        }
    }

    long double P(double t, int p)
    {
        int i, j, k;
        long double res;

        if (size < M)
            return 0;

        res = 0;
        for (i = p; i < M; i++){
            k = 1;
            for (j = 0; j < p; j++)
                k *= (i-j);
            res += k*alfa[i]*powl(t, i-p);
        }

        return res;
    }

    double pi(double a)
    {
        int c;

        c = a/2/M_PI;
        a -= 2*M_PI*c;
        if (a <= M_PI)
            return a;
        else
            return a-2*M_PI;
    }

    //Вычисляем амплитуду
    double computedMag()
    {
        if (size < M)
            return 0;

        double t = ti[size/2];
        double f = 1/2.0/M_PI*sqrt(-P(t, 2)/P(t, 0));
        return sqrt((pow(2*M_PI*f*P(t, 3), 2) + pow(P(t, 4), 2))/pow(2*M_PI*f, 8));
    }

    //Вычисляем частоту
    double computedFreq()
    {
        if (size < M)
            return 0;

        double t = ti[size/2];
        return 1/2.0/M_PI*sqrt(-P(t, 2)/P(t, 0));
    }

    //Вычисляем фазу
    double computedAng()
    {
        if (size < M)
            return 0;

        double a;
        double t = ti[size/2];
        double f = 1/2.0/M_PI*sqrt(-P(t, 2)/P(t, 0));
        double v = sqrt((pow(2*M_PI*f*P(t, 3), 2) + pow(P(t, 4), 2))/pow(2*M_PI*f, 8));

        a = pi(2*M_PI*f*t);

        return atan(-P(t, 4)/P(t, 3)/2/M_PI/f) - a;
    }

    double computedP()
    {
        if (size < M)
             return 0;

        double t = ti[size/2];
        return P(t, 0);
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
