#include "conmain.h"

#define M 11

#ifndef DEBUG
    #define DEBUG 0
#endif


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
        ang += step;
        return mag*sin(ang);
    }
};


class data {
    int size;
    double ti[M];
    double vi[M];
    double *alfa;
    double freq;

    // По коэффициентам интерполяции вычисляем функцию в точке t (p == 0)
    // либо ее производную порядка p.
    double P(double t, int p)
    {
        int i, j, k;
        double res;

        if (size < M)
            return 0;

        res = 0;
        for (i = p; i < M; i++){
            k = 1;
            for (j = 0; j < p; j++)
                k *= (i-j);
            res += k*alfa[i]*pow(t, i-p);
        }

        return res;
    }

    // Сворачиваем угол в интервал (-pi; pi)
    double pi(double a)
    {
        int c;

        c = a/2/M_PI;
        a -= 2*M_PI*c;

        if (a >= 0){
            if (a <= M_PI)
                return a;
            else
                return a-2*M_PI;
        }else{
            if (a >= -M_PI)
                return a;
            else
                return 2*M_PI+a;
        }
    }

public:
    // Конструктор, внес сюда всю инициализацию
    data(double samplingFrequency)
    {
        int i;

        freq = 1/samplingFrequency;
        alfa = (double *) calloc(M, sizeof(double));
        size = 0;
        for (i = 0; i < M; i++){
            ti[i] = 0;
            vi[i] = 0;
        }
    }

    //Добавляем новый элемент в массивы данных
    void add(double v)
    {
        int i;
        double tmpti[M];

        if (size < M){
            vi[size] = v;
            ti[size] = (size+1)*freq;
            size++;
            if (size < M)
                return;
        }else{
            for (i = 0; i < M-1; i++){
                vi[i] = vi[i+1];
                ti[i] = ti[i+1];
            }
            vi[M-1] = v;
            ti[M-1] += freq;
        }

        // Ход конем, для более точной интерполяции по оси OX берем не реальный
        // промежуток времени, а точки 0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10.
        // Далее все велечины (частота, амплитуда, фаза) перещитываем исходя из
        // сделанных тут допущений.
        for (i = 0; i < M; i++)
            tmpti[i] = i;

        free(alfa);
        //Расчет коэффициентов интерполяции методом Гаусса
        if ((alfa = gauss(M, tmpti, vi)) == NULL){
            fprintf(stderr, "[ERROR GAUSS] Ошибка в определении коэфициентов\n");
            exit(ERR_GAUSS);
        }

        if (DEBUG){
            for (i = 0; i < M; i++)
                fprintf((DEBUG == 1) ? stdout : stderr, "%s%11f%s", 
                        (i) ? ", " : "\nt >>> [", tmpti[i], (i == M-1) ? "]\n" : "");
            for (i = 0; i < M; i++)
                fprintf((DEBUG == 1) ? stdout : stderr, "%s%11f%s", 
                        (i) ? ", " : "v >>> [", vi[i], (i == M-1) ? "]\n" : "");
            for (i = 0; i < M; i++)
                fprintf((DEBUG == 1) ? stdout : stderr,
                        (fabs(alfa[i]) < 0.0001) ? "%s%.6e%s" : "%s%12f%s", 
                        (i) ? " " : "a >>> ", alfa[i], (i == M-1) ? "\ny = " : "");
            for (i = 0; i < M; i++)
                fprintf((DEBUG == 1) ? stdout : stderr,
                        (fabs(alfa[i]) < 0.0001) ? "%s%.12e*x**%d%s" : "%s%.12f*x**%d%s", 
                        (alfa[i] > 0) ? "+ " : "- ", fabs(alfa[i]), i, (i == M-1) ? "\n" : " ");
        }
    }

    //Вычисляем амплитуду
    double computedMag()
    {
        if (size < M)
            return 0;

        double t = 5;
        double f = freq*computedFreq(); // Не забываем про нашего коня
        return sqrt((pow(2*M_PI*f*P(t, 3), 2) + pow(P(t, 4), 2))/pow(2*M_PI*f, 8));
    }

    //Вычисляем частоту
    double computedFreq()
    {
        if (size < M)
            return 0;

        double t = 5;
        return 1/2.0/M_PI*sqrt(-P(t, 2)/P(t, 0))/freq; // Не забываем про нашего коня
    }

    //Вычисляем фазу
    double computedAng()
    {
        if (size < M)
            return 0;

        double a;
        double t = 5;
        double f = freq*computedFreq(); // Не забываем про нашего коня
        double v = computedMag();

        a = 2*M_PI*f*ti[5]/freq;
        // Если первая производная < 0, то улол лежит в диапазоне
        // [-pi, -pi/2) U (pi/2, pi] во второй и третьей четверти
        // (тк P1 = cos). То сдвигаем угол на pi в первый и третий
        // квадрант.
        if (P(t, 1) < 0)
            a += M_PI;

        return pi(atan(-P(t, 4)/P(t, 3)/2/M_PI/f) - a);
    }
};


int main(int argc, char **argv) {
    int i;
    const double samplingFrequency = 4000.0;
    const double modelMag = 380;
    const double modelAng = 0.25 * M_PI;
    const double modelFreq = 50.123456789;
 
    model mdl(modelMag, modelAng, modelFreq, samplingFrequency);
    data dat(samplingFrequency);

    for (i = 0; i < 1000; i++) {
        double samp = mdl.next();
        dat.add(samp);
        double computedMag = dat.computedMag();
        double computedAng = dat.computedAng();
        double computedFreq = dat.computedFreq();
        printf("samp %.9f mag %.9f ang %.9f freq %.9f\n",
                samp, computedMag, computedAng, computedFreq);
    }

    return 0;
}


// Решение СЛАУ методом Гаусса, определитель 
// Вандермонда составляем на основании вектора t
double * gauss(int m, double *t, double *f)
{
    int size, i, j, d, J;
    double max, ma[m][m+1], rotate[m+1], *X;

    // Выделили память под ответ
    X = (double *) calloc(m, sizeof(double));

    // Формируем матрицу из определителя Вандермонда, и узлов интерполяции
    for (i = 0; i < m; i++){
        for (j = 0; j < m; j++){
            ma[i][j] = pow(t[i], j);
        }
        ma[i][m] = f[i];
    }

    // Прямой ход метода, переводим матрицу к диагональному виду
    for (i = 0; i < m; i++){
        max = fabs(ma[i][i]);
        J = i;
        for (j = i+1; j < m; j++){
            if (fabs(ma[j][i]) > max){
                max = fabs(ma[j][i]);
                J = j;
            }
        }

        if (J != i){
            memcpy(rotate, &ma[i][i], sizeof(double)*(m+1-i));
            memcpy(&ma[i][i], &ma[J][i], sizeof(double)*(m+1-i));
            memcpy(&ma[J][i], rotate, sizeof(double)*(m+1-i));
        }

        if (ma[i][i] != 0){
            for (j = m; j >= i; j--)
                ma[i][j] /= ma[i][i];
        }else{
            fprintf(stderr, "Система не совместна\nma[%d]=%.2f ma[%d][%d] = %.2f\n",
                    i, ma[i][i], i, j, ma[i][j]);
            return NULL;
        }

        for (j = i+1; j < m; j++){
            for (d = m; d >= i; d--){
                ma[j][d] -= ma[j][i]*ma[i][d];
            }
        }
    }

    // Обратный ход метода
    for ( i = 0; i < m; i++ )
        X[i] = ma[i][m];

    for ( i = m - 2; i >= 0; i-- )
        for ( j = i + 1; j < m; j++ )
            X[i] -= X[j] * ma[i][j];

    return X;
}
