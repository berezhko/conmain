#include "conmain.h"


double * gauss(int M, double *t, double *f)
{
    int size, i, j, v, k, d, p, J;
    double max, ma[M][M+1], rotate[M+1], *X;

    X = (double *) malloc(M*sizeof(double));

    for (i = 0; i < M; i++){
        for (j = 0; j < M; j++){
            ma[i][j] = pow(t[i], j);
        }
        ma[i][M] = f[i];
    }

    for (i = 0; i < M; i++){
        max = fabs(ma[i][i]);
        J = i;
        for (j = i+1; j < M; j++){
            if (fabs(ma[j][i]) > max){
                max = fabs(ma[j][i]);
                J = j;
            }
        }

        if (J != i){
            memcpy(rotate, &ma[i][i], sizeof(double)*(M+1-i));
            memcpy(&ma[i][i], &ma[J][i], sizeof(double)*(M+1-i));
            memcpy(&ma[J][i], rotate, sizeof(double)*(M+1-i));
        }

        if (ma[i][i] != 0){
            for (j = M; j >= i; j--)
                ma[i][j] /= ma[i][i];
        }else{
            fprintf(stderr, "Система не совместна\nma[%d]=%.2f ma[%d][%d] = %.2f\n",
                    i, ma[i][i], i, j, ma[i][j]);
            return NULL;
        }

        for (k = i+1; k < M; k++){
            for (d = M; d >= i; d--){
                ma[k][d] -= ma[k][i]*ma[i][d];
            }
        }
    }

    for ( i = 0; i < M; i++ )
        X[i] = ma[i][M];

    for ( i = M - 2; i >= 0; i-- )
        for ( j = i + 1; j < M; j++ )
            X[i] -= X[j] * ma[i][j];

    return X;
}
