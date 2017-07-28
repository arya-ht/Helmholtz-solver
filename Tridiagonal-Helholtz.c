#include <stdio.h>
#include <stdlib.h>
#include <math.h>
/*
 *created by Arya HajiTaheri
 */
#define N 10000
#define Q 2
#define U0 0
#define UL 0
#define L 1

#define f(x)  exp(-Q*x)
void calculate();

int main(void) {
	calculate();
	return 0;
}

void calculate() {
	static float mat[N][N] = { { 0 } }, u[N], result[N], exact[N], error = 0.0;
	double h = (double)1 / N, lambda = 1;
	int i, j;
	mat[0][0] = 2 + lambda*lambda*h*h;
	mat[0][1] = -1;
	for (i = 1; i < N - 1; i++) {
		mat[i][i - 1] = -1;
		mat[i][i] = 2 + lambda*lambda*h*h;
		mat[i][i + 1] = -1;
	}
	mat[N - 1][N - 2] = -1;
	mat[N - 1][N - 1] = 2 + lambda*lambda*h*h;
	result[0] = pow(h, 2)*f(0 * h) + U0;
	for (i = 1; i < N - 1; i++) {
		result[i] = pow(h, 2)*f(i*h);
	}
	result[N - 1] = pow(h, 2)*f((N - 1)*h) + UL;


	double alpha[N], g[N], epsilon[N], top, bottom;
	int k;

	alpha[0] = mat[0][0];
	g[0] = result[0];
	for (j = 1; j < N; j++){
		alpha[j] = mat[j][j] - (mat[j][j - 1] / alpha[j - 1])*mat[j - 1][j];
		g[j] = result[j] - (mat[j][j - 1] / alpha[j - 1])*g[j - 1];
	}
	u[N - 1] = g[N - 1] / alpha[N - 1];
	for (k = N - 2; k > -1; k--){
		u[k] = (g[k] - mat[k][k + 1] * u[k + 1]) / alpha[k];
	}
	for (i = 0; i < N; i++) {
		exact[i] = U0*sinh(lambda*(L - i*h)) / sinh(lambda*L);
		exact[i] = exact[i] - ((-1 * pow(f(i*h), 2) / lambda)*((lambda*cosh(lambda*i*h) + Q*sinh(lambda*i*h) - lambda) / (lambda*lambda - Q*Q)));
		exact[i] = exact[i] + (sinh(lambda*i*h) / sinh(lambda * 1))*(UL + ((-1 * pow(f(L), 2) / (lambda*lambda - Q*Q))))*(lambda*cosh(lambda*L) + Q*sinh(lambda*L) - lambda);

	}
	for (j = 0; j < N; j++) {
		error = fabs(exact[j] - u[j]) / fabs(exact[j]);
		printf("U[%d] = %e\t\tExact[%d] = %e\t\tError: %lf\n", j + 1, u[j], j + 1, exact[j], error);
	}
}
