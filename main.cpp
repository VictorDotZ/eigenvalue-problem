#include <cmath>
#include <iostream>
#include <vector>

double euclideanNorm(const std::vector<double>& u)
{
	double result = 0;
	for (auto elem : u)
		result += elem * elem;

	return std::sqrt(result);
}

std::vector<double> pseudoMatrixToVectorMultiplication(const std::vector<double>& vec, const double h)
{
	size_t vecSize = vec.size();
	std::vector<double> result(vecSize);

	result[0] = (2.0 * vec[0] - 2.0 * vec[1]) / h / h;

	for (size_t i = 1; i < vecSize - 1; ++i)
		result[i] = (-vec[i - 1] + 2.0 * vec[i] - vec[i + 1]) / h / h;

	result[vecSize - 1] = (-vec[vecSize - 2] + vec[vecSize - 1]) / h / h;

	return result;
}

/* явный вид собственных чисел */
double eigenValue(const uint64_t n, const double h, const uint64_t N)
{
	double lambda = 2.0 * std::sin(M_PI * n / (2.0 * N - 1)) / h;
	return lambda * lambda;
}

/* явный вид компонент собственных векторов */
double eigenFunctionComponent(const double C, const uint64_t n, const uint64_t k, const uint64_t N)
{
	return C * std::cos(2.0 * M_PI * static_cast<double>(n) * static_cast<double>(k) / (2.0 * static_cast<double>(N) - 1.0));
}

std::vector<double> eigenFunction(const double C, const uint64_t n, const uint64_t N)
{
	std::vector<double> y_n(N);

	for (size_t k = 0; k < N; ++k)
		y_n[k] = eigenFunctionComponent(C, n, k, N);

	return y_n;
}

/* в моей задаче несимметрия только в нулевой компоненте, коэфф 0.5*/
double dotProduct(const std::vector<double>& u, const std::vector<double>& v, const double h)
{
	double sum = u[0] * v[0] * 0.5;

	for (size_t i = 1; i < u.size(); ++i)
		sum += u[i] * v[i];

	return h * sum;
}

double maxDotProductInaccuracy(const double C, const double h, const uint64_t N)
{
	double max = 0.0;

	for (size_t i = 0; i < N; ++i) {
		for (size_t j = i + 1; j < N; ++j) {
			double tmp = std::abs(dotProduct(eigenFunction(C, i, N), eigenFunction(C, j, N), h));
			max = (tmp > max) ? tmp : max;
		}
	}

	return max;
}

double eigenVectorInaccuracy(std::vector<double> vec, const double lambda, const double h)
{
	std::vector<double> matrixToVec = pseudoMatrixToVectorMultiplication(vec, h);

	double vecNorm = euclideanNorm(vec);

	for (auto& elem : vec)
		elem *= lambda;

	for (size_t i = 0; i < matrixToVec.size(); ++i)
		matrixToVec[i] -= vec[i];

	if ((std::abs(lambda) >= std::numeric_limits<double>::epsilon()))
		vecNorm *= std::abs(lambda);

	return euclideanNorm(matrixToVec) / vecNorm;
}

double maxEigenVectorInaccuracy(const double C, const double h, const uint64_t N)
{
	double max = 0.0;

	for (size_t i = 0; i < N; ++i) {
		double tmp = eigenVectorInaccuracy(eigenFunction(C, i, N), eigenValue(i, h, N), h);
		max = (tmp > max) ? tmp : max;
	}

	return max;
}

int main(int argc, char** argv)
{
	if (argc != 2)
		return -1;

	uint64_t N = std::stoull(argv[1]);

	/* значение h убирающее зависимость от N в константах C */
	double h = 1.0 / (static_cast<double>(N) - 0.5);

	/* константа для ортонормированности собственных векторов */
	double C = std::sqrt(2);

	std::cout << "max(e_i, e_j), i != j:\t\t\t\t" << std::scientific << maxDotProductInaccuracy(C, h, N) << std::endl;

	std::cout << "max ||A * e - lambda * e|| / ||lambda * e||\t" << std::scientific << maxEigenVectorInaccuracy(C, h, N)
	          << std::endl;

	return 0;
}