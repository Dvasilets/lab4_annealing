#include <iostream>
#include <sstream>
#include <algorithm>
#include <fstream>
#include <vector>
#include <random>

using namespace std;

typedef vector<vector<double> > Matrix;

double minus_temp(double t, double i) {
	return t * 0.1 / i;
}

int genint(int b = INT_MAX) {
	random_device rd;
	mt19937 mersenne(rd());
	uniform_int_distribution<> dis(0, b);
	return dis(mersenne);
}

double gendbl(double b = 1.) {
	random_device rd;
	mt19937 mersenne(rd());
	uniform_real_distribution<> dis(0, b);
	return dis(mersenne);
}

double trans(double dE, double T) {
	auto val = exp(-dE / T);
	if (val >= 0) return val;
	else return 1;
}

bool dotrans(double prob) {
	double val = gendbl();
	if (val <= prob) {
		return true;
	}
	else {
		return false;
	}
}

vector<int> gencand(const vector<int> & seq) {
	int i = genint(seq.size() - 1);
	int j = genint(seq.size() - 1);

	vector<int> mbseq(seq);

	if (i > j) {
		swap(i, j);
	}
	for (int lk = i, rk = j; lk < rk; ++lk, --rk) {
		swap(mbseq[lk], mbseq[rk]);
	}

	return move(mbseq);
}

double energy(const Matrix & matrix, const vector<int> & state) {
	double sum = 0.;
	int last = state.size() - 1, idx, next;
	for (int i = 0; i < last; ++i) {
		idx = state.at(i);
		next = state.at(i + 1);
		sum += matrix.at(idx).at(next);
	}
	idx = state.at(last);
	next = state.at(0);
	sum += matrix.at(idx).at(next);

	return sum;
}

double annealing(const Matrix & matrix, double temp, double end_temp) {
	vector<int> state;
	int matrixSize = matrix.size();
	for (int i = 0; i < matrixSize; ++i) {
		state.push_back(i);
	}
	shuffle(state.begin(), state.end(), mt19937());

	double now_energy = energy(matrix, state);
	double T = temp;

	vector<int> stcand;
	double mb_energy, p;

	const int n = 1e6;
	for (int i = 0; i < n; ++i) {
		stcand = gencand(state);
		mb_energy = energy(matrix, stcand);

		if (mb_energy < now_energy) {
			now_energy = mb_energy;
			state = stcand;
		}
		else {
			p = trans(mb_energy - now_energy, T);

			if (dotrans(p)) {
				now_energy = mb_energy;
				state = stcand;
			}
		}

		T = minus_temp(temp, i);
		if (T <= end_temp) {
			break;
		}
	}
	return now_energy;
}

Matrix read(const string & fileName)
{
	ifstream fin(fileName);
	if (fin.is_open() == false)
	{
		cout << "cant open the file" << endl;
		return Matrix();
	}
	Matrix matrix;
	int rows = 0, cols = 0, i = 0;
	string line;
	double num;
	while (getline(fin, line))
	{
		if (line.size() == 0)
		{
			continue;
		}
		cols = 0;
		istringstream iss(line);
		matrix.emplace_back();
		++rows;
		while (iss >> num)
		{
			++cols;
			matrix.at(i).push_back(num);
		}
		++i;
	}
	cout << rows << " " << cols << endl;;
	return move(matrix);
}

int main() {
	ios::sync_with_stdio(false);
	string num;
	cout << "Number: ";
	cin >> num;
	int ans = annealing(read("dataset" + num + ".txt"), 10, 0.0001);
	cout << ans << endl;
	system("pause");
	return 0;
}