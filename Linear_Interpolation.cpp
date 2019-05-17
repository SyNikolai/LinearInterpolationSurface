#include <iostream>
#include <fstream>
#include <vector>
#include <math.h>
#include <cmath>
#include <string>
using namespace std;

#define FILENAME1 "didymosl162_vertices.txt"
#define FILENAME2 "didymosl162_indices.txt"
#define COLS 3

template <typename T>

// Reading data files, containing either (x,y,z) FLOAT coordinates (FILENAME1), or INT indices (FILENAME2)
void data( string filename, vector<vector<T>> &df, int *p){
	fstream file;
	vector <T> rowVector(COLS);
	int row = 0;
	file.open(filename, ios::in);
	if (file.is_open()){
		while (file.good()){
			df.push_back(rowVector);
			for (int col = 0; col < COLS; col++){
			 file >> vec[row][col];
			}
			row++;
		}
	}
	else cout << "Error reading data file: " << filename << endl;
	file.close();
	*p = row;

}

//Point distance
double Distance(int n, int m){
	if (n == m)
		return 0;
	return sqrt( (vec[m][0] - vec[n][0])*(vec[m][0] - vec[n][0]) + (vec[m][1] - vec[n][1])*(vec[m][1] - vec[n][1]) + (vec[m][2] - vec[n][2])*(vec[m][2] - vec[n][2]) );
}

//Calculate the distances of all adjacent vertices
void find_adjacent_distances(vector <double> &adjacent_distances, int *L){
	int j =0;
	for (int i = 0; i < M; i++){
	adjacent_distances.push_back(1);
	adjacent_distances[j] = Distance(ind[i][0], ind[i][1]);
	j++;
	adjacent_distances.push_back(1);
	adjacent_distances[j] = Distance(ind[i][1], ind[i][2]);
	j++;
	adjacent_distances.push_back(1);
	adjacent_distances[j] =  Distance(ind[i][2], ind[i][0]);
	j++;
	}
	*L = j - 1;
}

//Minimum point distance
double minDistance(int i, int N){
	double minDist;
	if (i < N - 1) minDist = Distance(i, i+1);
	else minDist = Distance(i, i-1);
	for (int j = 0; j < N; j++){
		if (j==i) continue;
		double dist = Distance(i, j);
		if (dist < minDist)	minDist = dist;
	}
	return minDist;
}

//Max value of minimum distances
double max_minDistance(int N){
	double dmax = 0;
	for (int i = 0; i < N; i++){
		double d = minDistance(i, N);
		if (d > dmax){
			dmax = d;
		}
	}
	return dmax;
}

//Reverse function takes as input the distance and finds the two elements
void distance_to_coordinates(double distance, int N , int *max1, int *max2){
	for(int i = 0; i < N; i++){
		for (int j = 1; j < N; j++){
			if (Distance(i,j) == distance){
				*max1 = i;
				*max2 = j;
				i = N;
				break;
			}
		}
	}
}

//Find the middle point of the max minimum distance and add that value as the new last element of the vector : vec
void linear_interpolation(int i , int j , int &N){
	vector < double > rowVector(COLS);
	vec.push_back(rowVector);
	vec[N][0] = (vec[i][0] + vec[j][0])/2;
	vec[N][1] = (vec[i][1] + vec[j][1])/2;
	vec[N][2] = (vec[i][2] + vec[j][2])/2;
	N++;
}// after real N = N - 1

//Output in original file new interpolated data
void output(int N){
	fstream file(FILENAME1);
	if (file.is_open()){
		for (int i = 0; i < N; i++){
			file  << vec[i][0] << "   " << vec[i][1] << "   " << vec[i][2] << endl;
		}
		file.close();
	}
}

//Main function
int main(int argc, char** argv) {
	int max1=0 , max2=0;
	double dmax;
	data(FILENAME1, vec, &N);
	indexdata(FILENAME2, ind, &M);
	dmax = max_minDistance(N);
    distance_to_coordinates(dmax, N, &max1, &max2);

	//Interpolate max distance
	linear_interpolation(max1,max2,N);

	//Calculate adjacent distances
	find_adjacent_distances(adjacent_distances, &L);
	int n  = N;
	//Continue to interpolate with new points any time the distance between two vertices is greater than dmax/2
	for(int i =0; i < L; i++){
		if (adjacent_distances[i] > dmax/2){
		distance_to_coordinates(adjacent_distances[i], n, &max1, &max2);
		linear_interpolation(max1,max2,N);
		}
	}

	//Add new vertices' coordinates in the original data
	output(N);
	cout  << "New max min distance is   " <<  max_minDistance(N); //error checking
}