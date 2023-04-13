#include <filesystem>
#include <string>
#include <vector>
#include <iostream>
#include <fstream>
#include "Eigen/Dense"
#include "Eigen/Eigenvalues"

extern "C" {
#include "jacobi.h"
}
#include "ReadImage.h"
#include "image.h"

namespace fs = std::filesystem;

/* getDirFiles():
 * 	Get the list of files contained within a directory.
 * args:
 * 	@dirName: The path to the directory.
 * 	@dest: The output to store the file paths.
 * return:
 * 	void
 */
void getDirFiles(std::string dirName, std::vector<std::string>& dest) {
	for(auto const& dir_entry : fs::directory_iterator(dirName)) {
		dest.push_back(dir_entry.path().string());
	}
}

/* getImageSize():
 * 	Gets the height and width of an image.
 * args:
 * 	@fName: Name of the file to compute size for.
 * 	@rows: Location to store the number of rows.
 * 	@cols: Location to store the number of cols.
 * return:
 * 	void
 */
void getImageRC(std::string fName, int& rows, int& cols) {
	ImageType img;
	int maxVal;
	readImagePGM((char*)fName.c_str(), img);
	img.getImageInfo(rows, cols, maxVal);
}

/* getImagesCovm():
 * 	Gets the covariance matrix for a set of matrices.
 *	Uses the (A^T)A trick.
 * args:
 * 	@fPaths: The paths to files to add to matrix.
 *	@dest: The location to store the covariance matrix.
 * return:
 * 	void
 */
void getImagesCovm(std::vector<std::string> fPaths, Eigen::MatrixXf& dest) {
	// Get size for images
	int rows, cols;
	getImageRC(fPaths[0], rows, cols);

	// Initialize matrix
	Eigen::MatrixXf mat(fPaths.size(), rows * cols); // images x pixel values

	// Fill in matrix
	for(int i = 0; i < fPaths.size(); i++) {
		// Get image info
		ImageType img;
		readImagePGM((char*)fPaths[i].c_str(), img);
		
		// Add image to matrix
		for(int r = 0; r < rows; r++) {
			for(int c = 0; c < cols; c++) {
				int pixVal; 
				img.getPixelVal(r, c, pixVal);
				mat(i, (r * cols) + c) = pixVal;
			}
		}
	}

	// Get covariance matrix
	dest = (float(1)/mat.cols()) * mat * mat.transpose();
}

/* saveEigen():
 * 	Saves the eigenvalues and eigenvectors of a matrix
 * 	to evalue.txt and evector.txt, respectively.
 * args:
 * 	@source: The matrix to get the eigen values and vectors from.
 * return:
 * 	void
 */
void saveEigen(Eigen::MatrixXf source) {
	// Variables
	double** covm;
	double** evectors;
	double* evalues;
	int rows = source.rows() + 1;

	// Allocate memory
	covm = new double*[rows];
	evectors = new double*[rows];
	evalues = new double[rows];
	for(int i = 0; i < rows; i++) {
		covm[i] = new double[rows];
		evectors[i] = new double[rows];
	}

	// Eigen matrix to 2d array
	for(int i = 1; i < rows; i++) {
		for(int j = 1; j < rows; j++) {
			covm[i][j] = source(i-1, j-1);
		}
	}

	// Compute eigen values and vectors
	jacobi(covm, rows - 1, evalues, evectors);

	// Print out eigen values
	for(int i = 0; i < rows; i++) {
		std::cout << evalues[i] << std::endl;
	}

	// Open file to save eigen values
	std::ofstream outF("evalues.txt");
	if(!outF.is_open()) {
		std::cout << "Could not access file evalues.txt" << std::endl;
		exit(1);
	}
	
	// Save eigen values
	for(int i = 1; i < rows; i++) {
		outF << evalues[i] << std::endl;
	}
	outF.close();

	// Open file to save eigen vectors
	outF.open("evectors.txt");
	if(!outF.is_open()) {
		std::cout << "Could not access file evectors.txt" << std::endl;
		exit(1);
	}

	// Save eigen vectors
	for(int i = 1; i < rows; i++) {
                for(int j = 1; j < rows; j++) {
			outF << evectors[i][j] << " ";
                }
		outF << std::endl;
        }
	outF.close();

	// Deallociate memory
	for(int i = 0; i < rows; i++) {
		delete[] covm[i];
		delete[] evectors[i];
	}
	delete covm;
	delete evectors;
	delete evalues;
}

int main(int argc, char** argv) {
	// Test for proper command-line call
	if(argc != 2) {
		std::cout << "Usage:\n\t $ ./main <dir>" << std::endl;
		return 1;
	}

	// Variables
	std::vector<std::string> filePaths;

	// Get file paths
	std::cout << "Getting files in directory..." << std::endl;
	getDirFiles(argv[1], filePaths);

	// Compute covariance matrix
	std::cout << "Computing sample covariance..." << std::endl;
	Eigen::MatrixXf covm;
	getImagesCovm(filePaths, covm);

	// Save covariance matrix
	std::ofstream outF("covm.txt");
	for(int i = 0; i < covm.rows(); i++) {
		for(int j = 0; j < covm.cols(); j++) {
			outF << covm(i, j) << " ";
		}
		outF << std::endl;
	}

	// Save eigen values and vectors
	std::cout << "Computing and saving eigen values and eigen vectors..." << std::endl;
	saveEigen(covm);

	return 0;
}
