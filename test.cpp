#include <iostream>
#include <fstream>
#include <vector>

typedef unsigned char B;

int main(int argc, char **argv){
	if (argc < 2)
		return -1;
	std::ifstream file(argv[1], std::ios::binary);
	if (!file)
		return -1;
	std::vector<double> data(0x80000 / sizeof(double));
	file.read((char *)&data[0], data.size() * sizeof(double));
	for (size_t i = 0; i < data.size(); i++)
		std::cout << data[i] << std::endl;
	return 0;
	for (size_t i = data.size(); i--;){
		if (data[i] < 0){
			std::cout << "negative at " << i << std::endl;
		}
	}
	return 0;
}