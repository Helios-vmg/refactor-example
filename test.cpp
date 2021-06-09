#include <iostream>
#include <fstream>
#include <vector>

typedef unsigned char B;

int main(int argc, char **argv){
	if (argc < 2)
		return -1;
	int dimensions = 1;
	if (argc > 2)
		dimensions = atoi(argv[2]);
	std::ifstream file(argv[1], std::ios::binary);
	if (!file)
		return -1;
	file.seekg(0x50);
	std::vector<double> data(0x80000 / sizeof(double));
	file.read((char *)&data[0], data.size() * sizeof(double));
	for (size_t i = 0; i < data.size();){
		for (int j = 0; j < dimensions && i < data.size(); j++, i++){
			if (j)
				std::cout << '\t';
			std::cout << data[i];
		}
		std::cout << std::endl;
	}
	return 0;
	for (size_t i = data.size(); i--;){
		if (data[i] < 0){
			std::cout << "negative at " << i << std::endl;
		}
	}
	return 0;
}