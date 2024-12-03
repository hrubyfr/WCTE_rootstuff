#include <vector>
#include <iostream>

int main(int argc, char* argv[]){
	if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <number of vectors>" << std::endl;
        return 1;
    }

    int nvectors = std::stoi(argv[1]);

    std::vector<double> vec1(nvectors);
    std::vector<double> vec2 = std::vector<double> (nvectors);
    std::vector<std::vector<double>> TwoDvec1(nvectors);
    std::vector<std::vector<double>> TwoDvec2 = std::vector<std::vector<double>> (nvectors);
    for (int i = 0; i < vec1.size(); i++){
    	vec1[i] = i;
    	vec2[i]  = i + i;
    	for (double j = 0; j < 10.0; j++){

    		TwoDvec1[i].push_back(j);
    		TwoDvec2[i].push_back(j);
    	}
    	

    }
    for (int i = 0; i < vec1.size(); i++){
    	std::cout << vec1[i] << " " << vec2[i] << std::endl;
    	for(int j = 0; j < TwoDvec1[i].size(); j++){

    		std::cout << " " << TwoDvec1[i][j] << std::endl;
    	}
}
	return 0;
}