#include <iostream>
#include <TRandom3.h>
#include <TMath.h>

int generate_coords() {
    // Initialize random number generator
    TRandom3 *rand = new TRandom3();

    for (int i = 0; i < 100; ++i){
	// Generate uniform random number
    	Double_t uniform_rand = rand->Uniform();

    // Generate random number according to cos^2(theta) distribution
    	Double_t scaled_cos_theta_rand = TMath::ACos(TMath::Sqrt(rand->Uniform()));

    // Print the results
    	//std::cout << "Uniform Random Number: " << uniform_rand << std::endl;
    	std::cout << "Cos^2(theta) Random Number: " << scaled_cos_theta_rand << std::endl;

}

    

    // Cleanup
    delete rand;

    return 0;
}