#include <TRandom3.h>
#include <TF1.h>

int randomnumbers() {
    // Define the function dependence
    TF1* func = new TF1("func", "x^2", 0, 10); // Example: x^2 function from 0 to 10

    func -> SetNpx(10000);

    // Calculate the integral of the function over the specified range
    double integral = func->Integral(0, 10);

    // Create TRandom3 object for random number generation
    TRandom3 randGen;

    // Number of random points to generate
    int numPoints = 10;

    // Generate random numbers weighted by the functional dependence
    for (int i = 0; i < numPoints; ++i) {
        // Generate a random number between 0 and 1
        double randNumber = randGen.Uniform(0, 1);

        // Calculate x value corresponding to the random number
        double randomX = func->GetX(randNumber * integral);

        std::cout << "Random x: " << randomX << std::endl;
    }

    return 0;
}
