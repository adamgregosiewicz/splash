#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include <boost/math/distributions/normal.hpp>

typedef boost::multiprecision::uint512_t bigInt;
typedef boost::multiprecision::cpp_dec_float_100 bigFloat;

template <typename T>
class vector3d {
    public:
        vector3d(size_t d1 = 0, size_t d2 = 0, size_t d3 = 0, T const & t = T()) :
            d1(d1), d2(d2), d3(d3), data(d1 * d2 * d3, t)
        {}

        T & operator()(size_t i, size_t j, size_t k) {
            return data[i * d2 * d3 + j * d3 + k];
        }

        T const & operator()(size_t i, size_t j, size_t k) const {
            return data[i* d2 * d3 + j * d3 + k];
        }

    private:
        size_t d1, d2, d3;
        std::vector<T> data;
};

class Partitions {
    private:
        vector3d<bigInt> partitionsMatrix;

    public:
        Partitions(size_t maxNumber, size_t numPartsMax, size_t partSizeMax) {
            partitionsMatrix = vector3d<bigInt>(maxNumber + 1, numPartsMax + 1, partSizeMax + 1, 0);
        }

        bigInt numberOfPartitions(size_t number, size_t parts, size_t partSizeMax) {
        // Count the number of partitions of the number into parts parts with restrictions
        // on each part to be at least 1 and at most partSizeMax

            // already counted
            if (partitionsMatrix(number, parts, partSizeMax) != 0)
                return partitionsMatrix(number, parts, partSizeMax);
            
            // partition is impossible, the number is too large or too small
            if (partSizeMax * parts < number || number < parts) 
                return 0;

            // only one possibility: all ones or all partSizeMax
            if (partSizeMax * parts == number || number <= parts + 1)
                return (partitionsMatrix(number, parts, partSizeMax) = 1);
            
            if (parts == 1)
                return (partitionsMatrix(number, parts, partSizeMax) = 1);
            
            if (parts == 2) {
                if (partSizeMax * 2 >= number) {
                    // partition is possible
                    partSizeMax = std::min(partSizeMax, number - 1);
                    return (partitionsMatrix(number, parts, partSizeMax) = number / parts - (number - 1 - partSizeMax));
                } else {
                    return 0;
                }
            }

            // real counting, use a formula from the paper
            bigInt count = 0;
            size_t iterNum = number / parts;
            for (size_t i = 0; i < iterNum; ++i) {
                count += (partitionsMatrix(number - 1, parts - 1, partSizeMax) = numberOfPartitions(number - 1, parts - 1, partSizeMax));
                number -= parts;
                --partSizeMax;
            }
            
            return (partitionsMatrix(number, parts, partSizeMax) = count);
        }

        // Return number of partitions with parts bounded by partSizeMin and partSizeMax
        bigInt numberOfPartitions(size_t number, size_t parts, size_t partSizeMin, size_t partSizeMax) {
            return numberOfPartitions(number - parts * (partSizeMin - 1), parts, partSizeMax - partSizeMin + 1);
        }
};

// discrete normal distribution restricted to [ceil(min),floor(max)]
std::vector<bigFloat> calculateNormalDiscrete(double mean, double std, double min, double max) {

    int minCeil = ceil(min);
    int maxFloor = floor(max);

    std::vector<bigFloat> probabilities;
    probabilities = std::vector<bigFloat>(maxFloor + 1, 0);
    
    boost::math::normal normal(mean,std);
    double normalization = cdf(normal, max) - cdf(normal, min);

    probabilities[minCeil] = cdf(normal, minCeil + 0.5) - cdf(normal, min);
    probabilities[maxFloor] = cdf(normal, max) - cdf(normal, maxFloor - 0.5);

    for(int i = minCeil + 1; i < maxFloor; ++i) {
	    probabilities[i] = (cdf(normal, i + 0.5) - cdf(normal, i - 0.5)) / normalization;
    }

    return probabilities;
}

struct Parameters {
    int eMean;
    double eStd;
    double eMin;
    double eMax;
    int eMinDiscrete;
    int eMaxDiscrete;
    int numPartsMin;
    int numPartsMax;
    int partSizeMin;
    int partSizeMax;

    Parameters(char* argv[]) {
        eMean = atoi(argv[1]);
        eStd = atof(argv[2]);
        eMin = atof(argv[3]);
        eMax = atof(argv[4]);
        partSizeMin = atoi(argv[5]);
        partSizeMax = atoi(argv[6]);
        eMinDiscrete = (int)ceil(eMin);
        eMaxDiscrete = (int)floor(eMax);
        numPartsMin = (size_t)ceil(eMin / (double)partSizeMax);
        numPartsMax = (size_t)floor(eMax / (double)partSizeMin);
    };
};


int main(int argc, char* argv[]) {

    if(argc != 7) {
        std::cout << "Execution: ./partitions energyMean energyStd energyMin energyMax partSizeMin partSizeMax" << std::endl;
        return 0;
    }

    Parameters Params(argv);

    std::vector<bigFloat> probabilitiesNormalDiscrete;
    probabilitiesNormalDiscrete = calculateNormalDiscrete(Params.eMean, Params.eStd, Params.eMin, Params.eMax);

    std::vector<bigInt> sumsOfPartitions = std::vector<bigInt>(Params.eMaxDiscrete + 1, 0);

    Partitions P(Params.eMaxDiscrete, Params.numPartsMax, Params.partSizeMax);

    // calculate partitions
    for(int energy = Params.eMinDiscrete; energy <= Params.eMaxDiscrete; ++energy) {
		// if(energy % 1000 == 0) std::cout << energy << std::endl;
        for(int parts = (int)ceil(energy / (double)Params.partSizeMax); parts <= energy / Params.partSizeMin; ++parts) {
            sumsOfPartitions[energy] += P.numberOfPartitions(energy, parts, Params.partSizeMin, Params.partSizeMax);
		}
	}

    std::vector<bigFloat> probabilitiesPartitions(Params.eMaxDiscrete + 1, 0);

    // calculate distribution
    for(int energy = Params.eMinDiscrete; energy <= Params.eMaxDiscrete; ++energy) {
		// if(energy % 1000 == 0) std::cout << energy << std::endl;
        for(int parts = (int)ceil(energy / (double)Params.partSizeMax); parts <= energy / Params.partSizeMin; ++parts) {
			probabilitiesPartitions[parts] += (bigFloat)P.numberOfPartitions(energy, parts, Params.partSizeMin, Params.partSizeMax)
                                              / (bigFloat)sumsOfPartitions[energy]
                                              * probabilitiesNormalDiscrete[energy];
		}
	}

    // print distribution
    // std::cout << Params.numPartsMin << " " << Params.numPartsMax << std::endl;
	
    std::cout << "no,prob" << std::endl;
	for(int parts = Params.numPartsMin; parts <= Params.numPartsMax; ++parts)
		std::cout << parts << "," << std::setprecision(10) << probabilitiesPartitions[parts] << std::endl;

	return 0;
}
