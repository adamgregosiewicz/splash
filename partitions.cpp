#include <vector>
#include <algorithm>
#include <iostream>
#include <cmath>
//#include <gmpxx.h>
//#include <cln/integer.h>
//#include <cln/real.h>
//using namespace cln;
#include <boost/multiprecision/cpp_int.hpp>
//#include "uint256_t.h"

#include <boost/multiprecision/cpp_dec_float.hpp> 
#include <boost/math/constants/constants.hpp> 
#include <boost/math/distributions/normal.hpp>

using namespace boost::math;
using boost::multiprecision::cpp_dec_float_50;

#define nnmax 10260
#define mmmax 130
#define kkmax 270

#define mNmax 14566
#define mNmin 5502
#define kjmin 113
#define kjmax 288

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
        vector3d<bigInt> P;

    public:
        Partitions(size_t maxNumber, size_t maxPartsNumber, size_t maxPartSize) {
            P = vector3d<bigInt>(maxNumber + 1, maxPartsNumber + 1, maxPartSize + 1, 0);
        }

        bigInt numberOfPartitions(size_t number, size_t parts, size_t maxPart) {
        // Count the number of partitions of the number into parts parts with restrictions
        // on each part to be at least 1 and at most maxPart

            // already counted
            if (P(number, parts, maxPart) != 0)
                return P(number, parts, maxPart);
            
            // partition is impossible, the number is too large or too small
            if (maxPart * parts < number || number < parts) 
                return 0;

            // only one possibility: all ones or all maxPart
            if (maxPart * parts == number || number <= parts + 1)
                return (P(number, parts, maxPart) = 1);
            
            if (parts == 1)
                return (P(number, parts, maxPart) = 1);
            
            if (parts == 2) {
                if (maxPart * 2 >= number) {
                    // partition is possible
                    maxPart = std::min(maxPart, number - 1);
                    return (P(number, parts, maxPart) = number / parts - (number - 1 - maxPart));
                } else {
                    return 0;
                }
            }

            bigInt count = 0;
            size_t iterNum = number / parts;
            for (size_t i = 0; i < iterNum; ++i)
            //for (; iterNum--; number -= parts, --maxPart)
            {
                count += (P(number-1, parts-1, maxPart) = numberOfPartitions(number - 1, parts - 1, maxPart));
                number -= parts;
                --maxPart;
            }
            
            return (P(number,parts,maxPart) = count);
        }

        bigInt numberOfPartitions(size_t number, size_t parts, size_t minPart, size_t maxPart) {
            return numberOfPartitions(number - parts * (minPart - 1), parts, maxPart - minPart + 1);
        }
};

int width;
long long blockSize;
static vector3d<bigInt> memoize;
//static std::vector<std::vector<std::vector<bigInt>>> memoize;

long long blockCalc(int n, int m, int myMax){
    return myMax * blockSize + (n - m) * width + m - 2;
}

bigInt CountPartLenCap(int n, int m, int myMax) {
    
    //const int block = myMax * blockSize + (n - m) * width + m - 2;
    
    if (memoize(n,m,myMax) != 0) return memoize(n,m,myMax);    
    
    if (myMax * m < n || n < m) return 0;
    if (myMax * m == n || n <= m + 1) return (memoize(n,m,myMax) = 1);
    if (m < 2) return (memoize(n,m,myMax) = m);
    
    if (m == 2) {
        if (myMax * 2 >= n) {
            myMax = std::min(myMax, n - 1);
            return (memoize(n,m,myMax) = n / m - (n - 1 - myMax));
        } else {
            return 0;
        }
    }

    bigInt count = 0;
    int niter = n / m;
    for (; niter--; n -= m, --myMax) {
        count += (memoize(n-1,m-1,myMax) = CountPartLenCap(n - 1, m - 1, myMax));
    }
    
    return (memoize(n,m,myMax) = count);
}


int main() {
	
	int minn = mNmin;
    int maxn = mNmax;
    
    std::vector<bigFloat> normProb;
    normProb = std::vector<bigFloat>(15000, 0);
    
    bigFloat sum = 0.0;
    
    for(int i=minn; i<=maxn; i++) {
	    boost::math::normal norm(10000,2674.377);
	    normProb[i] = cdf(norm, i+0.5)-cdf(norm, i-0.5);
	    sum += normProb[i];
    }
    
    std::cout << std::setprecision(20) << normProb[6000] << std::endl;
    std::cout << std::setprecision(50) << sum << std::endl;
    
    //return 0;

    int n,nlow,m,low,up,myMax;
    //std::cin >> n >> m >> low >> up;
    
    nlow = mNmin;
    n = mNmax;
    m = 10;
    //myMax = 60;
    
    int kmax,kmin;
    kmin = kjmin;
    kmax = kjmax;
    myMax = kmax;
    
    //int mmax = 150;

    width = m;
    blockSize = m * (n - m + 1);
    //memoize = vector3d<bigInt>(nnmax,mmmax,kkmax, 0);
    memoize = vector3d<bigInt>(100,20,50, 0);
    
    std::vector<bigInt> sums;
    //sums = std::vector<bigInt>(15000, 0);
    //return 0;
    /*memoize.resize(nnmax+1);
    for(int i = 0; i <= nnmax; i++)
    {
        memoize[i].resize(mmmax+1);
        for(int j = 0; j <= mmmax; j++)
            memoize[i][j].resize(kkmax+1);
    }*/
    /*for(int i=0;i<=nnmax;i++)
        for(int j=0;j<=mmmax;j++)
            for(int k=0;k<=kkmax;k++)
                memoize[i][j][k] = 0;*/
    
    //std::cout << CountPartLenCap(n, 101, 50) << std::endl;
    //std::cout << CountPartLenCap(n, 100, 50) << std::endl;
    //return 0;

    std::cout << CountPartLenCap(80, 20, 60) << std::endl;

    Partitions P = Partitions(80, 20, 60);
    std::cout << P.numberOfPartitions(80, 20, 60) << std::endl;
    std::cout << P.numberOfPartitions(80, 15, 60) << std::endl;
    std::cout << P.numberOfPartitions(10, 3, 3, 60) << std::endl;
    return 0;
    
    /*mpz_int a=1;
    
    for(int i=1;i<10;i++) a *= i;
    std::cout << a << std::endl;*/
    //std::cout << CountPartLenCap(n, 101, 50) << std::endl;
    for(int i = nlow; i <= n; i++) {
		if(i%1000 == 0) std::cout << i << std::endl;
        for(int j = (int)ceil(i/(double)kmax); j <= i/kmin; j++){
			int nn = i-j*(kmin-1);
			int mmax = kmax-(kmin-1);
			sums[i] += CountPartLenCap(nn,j,mmax);
		}
	}
	
	int nimin = (int)ceil(mNmin/(double)kjmax);
	int nimax = mNmax/kjmin;

	std::cout << "policzone podzialy" << std::endl;
	int nnn = 6000;
	int uuu = 244;
	std::cout << CountPartLenCap(nnn-50*111, 50, 244-111) << std::endl;
	std::cout << sums[6000] << std::endl;
	
	std::cout << std::setprecision(50) << (bigFloat)CountPartLenCap(nnn-50*111,50,244-111)/(bigFloat)sums[6000] << std::endl;
	//return 0;

    std::vector<bigFloat> prob;
    prob = std::vector<bigFloat>(15000, 0);

    for(int i = nlow; i <= n; i++) {
		if(i%1000 == 0) std::cout << i << std::endl;
        for(int j = (int)ceil(i/(double)kmax); j <= i/kmin; j++){
			int nn = i-j*(kmin-1);
			int mmax = kmax-(kmin-1);
			prob[j] += (bigFloat)CountPartLenCap(nn,j,mmax)/(bigFloat)sums[i] * normProb[i] / sum;
		}
	}
	
	std::cout << nimin << " " << nimax << std::endl;
	for(int i = nimin; i <= nimax; i++)
		//std::cout << i << "," << std::setprecision(10) << prob[i] << std::endl;
		std::cout << std::setprecision(10) << prob[i] << std::endl;
	//return 0;
            
    /*std::cout << CountPartLenCap(10, 7, 10) << std::endl;
    std::cout << "czekam" << std::endl;
    std::cout << CountPartLenCap(13, 2, 8) << std::endl;
    std::cout << CountPartLenCap(20, 5, 10) << std::endl;
    std::cout << CountPartLenCap(200, 10, 50) << std::endl;*/
    std::cout << CountPartLenCap(1000, 101, 50) << std::endl;
    std::cout << CountPartLenCap(1000, 100, 50) << std::endl;
    return 0;
    while(1){
        std::cin >> n >> m >> low >> up;
        n -= m*(low-1);
        up -= (low-1);
        std::cout << CountPartLenCap(n, m, up) << std::endl;
    }
    
	//std::cout << CountPartLenCap(n, m, up) << std::endl;
	return 0;
}
