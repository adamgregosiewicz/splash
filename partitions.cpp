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

using namespace boost::multiprecision;

typedef uint512_t bigtype;

template <typename T>
class vector3d {
public:
    vector3d(size_t d1=0, size_t d2=0, size_t d3=0, T const & t=T()) :
        d1(d1), d2(d2), d3(d3), data(d1*d2*d3, t)
    {}

    T & operator()(size_t i, size_t j, size_t k) {
        return data[i*d2*d3 + j*d3 + k];
    }

    T const & operator()(size_t i, size_t j, size_t k) const {
        return data[i*d2*d3 + j*d3 + k];
    }

private:
    size_t d1,d2,d3;
    std::vector<T> data;
};

int width;
long long blockSize;
static vector3d<bigtype> memoize;
//static std::vector<std::vector<std::vector<bigtype>>> memoize;

long long blockCalc(int n, int m, int myMax){
    return myMax * blockSize + (n - m) * width + m - 2;
}

bigtype pStdCap(int n, int m, int myMax) {
    
    const int block = blockCalc(n,m,myMax);
    //std::cout << n << " " << m << " " << myMax << " " << block << " ";
    
    if (memoize(n,m,myMax) != 0) { //std::cout << memoize(n,m,myMax) << std::endl; 
    return memoize(n,m,myMax); }
    if (myMax * m < n || n < m) { //std::cout << 0 << std::endl; 
    return 0; }
    if (myMax * m == n || n <= m + 1) { //std::cout << 1 << std::endl; 
    return (memoize(n,m,myMax) = 1); }
    if (m < 2) { //std::cout << m << std::endl; 
    return (memoize(n,m,myMax) = m); }
    
    int niter = n / m;
    
    if (m == 2) {
        if (myMax * 2 >= n) {
            myMax = std::min(myMax, n - 1);
            //std::cout << (niter - (n - 1 - myMax)) << std::endl;
            return (memoize(n,m,myMax) = (niter - (n - 1 - myMax)));
        } else {
            //std::cout << 0 << std::endl;
            return 0;
        }
    }
    
    bigtype count = 0;
    //std::cout << "iter" << std::endl;
    for (; niter--; n -= m, --myMax) {
        count += (memoize(n-1,m-1,myMax) = pStdCap(n - 1, m - 1, myMax));
    }
    
    memoize(n,m,myMax) = count;
    //std::cout << count << std::endl;
    return count;
}

bigtype CountPartLenCap(int n, int m, int myMax) {
    
    const int block = blockCalc(n,m,myMax);
    //std::cout << n << " " << m << " " << myMax << " " << block << std::endl;
    
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
    
    memoize(n,m,myMax) = pStdCap(n, m, myMax);
    return memoize(n,m,myMax);
}


int main() {
	
	int minn = mNmin;
    int maxn = mNmax;
    
    std::vector<cpp_dec_float_100> normProb;
    normProb = std::vector<cpp_dec_float_100>(15000, 0);
    
    cpp_dec_float_100 sum = 0.0;
    
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
    memoize = vector3d<bigtype>(nnmax,mmmax,kkmax, 0);
    
    std::vector<bigtype> sums;
    sums = std::vector<bigtype>(15000, 0);
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
	
	std::cout << std::setprecision(50) << (cpp_dec_float_100)CountPartLenCap(nnn-50*111,50,244-111)/(cpp_dec_float_100)sums[6000] << std::endl;
	//return 0;

    std::vector<cpp_dec_float_100> prob;
    prob = std::vector<cpp_dec_float_100>(15000, 0);

    for(int i = nlow; i <= n; i++) {
		if(i%1000 == 0) std::cout << i << std::endl;
        for(int j = (int)ceil(i/(double)kmax); j <= i/kmin; j++){
			int nn = i-j*(kmin-1);
			int mmax = kmax-(kmin-1);
			prob[j] += (cpp_dec_float_100)CountPartLenCap(nn,j,mmax)/(cpp_dec_float_100)sums[i] * normProb[i] / sum;
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
