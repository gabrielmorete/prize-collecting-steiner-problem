#ifndef FRAC
#define FRAC

#include <cmath>

struct frac{
	long long a, b;
	frac() : a(0), b(1) {}
	frac(int _a, int _b){ 
		long long c = __gcd(_a, _b); 
		a = _a / c;
		b = _b / c;
		simpfy();
	}

	frac simpfy(){
		long long c = __gcd(abs(a), abs(b));
		a = ((a < 0) ^ (b < 0) ? -abs(a) : abs(a));
		a = a / c;
		b = abs(b) / c;
		return *this;
	}

	frac operator+(frac p){
		long long top = a * p.b + b * p.a;
		long long bot = b * p.b;
		if (top == 0) return frac(0, 1);
		return frac(top, bot).simpfy(); 
	}
	frac operator-(frac p){
		long long top = a * p.b - b * p.a;
		long long bot = b * p.b;
		if (top == 0) return frac(0, 1);
		return frac(top, bot).simpfy(); 
	}
	frac operator*(frac p){
		long long top = a * p.a;
		long long bot = b * p.b;
		return frac(top, bot).simpfy(); 
	}
	frac operator*(long long p){
		long long top = a * p;
		long long bot = b;
		return frac(top, bot).simpfy(); 
	}
	frac operator/(long long p){
		long long top = a;
		long long bot = b * p;
		return frac(top, bot).simpfy(); 
	}

	double val(){ return ((double) a) / ((double) b);}

	inline bool operator<(frac &p){ return a * p.b < p.a * b; }

	void print(){ cout<<a<<'/'<<b<<endl; }
};

#endif