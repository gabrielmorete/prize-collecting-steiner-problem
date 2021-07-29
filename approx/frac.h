/* Author : Gabriel Morete de Azevedo
   Simple library for rational numbers
*/

#ifndef FRAC
#define FRAC

#include <cmath>
#include <ostream>
#include <iostream>

struct frac{
	long long a, b;
	frac() : a(0), b(1) {}
	frac(long long _a, long long _b){ 
		a = _a;
		b = _b;
		simpfy();
	}

	frac(long long _a){ a = _a; b = 1; }
	// frac(double val){ *this = frac((long long) (val * 1e10), (long long) 1e10);}

	frac simpfy(){
		long long c = std::__gcd(abs(a), abs(b));
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

	inline void operator=(frac p){ a = p.a; b = p.b;}
	inline bool operator==(frac &p){ return (a == p.a) and (b == p.b); }
	inline bool operator==(frac p){ return (a == p.a) and (b == p.b); }
	inline bool operator!=(frac p){ return (a != p.a) or (b != p.b); }
	// inline bool operator!=(frac p){ return (a != p.a) or (b != p.b); }
	inline bool operator<(const frac p) const { return a * p.b < p.a * b; }
	inline bool operator>(const frac p) const { return a * p.b > p.a * b; }
	inline bool operator<=(const frac p) const { return a * p.b <= p.a * b; }
	inline bool operator>=(const frac p) const { return a * p.b >= p.a * b; }

	double val(){ return ((double) a) / ((double) b);}

	void print(){ std::cout<<a<<'/'<<b<<std::endl; }
};

 
std::basic_ostream<char>& operator<<(std::basic_ostream<char> &o, const frac &f){
	o << f.a << "/" << f.b;
	return o;
}
#endif