#pragma once
//Matrix 3x3 computations 
#include <cmath>
#include <iostream>

template<typename T> class Matrix33;

template <class T> bool operator!= (const Matrix33<T> & p1, const Matrix33<T> & p2) {
	return (p1[0] != p2[0] || p1[1] != p2[1] || p1[2] != p2[2]);
}

template <class T> const Matrix33<T> operator* (const Matrix33<T> & p, float factor) {
	return Matrix33<T>(p[0] * factor, p[1] * factor, p[2] * factor);
}

template <class T> const Matrix33<T> operator* (float factor, const Matrix33<T> & p) {
	return Matrix33<T>(p[0] * factor, p[1] * factor, p[2] * factor);
}
/* not yet implemented
template <class T> const Matrix33<T> operator* (const Matrix33<T> & p1, const Matrix33<T> & p2) {
	return Matrix33<T>(p1[0] * p2[0], p1[1] * p2[1], p1[2] * p2[2]);
}
*/
template <class T> const Matrix33<T> operator+ (const Matrix33<T> & p1, const Matrix33<T> & p2) {
	return Matrix33<T>(p1[0] + p2[0], p1[1] + p2[1], p1[2] + p2[2]);
}

template <class T> const Matrix33<T> operator- (const Matrix33<T> & p1, const Matrix33<T> & p2) {
	return Matrix33<T>(p1[0] - p2[0], p1[1] - p2[1], p1[2] - p2[2]);
}

template <class T> const Vec3D<T> operator- (const Matrix33<T> & p) {
	return Matrix33<T>(-p[0], -p[1], -p[2]);
}

template <class T> const Matrix33<T> operator/ (const Matrix33<T> & p, float divisor) {
	return Matrix33<T>(p[0] / divisor, p[1] / divisor, p[2] / divisor);
}

template <class T> bool operator== (const Matrix33<T> & p1, const Matrix33<T> & p2) {
	return (p1[0] == p2[0] && p1[1] == p2[1] && p1[2] == p2[2]);
}

template <class T> bool operator< (const Matrix33<T> & a, const Matrix33<T> & b) {
	return (a[0] < b[0] && a[1] < b[1] && a[2] < b[2]);
}

template <class T> bool operator>= (const Matrix33<T> & a, const Matrix33<T> & b) {
	return (a[0] >= b[0] || a[1] >= b[1] || a[2] >= b[2]);
}


/**
* Matrix with dimensions 3x3, with basics operators overloaded.
*/
template <typename T>
class Matrix33 {
public:
	inline Matrix33(void) {
		p[0] = p[1] = p[2] = T();
	}
	inline Matrix33(Vec3D<T> p0, Vec3D<T> p1, Vec3D<T> p2) {
		p[0] = p0;
		p[1] = p1;
		p[2] = p2;
	};
	inline Matrix33(const Matrix33 & m) {
		init(m[0], m[1], m[2]);
	}
	// ---------
	// Operators
	// ---------
	inline Vec3D<T>& operator[] (int Index) {
		return (p[Index]);
	};
	inline const Vec3D<T>& operator[] (int Index) const {
		return (p[Index]);
	};
	inline Matrix33& operator= (const Matrix33<T> & P) {
		p[0] = P[0];
		p[1] = P[1];
		p[2] = P[2];
		return (*this);
	};
	inline Matrix33& operator+= (const Matrix33<T> & P) {
		p[0] += P[0];
		p[1] += P[1];
		p[2] += P[2];
		return (*this);
	};
	inline Matrix33& operator-= (const Matrix33<T> & P) {
		p[0] -= P[0];
		p[1] -= P[1];
		p[2] -= P[2];
		return (*this);
	};/* not sure if we should implement this, as multiplication isn't commutative
	inline Matrix33& operator*= (const Vec3D<T> & P) {
		p[0] *= P[0];
		p[1] *= P[1];
		p[2] *= P[2];
		return (*this);
	};*/
	inline Matrix33& operator*= (T s) {
		p[0] *= s;
		p[1] *= s;
		p[2] *= s;
		return (*this);
	};
	inline Matrix33& operator/= (T s) {
		p[0] /= s;
		p[1] /= s;
		p[2] /= s;
		return (*this);
	};

	//---------------------------------------------------------------

	inline Matrix33 & init(Vec3D<T> x, Vec3D<T> y, Vec3D<T> z) {
		p[0] = x;
		p[1] = y;
		p[2] = z;
		return (*this);
	};
	inline T det() const {
		// add the down-diagonals, subtract the up-diagonals
		return p[0][0] * p[1][1] * p[2][2] + p[1][0] * p[2][1] * p[0][2] + p[2][0] * p[0][1] * p[1][2]
			 - p[0][2] * p[1][1] * p[2][0] - p[1][2] * p[2][1] * p[0][0] - p[2][2] * p[0][1] * p[1][0];
	}
	// matrix.solve(vector) solves: matrix * X = vector
	inline Vec3D<T> solve(Vec3D<T> constant) {
		T det = (*this).det();
		if (det == 0)
			return NULL;
		return Vec3D<T>(
			Matrix33<T>(constant, p[1], p[2]).det() / det,
			Matrix33<T>(p[0], constant, p[2]).det() / det,
			Matrix33<T>(p[0], p[1], constant).det() / det);
	}

	T * pointer()
	{
		return p;
	}

	const T * pointer() const
	{
		return p;
	}


	Vec3D<T> p[3];
};

template <class T> inline Matrix33<T> swap(Matrix33<T> & P, Matrix33<T> & Q) {
	Matrix33<T> tmp = P;
	P = Q;
	Q = tmp;
}

template <class T> std::ostream & operator<< (std::ostream & output, const Matrix33<T> & v) {
	output << v[0] << " | " << v[1] << " | " << v[2];
	return output;
}

template <class T> std::istream & operator>> (std::istream & input, Matrix33<T> & v) {
	input >> v[0] >> v[1] >> v[2];
	return input;
}

typedef Matrix33<float> Matrix33f;
typedef Matrix33<double> Matrix33d;
typedef Matrix33<int> Matrix33i;
