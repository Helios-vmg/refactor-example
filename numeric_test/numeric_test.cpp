#include <iostream>
#include <vector>
#include <cmath>

const double Lx = 200000;
const double Ly = 80000;
const int nx = 256;
const int ny = 256;
const int nyk = ny / 2 + 1;
const double pi = 3.1415926535897932384626433832795;
const double tau = pi * 2;

struct P{
	double x;
	double y;

	double sqnorm() const{
		return x * x + y * y;
	}
};

struct R{
	std::vector<P> k0;
	std::vector<P> k;
	std::vector<double> ksqu;
	std::vector<double> ninvksqu;

	R() = default;
	R(R &&other){
		*this = std::move(other);
	}
	R &operator=(R &&other){
		this->k = std::move(other.k);
		this->ksqu = std::move(other.ksqu);
		this->ninvksqu = std::move(other.ninvksqu);
		return *this;
	}
};

R method1(){
	R ret;
	ret.k0.resize(nx * nyk);
	ret.k.resize(nx * nyk);
	ret.ksqu.resize(nx * nyk);
	ret.ninvksqu.resize(nx * nyk);
	for (int i = 0; i < nx / 2 ; i++){
		for (int j = 0; j < nyk; j++){
			auto &p = ret.k[j + nyk * i];
			p.x = 2 * pi * i / Lx; 
		}
	}
	for (int i = nx / 2; i < nx ; i++){
		for (int j = 0; j < nyk; j++){
			auto &p = ret.k[j + nyk * i];
			int k_counter = i - nx/2;
			p.x =  2.* pi * (-i + 2.*k_counter) / Lx;
		}
	}
	for (int i = 0; i < nx; i++){
		for (int j = 0; j < nyk; j++){
			auto &p = ret.k[j + nyk * i];
			p.y = 2 * pi * j / Ly;
		}	
	}
	double dx = Lx / nx;
	double dy = Ly / ny;

	for (int i = 0; i < nx; i++){
		for (int j = 0; j < nyk; j++){
			auto &p = ret.k[j + nyk * i];
			auto &p0 = ret.k0[j + nyk * i];
			p0 = p;
			auto &KX = ret.k;
			auto &KY = ret.k;
			ret.ksqu[j + nyk * i] = ((sin(KX[j + nyk*i].x * dx/2)/(dx/2))*
				                     (sin(KX[j + nyk*i].x * dx/2)/(dx/2))+
				                     (sin(KY[j + nyk*i].y * dy/2)/(dy/2))*
				                     (sin(KY[j + nyk*i].y * dy/2)/(dy/2)));
			ret.ninvksqu[j + nyk * i] = -1 /ret.ksqu[j + nyk * i];
			KX[j+ nyk*i].x = sin(KX[j + nyk * i].x * dx) / dx;
			KY[j+ nyk*i].y = sin(KY[j + nyk * i].y * dy) / dy;
		}	
	}

	ret.ninvksqu[0] = -1;
	return ret;
}

R method2(){
	R ret;
	ret.k0.resize(nx * nyk);
	ret.k.resize(nx * nyk);
	ret.ksqu.resize(nx * nyk);
	ret.ninvksqu.resize(nx * nyk);

	const auto w = nx;
	const auto h = nyk;
	
	for (size_t i = 0; i < w / 2; i++){
		for (size_t j = 0; j < h; j++){
			auto &p = ret.k[j + nyk * i];
			p.x = tau * i / Lx;
		}
	}
	for (size_t i = w / 2; i < w; i++){
		double k_counter = (int)i - (int)w;
		for (size_t j = 0; j < h; j++){
			auto &p = ret.k[j + nyk * i];
			p.x = tau * k_counter / Lx;
		}
	}
	for (size_t i = 0; i < w; i++){
		auto k_counter = i - w;
		for (size_t j = 0; j < h; j++){
			auto &p = ret.k[j + nyk * i];
			p.y = tau * j / Ly;
		}
	}
	
	double dx = Lx / w;
	double dy = Ly / ((h - 1) * 2);
	double dx2 = dx / 2;
	double dy2 = dy / 2;
	for (size_t i = 0; i < w; i++){
		for (size_t j = 0; j < h; j++){
			auto &p = ret.k[j + nyk * i];
			auto &p0 = ret.k0[j + nyk * i];
			p0 = p;
			auto x2 = sin(p.x * dx2) / dx2;
			auto y2 = sin(p.y * dy2) / dy2;
			auto x = sin(p.x * dx) / dx;
			auto y = sin(p.y * dy) / dy;
			p = { x, y };
			auto sqnorm = P{x2, y2}.sqnorm();
			ret.ksqu[j + nyk * i] = sqnorm;
			ret.ninvksqu[j + nyk * i] = -1.0 / sqnorm;
		}
	}
	ret.ninvksqu[0] = -1;
	return ret;
}

struct MatrixDimensions{
	size_t cols;
	size_t rows;
};

template <typename T>
class Matrix{
	static_assert(std::is_pod<T>::value, "T must be a POD");

	struct Releaser{
		void operator()(T *p){
			if (p)
				free(p);
		}
	};

	std::unique_ptr<T, Releaser> matrix;
	size_t w = 0;
	size_t h = 0;

	static std::unique_ptr<T, Releaser> allocate(size_t w, size_t h){
		std::unique_ptr<T, Releaser> ret((T *)malloc(w * h * sizeof(T)));
		if (!ret)
			throw std::bad_alloc();
		return ret;
	}
public:
	Matrix() = default;
	Matrix(size_t w, size_t h, const T &init = {}): w(w), h(h){
		this->matrix = allocate(w, h);
		auto p = this->matrix.get();
		std::fill(p, p + w * h, init);
	}
	Matrix(const MatrixDimensions &d, const T &init = {}): Matrix(d.cols, d.rows, init){}
	Matrix(const Matrix &other){
		if (!other)
			return;
		this->w = other.w;
		this->h = other.h;
		this->matrix = allocate(this->w, this->h);
		memcpy(this->matrix.get(), other.matrix.get(), this->w * this->h * sizeof(T));
	}
	Matrix &operator=(const Matrix &other){
		if (&other == this)
			return *this;
		if (!other){
			this->matrix.reset();
			this->w = this->h = 0;
			return *this;
		}
		if (other.w != this->w || other.h != this->h){
			this->w = other.w;
			this->h = other.h;
			this->matrix = allocate(this->w, this->h);
		}
		memcpy(this->matrix.get(), other.matrix.get(), this->w * this->h * sizeof(T));
		return *this;
	}
	Matrix(Matrix &&other){
		*this = std::move(other);
	}
	Matrix &operator=(Matrix &&other){
		this->matrix = std::move(other.matrix);
		this->w = other.w;
		this->h = other.h;
		other.w = other.h = 0;
		return *this;
	}
	bool operator!() const{
		return !this->matrix;
	}
	T &get(size_t row, size_t col){
#ifdef CHECK_ACCESSES
		if (col >= this->w || row >= this->h)
			throw std::runtime_error("bad access");
#endif
		return this->matrix.get()[col * this->h + row];
	}
	const T &get(size_t row, size_t col) const{
#ifdef CHECK_ACCESSES
		if (col >= this->w || row >= this->h)
			throw std::runtime_error("bad access");
#endif
		return this->matrix.get()[col * this->h + row];
	}
	MatrixDimensions geom() const{
		return { w, h };
	}
	size_t cols() const{
		return this->w;
	}
	size_t rows() const{
		return this->h;
	}
	template <typename F>
	void for_each(const F &f){
		for (size_t i = 0; i < w; i++)
			for (size_t j = 0; j < h; j++)
				f(j, i, this->get(j, i));
	}
	template <typename F, typename T2>
	void for_each_with(const Matrix<T2> &other, const F &f){
		for (size_t i = 0; i < w; i++)
			for (size_t j = 0; j < h; j++)
				f(j, i, this->get(j, i), other.get(j, i));
	}
	template <typename F>
	void for_each_unordered(const F &f){
		auto p = this->matrix.get();
		for (size_t i = w * h; i--;)
			f(p[i]);
	}
	template <typename F>
	void for_each_unordered(const F &f) const{
		auto p = this->matrix.get();
		for (size_t i = w * h; i--;)
			f(p[i]);
	}
	template <typename T2>
	Matrix<T> &operator*=(const T2 &x){
		this->for_each_unordered([x](auto &p){
			p *= x;
		});
		return *this;
	}
	template <typename T2>
	Matrix<T> operator*(const T2 &x) const &{
		auto ret = *this;
		ret *= x;
		return ret;
	}
	template <typename T2>
	Matrix<T> operator*(const T2 &x) &&{
		*this *= x;
		return std::move(*this);
	}
	template <typename T2>
	void elementwise_multiplication(const Matrix<T2> &other){
		this->for_each([other](auto j, auto i, auto &p){
			p *= other.get(j, i);
		});
	}
	Matrix<T> &operator+=(const Matrix &other){
		this->for_each([&other](auto j, auto i, auto &p){
			p += other.get(j, i);
		});
		return *this;
	}
	Matrix<T> operator+(const Matrix &other) const{
		auto ret = *this;
		ret += other;
		return ret;
	}
	void nan_check() const;
};

class point2d{
public:
	double x;
	double y;

	double sqnorm() const{
		return x * x + y * y;
	}
	double norm() const{
		return sqrt(this->sqnorm());
	}
	point2d &operator*=(double x){
		this->x *= x;
		this->y *= x;
		return *this;
	}
	point2d &operator*=(const point2d &other){
		this->x *= other.x;
		this->y *= other.y;
		return *this;
	}
};

class FourierMesh2d{
public:
	Matrix<point2d> k;
	Matrix<double> ksqu;
	Matrix<double> ninvksqu;

	FourierMesh2d() = default;
	FourierMesh2d(double Lx, double Ly, size_t w, size_t h, int FourierMeshType);
	FourierMesh2d(const FourierMesh2d &) = delete;
	FourierMesh2d &operator=(const FourierMesh2d &) = delete;
	FourierMesh2d(FourierMesh2d &&other){
		*this = std::move(other);
	}
	FourierMesh2d &operator=(FourierMesh2d &&other){
		this->k = std::move(other.k);
		this->ksqu = std::move(other.ksqu);
		this->ninvksqu = std::move(other.ninvksqu);
		return *this;
	}
};

FourierMesh2d::FourierMesh2d(double Lx, double Ly, size_t w, size_t h, int FourierMeshType){
	this->k = Matrix<point2d>(w, h);

	for (size_t i = 0; i < w / 2; i++)
		for (size_t j = 0; j < h; j++)
			this->k.get(j, i).x = tau * i / Lx;
	for (size_t i = w / 2; i < w; i++){
		double k_counter = (int)i - (int)w;
		for (size_t j = 0; j < h; j++)
			this->k.get(j, i).x = tau * k_counter / Lx;
	}

	// Make ky. Because we are using special FFTs for purely real valued functions, the last dimension (y in this case) only need n/2 + 1 points due to symmetry about the imaginary axis.
	this->k.for_each([Ly](auto j, auto i, auto &p){
		p.y = tau * j / Ly;
	});

	double dx = Lx / w;
	double dy = Ly / ((h - 1) * 2);
	double dx2 = dx / 2;
	double dy2 = dy / 2;
	this->ksqu = Matrix<double>(w, h);
	this->ninvksqu = ksqu;
	if (FourierMeshType == 0){
		this->k.for_each([&](auto j, auto i, const auto &p){
			auto sqnorm = p.sqnorm();
			this->ksqu.get(j, i) = sqnorm;
			this->ninvksqu.get(j, i) = -1.0 / sqnorm;
		});
	}else if (FourierMeshType == 1){
		this->k.for_each([&](auto j, auto i, auto &p){
			//Use approximations for Kx, Ky, K^2
			auto x2 = sin(p.x * dx2) / dx2;
			auto y2 = sin(p.y * dy2) / dy2;
			auto x = sin(p.x * dx) / dx;
			auto y = sin(p.y * dy) / dy;
			p = { x, y };
			auto sqnorm = point2d{x2, y2}.sqnorm();
			this->ksqu.get(j, i) = sqnorm;
			this->ninvksqu.get(j, i) = -1.0 / sqnorm;
		});
	}

	//Account for the [0][0] point since there is a divide by 0 issue.
	this->ninvksqu.get(0, 0) = -1;
}

double div(double a, double b){
	if (!b && !a){
		return 1;
	}
	return a / b;
}

int main(){
	FourierMesh2d fourier_mesh(Lx, Ly, nx, nyk, 1);

	std::cout << fourier_mesh.k.get(0, 128).x << ", " << fourier_mesh.k.get(0, 128).y << std::endl;
	
	auto r1 = method1();
	auto r2 = method2();

	std::cout <<
		"i,j,"
		"r1,r2,r3,r4,"
		"r1.k0.x,r1.k0.y,r1.k.x,r1.k.y,r1.ksqu,r1.ninvksqu,"
		"r2.k0.x,r2.k0.y,r2.k.x,r2.k.y,r2.ksqu,r2.ninvksqu,"
		"\n";
	for (size_t i = 0; i < r1.k.size(); i++){
		auto p01 = r1.k0[i];
		auto p1 = r1.k[i];
		auto a1 = r1.ksqu[i];
		auto b1 = r1.ninvksqu[i];
		auto p02 = r2.k0[i];
		auto p2 = r2.k[i];
		auto a2 = r2.ksqu[i];
		auto b2 = r2.ninvksqu[i];
		auto x = i / nyk;
		auto y = i % nyk;
		std::cout
			<< x << ','
			<< y << ','
			<< div(p1.x, p2.x) << ','
			<< div(p1.y, p2.y) << ','
			<< div(a1, a2) << ','
			<< div(b1, b2) << ','
			<< p01.x << ','
			<< p01.y << ','
			<< p1.x << ','
			<< p1.y << ','
			<< a1 << ','
			<< b1 << ','
			<< p02.x << ','
			<< p02.y << ','
			<< p2.x << ','
			<< p2.y << ','
			<< a2 << ','
			<< b2 << std::endl;
	}
	return 0;
}
