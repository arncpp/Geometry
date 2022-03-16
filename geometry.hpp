#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace Geometry {
class Vector;
class IShape;
class Point;
class Segment;
class Ray;
class Line;
class Circle;
class Polygon;

// VECTOR CLASS

class Vector {
 private:
  int x_;
  int y_;

 public:
  // CONSTRUCTORS
  Vector();
  Vector(int x1, int y1);
  explicit Vector(const Point& p);
  Vector(const Point& p1, const Point& p2);

  // FUNCTIONS
  double GetLength() const;
  int GetX() const;
  int GetY() const;
  void SetX(int x);
  void SetY(int y);

  // OPERATORS
  Vector operator*(const double& val) const;
  int operator*(const Vector& vec) const;
  friend int operator^(const Vector& vec1, const Vector& vec2);
  Vector operator+(const Vector& vec) const;
  Vector operator-(const Vector& vec) const;
  Vector& operator+=(const Vector& vec);
  Vector& operator-=(const Vector& vec);
  Vector& operator*=(const double& val);
};

Vector::Vector() : x_(0), y_(0) {}
Vector::Vector(int x1, int y1) : x_(x1), y_(y1) {}

Vector Vector::operator*(const double& val) const {
  return Vector(x_ * static_cast<int>(val), y_ * static_cast<int>(val));
}
int Vector::operator*(const Vector& vec) const {
  return (x_ * vec.x_ + y_ * vec.y_);
}
Vector Vector::operator+(const Vector& vec) const {
  return Vector(x_ + vec.x_, y_ + vec.y_);
}
Vector& Vector::operator+=(const Vector& vec) {
  x_ += vec.x_;
  y_ += vec.y_;
  return *this;
}
Vector& Vector::operator*=(const double& val) {
  x_ *= static_cast<int>(val);
  y_ *= static_cast<int>(val);
  return *this;
}
Vector Vector::operator-(const Vector& vec) const {
  return Vector(x_ - vec.x_, y_ - vec.y_);
}
Vector& Vector::operator-=(const Vector& vec) {
  x_ -= vec.x_;
  y_ -= vec.y_;
  return *this;
}
inline double Vector::GetLength() const { return sqrt(x_ * x_ + y_ * y_); }
int Vector::GetX() const { return x_; }
int Vector::GetY() const { return y_; }
void Vector::SetY(int y) { y_ = y; }
void Vector::SetX(int x) { x_ = x; }

///////////////////////// ISHAPE CLASS /////////////////////////
class IShape {
 public:
  virtual ~IShape() = default;
  virtual IShape& Move(const Vector& v) = 0;
  virtual bool ContainsPoint(const Point& p) const = 0;
  virtual bool CrossesSegment(const Segment& s) const = 0;
  virtual IShape* Clone() const = 0;
  virtual std::string ToString() = 0;
};
///////////////////////// POINT CLASS /////////////////////////
class Point : public IShape {
 private:
  int x_;
  int y_;

 public:
  // constructors
  Point();
  Point(int x, int y);
  Point(Point const& p);

  // getters
  int GetX() const;
  int GetY() const;

  // setters
  void SetX(int x);
  void SetY(int y);

  // operators
  bool operator==(Point const& p) const;
  Vector operator-(const Point& p1) const;

  // functions
  double GetDistance(Point const& p) const;
  virtual Point& Move(const Vector& vec) override;
  virtual bool ContainsPoint(const Point& p) const override;
  virtual bool CrossesSegment(const Segment& s) const override;
  virtual Point* Clone() const override;
  virtual std::string ToString() override;
};
///////////////////////// LINE CLASS /////////////////////////
class Line : public IShape {
 private:
  int a_;
  int b_;
  int c_;

 public:
  Line();
  Line(int a, int b, int c);
  Line(const Line& l);
  Line(const Point& p1, const Point& p2);

  Vector BuildVectors() const;
  bool Parallelism(const Line& l) const;
  bool MatchLines(const Line& l) const;
  double GetDistance(const Line& l) const;
  Point Intersection(const Line& l) const;

  double GetA() const;
  double GetB() const;
  double GetC() const;

  virtual Line& Move(const Vector& vec) override;
  virtual bool ContainsPoint(const Point& p) const override;
  virtual bool CrossesSegment(const Segment& s) const override;
  virtual Line* Clone() const override;
  virtual std::string ToString() override;
};
///////////////////////// SEGMENT CLASS /////////////////////////
class Segment : public IShape {
 private:
  Point p1_;
  Point p2_;

 public:
  Segment();
  Segment(int x1, int y1, int x2, int y2);
  Segment(const Segment& s);
  Segment(const Point& pl, const Point& pr);

  Point GetP1() const;
  Point GetP2() const;

  virtual Segment& Move(const Vector& vec) override;
  virtual bool ContainsPoint(const Point& p) const override;
  virtual bool CrossesSegment(const Segment& s) const override;
  virtual Segment* Clone() const override;
  virtual std::string ToString() override;
};
///////////////////////// RAY CLASS /////////////////////////
class Ray : public IShape {
 private:
  Point p_;
  Vector vec_;

 public:
  Ray();
  Ray(const Point& p, const Vector& vec);
  Ray(const Point& p1, const Point& p2);
  Ray(int p_x, int p_y, int vec_x, int vec_y);

  virtual Ray& Move(const Vector& vec) override;
  virtual bool ContainsPoint(const Point& p) const override;
  virtual bool CrossesSegment(const Segment& s) const override;
  virtual Ray* Clone() const override;
  virtual std::string ToString() override;
};
///////////////////////// POLYGON CLASS /////////////////////////
class Polygon : public IShape {
 private:
  std::vector<Point> data_;

 public:
  Polygon(std::vector<Point> points);

  virtual Polygon& Move(const Vector& vec) override;
  virtual bool ContainsPoint(const Point& p) const override;
  virtual bool CrossesSegment(const Segment& s) const override;
  virtual Polygon* Clone() const override;
  virtual std::string ToString() override;
};
///////////////////////// CIRCLE CLASS /////////////////////////
class Circle : public IShape {
 private:
  Point mid_;
  uint64_t radius_;

 public:
  Circle();
  Circle(const Point& p, uint64_t r);
  Circle(const Circle& c);

  Point GetMid() const;
  uint64_t GetRadius() const;

  virtual Circle& Move(const Vector& vec) override;
  virtual bool ContainsPoint(const Point& p) const override;
  virtual bool CrossesSegment(const Segment& s) const override;
  virtual Circle* Clone() const override;
  virtual std::string ToString() override;
};

///////////////////////// CONSTRUCTORS /////////////////////////
// vector
Vector::Vector(const Point& p) : x_(p.GetX()), y_(p.GetY()) {}
Vector::Vector(const Point& p1, const Point& p2)
    : x_(p2.GetX() - p1.GetX()), y_(p2.GetY() - p1.GetY()) {}
int operator^(const Vector& vec1, const Vector& vec2) {
  return (vec1.GetX() * vec2.GetY() - vec1.GetY() * vec2.GetX());
}
// point
Point::Point() : x_(0), y_(0) {}
Point::Point(int x, int y) : x_(x), y_(y) {}
Point::Point(const Point& p) : x_(p.x_), y_(p.y_) {}

// line
Line::Line() : a_(0), b_(0), c_(0) {}
Line::Line(int a, int b, int c) : a_(a), b_(b), c_(c) {}
Line::Line(const Line& l) : a_(l.GetA()), b_(l.GetB()), c_(l.GetC()) {}
Line::Line(const Point& p1, const Point& p2)
    : a_(p2.GetY() - p1.GetY()),
      b_(p1.GetX() - p2.GetX()),
      c_(p2.GetX() * p1.GetY() - p1.GetX() * p2.GetY()) {}

// segment
Segment::Segment() {
  p1_.SetX(0);
  p1_.SetY(0);
  p2_.SetX(0);
  p2_.SetY(0);
}
Segment::Segment(int x1, int y1, int x2, int y2) {
  p1_.SetX(x1);
  p1_.SetY(y1);
  p2_.SetX(x2);
  p2_.SetY(y2);
}
Segment::Segment(const Segment& s) : p1_(s.GetP1()), p2_(s.GetP2()) {}
Segment::Segment(const Point& pl, const Point& pr) {
  p1_.SetX(pl.GetX());
  p1_.SetY(pl.GetY());
  p2_.SetX(pr.GetX());
  p2_.SetY(pr.GetY());
}

// ray
Ray::Ray() {
  p_.SetX(0);
  p_.SetY(0);
  vec_.SetX(0);
  vec_.SetY(0);
}
Ray::Ray(const Point& p, const Vector& vec) {
  p_ = p;
  vec_.SetX(vec.GetX());
  vec_.SetY(vec.GetY());
}
Ray::Ray(int p_x, int p_y, int vec_x, int vec_y) {
  p_.SetX(p_x);
  p_.SetY(p_y);
  vec_.SetX(vec_x);
  vec_.SetY(vec_y);
}
Ray::Ray(const Point& p1, const Point& p2) {
  p_ = p1;
  vec_.SetX(p2.GetX() - p1.GetX());
  vec_.SetY(p2.GetY() - p1.GetY());
}

// polygon
Polygon::Polygon(std::vector<Point> points) { data_ = std::move(points); }

// circle
Circle::Circle() {
  mid_.SetX(0);
  mid_.SetY(0);
  radius_ = 0;
}
Circle::Circle(const Point& p, uint64_t r) : mid_(p), radius_(r) {}
Circle::Circle(const Circle& c) : mid_(c.GetMid()), radius_(c.GetRadius()) {}

///////////////////////// GETTERS /////////////////////////
// point
int Point::GetX() const { return x_; }
int Point::GetY() const { return y_; }

// line
double Line::GetA() const { return a_; }
double Line::GetB() const { return b_; }
double Line::GetC() const { return c_; }

// segment
Point Segment::GetP1() const { return p1_; }
Point Segment::GetP2() const { return p2_; }

// circle
uint64_t Circle::GetRadius() const { return radius_; }
Point Circle::GetMid() const { return mid_; }

///////////////////////// SETTERS /////////////////////////
// point
void Point::SetX(int x) { x_ = x; }
void Point::SetY(int y) { y_ = y; }

///////////////////////// OPERATORS /////////////////////////
// point
bool Point::operator==(const Point& p) const {
  return (p.x_ == x_ && p.y_ == y_);
}
Vector Point::operator-(const Point& p1) const {
  return Vector(x_ - p1.GetX(), y_ - p1.GetY());
}

///////////////////////// FUNCTIONS /////////////////////////
/* other */
double Distance(const Segment& s, const Point& p) {
  Point p1 = s.GetP1();
  Point p2 = s.GetP2();
  Vector v1(p1, p2);
  Vector v2(p1, p);
  Vector v3(p2, p);
  Vector v4(p2, p1);
  double ds = 0;
  double dist1 = sqrt((p1.GetX() - p.GetX()) * (p1.GetX() - p.GetX()) +
                      (p1.GetY() - p.GetY()) * (p1.GetY() - p.GetY()));
  double dist2 = sqrt((p2.GetX() - p.GetX()) * (p2.GetX() - p.GetX()) +
                      (p2.GetY() - p.GetY()) * (p2.GetY() - p.GetY()));
  if (v2 * v1 <= 0 || v3 * v4 <= 0) {
    ds = std::min(dist1, dist2);
  } else {
    ds = ((p1.GetY() - p2.GetY()) * p.GetX() +
          (p1.GetX() - p2.GetX()) * p.GetY() +
          (p1.GetX() * p2.GetY() - p2.GetX() * p1.GetY())) /
         v1.GetLength();
  }
  return ds;
}
inline int Area(const Point& a, const Point& b, const Point& c) {
  return (b.GetX() - a.GetX()) * (c.GetY() - a.GetY()) -
         (b.GetY() - a.GetY()) * (c.GetX() - a.GetX());
}
bool Intersect1(int a, int b, int c, int d) {
  if (a > b) {
    std::swap(a, b);
  }
  if (c > d) {
    std::swap(c, d);
  }
  return std::max(a, c) <= std::min(b, d);
}
bool PointInCircle(const Point& p, const Circle& c) {
  bool res = false;
  int x = p.GetX();
  int y = p.GetY();
  Point mid(c.GetMid());
  int radius = c.GetRadius();
  double equation = (x - mid.GetX()) * (x - mid.GetX()) +
                    ((y - mid.GetY())) * (y - mid.GetY());
  if (equation < radius * radius) {
    res = true;
  }
  return res;
}
inline int Det(int a, int b, int c, int d) { return a * d - b * c; }
inline bool Between(int a, int b, double c) {
  return std::min(a, b) <= c && c <= std::max(a, b);
}
double Point::GetDistance(const Point& p) const {
  double x = p.x_ - x_;
  double y = p.y_ - y_;
  return sqrt(x * x + y * y);
}
Vector Line::BuildVectors() const { return Vector(b_, -a_); }
bool Line::Parallelism(const Line& l) const {
  return (a_ * l.GetB() - l.GetA() * b_ == 0);
}
double Line::GetDistance(const Line& l) const {
  double m = 0;
  double a = 0;
  double b = 0;
  double c1 = 0;
  double c2 = 0;

  if (a_ > l.GetA()) {
    m = a_ / l.GetA();
    a = l.GetA();
    b = l.GetB();
    c1 = c_ / m;
    c2 = l.GetC();
  } else {
    m = l.GetA() / a_;
    a = a_;
    b = b_;
    c2 = l.GetC() / m;
    c1 = c_;
  }
  double diff = c1 - c2;
  diff = (diff < 0) ? -diff : diff;
  return (diff / sqrt(a * a + b * b));
}
Point Line::Intersection(const Line& l) const {
  double x = 0;
  double y = 0;
  if (!Parallelism(l)) {
    x = (b_ * l.GetC() - l.GetB() * c_) / (a_ * l.GetB() - b_ * l.GetA());
    y = (c_ * l.GetA() - l.GetC() * a_) / (l.GetB() * a_ - b_ * l.GetA());
    return Point(x, y);
  }
  return Point();
}
bool Line::MatchLines(const Line& l) const {
  bool res = false;
  if (a_ * l.GetB() - l.GetA() * b_ == 0 &&
      a_ * l.GetC() - l.GetA() * c_ == 0 &&
      c_ * l.GetB() - l.GetC() * b_ == 0) {
    return !res;
  }
  return res;
}

/* move */
Point& Point::Move(const Vector& vec) {
  x_ += vec.GetX();
  y_ += vec.GetY();
  return *this;
}
Line& Line::Move(const Vector& vec) {
  int x0 = vec.GetX();
  int y0 = vec.GetY();
  int new_c = (c_ - a_ * x0 - b_ * y0);
  c_ = new_c;
  return *this;
}
Segment& Segment::Move(const Vector& vec) {
  p1_.Move(vec);
  p2_.Move(vec);
  return *this;
}
Ray& Ray::Move(const Vector& vec) {
  p_.Move(vec);
  return *this;
}
Polygon& Polygon::Move(const Vector& vec) {
  size_t size = data_.size();
  for (size_t i = 0; i < size; ++i) {
    int x = data_[i].GetX();
    data_[i].SetX(x + vec.GetX());
    int y = data_[i].GetY();
    data_[i].SetY(y + vec.GetY());
  }
  return *this;
}
Circle& Circle::Move(const Vector& vec) {
  mid_.Move(vec);
  return *this;
}

/* contains point */
bool Point::ContainsPoint(const Point& p) const {
  return (x_ == p.GetX() && y_ == p.GetY());
}
bool Line::ContainsPoint(const Point& p) const {
  return (a_ * p.GetX() + b_ * p.GetY() + c_ == 0);
}
bool Segment::ContainsPoint(const Point& p) const {
  Point p1 = p1_;
  Point p2 = p2_;
  if (p1 == p2) {
    return p == p1;
  }
  return (p.CrossesSegment(*this));
}
bool Ray::ContainsPoint(const Point& p) const {
  return ((p.GetX() - p_.GetX()) * (vec_.GetY()) ==
          (vec_.GetX()) * (p.GetY() - p_.GetY())) &&
         ((p.GetX() - p_.GetX()) * (vec_.GetX()) >= 0) &&
         ((p.GetY() - p_.GetY()) * (vec_.GetY()) >= 0);
}
bool Polygon::ContainsPoint(const Point& p) const {
  int size = data_.size();
  bool res = false;
  int i = 0;
  int j = 0;
  for (i = 0, j = size - 1; i < size; j = i++) {
    if (((data_[i].GetY() > p.GetY()) != (data_[j].GetY() > p.GetY())) &&
        (p.GetX() < (data_[j].GetX() - data_[i].GetX()) *
                            (p.GetY() - data_[i].GetY()) /
                            (data_[j].GetY() - data_[i].GetY()) +
                        data_[i].GetX())) {
      res = !res;
    }
    if (data_[i].GetX() == p.GetX() && data_[i].GetY() == p.GetY()) {
      return true;
    }
  }
  return res;
}
bool Circle::ContainsPoint(const Point& p) const {
  bool res = false;
  int x = p.GetX();
  int y = p.GetY();
  double equation = (x - mid_.GetX()) * (x - mid_.GetX()) +
                    ((y - mid_.GetY())) * (y - mid_.GetY());
  if (equation <= static_cast<double>(radius_ * radius_)) {
    res = true;
  }
  return res;
}

/* crosses segment */
bool Point::CrossesSegment(const Segment& s) const {
  Point p1 = s.GetP1();
  Point p2 = s.GetP2();
  if (p1 == p2) {
    return (p1 == *this);
  }
  Point a(p1.GetX() - x_, p1.GetY() - y_);
  Point b(p2.GetX() - x_, p2.GetY() - y_);
  return b.GetX() * a.GetY() == a.GetX() * b.GetY() &&
         (a.GetX() * b.GetX() + a.GetY() * b.GetY() <= 0);
}
bool Line::CrossesSegment(const Segment& s) const {
  Point p1 = s.GetP1();
  Point p2 = s.GetP2();
  int a2 = s.GetP1().GetY() - s.GetP2().GetY();
  int b2 = s.GetP2().GetX() - s.GetP1().GetX();
  int c2 = -a2 * s.GetP1().GetX() - b2 * s.GetP1().GetY();
  double zn = Det(a_, b_, a2, b2);
  if (zn != 0) {
    double x = -Det(c_, b_, c2, b2) * 1. / zn;
    double y = -Det(a_, c_, a2, c2) * 1. / zn;
    return Between(p1.GetX(), p2.GetX(), x) && Between(p1.GetY(), p2.GetY(), y);
  }
  return Det(a_, c_, a2, c2) == 0 && Det(b_, c_, b2, c2) == 0;
}
bool Segment::CrossesSegment(const Segment& s) const {
  int a1 = p1_.GetY() - p2_.GetY();
  int b1 = p2_.GetX() - p1_.GetX();
  int c1 = -a1 * p1_.GetX() - b1 * p1_.GetY();
  int a2 = s.GetP1().GetY() - s.GetP2().GetY();
  int b2 = s.GetP2().GetX() - s.GetP1().GetX();
  int c2 = -a2 * s.GetP1().GetX() - b2 * s.GetP1().GetY();
  double zn = Det(a1, b1, a2, b2);
  if (zn != 0) {
    double x = -Det(c1, b1, c2, b2) / zn;
    double y = -Det(a1, c1, a2, c2) / zn;
    return Between(p1_.GetX(), p2_.GetX(), x) &&
           Between(p1_.GetY(), p2_.GetY(), y) &&
           Between(s.GetP1().GetX(), s.GetP2().GetX(), x) &&
           Between(s.GetP1().GetY(), s.GetP2().GetY(), y);
  }
  return Det(a1, c1, a2, c2) == 0 && Det(b1, c1, b2, c2) == 0 &&
         Intersect1(p1_.GetX(), p2_.GetX(), s.GetP1().GetX(),
                    s.GetP2().GetX()) &&
         Intersect1(p1_.GetY(), p2_.GetY(), s.GetP1().GetY(), s.GetP2().GetY());
}

bool ProcessCrossesRay(const Vector& vec, const Point& p, const double kX,
                       const double kY) {
  if (vec.GetX() > p.GetX()) {
    if (vec.GetY() < p.GetY()) {
      if (kX < p.GetX() || kY > p.GetY()) {
        return false;
      }
    }
    if (vec.GetY() > p.GetY()) {
      if (kX < p.GetX() || kY < p.GetY()) {
        return false;
      }
    }
  } else {
    if (vec.GetY() < p.GetY()) {
      if (kX > p.GetX() || kY > p.GetY()) {
        return false;
      }
    }
    if (vec.GetY() > p.GetY()) {
      if (kX > p.GetX() || kY < p.GetY()) {
        return false;
      }
    }
  }
  return true;
}

bool Ray::CrossesSegment(const Segment& s) const {
  Point p1 = s.GetP1();
  Point p2 = s.GetP2();
  int a2 = p1.GetY() - p2.GetY();
  int b2 = p2.GetX() - p1.GetX();
  int c2 = -a2 * p1.GetX() - b2 * p1.GetY();
  Line ray(p_, Point(vec_.GetX(), vec_.GetY()));
  double zn =
      Det(static_cast<int>(ray.GetA()), static_cast<int>(ray.GetB()), a2, b2);
  if (zn == 0) {
    return Det(static_cast<int>(ray.GetA()), static_cast<int>(ray.GetC()), a2,
               c2) == 0 &&
           Det(static_cast<int>(ray.GetB()), static_cast<int>(ray.GetC()), b2,
               c2) == 0;
  }
  double x =
      -Det(static_cast<int>(ray.GetC()), static_cast<int>(ray.GetB()), c2, b2) /
      zn;
  double y =
      -Det(static_cast<int>(ray.GetA()), static_cast<int>(ray.GetC()), a2, c2) /
      zn;
  bool res = ProcessCrossesRay(vec_, p_, x, y);
  return res;
}
bool Polygon::CrossesSegment(const Segment& s) const {
  bool res = false;
  int size = data_.size();
  for (int i = 0; i < size - 1; ++i) {
    if (s.CrossesSegment(Segment(data_[i], data_[i + 1]))) {
      res = true;
    }
  }
  if (s.CrossesSegment(Segment(data_[0], data_[size - 1]))) {
    res = true;
  }
  return res;
}
bool Circle::CrossesSegment(const Segment& s) const {
  Point p1 = s.GetP1();
  Point p2 = s.GetP2();
  if (PointInCircle(p1, *this) && PointInCircle(p2, *this)) {
    return false;
  }
  return (!(Distance(s, mid_) > radius_));
}

/* clone */
Point* Point::Clone() const { return new Point(*this); }
Line* Line::Clone() const { return new Line(*this); }
Segment* Segment::Clone() const { return new Segment(*this); }
Ray* Ray::Clone() const { return new Ray(*this); }
Polygon* Polygon::Clone() const { return new Polygon(*this); }
Circle* Circle::Clone() const { return new Circle(*this); }

/* to string */
std::string Point::ToString() {
  return "Point(" + std::to_string(x_) + ", " + std::to_string(y_) + ")";
}
std::string Line::ToString() {
  return "Line(" + std::to_string(a_) + ", " + std::to_string(b_) + ", " +
         std::to_string(c_) + ")";
}
std::string Segment::ToString() {
  std::string p1 = p1_.ToString();
  std::string p2 = p2_.ToString();
  return "Segment(" + p1 + ", " + p2 + ")";
}
std::string Ray::ToString() {
  std::string p = p_.ToString();
  std::string vec = "Vector(" + std::to_string(vec_.GetX()) + ", " +
                    std::to_string(vec_.GetY()) + ")";
  return "Ray(" + p + ", " + vec + ")";
}
std::string Polygon::ToString() {
  int size = data_.size();
  std::string res = "Polygon(";
  for (int i = 0; i < size; ++i) {
    if (i != size - 1) {
      std::string temp = data_[i].ToString() + ", ";
      res += temp;
    } else {
      std::string temp = data_[i].ToString();
      res += temp;
    }
  }
  res += ")";
  return res;
}
std::string Circle::ToString() {
  std::string p = mid_.ToString();
  return "Circle(" + p + ", " + std::to_string(radius_) + ")";
}
}  // namespace Geometry

template <class SmartPtrT>
void Delete(const SmartPtrT& ptr) {}

template <class T>
void Delete(T* ptr) {
  delete ptr;
}

void CheckFunctions(const Geometry::IShape* shape_ptr,
                    const Geometry::Point& point_a,
                    const Geometry::Point& point_b) {
  std::cout << "Given shape "
            << (shape_ptr->ContainsPoint(point_a) ? "contains"
                                                  : "does not contain")
            << " point A\n";

  const auto kSegmentAb = Geometry::Segment(point_a, point_b);
  std::cout << "Given shape "
            << (shape_ptr->CrossesSegment(kSegmentAb) ? "crosses"
                                                      : "does not cross")
            << " segment AB\n";

  const auto kVectorAb = point_b - point_a;
  auto* const kClonedShapePtr =
      shape_ptr->Clone();  // may return either raw or smart pointer
  std::cout << kClonedShapePtr->Move(kVectorAb).ToString();

  Delete(kClonedShapePtr);
}
