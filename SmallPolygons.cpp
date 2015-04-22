#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <climits>
#include <cfloat>
#include <ctime>
#include <cassert>
#include <map>
#include <utility>
#include <set>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include <algorithm>
#include <functional>
#include <sstream>
#include <complex>
#include <stack>
#include <queue>
#include <numeric>
#include <list>
#include <iomanip>
#include <fstream>
#include <bitset>

using namespace std;

#define foreach(it, c) for (__typeof__((c).begin()) it=(c).begin(); it != (c).end(); ++it)
template <typename T> void print_container(ostream& os, const T& c) { const char* _s = " "; if (!c.empty()) { __typeof__(c.begin()) last = --c.end(); foreach (it, c) { os << *it; if (it != last) os << _s; } } }
template <typename T> ostream& operator<<(ostream& os, const vector<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const set<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const multiset<T>& c) { print_container(os, c); return os; }
template <typename T> ostream& operator<<(ostream& os, const deque<T>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const map<T, U>& c) { print_container(os, c); return os; }
template <typename T, typename U> ostream& operator<<(ostream& os, const pair<T, U>& p) { os << "(" << p.first << ", " << p.second << ")"; return os; }

template <typename T> void print(T a, int n, const string& split = " ") { for (int i = 0; i < n; i++) { cout << a[i]; if (i + 1 != n) cout << split; } cout << endl; }
template <typename T> void print2d(T a, int w, int h, int width = -1, int br = 0) { for (int i = 0; i < h; ++i) { for (int j = 0; j < w; ++j) { if (width != -1) cout.width(width); cout << a[i][j] << ' '; } cout << endl; } while (br--) cout << endl; }
template <typename T> void input(T& a, int n) { for (int i = 0; i < n; ++i) cin >> a[i]; }
#define dump(v) (cerr << #v << ": " << v << endl)

#define rep(i, n) for (int i = 0; i < (int)(n); ++i)
#define erep(i, n) for (int i = 0; i <= (int)(n); ++i)
#define all(a) (a).begin(), (a).end()
#define rall(a) (a).rbegin(), (a).rend()
#define clr(a, x) memset(a, x, sizeof(a))
#define sz(a) ((int)(a).size())
#define mp(a, b) make_pair(a, b)
#define ten(n) ((long long)(1e##n))

template <typename T, typename U> void upmin(T& a, const U& b) { a = min<T>(a, b); }
template <typename T, typename U> void upmax(T& a, const U& b) { a = max<T>(a, b); }
template <typename T> void uniq(T& a) { sort(a.begin(), a.end()); a.erase(unique(a.begin(), a.end()), a.end()); }
template <class T> string to_s(const T& a) { ostringstream os; os << a; return os.str(); }
template <class T> T to_T(const string& s) { istringstream is(s); T res; is >> res; return res; }
void fast_io() { cin.tie(0); ios::sync_with_stdio(false); }
bool in_rect(int x, int y, int w, int h) { return 0 <= x && x < w && 0 <= y && y < h; }

typedef long long ll;
typedef pair<int, int> pint;


struct Pos
{
    int x, y;
    Pos(int x, int y)
        : x(x), y(y)
    {
    }
    Pos()
        : x(0), y(0)
    {
    }

    bool operator==(const Pos& other) const
    {
        return x == other.x && y == other.y;
    }
    bool operator !=(const Pos& other) const
    {
        return x != other.x || y != other.y;
    }

    void operator+=(const Pos& other)
    {
        x += other.x;
        y += other.y;
    }
    void operator-=(const Pos& other)
    {
        x -= other.x;
        y -= other.y;
    }

    Pos operator+(const Pos& other) const
    {
        Pos res = *this;
        res += other;
        return res;
    }
    Pos operator-(const Pos& other) const
    {
        Pos res = *this;
        res -= other;
        return res;
    }
    Pos operator*(int a) const
    {
        return Pos(x * a, y * a);
    }

    bool operator<(const Pos& other) const
    {
        if (x != other.x)
            return x < other.x;
        else
            return y < other.y;
    }
};
Pos operator*(int a, const Pos& pos)
{
    return pos * a;
}
ostream& operator<<(ostream& os, const Pos& pos)
{
    os << "(" << pos.x << ", " << pos.y << ")";
    return os;
}

int norm2(const Pos& p)
{
    return p.x * p.x + p.y * p.y;
}
double norm(const Pos& p)
{
    return sqrt(norm2(p));
}
int dot(const Pos& a, const Pos& b)
{
    return a.x * b.x + a.y * b.y;
}
int cross(const Pos& a, const Pos& b)
{
    return a.x * b.y - a.y * b.x;
}
int dist2(const Pos& a, const Pos& b)
{
    return norm2(a - b);
}
double dist(const Pos& a, const Pos& b)
{
    return sqrt(dist2(a, b));
}

double dist(const Pos& a, const Pos& b, const Pos& p)
{
    const Pos vec = b - a;
    if (dot(vec, p - a) <= 0)
        return dist(p, a);
    if (dot(vec, p - b) >= 0)
        return dist(p, b);
    return abs(-vec.y*p.x+vec.x*p.y+a.x*b.y-a.y*b.x) / norm(vec);
}
double dist(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2)
{
    return min(dist(p1, p2, q1), dist(p1, p2, q2));
}
bool eq(double a, double b)
{
    return abs(a - b) < 1e-9;
}

bool intersect(const Pos& p1, const Pos& p2, const Pos& q1, const Pos& q2)
{
    //do edges "this" and "other" intersect?
    if (min(p1.x,p2.x) > max(q1.x,q2.x)) return false;
    if (max(p1.x,p2.x) < min(q1.x,q2.x)) return false;
    if (min(p1.y,p2.y) > max(q1.y,q2.y)) return false;
    if (max(p1.y,p2.y) < min(q1.y,q2.y)) return false;


    const Pos p_vec = p2 - p1;
    const Pos q_vec = q2 - q1;
    const int den = q_vec.y*p_vec.x-q_vec.x*p_vec.y;
    const int num1 = q_vec.x*(p1.y-q1.y)-q_vec.y*(p1.x-q1.x);
    const int num2 = p_vec.x*(p1.y-q1.y)-p_vec.y*(p1.x-q1.x);

    //parallel edges
    if (den==0)
    {
        if (dist(q1, q2, p1, p2) > 0)
            return false;

        //on the same line - "not intersect" only if one of the vertices is common,
        //and the other doesn't belong to the line
        const double sum_len = dist(p1, p2) + dist(q1, q2);
        if ((p1==q1 && eq(dist(p2, q2) , sum_len)) || 
            (p1==q2 && eq(dist(p2, q1) , sum_len)) ||
            (p2==q1 && eq(dist(p1, q2) , sum_len)) ||
            (p2==q2 && eq(dist(p1, q1) , sum_len)))
            return false;

        return true;
    }

    //common vertices
    if (p1==q1 || p1==q2 || p2==q1 || p2==q2)
        return false;

    double u1 = num1*1./den;
    double u2 = num2*1./den;
    if (u1<0 || u1>1 || u2<0 || u2>1)
        return false;
    return true;
}

typedef vector<Pos> Poly;

bool is_valid_poly(const Poly& poly)
{
    const int n = poly.size();

    if (n < 3)
        return false;

    rep(j, n) rep(i, j)
    {
        if (intersect(poly[i], poly[(i + 1) % n], poly[j], poly[(j + 1) % n]))
            return false;
    }

    return true;
}
int area2(const Poly& poly)
{
    int s = 0;
    rep(i, poly.size())
        s += cross(poly[i], poly[(i + 1) % poly.size()]);
    return abs(s);
}
double area(const Poly& poly)
{
    return area2(poly) / 2.0;
}

bool intersect(const Poly& poly, const Pos& p1, const Pos& p2)
{
    rep(i, poly.size())
        if (intersect(poly[i], poly[(i + 1) % poly.size()], p1, p2))
            return true;
    return false;
}
bool intersect(const Poly& a, const Poly& b)
{
    rep(i, a.size()) rep(j, b.size())
    {
        if (intersect(a[i], a[(i + 1) % a.size()], b[j], b[(j + 1) % b.size()]))
            return true;
    }
    return false;
}


Poly build_poly(const vector<Pos>& points)
{
    const int n = points.size();

    vector<Pos> rem = points;
    Poly poly;
    {
        int min_dist = ten(9);
        Pos a, b;
        rep(j, n) rep(i, j)
        {
            int d = dist2(points[i], points[j]);
            if (d < min_dist)
            {
                min_dist = d;
                a = points[i];
                b = points[j];
            }
        }

        poly.push_back(a);
        poly.push_back(b);
        rem.erase(find(all(rem), a));
        rem.erase(find(all(rem), b));
    }

    while (!rem.empty())
    {
//         dump(poly.size());

        vector<tuple<int, int, int>> area_i_j;
        rep(i, poly.size())
        {
            rep(j, rem.size())
            {
                int add_area = area2({poly[i], poly[(i + 1) % poly.size()], rem[j]});
                area_i_j.push_back(make_tuple(add_area, i, j));
            }
        }
        sort(all(area_i_j));

        int best_i = -1, best_j;
        for (auto& it : area_i_j)
        {
            int add_area, i, j;
            tie(add_area, i, j) = it;
            if (!intersect(poly, poly[i], rem[j]) && !intersect(poly, poly[(i + 1) % poly.size()], rem[j]))
            {
                best_i = i;
                best_j = j;
                break;
            }
        }
        if (best_i == -1)
            return poly;
        assert(best_i != -1);

        poly.insert(poly.begin() + best_i + 1, rem[best_j]);
        rem.erase(rem.begin() + best_j);
//         if (!is_valid_poly(poly))
//             return poly;
        assert(is_valid_poly(poly));
    }

    return poly;
}
vector<Poly> solve(const vector<Pos>& points, const int max_polys)
{
    return {build_poly(points)};
}

class SmallPolygons
{
public:
    vector<string> choosePolygons(vector<int> _points, int max_polys)
    {
        vector<Pos> points;
        for (int i = 0; i < _points.size(); i += 2)
            points.push_back(Pos(_points[i], _points[i + 1]));

        vector<Poly> polys = solve(points, max_polys);

        map<Pos, int> index;
        rep(i, points.size())
            index[points[i]] = i;
        vector<string> res;
        for (auto& poly : polys)
        {
            string s;
            s += to_s(index[poly[0]]);
            for (int i = 1; i < poly.size(); ++i)
                s += " " + to_s(index[poly[i]]);
            res.push_back(s);
        }
        return res;
    }
};

#ifdef LOCAL
int main(int argc, char** argv)
{
    int n;
    cin >> n;
    vector<int> points(n);
    input(points, n);
    int max_polys;
    cin >> max_polys;

    vector<string> res = SmallPolygons().choosePolygons(points, max_polys);
    cout << res.size() << endl;
    for (auto& s : res)
        cout << s << endl;
    cout.flush();
}
#endif
