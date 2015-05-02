#ifndef LOCAL
#define NDEBUG
#endif

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
typedef unsigned long long ull;


int getms_calls = 0;
#ifdef _MSC_VER
#include <Windows.h>
#else
#include <sys/time.h>
#endif
class Timer
{
    typedef double time_type;
    typedef unsigned int skip_type;

private:
    time_type start_time;
    time_type elapsed;

#ifdef _MSC_VER
    time_type get_ms() { return (time_type)GetTickCount64() / 1000; }
#else
    time_type get_ms() { ++getms_calls; struct timeval t; gettimeofday(&t, NULL); return (time_type)t.tv_sec * 1000 + (time_type)t.tv_usec / 1000; }
//     time_type get_ms() { ++getms_calls; return 0; }
#endif

public:
    Timer() {}

    void start() { start_time = get_ms(); }
    time_type get_elapsed() { return elapsed = get_ms() - start_time; }
};

class Random
{
private:
    unsigned int  x, y, z, w;
public:
    Random(unsigned int x
             , unsigned int y
             , unsigned int z
             , unsigned int w)
        : x(x), y(y), z(z), w(w) { }
    Random() 
        : x(123456789), y(362436069), z(521288629), w(88675123) { }
    Random(unsigned int seed)
        : x(123456789), y(362436069), z(521288629), w(seed) { }

    unsigned int next()
    {
        unsigned int t = x ^ (x << 11);
        x = y;
        y = z;
        z = w;
        return w = (w ^ (w >> 19)) ^ (t ^ (t >> 8));
    }

    int next_int() { return next(); }

    // [0, upper)
    int next_int(int upper) { return next() % upper; }

    // [low, high]
    int next_int(int low, int high) { return next_int(high - low + 1) + low; }

    double next_double(double upper) { return upper * next() / UINT_MAX; }
    double next_double(double low, double high) { return next_double(high - low) + low; }

    template <typename T>
    int select(const vector<T>& ratio)
    {
        T sum = accumulate(ratio.begin(), ratio.end(), (T)0);
        T v = next_double(sum) + (T)1e-6;
        for (int i = 0; i < (int)ratio.size(); ++i)
        {
            v -= ratio[i];
            if (v <= 0)
                return i;
        }
        return 0;
    }
};
Random g_rand;


#ifdef LOCAL
const double G_TL = 10.0 * 1000.0;
#else
const double G_TL = 9.6 * 1000.0;
#endif
Timer g_timer;



struct Point
{
    int x, y;
    Point(int x, int y)
        : x(x), y(y)
    {
    }
    Point()
        : x(0), y(0)
    {
    }

    bool operator==(const Point& other) const
    {
        return x == other.x && y == other.y;
    }
    bool operator !=(const Point& other) const
    {
        return x != other.x || y != other.y;
    }

    void operator+=(const Point& other)
    {
        x += other.x;
        y += other.y;
    }
    void operator-=(const Point& other)
    {
        x -= other.x;
        y -= other.y;
    }

    Point operator+(const Point& other) const
    {
        Point res = *this;
        res += other;
        return res;
    }
    Point operator-(const Point& other) const
    {
        Point res = *this;
        res -= other;
        return res;
    }
    Point operator*(int a) const
    {
        return Point(x * a, y * a);
    }

    bool operator<(const Point& other) const
    {
        if (x != other.x)
            return x < other.x;
        else
            return y < other.y;
    }
};
Point operator*(int a, const Point& pos)
{
    return pos * a;
}
ostream& operator<<(ostream& os, const Point& pos)
{
    os << "(" << pos.x << ", " << pos.y << ")";
    return os;
}

int norm2(const Point& p)
{
    return p.x * p.x + p.y * p.y;
}
double norm(const Point& p)
{
    return sqrt(norm2(p));
}
int dot(const Point& a, const Point& b)
{
    return a.x * b.x + a.y * b.y;
}
int cross(const Point& a, const Point& b)
{
    return a.x * b.y - a.y * b.x;
}
int dist2(const Point& a, const Point& b)
{
    return norm2(a - b);
}
double dist(const Point& a, const Point& b)
{
    return sqrt(dist2(a, b));
}

double dist(const Point& a, const Point& b, const Point& p)
{
    const Point vec = b - a;
    if (dot(vec, p - a) <= 0)
        return dist(p, a);
    if (dot(vec, p - b) >= 0)
        return dist(p, b);
    return abs(-vec.y*p.x+vec.x*p.y+a.x*b.y-a.y*b.x) / norm(vec);
}
double dist(const Point& p1, const Point& p2, const Point& q1, const Point& q2)
{
    return min(dist(p1, p2, q1), dist(p1, p2, q2));
}
bool eq(double a, double b)
{
    return abs(a - b) < 1e-9;
}

bool intersect(const Point& p1, const Point& p2, const Point& q1, const Point& q2)
{
    //do edges "this" and "other" intersect?
    if (min(p1.x,p2.x) > max(q1.x,q2.x)) return false;
    if (max(p1.x,p2.x) < min(q1.x,q2.x)) return false;
    if (min(p1.y,p2.y) > max(q1.y,q2.y)) return false;
    if (max(p1.y,p2.y) < min(q1.y,q2.y)) return false;


    const Point p_vec = p2 - p1;
    const Point q_vec = q2 - q1;
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

typedef vector<Point> Poly;

bool is_valid_poly(const Poly& poly)
{
    const int n = poly.size();

    if (n < 3)
        return false;

    rep(j, n) rep(i, j)
    {
        if (intersect(poly[i], poly[(i + 1) % n], poly[j], poly[(j + 1) % n]))
        {
            return false;
        }
    }

    return true;
}
int signed_area2(const Poly& poly)
{
    int s = 0;
    rep(i, poly.size())
        s += cross(poly[i], poly[(i + 1) % poly.size()]);
    return s;
}
int area2(const Poly& poly)
{
    return abs(signed_area2(poly));
}
double area(const Poly& poly)
{
    return area2(poly) / 2.0;
}
int area2(const vector<Poly>& polys)
{
    int sum = 0;
    for (auto& poly : polys)
        sum += area2(poly);
    return sum;
}

int triangle_area2(Point a, Point b, const Point& c)
{
    a -= c;
    b -= c;
    return abs(cross(a, b));
}

Poly to_counter_clockwise(const Poly& poly)
{
    if (signed_area2(poly) > 0)
        return poly;
    else
        return Poly(rall(poly));
}

bool intersect(const Poly& poly, const Point& p1, const Point& p2)
{
    rep(i, poly.size())
        if (intersect(poly[i], poly[(i + 1) % poly.size()], p1, p2))
            return true;
    return false;
}
bool intersect(const Poly& poly, int start_i, const Point& p1, const Point& p2)
{
    const int n = poly.size();
    vector<bool> f(n);
    for (int i = (start_i + 1) % n, j = start_i; !f[i] || !f[j]; ++i %= n, j = (j - 1 + n) % n)
    {
        if (intersect(poly[i], poly[(i + 1) % n], p1, p2) || intersect(poly[j], poly[(j + 1) % n], p1, p2))
            return true;
        f[i] = f[j] = true;
    }
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
bool intersect_polys(const vector<Poly>& polys)
{
    rep(j, polys.size()) rep(i, j)
        if (intersect(polys[i], polys[j]))
            return true;
    return false;
}

enum ContainRes
{
    OUT,
    ON,
    IN
};
ContainRes contain(const Poly& poly, const Point& p)
{
    bool in = false;
    rep(i, poly.size())
    {
        Point a = poly[i] - p;
        Point b = poly[(i + 1) % poly.size()] - p;
        if (a.y > b.y)
            swap(a, b);

        if (a.y <= 0 && 0 < b.y && cross(a, b) < 0)
            in = !in;
        if (cross(a, b) == 0 && dot(a, b) <= 0)
            return ON;
    }
    return in ? IN : OUT;
}

Poly init_nearest_pair(vector<Point>& rem)
{
    Poly poly;
    int min_dist = ten(9);
    Point a, b;
    rep(j, rem.size()) rep(i, j)
    {
        int d = dist2(rem[i], rem[j]);
        if (d < min_dist)
        {
            min_dist = d;
            a = rem[i];
            b = rem[j];
        }
    }

    poly.push_back(a);
    poly.push_back(b);
    rem.erase(find(all(rem), a));
    rem.erase(find(all(rem), b));

    return poly;
}

Poly init_rand_triangle(vector<Point>& ps)
{
    assert(ps.size() >= 3);

    rep(try_i, 10)
    {
        auto copy_ps = ps;

        int rand_i = g_rand.next_int(copy_ps.size());
        Point rand_p = copy_ps[rand_i];
        copy_ps.erase(copy_ps.begin() + rand_i);

        vector<pair<int, Point>> near;
        for (auto& p : copy_ps)
        {
            near.push_back(make_pair(dist2(rand_p, p), p));
            sort(all(near));
            if (near.size() > 2)
                near.pop_back();
        }

        copy_ps.erase(find(all(copy_ps), near[0].second));
        copy_ps.erase(find(all(copy_ps), near[1].second));

        Poly triangle = to_counter_clockwise({rand_p, near[0].second, near[1].second});
        if (triangle_area2(triangle[0], triangle[1], triangle[2]) > 0)
        {
            ps = copy_ps;
            return triangle;
        }
    }
    return {};
}

Poly build_poly(const Poly& init_poly, vector<Point> rem)
{
    Poly poly = init_poly;

    vector<Point> ps(all(init_poly));
    ps.insert(ps.end(), all(rem));
    static int ps_index[1024][1024];
    rep(i, ps.size())
        ps_index[ps[i].y][ps[i].x] = i;

    static bool is_remain[1024][1024];
    for (auto& p : rem)
        is_remain[p.y][p.x] = true;


    const auto encode = [&](const Point& a, const Point& b, const Point& p)
    {
        ull res = triangle_area2(a, b, p);
        res <<= 11;
        res |= ps_index[a.y][a.x];
        res <<= 11;
        res |= ps_index[b.y][b.x];
        res <<= 11;
        res |= ps_index[p.y][p.x];
        return res;
    };
    const auto decode = [&](ull e)
    {
        Point p = ps[e & ((1 << 11) - 1)];
        e >>= 11;
        Point b = ps[e & ((1 << 11) - 1)];
        e >>= 11;
        Point a = ps[e & ((1 << 11) - 1)];
        e >>= 11;
        int area2 = e;
        return make_tuple(area2, a, b, p);
    };

    priority_queue<ull, vector<ull>, greater<ull>> cand;
    rep(i, poly.size()) for (auto& p : rem)
        cand.push(encode(poly[i], poly[(i + 1) % poly.size()], p));

    while (!rem.empty())
    {
        static int poly_index[1024][1024];
        rep(i, poly.size())
            poly_index[poly[i].y][poly[i].x] = i;


        pair<Point, Point> best_pair;
        best_pair.first.x = -1;
        while (!cand.empty())
        {
            ull it_e = cand.top();
            cand.pop();

            auto it = decode(it_e);
            const Point& a = get<1>(it);
            const Point& b = get<2>(it);
            const Point& p = get<3>(it);

            const int poly_i = poly_index[a.y][a.x];

            if (b == poly[(poly_i + 1) % poly.size()] && is_remain[p.y][p.x]
                && !intersect(poly, poly_i, a, p) && !intersect(poly, poly_i, b, p))
            {
                best_pair = make_pair(a, p);
                break;
            }
        }
        if (best_pair.first.x == -1)
            return {};

        const int poly_i = poly_index[best_pair.first.y][best_pair.first.x];

        poly.insert(poly.begin() + poly_i + 1, best_pair.second);
        rem.erase(find(all(rem), best_pair.second));
        is_remain[best_pair.second.y][best_pair.second.x] = false;

        for (auto& p : rem)
        {
            cand.push(encode(poly[poly_i], poly[(poly_i + 1) % poly.size()], p));
            cand.push(encode(poly[(poly_i + 1) % poly.size()], poly[(poly_i + 2) % poly.size()], p));
        }
    }

    assert(is_valid_poly(poly));

    return poly;
}

vector<Poly> build_polys(const vector<Poly>& init_polys, vector<Point> rem)
{
    for (auto& poly : init_polys)
        assert(is_valid_poly(poly));
    assert(!intersect_polys(init_polys));

    vector<Poly> polys = init_polys;

    vector<Point> ps;
    for (auto& poly : polys)
        ps.insert(ps.end(), all(poly));
    ps.insert(ps.end(), all(rem));
    static int ps_index[1024][1024];
    rep(i, ps.size())
        ps_index[ps[i].y][ps[i].x] = i;

    static bool is_remain[1024][1024];
    for (auto& p : rem)
        is_remain[p.y][p.x] = true;


    const auto encode = [&](const Point& a, const Point& b, const Point& p)
    {
        ull res = triangle_area2(a, b, p);
        res <<= 11;
        res |= ps_index[a.y][a.x];
        res <<= 11;
        res |= ps_index[b.y][b.x];
        res <<= 11;
        res |= ps_index[p.y][p.x];
        return res;
    };
    const auto decode = [&](ull e)
    {
        Point p = ps[e & ((1 << 11) - 1)];
        e >>= 11;
        Point b = ps[e & ((1 << 11) - 1)];
        e >>= 11;
        Point a = ps[e & ((1 << 11) - 1)];
        e >>= 11;
        int area2 = e;
        return make_tuple(area2, a, b, p);
    };

    priority_queue<ull, vector<ull>, greater<ull>> cand;
    for (Poly& poly : polys) rep(i, poly.size()) for (auto& p : rem)
        cand.push(encode(poly[i], poly[(i + 1) % poly.size()], p));

    while (!rem.empty())
    {
        static pint poly_index[1024][1024];
        rep(polys_i, polys.size())
        {
            Poly& poly = polys[polys_i];
            rep(i, poly.size())
                poly_index[poly[i].y][poly[i].x] = pint(polys_i, i);
        }


        int min_add_area = ten(9);
        pair<Point, Point> best_pair;
        best_pair.first.x = -1;
        while (!cand.empty())
        {
            ull it_e = cand.top();
            cand.pop();

            auto it = decode(it_e);
            const Point& a = get<1>(it);
            const Point& b = get<2>(it);
            const Point& p = get<3>(it);

            int polys_i, poly_i;
            tie(polys_i, poly_i) = poly_index[a.y][a.x];

            Poly& poly = polys[polys_i];

            if (b == poly[(poly_i + 1) % poly.size()] && is_remain[p.y][p.x]
                && !intersect(poly, poly_i, a, p) && !intersect(poly, poly_i, b, p))
            {
                bool is = false;
                rep(i, polys.size())
                {
                    if (i != polys_i && (intersect(polys[i], a, p) || intersect(polys[i], b, p)))
                    {
                        is = true;
                        break;
                    }
                }

                if (!is)
                {
                    best_pair = make_pair(a, p);
                    break;
                }
            }
        }
        if (best_pair.first.x == -1)
            return {};

        const int polys_i = poly_index[best_pair.first.y][best_pair.first.x].first;
        const int poly_i = poly_index[best_pair.first.y][best_pair.first.x].second;
        Poly& poly = polys[polys_i];

        poly.insert(poly.begin() + poly_i + 1, best_pair.second);
        rem.erase(find(all(rem), best_pair.second));
        is_remain[best_pair.second.y][best_pair.second.x] = false;

        for (auto& p : rem)
        {
            cand.push(encode(poly[poly_i], poly[(poly_i + 1) % poly.size()], p));
            cand.push(encode(poly[(poly_i + 1) % poly.size()], poly[(poly_i + 2) % poly.size()], p));
        }
    }

    for (auto& poly : polys)
        assert(is_valid_poly(poly));
    assert(!intersect_polys(polys));

    return polys;
}

pair<Poly, vector<Point>> remove_points(Poly poly, int begin, int num)
{
    assert(num <= poly.size());

    rotate(poly.begin(), poly.begin() + begin, poly.end());

    vector<Point> removed_points(poly.begin(), poly.begin() + num);
    Poly remain_poly(poly.begin() + num, poly.end());
    assert(remain_poly.size() >= 3);

    rep(i, (int)remain_poly.size() - 1)
    {
        if (intersect(remain_poly[0], remain_poly.back(), remain_poly[i], remain_poly[i + 1]))
        {
            return make_pair(Poly(), vector<Point>());
        }
    }

    for (auto& p : removed_points)
    {
        if (contain(remain_poly, p) != OUT)
        {
            return make_pair(Poly(), vector<Point>());
        }
    }

    return make_pair(remain_poly, removed_points);
}

vector<Poly> divide_poly(const vector<Poly>& init_polys)
{
    // とりあえずpolys[0]をばらまくぞこの野郎
    int polys_i = 0;
    if (init_polys[polys_i].size() < 3 + 3)
        return {};

    int remove_begin = g_rand.next_int(init_polys[polys_i].size());
    int remove_num = g_rand.next_int(3, min<int>(50, (int)init_polys[polys_i].size() - 3));
    Poly div_poly;
    vector<Point> remain_points;
    tie(div_poly, remain_points) = remove_points(init_polys[polys_i], remove_begin, remove_num);
    if (div_poly.empty())
        return {};

    auto remain_polys = init_polys;
    remain_polys[polys_i] = div_poly;

    int best_area2 = ten(9);
    vector<Poly> best_polys;
    rep(try_i, 10)
    {
        vector<Poly> polys = remain_polys;
        vector<Point> rem = remain_points;
        Poly new_poly = init_rand_triangle(rem);
        if (new_poly.empty())
            continue;
        assert(new_poly.size() == 3);
        polys.push_back(new_poly);

        //TODO: 速度的にやばそう
        if (intersect_polys(polys))
            continue;

        polys = build_polys(polys, rem);

        int ar2 = area2(polys);
        if (ar2 < best_area2)
        {
            best_area2 = ar2;
            best_polys = polys;
        }
    }
    return best_polys;
}

Poly rebuild_poly(Poly poly)
{
    if (poly.size() <= 3)
        return {};

    int remove_begin = g_rand.next_int(poly.size());
    int remove_num = g_rand.next_int(1, min<int>(50, (int)poly.size() - 3));

    Poly remain_poly;
    vector<Point> remain_points;
    tie(remain_poly, remain_points) = remove_points(poly, remove_begin, remove_num);
    if (remain_poly.empty())
    {
        return {};
    }

    return build_poly(remain_poly, remain_points);
}


vector<vector<Point>> separate_points(const vector<Point>& points, int max_polys)
{
    vector<vector<Point>> separated;

    queue<tuple<int, int, int, int, vector<Point>>> q;
    q.push(make_tuple(0, 0, 700, 700, points));
    int rem_sepa = max_polys - 1;
    while (!q.empty() && rem_sepa > 0)
    {
        int x1, y1, x2, y2;
        vector<Point> ps;
        tie(x1, y1, x2, y2, ps) = q.front();
        q.pop();

        if (y2 - y1 >= x2 - x1)
        {
            int split_y = (y1 + y2) / 2;
            vector<Point> a, b;
            for (auto& p : ps)
                (p.y < split_y ? a : b).push_back(p);

            if (a.size() >= 3 && b.size() >= 3)
            {
                q.push(make_tuple(x1, y1, x2, split_y, a));
                q.push(make_tuple(x1, split_y, x2, y2, b));
                --rem_sepa;
            }
            else
                separated.push_back(ps);
        }
        else
        {
            int split_x = (x1 + x2) / 2;
            vector<Point> a, b;
            for (auto& p : ps)
                (p.x < split_x ? a : b).push_back(p);

            if (a.size() >= 3 && b.size() >= 3)
            {
                q.push(make_tuple(x1, y1, split_x, y2, a));
                q.push(make_tuple(split_x, y1, x2, y2, b));
                --rem_sepa;
            }
            else
                separated.push_back(ps);
        }
    }
    while (!q.empty())
    {
        separated.push_back(get<4>(q.front()));
        q.pop();
    }

    return separated;
}
Poly build_poly(const vector<Point>& points)
{
    Poly poly;
    while (poly.empty())
    {
        vector<Point> rem = points;
        Poly init = init_rand_triangle(rem);
        poly = build_poly(init, rem);
    }
    return poly;
}
Poly improve_poly(const Poly& init_poly, const double tl)
{
    assert(is_valid_poly(init_poly));

    int ar2 = area2(init_poly);
    Poly poly = init_poly;

    int succ = 0;
    int loops = 0;
    for (int loop_i = 0;
//          loop_i < ten(3) * 3;
            ;
         ++loop_i, ++loops)
    {
        if (loop_i % 1000 == 0 && g_timer.get_elapsed() > tl)
            break;

        Poly next = rebuild_poly(poly);
        if (!next.empty())
        {
            ++succ;
            int nar2 = area2(next);
            if (nar2 < ar2)
            {
                ar2 = nar2;
                poly = next;
            }
        }
    }
    return poly;
}
vector<Poly> improve_polys(const Poly& init_poly, const int max_polys, const double tl)
{
    assert(is_valid_poly(init_poly));
    assert(max_polys >= 2);

    vector<vector<Poly>> best_polys(max_polys);
    vector<int> best_area2(max_polys, ten(9));
    best_polys[0] = {init_poly};
    best_area2[0] = area2(init_poly);

    int loop_i = 0, poly_iter = 0;
    int updates = 0;
    for (;
        ;
        ++loop_i)
    {
        if (loop_i % 100 == 0 && g_timer.get_elapsed() > tl)
            break;

        auto polys = divide_poly(best_polys[poly_iter]);
        if (!polys.empty())
        {
            ++updates;
            int ar2 = area2(polys);
            if (ar2 < best_area2[poly_iter] && ar2 < best_area2[poly_iter + 1])
            {
                best_polys[poly_iter + 1] = polys;
                best_area2[poly_iter + 1] = ar2;
            }
        }

        ++poly_iter %= max_polys - 1;
        while (best_area2[poly_iter] == ten(9))
            ++poly_iter %= max_polys - 1;
    }
    dump(loop_i);
    dump(updates);

    rep(i, max_polys)
        fprintf(stderr, "%2d: %f\n", i, best_area2[i] / 2.0);

    for (int i = max_polys - 1; i >= 0; --i)
        if (best_area2[i] != ten(9))
            return best_polys[i];
    return {};
}
vector<Poly> solve(const vector<Point>& points, const int max_polys)
{
    Poly poly = build_poly(points);

    poly = improve_poly(poly, G_TL * (points.size() < 400 ? 0.3 : (points.size() < 1000 ? 0.5 : 0.8)));
    dump(g_timer.get_elapsed());
//     return {poly};

    vector<Poly> polys = improve_polys(poly, max_polys, G_TL);
    dump(g_timer.get_elapsed());

    dump(getms_calls);
    return polys;
}

class SmallPolygons
{
public:
    vector<string> choosePolygons(const vector<int>& _points, int max_polys)
    {
        g_timer.start();

        vector<Point> points;
        for (int i = 0; i < _points.size(); i += 2)
            points.push_back(Point(_points[i], _points[i + 1]));

        vector<Poly> polys = solve(points, max_polys);

        map<Point, int> index;
        rep(i, points.size())
            index[points[i]] = i;
        vector<string> res;
        for (auto& poly : polys)
        {
            assert(poly.size() >= 3);
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
