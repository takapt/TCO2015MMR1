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
    time_type get_ms() { struct timeval t; gettimeofday(&t, NULL); return (time_type)t.tv_sec * 1000 + (time_type)t.tv_usec / 1000; }
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
const double G_TLE = 6.0 * 1000;
#else
const double G_TLE = 9.6 * 1000;
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
            return false;
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

    int rand_i = g_rand.next_int(ps.size());
    Point rand_p = ps[rand_i];
    ps.erase(ps.begin() + rand_i);

    vector<pair<int, Point>> near;
    for (auto& p : ps)
    {
        near.push_back(make_pair(dist2(rand_p, p), p));
        sort(all(near));
        if (near.size() > 2)
            near.pop_back();
    }

    ps.erase(find(all(ps), near[0].second));
    ps.erase(find(all(ps), near[1].second));

    Poly triangle = to_counter_clockwise({rand_p, near[0].second, near[1].second});
    return triangle;
}

Poly build_poly(const vector<Point>& init_poly, vector<Point> rem)
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
        if (g_timer.get_elapsed() > G_TLE)
            return {};

        static int poly_index[1024][1024];
        rep(i, poly.size())
            poly_index[poly[i].y][poly[i].x] = i;


        int min_add_area = ten(9);
        pair<Point, Point> best_pair;
        vector<ull> to_remove;
        while (!cand.empty())
        {
            ull it_e = cand.top();
            cand.pop();

            auto it = decode(it_e);
            const Point& a = get<1>(it);
            const Point& b = get<2>(it);
            const Point& p = get<3>(it);

            const int poly_i = poly_index[a.y][a.x];

            int add_area = get<0>(it);
            if (b == poly[(poly_i + 1) % poly.size()] && is_remain[p.y][p.x]
//                 && !intersect(poly, a, p) && !intersect(poly, b, p))
                && !intersect(poly, poly_i, a, p) && !intersect(poly, poly_i, b, p))
            {
                if (add_area < min_add_area)
                {
                    min_add_area = add_area;
                    best_pair = make_pair(a, p);
                    break;
                }
            }
            else
                to_remove.push_back(it_e);
        }
        if (min_add_area == ten(9))
            return {};
        assert(min_add_area != ten(9));

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

Poly rebuild_poly(Poly poly)
{
    if (poly.size() <= 3)
        return {};

    int remove_begin = g_rand.next_int(poly.size());
    int remove_num = g_rand.next_int(1, min<int>(100, (int)poly.size() - 3));
//     if (g_timer.get_elapsed() > G_TLE * 0.5)
//         remove_num = g_rand.next_int(1, min<int>(10, (int)poly.size() - 3));

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
vector<Poly> solve(const vector<Point>& points, const int max_polys)
{
    auto separated = separate_points(points, max_polys);
//     auto separated = separate_points(points, 1);

    vector<Poly> polys(separated.size());
    vector<Poly> best_polys(separated.size());
    rep(i, separated.size())
    {
        Poly poly;
        while (poly.empty())
        {
            vector<Point> rem = separated[i];
            Poly init = init_rand_triangle(rem);
            poly = build_poly(init, rem);
        }
        best_polys[i] = polys[i] = poly;
    }
    dump(g_timer.get_elapsed());

    int loops = 0;
    for (int loop_i = 0; loop_i < ten(9); ++loop_i)
    {
        if (loop_i && loop_i % 10000 == 0)
        {
            rep(i, separated.size())
            {
                Poly poly;
                while (poly.empty())
                {
                    vector<Point> rem = separated[i];
                    Poly init = init_rand_triangle(rem);
                    poly = build_poly(init, rem);
                }
                polys[i] = poly;
            }
        }
        rep(i, separated.size())
        {
            if (g_timer.get_elapsed() > G_TLE)
                goto END;

            Poly next = rebuild_poly(polys[i]);
            if (!next.empty())
            {
                if (area2(next) < area2(polys[i]))
                {
                    polys[i] = next;
                    if (area2(next) < area2(best_polys[i]))
                    {
//                         fprintf(stderr, "%6d: %5.1f -> %5.1f\n", loop_i, area(best_polys[i]), area(next));
                        best_polys[i] = next;
                    }
                }
            }
        }
        ++loops;
    }
END:;
    dump(g_timer.get_elapsed());
    dump(loops);
    return best_polys;
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
