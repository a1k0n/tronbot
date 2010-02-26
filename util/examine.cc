#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#include <assert.h>
#include <map>
#include <vector>

#include "../cpp/artictbl.h"

#define TIMEOUT_USEC 950000
#define FIRSTMOVE_USEC 2950000
#define DEPTH_INITIAL 1
#define DEPTH_MAX 100
#define DRAW_PENALTY 0 // (10*itr) // -500
#define VERBOSE 1

// {{{ position
struct position {
  static const char dx[4], dy[4];

  int x,y;
  position() {}
  position(int _x, int _y): x(_x), y(_y) {}
  // direction encoding:
  //  1
  // 4 2
  //  3
  position next(int move) const { return position(x+dx[move], y+dy[move]); }
  position prev(int move) const { return position(x-dx[move], y-dy[move]); }
};

static bool operator==(const position &a, const position &b) { return a.x == b.x && a.y == b.y; }
//static bool operator<(const position &a, const position &b) { return a.x == b.x ? a.y < b.y : a.x < b.x; }

// note: the canonical order (above) is changed internally here in order to
// attain more symmetric play; this is mainly a failing of the evaluation
// function but it helps, e.g. when playing player 1 in joust
// so instead it's  1
//                 4 3
//                  2
const char position::dx[4]={ 0, 0, 1,-1};
const char position::dy[4]={-1, 1, 0, 0};
const int move_permute[4]={1,3,2,4};

static inline int _min(int a, int b) { return a<b ? a : b; }
static inline int _max(int a, int b) { return a>b ? a : b; }
// }}}

// {{{ Map
template <class T> struct Map {
  T *map;
  int width, height;
  Map() { map = NULL; }
  Map(int w, int h) { resize(w,h); }
  void resize(int w, int h) {
    width = w; height = h;
    map = new T[w*h];
    clear();
  }
  void clear(void) { memset(map, 0, width*height*sizeof(T)); }
  Map(const Map &m) { abort(); } // this shouldn't happen
  ~Map() { if(map) delete[] map; }
  T& operator()(position p) { return map[p.x + p.y*width]; }
  T& operator()(int x, int y) { return map[x + y*width]; }
  T& operator()(int idx) { return map[idx]; }
  T& M(position p) { return map[p.x + p.y*width]; }
  T& M(int x, int y) { return map[x + y*width]; }
  int idx(position p) { return p.x + p.y*width; }

  void dump(void) {
    for(int j=0;j<height;j++) {
      for(int i=0;i<width;i++) {
        int n = map[i+j*width];
        if(n == 0 || n == INT_MAX) fprintf(stderr, "  ");
        else fprintf(stderr, "%2d", n);
      }
      fprintf(stderr, "\n");
    }
  }
};
// }}}

// {{{ gamestate
struct gamestate {
  position p[2]; // position in current state
  int m[2]; // last move made

  // affects map using m[0], m[1]
  gamestate move(Map<char> M) {
    gamestate s = *this;
    M(p[0]) = 1;
    M(p[1]) = 1;
    s.p[0] = p[0].next(m[0]);
    s.p[1] = p[1].next(m[1]);
    return s;
  }
  // undoes effect on map
  void unmove(Map<char> M) {
    M(p[0]) = 0;
    M(p[1]) = 0;
  }
};

// }}}

static Map<char> M;
static Map<int> dp0, dp1;
static Map<int> low, num, articd; // for articulation point finding
static gamestate curstate;
static char _killer[DEPTH_MAX*2];
static int _maxitr=0;
static bool firstmove = true;

// {{{ imported map update garbage from original code
bool map_update()
{
  int x, y, c;
  int map_width, map_height;
  int num_items = fscanf(stdin, "%d %d\n", &map_width, &map_height);
  if (feof(stdin) || num_items < 2) {
    return false;
  }
  if(!M.map) {
    M.resize(map_width, map_height);
    dp0.resize(map_width, map_height);
    dp1.resize(map_width, map_height);
    num.resize(map_width, map_height);
    low.resize(map_width, map_height);
    articd.resize(map_width, map_height);
  }
  x = 0;
  y = 0;
  while (y < M.height && (c = fgetc(stdin)) != EOF) {
    switch (c) {
    case '\r':
      break;
    case '\n':
      if (x != M.width) {
	fprintf(stderr, "x != width in Board_ReadFromStream\n");
	return false;
      }
      ++y;
      x = 0;
      break;
    case '#':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return false;
      }
      M(x,y) = 1;
      ++x;
      break;
    case ' ':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return false;
      }
      M(x,y) = 0;
      ++x;
      break;
    case '1':
    case '2':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return false;
      }
      {
        position p(x,y);
        M(p) = 0;
        curstate.p[c - '1'] = p;
        curstate.m[c - '1'] = 0;
        ++x;
      }
      break;
    default:
      fprintf(stderr, "unexpected character %d in Board_ReadFromStream", c);
      return false;
    }
  }
  return true;
}
// }}}

// {{{ basic geometric stuff
static inline int color(position x) { return (x.x ^ x.y)&1; } // convention: 1=red, 0=black
static inline int color(int x, int y) { return (x ^ y)&1; } // convention: 1=red, 0=black

struct colorcount {
  int red, black, edges;
  colorcount() {}
  colorcount(int r, int b, int e): red(r), black(b), edges(e) {}
  int& operator()(const position &x) { return color(x) ? red : black; }
};

static colorcount operator+(const colorcount &a, const colorcount &b) { return colorcount(a.red+b.red, a.black+b.black, a.edges+b.edges); }

// number of fillable squares in area when starting on 'startcolor' (assuming starting point is not included)
int num_fillable(const colorcount &c, int startcolor) {
  if(startcolor) { // start on red?  then moves are black-red-black-red-black (2 red, 3 black: 5; 3 red 3 black: 6; 4 red 3 black
    return 2*_min(c.red-1, c.black) +
      (c.black >= c.red ? 1 : 0);
  } else { // moves are red-black-red-black-red
    return 2*_min(c.red, c.black-1) +
      (c.red >= c.black ? 1 : 0);
  }
}

static int degree(position x) {
  int idx = x.x+x.y*M.width;
  return 4 - M(idx-1) - M(idx+1) - M(idx-M.width) - M(idx+M.width);
}

static int degree(int idx) {
  return 4 - M(idx-1) - M(idx+1) - M(idx-M.width) - M(idx+M.width);
}

// return bitmask of neighbors, for table lookups
static int neighbors(position s) {
  return (M(s.x-1, s.y-1) |
          (M(s.x  , s.y-1)<<1) |
          (M(s.x+1, s.y-1)<<2) |
          (M(s.x+1, s.y  )<<3) |
          (M(s.x+1, s.y+1)<<4) |
          (M(s.x  , s.y+1)<<5) |
          (M(s.x-1, s.y+1)<<6) |
          (M(s.x-1, s.y  )<<7));
}

static int potential_articulation(position s) { return _potential_articulation[neighbors(s)]; }
// }}}

// {{{ connected components algorithm

struct Components {
  Map<int> c;
  std::vector<int> cedges, red, black;

  Components(Map<char> &M): c(M.width, M.height) { recalc(); }

  void recalc(void) {
    static std::vector<int> equiv;
    equiv.clear(); equiv.push_back(0);
    cedges.clear(); red.clear(); black.clear();
    int nextclass = 1;
    int mapbottom = M.width*(M.height-1)-1;
    for(int idx=M.width+1;idx<mapbottom;idx++) {
      if(M(idx)) continue; // wall
      int cup   = equiv[c(idx-M.width)],
          cleft = equiv[c(idx-1)];
      if(cup == 0 && cleft == 0) { // new component
        equiv.push_back(nextclass);
        c(idx) = nextclass++;
      } else if(cup == cleft) { // existing component
        c(idx) = cup;
      } else { // join components
        // deprecate the higher-numbered component in favor of the lower
        if(cleft == 0 || (cup != 0 && cup < cleft)) {
          c(idx) = cup;
          if(cleft != 0) _merge(equiv, cleft, cup);
        } else {
          c(idx) = cleft;
          if(cup != 0) _merge(equiv, cup, cleft);
        }
      }
    }
    cedges.resize(nextclass, 0);
    red.resize(nextclass, 0);
    black.resize(nextclass, 0);
    // now make another pass to translate equivalences and compute connected area
    for(int j=1,idx=M.width+1;j<M.height-1;j++,idx+=2) {
      for(int i=1;i<M.width-1;i++,idx++) {
        int e = equiv[c(idx)];
        c(idx) = e;
        cedges[e] += degree(idx);
        if(color(i,j)) red[e] ++; else black[e] ++;
      }
    }
  }

  void remove(position s) {
    c(s) = 0;
    if(potential_articulation(s)) {
      recalc();
    } else {
      cedges[c(s)] -= 2*degree(s);
      if(color(s)) red[c(s)] --; else black[c(s)] --;
    }
  }
  void add(position s) {
    for(int m=0;m<4;m++) {
      position r = s.next(m);
      if(M(r)) continue;
      if(c(s) != 0 && c(s) != c(r)) { recalc(); return; }
      c(s) = c(r);
    }
    cedges[c(s)] += 2*degree(s);
    if(color(s)) red[c(s)] ++; else black[c(s)] ++;
  }

  void dump() {
    for(size_t i=0;i<red.size();i++) {
      if(red[i])
        fprintf(stderr, "area %d: %d red %d black nodes, %d edges\n", (int)i, red[i], black[i], cedges[i]);
    }
    c.dump();
  }
  int component(const position &p) { return c(p); }
  int connectedarea(int component) { return red[component]+black[component]; }
  int connectedarea(const position &p) { return red[c(p)]+black[c(p)]; }
  // number of fillable squares in area when starting on 'startcolor' (assuming starting point is not included)
  int fillablearea(int component, int startcolor) {
    return num_fillable(colorcount(red[component], black[component], 0), startcolor);
  }
  // number of fillable squares starting from p (not including p)
  int fillablearea(const position &p) { return fillablearea(c(p), color(p)); }
  int connectedvalue(int component) { return cedges[component]; }
  int connectedvalue(const position &p) { return cedges[c(p)]; }
private:
#if 0
  int _find_equiv(std::map<int,int> &equiv, int c) {
    while(true) {
      std::map<int,int>::iterator e = equiv.find(c);
      if(e == equiv.end()) break;
      if(c < e->second)
        c = e->second;
      else
        break;
    }
    return c;
  }
#endif
  void _merge(std::vector<int> &equiv, int o, int n) {
    for(size_t k=0;k<equiv.size();k++)
      if(equiv[k] == o) equiv[k] = n;
  }
};

// }}}

// {{{ run timing
long _get_time()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_usec + tv.tv_sec*1000000;
}

static long _timer, _timeout;
static volatile bool _timed_out = false;
static int _ab_runs=0;
static int _spacefill_runs=0;

static void _alrm_handler(int sig) { _timed_out = true; }

static void reset_timer(long t)
{
  _timer = _get_time();
  itimerval timer;
  memset(&timer, 0, sizeof(timer));
  timer.it_value.tv_sec = t/1000000;
  timer.it_value.tv_usec = t%1000000;
  setitimer(ITIMER_REAL, &timer, NULL);
  _timed_out = false;
  _ab_runs = 0;
  _spacefill_runs = 0;
  _timeout = t;
}

static long elapsed_time() { return _get_time() - _timer; }
// }}}

// {{{ Dijkstra's
static void dijkstra(Map<int> &d, const position &s, Components &cp, int component)
{
  static std::vector<std::vector<position> > Q;
  static Map<int> loc;
  size_t min_dist=0;
  int siz = M.width*M.height;
  for(int idx=0;idx<siz;idx++)
    d(idx) = INT_MAX;
  if(!loc.map) loc.resize(d.width, d.height);

  Q.clear(); Q.push_back(std::vector<position>());
  Q[0].push_back(s);
  d(s) = 0;
  loc(s) = 0;
  while(min_dist != Q.size()) {
    position u = Q[min_dist].back();
    Q[min_dist].pop_back();
    for(int m=0;m<4;m++) {
      position v = u.next(m);
      if(M(v)) continue;
      int alt = 1 + d(u);
      int dist = d(v);
      if(dist == INT_MAX) {
        while(alt >= (int)Q.size())
          Q.push_back(std::vector<position>());
        int newloc = Q[alt].size();
        Q[alt].push_back(v);
        d(v) = alt;
        loc(v) = newloc;
      } else if(alt < dist) {
        // move last element to this one's spot in the pqueue
        int moveelem = Q[dist].size()-1;
        Q[dist][loc(v)] = Q[dist][moveelem];
        loc(Q[dist][moveelem]) = loc(v);
        Q[dist].pop_back();

        d(v) = alt;
        loc(v) = Q[alt].size()-1;
        Q[alt].push_back(v);
      }
    }
    while(min_dist < Q.size() && Q[min_dist].empty()) min_dist++;
  }
}

// }}}

// {{{ space-filling

static int floodfill(Components &ca, position s, bool fixup=true)
{
  // flood fill heuristic: choose to remove as few edges from the graph as
  // possible (in other words, move onto the square with the lowest degree)
  int bestv=0;
  position b = s;
  for(int m=0;m<4;m++) {
    position p = s.next(m);
    if(M(p)) continue;
    int v = ca.connectedvalue(p) + ca.fillablearea(p) - 2*degree(p) -
      4*potential_articulation(p);
    if(v > bestv) { bestv = v; b = p; }
  }
  if(bestv == 0)
    return 0;
  M(b) = 1; ca.remove(b);
  int a = 1+floodfill(ca, b);
  M(b) = 0; if(fixup) ca.add(b);
  return a;
}

// returns spaces unused (wasted); idea is to minimize waste
static int _spacefill(int &move, Components &ca, position p, int itr) {
  int bestv = 0;
  int spacesleft = ca.fillablearea(p);
  if(degree(p) == 0) { move=1; return 0; }
  if(_timed_out) {
    return 0;
  }
  if(itr == 0)
    return floodfill(ca, p);
  for(int m=0;m<4 && !_timed_out;m++) {
    position r = p.next(m);
    if(M(r)) continue;
    M(r) = 1; ca.remove(r);
    int _m, v = 1+_spacefill(_m, ca, r, itr-1);
    M(r) = 0; ca.add(r);
    if(v > bestv) { bestv = v; move = m; }
    if(v == spacesleft) break; // we solved it!
    if(itr == 0) break; // we can only use the first-chosen solution
  }
  return bestv;
}

// space-filling iterative deepening search
static int next_move_spacefill(Components &ca)
{
  int itr;
  int area = ca.fillablearea(curstate.p[0]);
  int bestv = 0, bestm = 1;
  for(itr=DEPTH_INITIAL;itr<DEPTH_MAX && !_timed_out;itr++) {
    int m;
    _maxitr = itr;
    int v = _spacefill(m, ca, curstate.p[0], itr);
    if(v > bestv) { bestv = v; bestm = m; }
    if(v <= itr) break; // we can't possibly search any deeper
#if VERBOSE >= 1
    struct timeval tv;
    gettimeofday(&tv, NULL);
    //M.dump();
    fprintf(stderr, "%d.%06d: area=%d/%d waste=%d (m=%d) @depth %d _spacefill_runs=%d\n", (int) tv.tv_sec, (int) tv.tv_usec, v, area, area-v, m, itr, _spacefill_runs);
#endif
    if(v >= area) break; // solved!
  }
  return bestm;
}
// }}}

// {{{ heuristic board evaluation

static int _art_counter=0;
static void reset_articulations()
{
  _art_counter=0;
  low.clear();
  num.clear();
  articd.clear();
}

// calculate articulation vertices within our voronoi region
// algorithm taken from http://www.eecs.wsu.edu/~holder/courses/CptS223/spr08/slides/graphapps.pdf
// DFS traversal of graph
static int calc_articulations(Map<int> *dp0, Map<int> *dp1, const position &v, int parent=-1)
{
  int nodenum = ++_art_counter;
  low(v) = num(v) = nodenum; // rule 1
  int children=0;
  int count = 0;
  for(int m=0;m<4;m++) {
    position w = v.next(m);
    if(M(w)) continue;
    if(dp0 && (*dp0)(w) >= (*dp1)(w)) continue; // filter out nodes not in our voronoi region
    if(!num(w)) { // forward edge
      children++;
      count += calc_articulations(dp0, dp1, w, nodenum);
      if(low(w) >= nodenum && parent != -1) {
        if(!articd(v)) count ++;
        articd(v) = 1;
      }
      if(low(w) < low(v)) low(v) = low(w);   // rule 3
    } else {
      if(num(w) < nodenum) { // back edge
        if(num(w) < low(v)) low(v) = num(w); // rule 2
      }
    }
  }
  if(parent == -1 && children > 1) {
    if(!articd(v)) count ++;
    articd(v) = 1;
  }
  return count;
}

// returns the maximum "weight" of connected reachable components: we find the
// "region" bounded by all articulation points, traverse each adjacent region
// recursively, and return the maximum traversable area
static colorcount _explore_space(Map<int> *dp0, Map<int> *dp1, std::vector<position> &exits, const position &v)
{
  colorcount c(0,0,0);
  c(v) ++;
  num(v) = 0;
  if(articd(v)) {
    // we're an articulation vertex; nothing to do but populate the exits
    for(int m=0;m<4;m++) {
      position w = v.next(m);
      if(M(w)) continue;
      c.edges++;
      if(dp0 && (*dp0)(w) >= (*dp1)(w)) { continue; }
      if(!num(w)) continue; // use 'num' from articulation vertex pass to mark nodes used
      exits.push_back(w);
    }
  } else {
    // this is a non-articulation vertex
    for(int m=0;m<4;m++) {
      position w = v.next(m);
      if(M(w)) continue;
      c.edges++;

      // filter out nodes not in our voronoi region
      if(dp0 && (*dp0)(w) >= (*dp1)(w)) { continue; }

      if(!num(w)) continue; // use 'num' from articulation vertex pass to mark nodes used
      if(articd(w)) { // is this vertex articulated?  then add it as an exit and don't traverse it yet
        num(w) = 0; // ensure only one copy gets pushed in here
        exits.push_back(w);
      } else {
        c = c + _explore_space(dp0,dp1,exits,w);
      }
    }
  }
  return c;
}

// this assumes the space is separated into a DAG of chambers
// if cycles or bidirectional openings really do exist, then we just get a bad estimate :/
static colorcount max_articulated_space(Map<int> *dp0, Map<int> *dp1, const position &v)
{
  std::vector<position> exits;
  colorcount space = _explore_space(dp0,dp1,exits,v);
  //fprintf(stderr, "space@%d,%d = (%d,%d,%d) exits: ", v.x,v.y, space.red, space.black, space.edges);
  //for(size_t i=0;i<exits.size();i++) fprintf(stderr, "%d,%d ", exits[i].x, exits[i].y);
  //fprintf(stderr, "\n");
  colorcount maxchild(0,0,0);
  int maxsteps=0;
  int entrancecolor = color(v);
  int localsteps[2] = {
    num_fillable(colorcount(space.red, space.black+1, 0), entrancecolor),
    num_fillable(colorcount(space.red+1, space.black, 0), entrancecolor)};
  for(size_t i=0;i<exits.size();i++) {
    int exitcolor = color(exits[i]);
    // space includes our entrance but not our exit node
    colorcount child = max_articulated_space(dp0,dp1,exits[i]);
    // child includes our exit node
    int steps = localsteps[exitcolor] + num_fillable(child, exitcolor);
    // now we need to figure out how to connect spaces via colored articulation vertices
    // exits[i] gets counted in the child space
    if(steps > maxsteps) {
      //fprintf(stderr, "space@%d,%d exit #%d steps=%d; new max\n", v.x, v.y, i, steps);
      maxsteps=steps; maxchild = child;
    }
  }
  return space+maxchild;
}

int _evaluate_territory(int *indicators, const gamestate &s, Components &cp, int comp, bool vis)
{
  dijkstra(dp0, s.p[0], cp, comp);
  dijkstra(dp1, s.p[1], cp, comp);
  reset_articulations();
  M(s.p[0])=0; M(s.p[1])=0;
  int a0 = calc_articulations(&dp0, &dp1, s.p[0]),
      a1 = calc_articulations(&dp1, &dp0, s.p[1]);
  colorcount ccount0 = max_articulated_space(&dp0, &dp1, s.p[0]),
             ccount1 = max_articulated_space(&dp1, &dp0, s.p[1]);
  int nc0_ = num_fillable(ccount0, color(s.p[0])),
      nc1_ = num_fillable(ccount1, color(s.p[1]));
  int nodecount = 0;
  memset(indicators, 0, 8*sizeof(int));
//  indicators[6] = INT_MAX;
//  // find distance
//  for(int m=1;m<=4;m++) { int d = dp0(s.p[1].next(m)); if(d < indicators[6]) indicators[6] = d+1; }

  indicators[0] = nc0_;
  indicators[1] = ccount0.edges;
  indicators[2] = a0;
  indicators[3] = nc1_;
  indicators[4] = ccount1.edges;
  indicators[5] = a1;

  //indicators[7] = num_paths(n, dp0, s.p[0], s.p[1]);
  return nodecount;
}

// }}}

void examine_board() {
  static int connected = false;
  static int movecount = 0, movestoclose = 0;
  static gamestate nextstate;
  static int v0, v1;
  Components cp(M);
  int comp;
  M(curstate.p[0]) = 1;
  M(curstate.p[1]) = 1;
  if((comp = cp.component(curstate.p[0])) == cp.component(curstate.p[1])) {
    if(!connected) {
      if(movecount == 0) {
        fprintf(stderr, "one of the players is stupid; giving up\n");
        exit(0);
      }
      movestoclose = 0;
      // here's where they first separated; count the maximum possible moves and emit a score
      assert(M(nextstate.p[0]) == 0);
      assert(M(nextstate.p[1]) == 0);
      Components cp2(M);
      M(nextstate.p[0]) = 1;
      M(nextstate.p[1]) = 1;
      assert(cp2.component(nextstate.p[0]) != cp2.component(nextstate.p[1]));
      int m;
      reset_timer(TIMEOUT_USEC);
      v0 = _spacefill(m, cp2, nextstate.p[0], 5);
      v1 = _spacefill(m, cp2, nextstate.p[1], 5);
      M(nextstate.p[0]) = 0;
      M(nextstate.p[1]) = 0;
      //fprintf(stderr, "v0 = %d v1 = %d movecount = %d at separation point\n", v0, v1, movecount);
    }
    int indicators[6];
    _evaluate_territory(indicators, curstate, cp, comp, true);
    printf("%d %d %d ", v0, v1, movestoclose);
    for(int i=0;i<6;i++) printf("%d ", indicators[i]);
    printf("\n");
    connected = true;
  } else {
    // ignore endgame!
  }
  movecount++; movestoclose++;
  nextstate = curstate; // we're moving backwards through time so the current state becomes the next one
}

int main() {
  memset(_killer, 1, sizeof(_killer));
  while (map_update()) {
    examine_board();
  }
  return 0;
}

// vim: sw=2:ts=8:et:foldmethod=marker
