#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#include <assert.h>
#include <signal.h>
#include <map>
#include <vector>

#include "artictbl.h"

#define TIMEOUT_USEC 990000
#define FIRSTMOVE_USEC 2950000
#define DEPTH_INITIAL 1
#define DEPTH_MAX 100
#define DRAW_PENALTY 0 // -itr // -500
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

bool operator==(const position &a, const position &b) { return a.x == b.x && a.y == b.y; }
bool operator<(const position &a, const position &b) { return a.x == b.x ? a.y < b.y : a.x < b.x; }

// note: the canonical order (above) is changed internally here in order to
// attain more symmetric play; this is mainly a failing of the evaluation
// function but it helps, e.g. when playing player 1 in joust
// so instead it's  1
//                 4 3
//                  2
const char position::dx[4]={ 0, 0, 1,-1};
const char position::dy[4]={-1, 1, 0, 0};
const int move_permute[4]={1,3,2,4};
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
  T& M(position p) { return map[p.x + p.y*width]; }
  T& M(int x, int y) { return map[x + y*width]; }

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
static char _killer[DEPTH_MAX*2+1];
static int _maxitr=0;

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
  for(int i=0;i<M.width;i++) { M(i,0) = 1; M(i,M.height-1)=1; }
  for(int j=0;j<M.height;j++) { M(0,j) = 1; M(M.width-1,j)=1; }
  return true;
}
// }}}

// {{{ basic geometric stuff
int runout(position p, int dir) {
  int r = 0;
  while(!M(p)) { r++; p = p.next(dir); }
  return r;
}

int degree(position x) {
  return 4 - M(x.next(1)) - M(x.next(2)) - M(x.next(3)) - M(x.next(0));
}

// return bitmask of neighbors, for table lookups
int neighbors(position s) {
  return (M(s.x-1, s.y-1) |
          (M(s.x  , s.y-1)<<1) |
          (M(s.x+1, s.y-1)<<2) |
          (M(s.x+1, s.y  )<<3) |
          (M(s.x+1, s.y+1)<<4) |
          (M(s.x  , s.y+1)<<5) |
          (M(s.x-1, s.y+1)<<6) |
          (M(s.x-1, s.y  )<<7));
}

int potential_articulation(position s) { return _potential_articulation[neighbors(s)]; }

int turn(int d, int n) { return 1+(d-1+n)&3; }
// }}}

// {{{ connected components algorithm

struct Components {
  Map<int> c;
  std::map<int,int> cedges, csize;

  Components(Map<char> &M): c(M.width, M.height) { recalc(); }

  void recalc(void) {
    static std::vector<int> equiv;
    equiv.clear(); equiv.push_back(0);
    cedges.clear(); csize.clear();
    int nextclass = 1;
    for(int j=1;j<M.height-1;j++) {
      for(int i=1;i<M.width-1;i++) {
        if(M(i,j)) continue; // wall
        int cup   = equiv[c(i, j-1)],
            cleft = equiv[c(i-1, j)];
        if(cup == 0 && cleft == 0) { // new component
          equiv.push_back(nextclass);
          c(i,j) = nextclass++;
        } else if(cup == cleft) { // existing component
          c(i,j) = cup;
        } else { // join components
          // deprecate the higher-numbered component in favor of the lower
          if(cleft == 0 || (cup != 0 && cup < cleft)) {
            c(i,j) = cup;
            if(cleft != 0) _merge(equiv, cleft, cup);
          } else {
            c(i,j) = cleft;
            if(cup != 0) _merge(equiv, cup, cleft);
          }
        }
      }
    }
    // now make another pass to compute connected area
    for(int j=1;j<M.height-1;j++) {
      for(int i=1;i<M.width-1;i++) {
        c(i,j) = equiv[c(i,j)];
        cedges[c(i,j)] += degree(position(i,j));
        csize[c(i,j)] ++;
      }
    }
  }

  void remove(position s) {
    c(s) = 0;
    if(potential_articulation(s)) {
      recalc();
    } else {
      cedges[c(s)] -= 2*degree(s);
      csize[c(s)] --;
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
    csize[c(s)] ++;
  }

  void dump() {
    std::map<int,int>::iterator i;
    for(i=csize.begin();i!=csize.end();i++) {
      fprintf(stderr, "area %d: %d nodes\n", i->first, i->second);
    }
    for(i=cedges.begin();i!=cedges.end();i++) {
      fprintf(stderr, "area %d: %d edges\n", i->first, i->second);
    }
    c.dump();
  }
  int component(const position &p) { return c(p); }
  int connectedarea(int component) { return csize[component]; }
  int connectedarea(const position &p) { return csize[c(p)]; }
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

void _alrm_handler(int sig) { _timed_out = true; }

void reset_timer(long t)
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

void stop_timer(void)
{
  itimerval timer = { { 0, 0 }, { 0, 0 } };
  setitimer(ITIMER_REAL, &timer, NULL);
}


long elapsed_time() { return _get_time() - _timer; }
// }}}

// {{{ Dijkstra's
void dijkstra(Map<int> &d, const position &s, Components &cp, int component)
{
  static std::vector<std::vector<position> > Q;
  static Map<int> loc;
  size_t min_dist=0;
  int i,j;
  for(j=0;j<M.height;j++)
    for(i=0;i<M.width;i++)
      d(i,j) = INT_MAX;
  if(!loc.map) loc.resize(d.width, d.height);

  Q.clear(); Q.push_back(std::vector<position>());
  Q[0].push_back(s);
  d(s) = 0;
  loc(s) = 0;
  while(min_dist != Q.size()) {
    position u = *(Q[min_dist].begin());
    Q[min_dist].erase(Q[min_dist].begin());
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

int xgradient(Map<int> &d, const position &s)
{
  int g=0,n=-1;
  if(!M(s.x-1,s.y)) { n++; g+=d(s.x-1,s.y) - d(s); }
  if(!M(s.x+1,s.y)) { n++; g+=d(s) - d(s.x+1,s.y); }
  return g>>n;
}

int ygradient(Map<int> &d, const position &s)
{
  int g=0,n=-1;
  if(!M(s.x,s.y-1)) { n++; g+=d(s.x,s.y)-1 - d(s); }
  if(!M(s.x,s.y+1)) { n++; g+=d(s) - d(s.x,s.y+1); }
  return g>>n;
}

// }}}

// {{{ space-filling

int floodfill(Components &ca, position s, bool fixup=true)
{
  // flood fill heuristic: choose to remove as few edges from the graph as
  // possible (in other words, move onto the square with the lowest degree)
  int bestv=0;
  position b = s;
  for(int m=0;m<4;m++) {
    position p = s.next(m);
    if(M(p)) continue;
    int v = ca.connectedvalue(p) + ca.connectedarea(p) - 1 - 2*degree(p) -
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
int _spacefill(int &move, Components &ca, position p, int itr) {
  int bestv = 0;
  int spacesleft = ca.connectedarea(p)-1;
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
int next_move_spacefill(Components &ca)
{
  int itr;
  int area = ca.connectedarea(curstate.p[0])-1;
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
void reset_articulations()
{
  _art_counter=0;
  low.clear();
  num.clear();
  articd.clear();
}

// calculate articulation vertices within our voronoi region
// algorithm taken from http://www.eecs.wsu.edu/~holder/courses/CptS223/spr08/slides/graphapps.pdf
// DFS traversal of graph
void calc_articulations(Map<int> &dp0, Map<int> &dp1, const position &v, int parent=-1)
{
  int nodenum = ++_art_counter;
  low(v) = num(v) = nodenum; // rule 1
  int children=0;
  for(int m=0;m<4;m++) {
    position w = v.next(m);
    if(M(w)) continue;
    if(dp0(w) >= dp1(w)) continue; // filter out nodes not in our voronoi region
    if(!num(w)) { // forward edge
      children++;
      calc_articulations(dp0, dp1, w, nodenum);
      if(low(w) >= nodenum && parent != -1)
        articd(v) = 1;
      if(low(w) < low(v)) low(v) = low(w);   // rule 3
    } else {
      if(num(w) < nodenum) { // back edge
        if(num(w) < low(v)) low(v) = num(w); // rule 2
      }
    }
  }
  if(parent == -1 && children > 1) {
    articd(v) = 1;
  }
}

// returns the maximum "weight" of connected reachable components: we find the
// "region" bounded by all articulation points, traverse each adjacent region
// recursively, and return the maximum traversable area
int _explore_space(Map<int> &dp0, Map<int> &dp1, std::vector<position> &exits, const position &v)
{
  int nodecount=1, edgecount=0, childcount=0;
  num(v) = 0;
  if(articd(v)) {
    // we're an articulation vertex; nothing to do but populate the exits
    for(int m=0;m<4;m++) {
      position w = v.next(m);
      if(M(w)) continue;
      edgecount++;
      if(dp0(w) >= dp1(w)) { continue; }
      if(!num(w)) continue; // use 'num' from articulation vertex pass to mark nodes used
      exits.push_back(w);
    }
  } else {
    // this is a non-articulation vertex
    for(int m=0;m<4;m++) {
      position w = v.next(m);
      if(M(w)) continue;
      edgecount++;

      // filter out nodes not in our voronoi region
      if(dp0(w) >= dp1(w)) { continue; }

      if(!num(w)) continue; // use 'num' from articulation vertex pass to mark nodes used
      if(articd(w)) { // is this vertex articulated?  then add it as an exit and don't traverse it yet
        num(w) = 0; // ensure only one copy gets pushed in here
        exits.push_back(w);
      } else {
        childcount += _explore_space(dp0,dp1,exits,w);
      }
    }
  }
  return 51*nodecount+170*edgecount+8*potential_articulation(v)+childcount;
  //return nodecount+edgecount+childcount;
  //return edgecount+childcount;
}

int max_articulated_space(Map<int> &dp0, Map<int> &dp1, const position &v)
{
  std::vector<position> exits;
  int space = _explore_space(dp0,dp1,exits,v);
  //fprintf(stderr, "space@%d,%d = %d exits: ", v.x,v.y, space);
  //for(size_t i=0;i<exits.size();i++) fprintf(stderr, "%d,%d ", exits[i].x, exits[i].y);
  //fprintf(stderr, "\n");
  int maxchild = 0;
  for(size_t i=0;i<exits.size();i++) {
    int child = max_articulated_space(dp0,dp1,exits[i]);
    if(child > maxchild) maxchild = child;
  }
  return space+maxchild;
}

int _evaluate_territory(const gamestate &s, Components &cp, int comp, bool vis)
{
  dijkstra(dp0, s.p[0], cp, comp);
  dijkstra(dp1, s.p[1], cp, comp);
  reset_articulations();
  M(s.p[0])=0; M(s.p[1])=0;
  calc_articulations(dp0, dp1, s.p[0]);
  calc_articulations(dp1, dp0, s.p[1]);
  int nc0_ = max_articulated_space(dp0, dp1, s.p[0]),
      nc1_ = max_articulated_space(dp1, dp0, s.p[1]);
#if VERBOSE >= 2
  int nc0=0, nc1=0;
  for(int j=0;j<M.height;j++)
    for(int i=0;i<M.width;i++) {
      position p(i,j);
      int diff = dp0(i,j) - dp1(i,j);
      // if the opponent's distance is shorter than ours, then this is "their"
      // node, otherwise it's ours
      //if(diff>0) { nc1 += degree(p); }//if(vis) fprintf(stderr, "nc1:(%d,%d)\n", p.x,p.y); }
      //else if(diff<0) { nc0 += degree(p); }//if(vis) fprintf(stderr, "nc0:(%d,%d)\n", p.x,p.y); }
      if(diff>0) { nc1 += 51 + 170*degree(p) + 8*potential_articulation(p); }
      else if(diff<0) { nc0 += 51 + 170*degree(p) + 8*potential_articulation(p); }
    }
#endif
  M(s.p[0])=1; M(s.p[1])=1;
  int nodecount = nc0_ - nc1_;
#if VERBOSE >= 2
  if(vis) {
    for(int j=0;j<M.height;j++) {
      for(int i=0;i<M.width;i++) {
        if(dp0(i,j) == INT_MAX) fprintf(stderr,M(i,j) ? " #" : "  ");
        else fprintf(stderr,"%2d", dp0(i,j));
      }
      fprintf(stderr," ");
      for(int i=0;i<M.width;i++) {
        if(dp1(i,j) == INT_MAX) fprintf(stderr,M(i,j) ? " #" : "  ");
        else fprintf(stderr,"%2d", dp1(i,j));
      }
      fprintf(stderr," ");
      for(int i=0;i<M.width;i++) {
        int d = dp1(i,j)-dp0(i,j);
        if(articd(i,j))
          fprintf(stderr,"-");
        else if(d == INT_MAX || d == -INT_MAX)
          fprintf(stderr,"#");
        else if(d == 0) fprintf(stderr,".");
        else {
          d = d<0 ? 2 : d>0 ? 1 : 0;
          fprintf(stderr,"%d", d);
        }
      }
      fprintf(stderr,"\n");
    }
    fprintf(stderr, "nodecount: %d (0: %d/%d, 1: %d/%d)\n", nodecount, nc0_,nc0, nc1_,nc1);
#if 0
    for(int j=0;j<M.height;j++) {
      for(int i=0;i<M.width;i++) {
        if(num(i,j) == 0) fprintf(stderr,"  %c", M(i,j) ? '#' : '.');
        else fprintf(stderr,"%3d", num(i,j));
      }
      fprintf(stderr," ");
      for(int i=0;i<M.width;i++) {
        if(low(i,j) == 0) fprintf(stderr,"  %c", M(i,j) ? '#' : '.');
        else fprintf(stderr,"%3d", low(i,j));
      }
      fprintf(stderr," ");
      for(int i=0;i<M.width;i++) {
        int d = num(i,j)-low(i,j);
        if(num(i,j) == 0)
          fprintf(stderr, " #");
        else if(d <= 0)
          fprintf(stderr," *");
        else fprintf(stderr," .");
      }
      fprintf(stderr,"\n");
    }
#endif
  }
#endif
  return nodecount;
}

static int evaluations=0;
int _evaluate_board(gamestate s, int player, bool vis=false)
{
  assert(player == 0); // we're always searching an even number of plies

  // remove players from the board when evaluating connected components,
  // because if a player is separating components he still gets to choose which
  // one to move into.
  M(s.p[0]) = 0; M(s.p[1]) = 0;
  Components cp(M); // pre-move components
  M(s.p[0]) = 1; M(s.p[1]) = 1;

  if(s.p[0] == s.p[1])
    return 0; // crash!

  evaluations++;
#if VERBOSE >= 2
  if(vis) {
    fprintf(stderr, "evaluating board: \n");
    M(s.p[0]) = 2; M(s.p[1]) = 3; M.dump();
    M(s.p[0]) = 1; M(s.p[1]) = 1;
  }
#endif
  int comp;
  // greedily follow the maximum territory gain strategy until we partition
  // space or crash
  if((comp = cp.component(s.p[0])) == cp.component(s.p[1])) {
    int v = _evaluate_territory(s, cp, comp, vis);
    return v;
  }

  // since each bot is in a separate component by definition here, it's OK to
  // destructively update cp for floodfill()
#if VERBOSE >= 2
  int cc0 = cp.connectedarea(s.p[0]);
  int cc1 = cp.connectedarea(s.p[1]);
#endif
  int _m;
  //int ff0 = floodfill(cp, s.p[0], false);
  //int ff1 = floodfill(cp, s.p[1], false);
  int ff0 = _spacefill(_m, cp, s.p[0], 1);
  int ff1 = _spacefill(_m, cp, s.p[1], 1);
  int v = 10000*(ff0-ff1);
  if(player == 1) v = -v;
#if VERBOSE >= 2
  if(vis) {
    fprintf(stderr, "player=%d connectedarea value: %d (0:%d/%d 1:%d/%d)\n", player, v, ff0,cc0, ff1,cc1);
  }
#endif
  return v;
}
// }}}

// {{{ alpha-beta iterative deepening search

// do an iterative-deepening search on all moves and see if we can find a move
// sequence that cuts off our opponent
int _alphabeta(char *moves, gamestate s, int player, int a, int b, int itr)
{
  // base cases: no more moves?  draws?
  *moves=1; // set default move
  _ab_runs++;
  if(s.p[0] == s.p[1]) { return DRAW_PENALTY; } // crash!  draw!
  if(degree(s.p[player]) == 0) {
    if(degree(s.p[player^1]) == 0) { // both boxed in; draw
      return DRAW_PENALTY;
    }
    return -INT_MAX;
  }
  if(degree(s.p[player^1]) == 0) {
    // choose any move
    for(int m=0;m<4;m++) if(!M(s.p[player].next(m))) break;
    return INT_MAX;
  }

  if(_timed_out) {
#if VERBOSE >= 1
    fprintf(stderr, "timeout; a=%d b=%d itr=%d\n", a,b,itr);
#endif
    return a;
  }

  // last iteration?
  if(itr == 0) {
#if VERBOSE >= 3
    int v = _evaluate_board(s, player, true);
    fprintf(stderr, "_alphabeta(itr=%d [%d,%d,%d]|[%d,%d,%d] p=%d a=%d b=%d) -> %d\n",
            itr, s.p[0].x, s.p[0].y, s.m[0],
            s.p[1].x, s.p[1].y, s.m[1], player, a,b,v);
#else
    int v = _evaluate_board(s, player);
#endif
    return v;
  }
#if VERBOSE >= 3
  fprintf(stderr, "_alphabeta(itr=%d [%d,%d,%d]|[%d,%d,%d] p=%d a=%d b=%d)\n",
          itr, s.p[0].x, s.p[0].y, s.m[0],
          s.p[1].x, s.p[1].y, s.m[1], player, a,b);
#endif

  // periodically check timeout.  if we do time out, give up, we can't do any
  // more work; whatever we found so far will have to do
  int kill = _killer[_maxitr-itr];
  char bestmoves[DEPTH_MAX*2+1];
  memset(bestmoves, 0, itr);
  for(int _m=-1;_m<4 && !_timed_out;_m++) {
    // convoluted logic: do "killer heuristic" move first
    if(_m == kill) continue;
    int m = _m == -1 ? kill : _m;
    if(M(s.p[player].next(m))) // impossible move?
      continue;
    gamestate r = s;
    r.m[player] = m;
    // after both players 0 and 1 make their moves, the game state updates
    if(player == 1) {
      r.p[0] = s.p[0].next(r.m[0]);
      r.p[1] = s.p[1].next(r.m[1]);
      M(r.p[0]) = 1;
      M(r.p[1]) = 1;
    }
    int a_ = -_alphabeta(moves+1, r, player^1, -b, -a, itr-1);
    if(a_ > a) {
      a = a_;
      bestmoves[0] = m;
      _killer[_maxitr-itr] = m;
      memcpy(bestmoves+1, moves+1, itr-1);
    }
    // undo game state update
    if(player == 1) {
      M(r.p[0]) = 0;
      M(r.p[1]) = 0;
      r.p[0] = s.p[0];
      r.p[1] = s.p[1];
    }

    if(_timed_out) // a_ is garbage if we timed out
      return -INT_MAX;

    if(a >= b) // beta cut-off
      break;
  }
  memcpy(moves, bestmoves, itr);
  return a;
}

int next_move_alphabeta()
{
  int itr;
  int lastv = -INT_MAX, lastm = 1;
  evaluations=0;
  char moves[DEPTH_MAX*2+1];
  memset(moves, 0, sizeof(moves));
  for(itr=DEPTH_INITIAL;itr<DEPTH_MAX && !_timed_out;itr++) {
    _maxitr = itr*2;
    int v = _alphabeta(moves, curstate, 0, -INT_MAX, INT_MAX, itr*2);
#if VERBOSE >= 1
    struct timeval tv;
    gettimeofday(&tv, NULL);
    //M.dump();
    fprintf(stderr, "%d.%06d: v=%d m=[", (int) tv.tv_sec, (int) tv.tv_usec, v);
    for(int i=0;i<(itr < 10 ? itr*2 : 20);i++) fprintf(stderr, "%d", move_permute[(int)moves[i]]);
    fprintf(stderr, "] @depth %d _ab_runs=%d\n",
            itr*2, _ab_runs);
#endif
    if(v == INT_MAX) // our opponent cannot move, so we win
      return moves[0];
    if(v == -INT_MAX) {
      // deeper searching is apparently impossible (either because there are no
      // more moves for us or because we don't have any search time left)
      break;
    }
    lastv = v;
    lastm = moves[0];
    memcpy(_killer, moves, itr*2);
  }
#if VERBOSE >= 1
  long e = elapsed_time();
  float rate = (float)evaluations*1000000.0/(float)e;
  fprintf(stderr, "%d evals in %ld us; %0.1f evals/sec; lastv=%d move=%d\n", evaluations, e, rate, lastv, move_permute[lastm]);
  if(e > TIMEOUT_USEC*11/10) {
    fprintf(stderr, "10%% timeout violation: %ld us\n", e);
  }
#endif
  memmove(_killer, _killer+2, sizeof(_killer)-2); // shift our best-move tree forward to accelerate next move's search
  return lastm;
}
// }}}

int next_move() {
  Components cp(M);
#if VERBOSE >= 2
  cp.dump();
  _evaluate_board(curstate, 0, true);
#endif
  M(curstate.p[0]) = 1;
  M(curstate.p[1]) = 1;
  if(degree(curstate.p[0]) == 1) {
    // only one possible move we can make, so make it and don't waste any time
    for(int m=0;m<4;m++)
      if(!M(curstate.p[0].next(m)))
        return m;
  }
  if(cp.component(curstate.p[0]) == cp.component(curstate.p[1])) {
    // start-midgame: try to cut off our opponent
    return next_move_alphabeta();
  } else {
    // endgame: use up space as efficiently as we can, and hope we have more
    // left than they do.
    return next_move_spacefill(cp);
  }
}

int main(int argc, char **argv) {
  memset(_killer, 0, sizeof(_killer));
  bool firstmove = true;
  signal(SIGALRM, _alrm_handler);
  setlinebuf(stdout);
  while (map_update()) {
    if(argc>1 && atoi(argv[1])) {
      position p = curstate.p[0];
      curstate.p[0] = curstate.p[1];
      curstate.p[1] = p;
    }
    if(argc>2 && atoi(argv[2])) {} else {
      reset_timer(firstmove ? FIRSTMOVE_USEC : TIMEOUT_USEC);
    }
    firstmove=false;
    printf("%d\n", move_permute[next_move()]);
  }
//#if VERBOSE >= 1
//  fprintf(stderr, "%d evaluations\n", evaluations);
//#endif
  return 0;
}

// vim: sw=2:ts=8:et:foldmethod=marker
