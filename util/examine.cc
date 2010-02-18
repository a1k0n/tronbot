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
  static const char dx[5], dy[5];

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

// note: the canonical order (above) is changed internally here in order to
// attain more symmetric play; this is mainly a failing of the evaluation
// function but it helps, e.g. when playing player 1 in joust
// so instead it's  1
//                 4 3
//                  2
const char position::dx[5]={0, 0, 0, 1,-1};
const char position::dy[5]={0,-1, 1, 0, 0};
const int move_permute[5]={0,1,3,2,4};
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
    memset(map, 0, w*h*sizeof(T));
  }
  Map(const Map &m) { abort(); } // this shouldn't happen
  ~Map() { if(map) delete[] map; }
  T& operator()(position p) { return map[p.x + p.y*width]; }
  T& operator()(int x, int y) { return map[x + y*width]; }
  T& M(position p) { return map[p.x + p.y*width]; }
  T& M(int x, int y) { return map[x + y*width]; }

  position argmin(void) {
    position p(0,0);
    for(int j=0;j<height;j++) {
      for(int i=0;i<width;i++) {
        if(M(i,j) < M(p))
          p = position(i,j);
      }
    }
    return p;
  }

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
int runout(position p, int dir) {
  int r = 0;
  while(!M(p)) { r++; p = p.next(dir); }
  return r;
}

int degree(position x) {
  return 4 - M(x.next(1)) - M(x.next(2)) - M(x.next(3)) - M(x.next(4));
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
    for(int m=1;m<=4;m++) {
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
static bool _timed_out = false;
static int _ab_runs=0;
void reset_timer(long t)
{
  _timer = _get_time();
  _timed_out = false;
  _ab_runs = 0;
  _timeout = t;
}

long elapsed_time() { return _get_time() - _timer; }
bool timeout() { _timed_out = elapsed_time() > _timeout; return _timed_out; }
// }}}

// {{{ Dijkstra's
void dijkstra(Map<int> &d, const position &s, Components &cp, int component)
{
  static std::vector<position> Q;
  int i,j;
  for(j=0;j<M.height;j++)
    for(i=0;i<M.width;i++) {
      d(i,j) = INT_MAX;
      if(cp.c(i,j) != component) continue;
      if(M(i,j)) continue; // the player and his opponent are considered walls here
      Q.push_back(position(i,j));
    }
  Q.push_back(s);
  d(s) = 0;
  while(!Q.empty()) {
    position u(0,0);
    int min_d = INT_MAX, min_i=0;
    for(i=0;i<(int)Q.size();i++)
      if(d(Q[i]) < min_d) { u = Q[i]; min_d = d(u); min_i = i; }
    Q[min_i] = Q[Q.size()-1]; Q.pop_back();
    if(min_d == INT_MAX) continue; // ??
    for(int m=1;m<=4;m++) {
      position v = u.next(m);
      if(M(v)) continue;
      int alt = 1 + d(u);
      if(alt < d(v))
        d(v) = alt;
    }
  }
}

// return number of edge-independent paths from source to sink
// uses dijkstra's output to compute it
int num_paths(Map<int> &_n, Map<int> &d, const position &s, const position &sink)
{
  if(_n(s)) return _n(s);
  int cur = d(s);
  int n = 0;
  for(int m=1;m<=4;m++) {
    position t = s.next(m);
    if(t == sink) { n++; continue; }
    if(M(t)) continue;
    if(d(t) <= cur) continue; // must be strictly increasing distance from source
    n += num_paths(_n, d, t, sink);
  }
  _n(s) = n;
  return n;
}


// }}}

// {{{ space-filling

int floodfill(Components &ca, position s, bool fixup=true)
{
  // flood fill heuristic: choose to remove as few edges from the graph as
  // possible (in other words, move onto the square with the lowest degree)
  int bestv=0;
  position b = s;
  for(int m=1;m<=4;m++) {
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
int _spacefill_runs=0;
int _spacefill(int &move, Components &ca, position p, int itr) {
  int bestv = 0;
  int spacesleft = ca.connectedarea(p);
  if(degree(p) == 0) { move=1; return 0; }
  if(_timed_out || ((_spacefill_runs++)&63) == 0 && timeout()) {
    return 0;
  }
  if(itr == 0)
    return floodfill(ca, p);
  for(int m=1;m<=4 && !_timed_out;m++) {
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
  int area = ca.connectedarea(curstate.p[0]);
  int bestv = 0, bestm = 1;
  reset_timer(firstmove ? FIRSTMOVE_USEC : TIMEOUT_USEC);
  firstmove=false;
  for(itr=DEPTH_INITIAL;itr<DEPTH_MAX && !timeout();itr++) {
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
  fprintf(stderr, "moving %d\n", bestm);
  return bestm;
}
// }}}

// {{{ heuristic board evaluation

int _evaluate_territory(int *indicators, const gamestate &s, Components &cp, int comp, bool vis)
{
  dijkstra(dp0, s.p[0], cp, comp);
  dijkstra(dp1, s.p[1], cp, comp);
  int nodecount = 0;
  memset(indicators, 0, 8*sizeof(int));
  indicators[6] = INT_MAX;
  indicators[7] = INT_MAX;
  for(int m=1;m<=4;m++) { int d = dp0(s.p[1].next(m)); if(d < indicators[6]) indicators[6] = d+1; }
  for(int j=0;j<M.height;j++)
    for(int i=0;i<M.width;i++) {
      position p(i,j);
      int diff = dp0(i,j) - dp1(i,j);
      // if the opponent's distance is shorter than ours, then this is "their" node
      if(diff>0) {
        nodecount -= degree(p) - potential_articulation(p);
        indicators[3] ++;
        indicators[4] += degree(p);
        indicators[5] += potential_articulation(p);
      }
      // otherwise it's ours
      if(diff<0) {
        nodecount += degree(p) - potential_articulation(p);
        indicators[0] ++;
        indicators[1] += degree(p);
        indicators[2] += potential_articulation(p);
      }
    }
#if VERBOSE >= 2
  if(vis) {
    dp0.dump();
    dp1.dump();
    fprintf(stderr, "player=%d nodecount: %d\n", player, nodecount);
  }
#endif
  Map<int> n(M.width, M.height);
  indicators[7] = num_paths(n, dp0, s.p[0], s.p[1]);
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
  if((comp = cp.component(s.p[0])) == cp.component(s.p[1])) {
    int v = _evaluate_territory(NULL, s, cp, comp, vis);
    return v;
  }

  // since each bot is in a separate component by definition here, it's OK to
  // destructively update cp for floodfill()
  int v = 1000*(floodfill(cp, s.p[0], false) -
                floodfill(cp, s.p[1], false)); // assume everyone else's floodfill is as bad as ours?
//                   cp.connectedarea(s.p[1]));
  if(player == 1) v = -v;
#if VERBOSE >= 2
  if(vis) {
    fprintf(stderr, "player=%d connectedarea value: %d\n", player, v);
  }
#endif
  return v;
}
// }}}

// {{{ alpha-beta iterative deepening search

// do an iterative-deepening search on all moves and see if we can find a move
// sequence that cuts off our opponent
int _alphabeta(int &move, gamestate s, int player, int a, int b, int itr)
{
  // base cases: no more moves?  draws?
  if(s.p[0] == s.p[1]) { return (player == 1 ? 1 : -1) * DRAW_PENALTY; } // crash!  draw!
  if(degree(s.p[player]) == 0) {
    if(degree(s.p[player^1]) == 0) { // both boxed in; draw
      return (player == 1 ? 1 : -1) * DRAW_PENALTY;
    }
    return -INT_MAX;
  }
  if(degree(s.p[player^1]) == 0) {
    // choose any move
    for(int m=1;m<=4;m++) if(!M(s.p[player].next(m))) break;
    return INT_MAX;
  }

  if(_timed_out || ((_ab_runs++)&31) == 0 && timeout()) {
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
  for(int _m=0;_m<=4 && !_timed_out;_m++) {
    // convoluted logic: do "killer heuristic" move first
    if(_m == kill) continue;
    int m = _m == 0 ? kill : _m;
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
    int m_; // next move; discard
    int a_ = -_alphabeta(m_, r, player^1, -b, -a, itr-1);
    if(a_ > a) {
      a = a_;
      move = m;
      _killer[_maxitr-itr] = m;
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
  return a;
}

int next_move_alphabeta()
{
  int itr;
  int lastv = -INT_MAX, lastm = 1;
  reset_timer(firstmove ? FIRSTMOVE_USEC : TIMEOUT_USEC);
  firstmove=false;
  evaluations=0;
  for(itr=DEPTH_INITIAL;itr<DEPTH_MAX && !timeout();itr++) {
    int m;
    _maxitr = itr*2;
    int v = _alphabeta(m, curstate, 0, -INT_MAX, INT_MAX, itr*2);
#if 0
    if(v >= 5000) {
#if VERBOSE >= 1
      struct timeval tv;
      gettimeofday(&tv, NULL);
      fprintf(stderr, "%d.%06d: v=%d best=%d (m=%d) -> found compelling move\n", (int) tv.tv_sec, (int) tv.tv_usec, v, lastv, m);
#endif
      return m;
    }
#endif
#if VERBOSE >= 1
    struct timeval tv;
    gettimeofday(&tv, NULL);
    //M.dump();
    fprintf(stderr, "%d.%06d: v=%d (m=%d) @depth %d _ab_runs=%d\n", (int) tv.tv_sec, (int) tv.tv_usec, v, m, itr*2, _ab_runs);
#endif
    if(v == INT_MAX) // our opponent cannot move, so we win
      return m;
    if(v == -INT_MAX) {
      // deeper searching is apparently impossible (either because there are no
      // more moves for us or because we don't have any search time left)
      break;
    }
    lastv = v;
    lastm = m;
  }
#if VERBOSE >= 1
  long e = elapsed_time();
  float rate = (float)evaluations*1000000.0/(float)e;
  fprintf(stderr, "%d evals in %ld us; %0.1f evals/sec; lastv=%d move=%d\n", evaluations, e, rate, lastv, lastm);
  if(e > TIMEOUT_USEC*11/10) {
    fprintf(stderr, "10%% timeout violation: %ld us\n", e);
  }
#endif
  memmove(_killer, _killer+2, sizeof(_killer)-2); // shift our best-move tree forward to accelerate next move's search
  return lastm;
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
    int indicators[8];
    _evaluate_territory(indicators, curstate, cp, comp, true);
    printf("%d %d %d ", v0, v1, movestoclose);
    for(int i=0;i<8;i++) printf("%d ", indicators[i]);
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
