#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <limits.h>
#include <map>
#include <vector>

#include "artictbl.h"

#define TIMEOUT_USEC 950000
#define FIRSTMOVE_USEC 2950000
#define DEPTH_INITIAL 1
#define DEPTH_MAX 100
#define DRAW_PENALTY (10*itr) // -500
#define VERBOSE 0

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
  position next(int move) { return position(x+dx[move], y+dy[move]); }
  position prev(int move) { return position(x-dx[move], y-dy[move]); }
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
void dijkstra(Map<int> &d, position s, Components &cp, int component)
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

// }}}

// {{{ space-filling

int floodfill(Components &ca, position s)
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
  M(b) = 0; ca.add(b);
  return a;
}

int next_move_spacefill()
{
  Components ca(M);

  // flood fill heuristic: choose to remove as few edges from the graph as
  // possible (in other words, move onto the square with the lowest degree)
  int bestm=1, bestv=0;
  for(int m=1;m<=4;m++) {
    position p = curstate.p[0].next(m);
    if(M(p)) continue;
    int v = ca.connectedvalue(p) + ca.connectedarea(p) - 1 - 2*degree(p) -
      4*potential_articulation(p);
    if(v > bestv) { bestv = v; bestm = m; }
#if VERBOSE >= 1
    fprintf(stderr, "move %d: edges=%d, nodes=%d, degree=%d, v=%d\n", m,
            ca.connectedvalue(p), ca.connectedarea(p), degree(p), v);
#endif
  }

  return bestm;
}
// }}}

// {{{ heuristic board evaluation
static int evaluations=0;
int _evaluate_board(gamestate s, int player, bool vis=false)
{
  // remove players from the board when evaluating connected components,
  // because if a player is separating components he still gets to choose which
  // one to move into.
  M(s.p[0]) = 0; M(s.p[1]) = 0;
  Components cp(M); // pre-move components
  M(s.p[0]) = 1; M(s.p[1]) = 1;

  evaluations++;
#if VERBOSE >= 2
  if(vis) {
    fprintf(stderr, "evaluating board: \n");
    M(s.p[player]) = 2; M(s.p[player^1]) = 3; M.dump();
    M(s.p[0]) = 1; M(s.p[1]) = 1;
  }
#endif
  int comp;
  if((comp = cp.component(s.p[0])) == cp.component(s.p[1])) {
    dijkstra(dp0, s.p[player], cp, comp);
    dijkstra(dp1, s.p[player^1], cp, comp);
#if VERBOSE >= 2
    Map<int> vor(M.width, M.height);
#endif
    int nodecount = 0;
    for(int j=0;j<M.height;j++)
      for(int i=0;i<M.width;i++) {
        position p(i,j);
        int diff = dp0(i,j) - dp1(i,j);
        // if the opponent's distance is shorter than ours, then this is "their" node
        if(diff>0) { nodecount -= degree(p) - potential_articulation(p); }
        // otherwise it's ours
        if(diff<0) { nodecount += degree(p) - potential_articulation(p); }
#if VERBOSE >= 2
        vor(i,j) = diff > 0 ? 2 : diff < 0 ? 1 : 0;
#endif
      }
#if VERBOSE >= 2
    if(vis) {
      dp0.dump();
      dp1.dump();
      vor.dump();
      fprintf(stderr, "player=%d nodecount: %d\n", player, nodecount);
    }
#endif
    return nodecount;
  } else {
    // since each bot is in a separate component by definition here, it's OK to
    // destructively update cp for floodfill()
    int v = 1000000*(floodfill(cp, s.p[0]) -
                     floodfill(cp, s.p[1])); // assume everyone else's floodfill is as bad as ours?
//                   cp.connectedarea(s.p[1]));
    if(player == 1) v = -v;
#if VERBOSE >= 2
    if(vis) {
      fprintf(stderr, "player=%d connectedarea value: %d\n", player, v);
    }
#endif
    return v;
  }
}
// }}}

// {{{ alpha-beta iterative deepening search

// do an iterative-deepening search on all moves and see if we can find a move
// sequence that cuts off our opponent
static char _killer[DEPTH_MAX*2];
static int _maxitr=0;
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

  if(_timed_out || ((_ab_runs++)&127) == 0 && timeout()) {
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

static bool firstmove = true;
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
    for(int m=1;m<=4;m++)
      if(!M(curstate.p[0].next(m)))
        return m;
  }
  if(cp.component(curstate.p[0]) == cp.component(curstate.p[1])) {
    // start-midgame: try to cut off our opponent
    return next_move_alphabeta();
  } else {
    // endgame: use up space as efficiently as we can, and hope we have more
    // left than they do.
    return next_move_spacefill();
  }
}

int main() {
  memset(_killer, 1, sizeof(_killer));
  while (map_update()) {
    printf("%d\n", move_permute[next_move()]);
    fflush(stdout);
  }
//#if VERBOSE >= 1
//  fprintf(stderr, "%d evaluations\n", evaluations);
//#endif
  return 0;
}

// vim: sw=2:ts=8:et:foldmethod=marker
