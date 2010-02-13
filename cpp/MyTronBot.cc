#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/time.h>
#include <map>

// we'll allot .9 seconds to searching
#define TIMEOUT_USEC 9000000

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

const char position::dx[5]={0,0,1,0,-1};
const char position::dy[5]={0,-1,0,1,0};
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
  ~Map() { if(map) delete[] map; }
  T& operator()(position p) { return map[p.x + p.y*width]; }
  T& operator()(int x, int y) { return map[x + y*width]; }
  T& M(position p) { return map[p.x + p.y*width]; }
  T& M(int x, int y) { return map[x + y*width]; }

  void dump(void) {
    for(int j=0;j<height;j++) {
      for(int i=0;i<height;i++) {
        int n = map[i+j*width];
        if(n == 0) fprintf(stderr, "  ");
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

Map<char> M;
gamestate curstate;

// {{{ imported map update garbage from original code
void map_update()
{
  int x, y, c;
  int map_width, map_height;
  int num_items = fscanf(stdin, "%d %d\n", &map_width, &map_height);
  if (feof(stdin) || num_items < 2) {
    exit(0); // End of stream means end of game. Just exit.
  }
  if(!M.map) M.resize(map_width, map_height);
  x = 0;
  y = 0;
  while (y < M.height && (c = fgetc(stdin)) != EOF) {
    switch (c) {
    case '\r':
      break;
    case '\n':
      if (x != M.width) {
	fprintf(stderr, "x != width in Board_ReadFromStream\n");
	return;
      }
      ++y;
      x = 0;
      break;
    case '#':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      M(x,y) = 1;
      ++x;
      break;
    case ' ':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
      }
      M(x,y) = 0;
      ++x;
      break;
    case '1':
    case '2':
      if (x >= M.width) {
	fprintf(stderr, "x >= width in Board_ReadFromStream\n");
	return;
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
      return;
    }
  }
}
// }}}

// {{{ connected components algorithm

struct Components {
  Map<int> c;
  std::map<int,int> csize;

  Components(Map<char> &M): c(M.width, M.height)
  {
    std::map<int,int> equiv;
    int nextclass = 1;
    for(int j=1;j<M.height-1;j++) {
      for(int i=1;i<M.width-1;i++) {
        if(M(i,j)) continue; // wall
        int cup   = c(i, j-1),
            cleft = c(i-1, j);
        if(cup == 0 && cleft == 0) { // new component
          c(i,j) = nextclass++;
        } else if(cup == cleft) { // existing component
          c(i,j) = cup;
        } else { // join components
          // deprecate the higher-numbered component in favor of the lower
          if(cleft == 0 || (cup != 0 && cup < cleft)) {
            c(i,j) = cup;
            if(cleft != 0) equiv[cleft] = cup;
          } else {
            c(i,j) = cleft;
            if(cup != 0) equiv[cup] = cleft;
          }
        }
      }
    }
    // now make another pass to patch up the equivalences
    for(int j=1;j<M.height-1;j++) {
      for(int i=1;i<M.width-1;i++) {
        while(true) {
          std::map<int,int>::iterator e = equiv.find(c(i,j));
          if(e == equiv.end()) break;
          c(i,j) = e->second;
        }
        csize[c(i,j)] ++;
      }
    }
  }
  void dump() {
    std::map<int,int>::iterator i;
    for(i=csize.begin();i!=csize.end();i++) {
      fprintf(stderr, "area %d: %d vertices\n", i->first, i->second);
    }
    c.dump();
  }
  int component(const position &p) { return c(p); }
  int connectedarea(int component) { return csize[component]; }
  int connectedarea(const position &p) { return csize[c(p)]; }
};

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

int turn(int d, int n) { return 1+(d-1+n)&3; }
// }}}

// {{{ run timing
long _get_time()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_usec + tv.tv_sec*1000000;
}

long _timer;
void reset_timer(void) { _timer = _get_time(); }

long elapsed_time()
{
  long t = _get_time();
  return t - _timer;
}

bool timeout() { return elapsed_time() > TIMEOUT_USEC; }
// }}}

// {{{ alpha-beta iterative deepening search

int _evaluate_board(gamestate s, int player)
{
  Components cp(M);
  if(cp.component(curstate.p[0]) == cp.component(curstate.p[1])) {
    // i have no idea how to evaluate this.  let's start with a small positive manhattan distance.
    // we need to come up with something better to figure out who controls the board
    return M.width + M.height - abs(curstate.p[0].x - curstate.p[1].x) - abs(curstate.p[0].y - curstate.p[1].y);
  } else {
    return 100*(cp.connectedarea(curstate.p[player]) - cp.connectedarea(curstate.p[player^1]));
  }
}

// do an iterative-deepening search on all moves and see if we can find a move
// sequence that cuts off our opponent
int _alphabeta(int &move, gamestate s, int player, int a, int b, int itr)
{
  if(s.p[0] == s.p[1]) { return 0; } // crash!  draw!
  if(itr == 0) { return _evaluate_board(s, player); }

  // at higher levels of the tree, periodically check timeout.  if we do time
  // out, give up, we can't do any more work; whatever we found so far will
  // have to do
  if(itr > 3 && timeout()) return a;

  for(int m=1;m<=4;m++) {
    if(M(s.p[player].next(m))) // impossible move?
      continue;
    gamestate r = s;
    r.m[player] = m;
    // after both players 0 and 1 make their moves, the game state updates
    if(player == 1) {
      M(s.p[0]) = 1;
      M(s.p[1]) = 1;
      r.p[0] = s.p[0].next(r.m[0]);
      r.p[1] = s.p[1].next(r.m[1]);
    }
    int m_; // next move; discard
    int a_ = -_alphabeta(m_, r, player^1, -b, -a, itr-1);
    if(a_ > a) {
      a = a_;
      move = m;
    }
    if(a >= b) // beta cut-off
      break;

    // undo game state update
    if(player == 1) {
      r.p[0] = s.p[0];
      r.p[1] = s.p[1];
      M(r.p[0]) = 0;
      M(r.p[1]) = 0;
    }
  }
  return a;
}

int next_move_alphabeta()
{
  int itr = 2;
  int bestv = -1000000, bestm=1;
  reset_timer();
  while(!timeout()) {
    int m;
    int v = _alphabeta(m, curstate, 0, -10000000, 10000000, itr);
    if(v >= 200) return m;
    if(v > bestv) { bestv = v; bestm = m;}
    itr++;
    struct timeval tv;
    gettimeofday(&tv, NULL);
    fprintf(stderr, "%d.%06d: best=%d deepening to %d\n", (int) tv.tv_sec, (int) tv.tv_usec, bestv, itr);
  }
  return bestm;
}
// }}}

// {{{ space-filling

// {{{ crappy greedy space filler that doesn't work as well as just wall-following
struct SpaceFiller
{
  int area, move;

  static SpaceFiller fill(position p, int max_area) {
    SpaceFiller f;
    f.area = _fill_area(f.move, p, max_area, 0);
    return f;
  }

#if 0
  static SpaceFiller greedy(position p) {
    SpaceFiller f;
    f.area = 0;
    for(int m=1;m<=4;m++) {
      int nm;
      int a = _fill_greedy(nm, p, 1);
      if(a > f.area) { f.area = a; f.move = nm; }
    return f;
  }
#endif

private:
#if 0
  int _fill_greedy(int &move, position p0, int majordir) {
    // try to fill in row-major or column-major order, looking ahead one row
    // (or column) for the largest opening
    
  }
  int _find_next_greedy(position p, int dir, int majordir) {
    //  0123456789ab
    //  >----------- ->    ##########...
    //---   ------       ..... 1##### 
    // travel in dir, finding all 0->1 transitions in the next column over (by majordir)
    // for each 0->1 transition, find the run length going -dir
    // optimize for this_runlen + next_runlen
    int ow = M(p.next(majordir));
    int n = 0;
    int besta=0, bestn=0;
    while(!M(p.next(dir))) {
      p = p.next(dir);
      int nw = M(p.next(majordir));
      if(nw && nw != ow) {
        int a = n+runout(p.prev(dir).next(majordir), turn(dir, 2));
        turnarounds.push_back(n);
      }
      n++;
    }
  }
#endif

  // this is exponential and sucks
  static int _fill_area(int &move, position p0, int max_area, int lastmove) {
//    printf("_fill_area([%d,%d], area=%d)\n", p0.x, p0.y, max_area);
    if(max_area <= 1) { move = 0; return 1; }
    int besta = 0, bestm = 0;
    for(int m=1;m<=4;m++) {
      int r = 0;
      position p = p0;
      if(m == lastmove) continue;
      if(M(p.next(m))) continue;
      while(!M(p)) { r++; M(p) = 1; p = p.next(m); } // do!
      while(r>0) {
        // undo! (we just ran into a wall, so undo that first...)
        r--;
        p = p.prev(m);
        M(p) = 0;
        if(r>0) {
          // now try longest to shortest...
          int _m;
//          printf("%d: ", m);
          int a = r + _fill_area(_m, p, max_area - r, m);
          if(a >= max_area) { move = m; return a; }
          if(a > besta) { besta = a; bestm = m; }
        }
      }
    }
    move = bestm;
    return besta;
  }
};
// }}}

const int degreescore[] = {0, 4, 3, 2, 1};
int next_move_spacefill()
{
  M(curstate.p[0]) = 1;
  M(curstate.p[1]) = 1;

  Components ca(M);

  int bestm=1, bestv=0;
  for(int m=1;m<=4;m++) {
    position p = curstate.p[0].next(m);
    if(M(p)) continue;
    int v = ca.connectedarea(p) + degreescore[degree(p)];
    if(v > bestv) { bestv = v; bestm = m; }
    fprintf(stderr, "move %d: ca=%d, degree=%d, v=%d\n", m, ca.connectedarea(p), degree(p), v);
  }
  

//  SpaceFiller s(curstate.p[0], c.connectedarea(curstate.p[0])/2);
//  return SpaceFiller::greedy(curstate.p[0]).move;
  return bestm;
}
// }}}

int next_move() {
  Components cp(M);
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
  setlinebuf(stdout);
  while (true) {
    map_update();
    printf("%d\n", next_move());
  }
  return 0;
}

// vim: sw=2:ts=8:et:foldmethod=marker
